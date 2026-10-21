"""bundle.py - cut the review bundle of every candidate from the trio's alignments, in the contract's layout.

The analysis cuts its bundles where the whole CRAMs are (bundle_candidates.py: local copies or HTTP range
reads of the public 1000 Genomes files; trio_post.py: on the cluster, before the CRAMs leave scratch). This
is the same cut for alignments at hand - the demo's BAMs, a bench trio - so that a run of one trio can be
reviewed without the cohort machinery:

  <out>/<trio>/<chrom>_<start>_<end>/     one directory per candidate (the table's own window)
    child.bam father.bam mother.bam (+ .bai)   (or .cram + .crai) the reads over the candidate's windows (own region and every
                                           junction partner, each padded by --pad, merged where they touch)
                                           and their mates wherever they map
    <member>.json                          what was cut (sample, source, windows, reads, bytes)
    candidate.json                         the candidate (type, tier, rank, its windows) and who shares the bundle
  <out>/bundle.tsv                         the index: candidate, trio, sv_type, review_tier, rank, member, sample,
                                           region, reads, bytes, source_url, status, path

The other end of a tier-1 pair has the same windows: its directory is a symlink to the first's (status
``shared=<dir>``), as the analysis does it. Members are written as BAM (what the viewer reads) or, with
``--cram``, as CRAM with the reference embedded (what the analysis' bundles are; ``svreview site`` converts
back); a CRAM source needs its reference (``--ref``) unless the reference is embedded.
"""

import collections
import datetime
import json
import os

import pysam

from . import contract

ROLES = contract.ROLES


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _open(path, ref=None):
    mode = "rc" if path.endswith(".cram") else "rb"
    kw = {}
    if ref:
        kw["reference_filename"] = ref
    return pysam.AlignmentFile(path, mode, **kw)


def contig_lengths(path, ref=None):
    with _open(path, ref) as f:
        return collections.OrderedDict(zip(f.references, f.lengths))


def cut_member(src, out, windows, ref=None, threads=1, dk_names=None, cram=False):
    """The reads overlapping *windows* (0-based half-open) and their mates wherever they map, as a sorted,
    indexed BAM (or, with *cram*, a CRAM with the reference embedded: needs *ref*); returns the record
    count. A mate outside every window is fetched at the mate's position, which is what
    ``samtools view -P`` does. Alignments of the reads in *dk_names* (``{(name, read_number)}``) get ``dk:i:1``."""
    tmp = out + ".unsorted.bam"
    n = 0

    def write(o, r):
        if dk_names is not None and (r.query_name, 1 if r.is_read1 else (2 if r.is_read2 else 0)) in dk_names:
            r.set_tag("dk", 1, "i")
        o.write(r)

    with _open(src, ref) as f:
        with pysam.AlignmentFile(tmp, "wb", template=f) as o:
            seen = set()
            mates = []
            for chrom, start, end in windows:
                if chrom not in f.references:
                    continue
                for r in f.fetch(chrom, start, end):
                    key = (r.query_name, r.flag & 0xC0, r.reference_id, r.reference_start, r.flag & 0x900, r.cigarstring)
                    if key in seen:
                        continue
                    seen.add(key)
                    write(o, r)
                    n += 1
                    if (r.is_paired and not r.mate_is_unmapped and r.next_reference_id >= 0
                            and not (r.is_secondary or r.is_supplementary)):
                        mchrom = f.get_reference_name(r.next_reference_id)
                        inside = any(c == mchrom and s <= r.next_reference_start < e for c, s, e in windows)
                        if not inside:
                            mates.append((r.query_name, r.flag & 0xC0, mchrom, r.next_reference_start))
            # the mates outside the windows: one fetch at each mate's position, matched by name and read number
            wanted = collections.defaultdict(set)
            for name, which, mchrom, mpos in mates:
                wanted[(mchrom, mpos)].add((name, which ^ 0xC0 if which in (0x40, 0x80) else which))
            for (mchrom, mpos), names in wanted.items():
                for r in f.fetch(mchrom, mpos, mpos + 1):
                    if r.reference_start != mpos or (r.query_name, r.flag & 0xC0) not in names:
                        continue
                    key = (r.query_name, r.flag & 0xC0, r.reference_id, r.reference_start, r.flag & 0x900, r.cigarstring)
                    if key in seen:
                        continue
                    seen.add(key)
                    write(o, r)
                    n += 1
    if cram:
        if not ref:
            raise contract.ContractError("a CRAM bundle needs the reference (--ref) to embed")
        sorted_bam = out + ".sorted.bam"
        pysam.sort("-@", str(max(1, threads)), "-o", sorted_bam, tmp)
        os.remove(tmp)
        pysam.view("-C", "-T", ref, "--output-fmt-option", "embed_ref=1", "-o", out, sorted_bam, catch_stdout=False)
        os.remove(sorted_bam)
    else:
        pysam.sort("-@", str(max(1, threads)), "-o", out, tmp)
        os.remove(tmp)
    pysam.index(out)
    return n


def bundle_candidates(candidates, members, out_dir, pad=10000, ref=None, samples=None, threads=1, log=print, cram=False):
    """Cut one bundle per candidate of *candidates* (``[Candidate]``) from *members*
    (``{role: path}``); *samples* names each member (``{role: sample}``; default: the file's SM or the role).
    With *cram*, the members are written as CRAMs with the reference embedded (as the analysis' bundles are).
    Returns the number of bundles written and the index path."""
    out_dir = os.path.abspath(out_dir)
    lengths = None
    names = dict(samples or {})
    for role in ROLES:
        if role not in members:
            raise contract.ContractError("no %s alignment given" % role)
        if lengths is None:
            lengths = contig_lengths(members[role], ref)
        if role not in names:
            with _open(members[role], ref) as f:
                rg = (f.header.to_dict().get("RG") or [{}])[0]
                names[role] = rg.get("SM", role)
    os.makedirs(out_dir, exist_ok=True)
    primaries = collections.OrderedDict()
    plan = []
    for cand in candidates:
        windows = contract.candidate_windows(cand, pad, lengths)
        windows = [w for w in windows if w[0] in lengths]
        d = os.path.join(out_dir, cand.trio, cand.name)
        key = (cand.trio, tuple(windows))
        shared_with = primaries[key]["name"] if key in primaries else None
        entry = collections.OrderedDict(cand=cand, windows=windows, dir=d, name=cand.name, shared_with=shared_with,
                                        shared_by=[])
        if shared_with:
            primaries[key]["shared_by"].append(collections.OrderedDict(
                candidate=cand.locus, name=cand.name, sv_type=cand.sv_type,
                review_tier=cand.review_tier if cand.review_tier is not None else contract.MISSING,
                rank=cand.rank if cand.rank is not None else contract.MISSING))
        else:
            primaries[key] = entry
        plan.append(entry)
    n_written = 0
    for entry in plan:
        cand = entry["cand"]
        if entry["shared_with"]:
            os.makedirs(os.path.dirname(entry["dir"]), exist_ok=True)
            if not os.path.lexists(entry["dir"]):
                os.symlink(entry["shared_with"], entry["dir"])
            log("%s %s: shared with %s (the same windows)" % (cand.trio, cand.name, entry["shared_with"]))
            continue
        os.makedirs(entry["dir"], exist_ok=True)
        status = collections.OrderedDict()
        parts = []
        for role in ROLES:
            out = os.path.join(entry["dir"], role + (".cram" if cram else ".bam"))
            for stale in (os.path.join(entry["dir"], role + ext) for ext in (".bam", ".bam.bai", ".cram", ".cram.crai")):
                if os.path.exists(stale):
                    os.remove(stale)
            n = cut_member(members[role], out, entry["windows"], ref, threads, cram=cram)
            info = collections.OrderedDict()
            info["member"], info["sample"], info["source"], info["source_url"] = role, names[role], "local", members[role]
            info["regions"] = [contract.region_str(c, s + 1, e) for c, s, e in entry["windows"]]
            info["reads"], info["bytes"] = n, os.path.getsize(out)
            info["fetched"] = now()
            with open(os.path.join(entry["dir"], role + ".json"), "w") as fh:
                json.dump(info, fh, indent=1)
            status[role] = "ok"
            parts.append("%s %d reads" % (role, n))
        cj = collections.OrderedDict()
        cj["trio"], cj["candidate"] = cand.trio, cand.locus
        cj["chrom"], cj["start"], cj["end"] = cand.chrom, cand.start, cand.end
        cj["sv_type"] = cand.sv_type
        cj["review_tier"] = cand.review_tier if cand.review_tier is not None else contract.MISSING
        cj["rank"] = cand.rank if cand.rank is not None else contract.MISSING
        cj["windows"] = [contract.region_str(c, s + 1, e) for c, s, e in entry["windows"]]
        cj["skipped_windows"] = []
        cj["pad"] = pad
        cj["members"] = collections.OrderedDict((r, names[r]) for r in ROLES)
        cj["status"] = status
        cj["shared_by"] = entry["shared_by"]
        cj["written"] = now()
        with open(os.path.join(entry["dir"], "candidate.json"), "w") as fh:
            json.dump(cj, fh, indent=1)
        n_written += 1
        log("%s %s %s tier %s rank %s: %s; windows %s" % (
            cand.trio, cand.name, cand.sv_type, cj["review_tier"], cj["rank"], "; ".join(parts), ",".join(cj["windows"])))
    index = write_index(out_dir)
    return n_written, index


def write_index(out_dir):
    """bundle.tsv rebuilt from what is on disk (bundle_candidates.py's layout and order)."""
    rows = []
    for trio in sorted(os.listdir(out_dir)):
        tdir = os.path.join(out_dir, trio)
        if trio.startswith(".") or not os.path.isdir(tdir):
            continue
        for name in sorted(os.listdir(tdir)):
            d = os.path.join(tdir, name)
            if os.path.islink(d):
                target = os.path.basename(os.readlink(d))
                cj = os.path.join(tdir, target, "candidate.json")
                if not os.path.isfile(cj):
                    continue
                with open(cj) as fh:
                    info = json.load(fh)
                meta = dict(info)
                for s in info.get("shared_by", []):
                    if s.get("name") == name:
                        meta.update(s)
                rows += _member_rows(out_dir, trio, meta.get("candidate", name), target, meta, "shared=%s " % target)
            elif os.path.isdir(d) and os.path.isfile(os.path.join(d, "candidate.json")):
                with open(os.path.join(d, "candidate.json")) as fh:
                    info = json.load(fh)
                rows += _member_rows(out_dir, trio, info.get("candidate", name), name, info, "")
    rows.sort(key=lambda r: (contract.to_int(r["review_tier"], 9), contract.to_int(r["rank"], 10 ** 9), r["trio"],
                             r["candidate"], ROLES.index(r["member"])))
    path = os.path.join(out_dir, "bundle.tsv")
    contract.write_tsv(path, list(contract.BUNDLE_INDEX_COLUMNS), rows)
    return path


def _member_rows(out_dir, trio, candidate, target, meta, status_prefix):
    rows = []
    tdir = os.path.join(out_dir, trio, target)
    for role in ROLES:
        row = collections.OrderedDict((c, contract.MISSING) for c in contract.BUNDLE_INDEX_COLUMNS)
        row["candidate"], row["trio"], row["member"] = candidate, trio, role
        row["sv_type"] = meta.get("sv_type", contract.MISSING)
        row["review_tier"], row["rank"] = meta.get("review_tier", contract.MISSING), meta.get("rank", contract.MISSING)
        members = contract.bundle_members(tdir)
        ext = os.path.splitext(members[role])[1] if role in members else ".bam"
        row["path"] = "%s/%s/%s%s" % (trio, target, role, ext)
        marker = os.path.join(tdir, role + ".json")
        if os.path.isfile(marker) and role in members:
            with open(marker) as fh:
                info = json.load(fh)
            row["sample"], row["region"] = info.get("sample", contract.MISSING), ",".join(info.get("regions", []))
            row["reads"], row["bytes"] = info.get("reads", contract.MISSING), info.get("bytes", contract.MISSING)
            row["source_url"] = info.get("source_url", contract.MISSING)
            row["status"] = status_prefix + "ok"
        elif os.path.isfile(os.path.join(tdir, role + ".failed")):
            with open(os.path.join(tdir, role + ".failed")) as fh:
                row["status"] = status_prefix + "failed: " + fh.read().strip().split("\t", 1)[-1]
        else:
            row["status"] = status_prefix + "absent"
        if row["sample"] == contract.MISSING:
            row["sample"] = meta.get("members", {}).get(role, contract.MISSING)
        rows.append(row)
    return rows
