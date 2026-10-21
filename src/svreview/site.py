"""site.py - the viewer's data from a review set: one manifest, mini-BAMs, the informative reads tagged.

The browser reads BAM, not CRAM (README: "CRAM or BAM"), so every member of every candidate becomes a BAM
here, where the reference is at hand: a copy when the bundle's member already is one, a conversion when it
is a CRAM (its own embedded reference, or the trio's / ``--ref``), a cut over the candidate's windows when the
set names per-trio alignments instead of bundles. The child's BAM gets a ``dk:i:1`` tag on every alignment
of a read the run's informative BAM holds, so the viewer can light the k-mer evidence up without a second
file. The manifest joins it all: the candidate rows, the windows, the members, the run's region metrics and
the non-human fraction, and where everything is, relative to the manifest::

  <out>/manifest.json
  <out>/<trio>/<chrom>_<start>_<end>/{child,father,mother}.bam (+ .bai)
  <out>/candidates.tsv  (and trios.tsv)   the tables, as given (the verdict merge reads the first back)
  <out>/verdicts.tsv                      written by `svreview serve` as verdicts come in (not by this)

With ``--with-viewer`` the viewer's own files (index.html, js/, css/) are copied beside the manifest, so the
directory is a complete static site: `svreview serve <out>` or any web server.
"""

import collections
import datetime
import json
import os
import shutil

import pysam

from . import CONTRACT_VERSION, __version__, bundle, contract, review_set

ROLES = contract.ROLES
VIEWER_FILES = ("index.html", "js", "css")


def convert_member(src, dst, ref=None, dk_names=None):
    """*src* (BAM or CRAM) as a BAM at *dst*, indexed; alignments of the reads in *dk_names* (a set of
    ``(name, read_number)``) get ``dk:i:1``. Returns ``(records, tagged)``."""
    n = tagged = 0
    with review_set.open_alignments(src, ref) as f:
        with pysam.AlignmentFile(dst + ".tmp", "wb", template=f) as o:
            for r in f:
                if dk_names is not None:
                    which = 1 if r.is_read1 else (2 if r.is_read2 else 0)
                    if (r.query_name, which) in dk_names:
                        r.set_tag("dk", 1, "i")
                        tagged += 1
                o.write(r)
                n += 1
    os.replace(dst + ".tmp", dst)
    pysam.index(dst)
    return n, tagged


def informative_names(informative_bam, windows):
    """``{(name, read_number)}`` of the informative BAM's reads placed in *windows* (0-based half-open), the
    unmapped ones included (they are placed at their mates)."""
    names = set()
    with pysam.AlignmentFile(informative_bam, "rb") as f:
        for chrom, start, end in windows:
            if chrom not in f.references:
                continue
            for r in f.fetch(chrom, start, end):
                names.add((r.query_name, 1 if r.is_read1 else (2 if r.is_read2 else 0)))
    return names


def _rel(path, base):
    return os.path.relpath(path, base).replace(os.sep, "/")


def _sample_of(path, role, ref=None):
    try:
        with review_set.open_alignments(path, ref) as f:
            rg = (f.header.to_dict().get("RG") or [{}])[0]
            return rg.get("SM", role)
    except contract.ContractError:
        return role


def _count_tagged(path):
    with pysam.AlignmentFile(path, "rb") as f:
        return sum(1 for r in f.fetch(until_eof=True) if r.has_tag("dk"))


def candidate_entry(cand):
    entry = collections.OrderedDict()
    entry["id"] = cand.id
    entry["trio"] = cand.trio
    entry["chrom"], entry["start"], entry["end"] = cand.chrom, cand.start, cand.end
    entry["class"], entry["sv_type"] = cand.cls, cand.sv_type
    entry["evidence"] = cand.evidence
    entry["junctions"] = [collections.OrderedDict(chrom=j.chrom, start=j.start, end=j.end, strands=j.strands,
                                                  sv_type=j.sv_type, reads=j.reads) for j in cand.junctions]
    entry["locus_flags"], entry["partner_flags"], entry["flags"] = cand.locus_flags, cand.partner_flags, cand.flags
    entry["review_tier"], entry["rank"], entry["both_ends"] = cand.review_tier, cand.rank, cand.both_ends
    entry["nhf_reads"], entry["nhf"] = cand.nhf_reads, cand.nhf
    entry["row"] = cand.row
    return entry


def build_site_from_set(rs, out, ref=None, title=None, with_viewer=False, viewer_dir=None, threads=1, log=print):
    """Write the site's data under *out* from a ReviewSet (opened from a directory, or ad hoc). Returns the
    manifest (a dict)."""
    out = os.path.abspath(out)
    os.makedirs(out, exist_ok=True)
    shutil.copyfile(rs.candidates_path, os.path.join(out, "candidates.tsv"))
    if rs.path and os.path.isfile(os.path.join(rs.path, "trios.tsv")):
        shutil.copyfile(os.path.join(rs.path, "trios.tsv"), os.path.join(out, "trios.tsv"))
    manifest = collections.OrderedDict()
    manifest["svreview"] = __version__
    manifest["contract"] = CONTRACT_VERSION
    manifest["format"] = review_set.FORMAT_VERSION
    manifest["title"] = title or rs.title or os.path.basename(rs.path or os.path.dirname(rs.candidates_path)) or "svreview"
    manifest["created"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    manifest["source"] = collections.OrderedDict(set=rs.path, candidates=rs.candidates_path, bundles=rs.bundles_dir)
    manifest["columns"] = rs.header
    manifest["trios"] = collections.OrderedDict()
    manifest["candidates"] = []
    converted = {}          # (trio, windows) -> {role: (dst, records, tagged, sample, bytes)} or the first's key
    for cand in rs.candidates:
        entry = candidate_entry(cand)
        try:
            m = rs.metrics_of(cand)
        except contract.ContractError as e:
            log("%s: metrics not read: %s" % (cand.id, e))
            m = None
        trio = rs.trio_of(cand)
        if cand.trio not in manifest["trios"]:
            t = collections.OrderedDict(members=collections.OrderedDict(trio.samples) if trio else collections.OrderedDict())
            t["child_sex"] = trio.child_sex if trio else "0"
            t["nhf"] = contract.nhf_block(m) if m else collections.OrderedDict(status="none", reason="no metrics given")
            t["informative_bam"] = bool(rs.informative_of(cand))
            manifest["trios"][cand.trio] = t
        regions = []
        for chrom, start, end in [(cand.chrom, cand.start, cand.end)] + cand.partners():
            r = collections.OrderedDict(chrom=chrom, start=start, end=end,
                                        own=(chrom, start, end) == (cand.chrom, cand.start, cand.end))
            mr = contract.metrics_region(m, chrom, start, end) if m else None
            if mr:
                for k, v in mr.items():
                    if k not in ("chrom", "start", "end"):
                        r[k] = v
            r["nhf_screen"] = contract.nhf_region(m, chrom, start, end) if m else None
            regions.append(r)
        entry["regions"] = regions
        kind, src = rs.source_of(cand)
        cand_ref = rs.reference_of(cand) or ref
        if kind is None:
            entry["source"], entry["bundle"] = None, None
            entry["windows"] = [collections.OrderedDict(chrom=c, start=s, end=e)
                                for c, s, e in contract.candidate_windows(cand, rs.pad)]
            entry["members"] = collections.OrderedDict()
            entry["missing"] = "no bundle under %s and no per-trio alignments" % (rs.bundles_dir or "bundles/")
            manifest["candidates"].append(entry)
            log("%s: nothing to draw (no bundle, no per-trio alignments)" % cand.id)
            continue
        try:
            windows, cj = rs.windows_of(cand)
        except contract.ContractError as e:
            raise contract.ContractError("%s: %s" % (cand.id, e))
        entry["source"] = kind
        entry["windows"] = [collections.OrderedDict(chrom=c, start=s, end=e) for c, s, e in windows]
        entry["bundle"] = _rel(src, rs.bundles_dir) if kind == "bundle" else "trio:" + cand.trio
        entry["shared_with"] = None
        if cj:
            entry["candidate_json"] = collections.OrderedDict((k, v) for k, v in cj.items() if k != "windows0")
            if cj.get("candidate") and cj.get("candidate") != cand.locus:
                entry["shared_with"] = cj.get("candidate")
        key = (cand.trio, tuple(windows))
        if key in converted and isinstance(converted[key], tuple):
            key = converted[key]
        if key not in converted:
            dest = os.path.join(out, cand.trio, os.path.basename(src) if kind == "bundle" else cand.name)
            os.makedirs(dest, exist_ok=True)
            inf = rs.informative_of(cand)
            dk = informative_names(inf, windows) if inf else None
            got = collections.OrderedDict()
            members = contract.bundle_members(src) if kind == "bundle" else src
            for role in ROLES:
                if role not in members or not members[role]:
                    got[role] = None
                    continue
                dst = os.path.join(dest, role + ".bam")
                names = dk if role == "child" else None
                if kind == "bundle":
                    n, tagged = convert_member(members[role], dst, cand_ref, names)
                else:
                    n = bundle.cut_member(members[role], dst, windows, cand_ref, threads, dk_names=names)
                    tagged = _count_tagged(dst) if names else 0
                sample = (trio.samples[role] if trio else None) or _sample_of(members[role], role, cand_ref)
                got[role] = (dst, n, tagged, sample, os.path.getsize(dst))
            converted[key] = got
            log("%s: %s" % (cand.id, "; ".join("%s %d reads%s" % (r, v[1], (" (%d informative)" % v[2]) if r == "child" and inf else "")
                                                for r, v in got.items() if v)))
        else:
            first = next(e for e in manifest["candidates"] if e.get("members") and e["trio"] == cand.trio
                         and e["members"].get("child") and e["members"]["child"]["bam"] == _rel(converted[key]["child"][0], out))
            if entry["shared_with"] is None and first["id"] != cand.id:
                entry["shared_with"] = first["chrom"] + ":%d-%d" % (first["start"], first["end"])
            log("%s: shares the windows of %s" % (cand.id, first["id"]))
        entry["members"] = collections.OrderedDict()
        for role, v in converted[key].items():
            if v is None:
                entry["members"][role] = None
                continue
            dst, n, tagged, sample, size = v
            entry["members"][role] = collections.OrderedDict(sample=sample, bam=_rel(dst, out), reads=n, bytes=size,
                                                             informative=(tagged if role == "child" else None))
            manifest["trios"][cand.trio]["members"].setdefault(role, sample)
        manifest["candidates"].append(entry)
    manifest["summary"] = collections.OrderedDict(
        candidates=len(rs.candidates), with_bundle=sum(1 for e in manifest["candidates"] if e.get("bundle")),
        tiers=collections.OrderedDict(sorted(collections.Counter(
            str(e["review_tier"]) for e in manifest["candidates"]).items())))
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    if with_viewer:
        copy_viewer(out, viewer_dir)
    log("manifest: %s (%d candidates, %d with reads)" % (os.path.join(out, "manifest.json"), len(rs.candidates),
                                                        manifest["summary"]["with_bundle"]))
    return manifest


def build_site(candidates_path, bundles, out, metrics=None, informative=None, ref=None, pad=10000, title=None,
               with_viewer=False, viewer_dir=None, trios=None, log=print):
    """The ad hoc form: a candidate table, a bundle directory, metrics and informative BAMs per trio
    (``{trio: path}``; a trio named ``*`` applies to every trio), an optional trios.tsv."""
    rs = review_set.ReviewSet.adhoc(candidates_path, bundles, metrics=metrics, informative=informative, pad=pad,
                                    title=title, trios_path=trios)
    return build_site_from_set(rs, out, ref=ref, title=title, with_viewer=with_viewer, viewer_dir=viewer_dir, log=log)


def default_viewer_dir():
    """The viewer shipped with the package (docs/ of the repository, or the copy inside the package)."""
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (os.path.join(here, "viewer"), os.path.join(here, os.pardir, os.pardir, "docs")):
        if os.path.isfile(os.path.join(d, "index.html")):
            return os.path.normpath(d)
    raise contract.ContractError("the viewer's files were not found beside the package; give --viewer-dir")


def copy_viewer(out, viewer_dir=None):
    src = viewer_dir or default_viewer_dir()
    for name in VIEWER_FILES:
        s, d = os.path.join(src, name), os.path.join(out, name)
        if os.path.isdir(s):
            if os.path.isdir(d):
                shutil.rmtree(d)
            shutil.copytree(s, d)
        elif os.path.isfile(s):
            shutil.copyfile(s, d)
    return out
