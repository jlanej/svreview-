"""site.py - the viewer's data from the analysis' files: one manifest, mini-BAMs, the informative reads tagged.

The browser reads BAM, not CRAM (README: "CRAM or BAM"), so every member of every bundle becomes a BAM here,
where the reference is at hand - a copy when it already is one, a conversion when it is a CRAM (its own
embedded reference, or ``--ref``). The child's BAM gets a ``dk:i:1`` tag on every alignment of a read the
run's informative BAM holds (``--informative``), so the viewer can light the k-mer evidence up without a
second file. The manifest joins it all: the candidate rows, the windows that were cut, the members, the
run's region metrics and the non-human fraction, and where everything is, relative to the manifest::

  <out>/manifest.json
  <out>/<trio>/<chrom>_<start>_<end>/{child,father,mother}.bam (+ .bai)
  <out>/candidates.tsv                       the table, as given (the verdict merge reads it back)
  <out>/verdicts.tsv                         written by `svreview serve` as verdicts come in (not by this)

With ``--with-viewer`` the viewer's own files (index.html, js/, css/) are copied beside the manifest, so the
directory is a complete static site: `svreview serve <out>` or any web server.
"""

import collections
import datetime
import json
import os
import shutil

import pysam

from . import CONTRACT_VERSION, __version__, contract

ROLES = contract.ROLES
VIEWER_FILES = ("index.html", "js", "css")


def _open_alignments(path, ref=None):
    if path.endswith(".cram"):
        try:
            return pysam.AlignmentFile(path, "rc", reference_filename=ref) if ref else pysam.AlignmentFile(path, "rc")
        except (OSError, ValueError) as e:
            raise contract.ContractError("%s: cannot open the CRAM%s: %s" % (
                path, "" if ref else " (its reference is not embedded: give --ref)", e))
    return pysam.AlignmentFile(path, "rb")


def convert_member(src, dst, ref=None, dk_names=None):
    """*src* (BAM or CRAM) as a BAM at *dst*, indexed; alignments of the reads in *dk_names* (a set of
    ``(name, read_number)``) get ``dk:i:1``. Returns ``(records, tagged)``."""
    n = tagged = 0
    with _open_alignments(src, ref) as f:
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
        with _open_alignments(path, ref) as f:
            rg = (f.header.to_dict().get("RG") or [{}])[0]
            return rg.get("SM", role)
    except contract.ContractError:
        return role


def _lengths_of(path, ref=None):
    with _open_alignments(path, ref) as f:
        return collections.OrderedDict(zip(f.references, f.lengths))


def build_site(candidates_path, bundles, out, metrics=None, informative=None, ref=None, pad=10000, title=None,
               with_viewer=False, viewer_dir=None, log=print):
    """Write the site's data under *out*. *metrics* and *informative* map a trio to its ``metrics.json`` /
    informative BAM (a trio named ``*`` applies to every trio). Returns the manifest (a dict)."""
    header, cands = contract.read_candidates(candidates_path)
    metrics = dict(metrics or {})
    informative = dict(informative or {})
    out = os.path.abspath(out)
    os.makedirs(out, exist_ok=True)
    shutil.copyfile(candidates_path, os.path.join(out, "candidates.tsv"))
    loaded_metrics = {}

    def metrics_of(trio):
        path = metrics.get(trio) or metrics.get("*")
        if not path:
            return None
        if path not in loaded_metrics:
            loaded_metrics[path] = contract.read_metrics(path)
        return loaded_metrics[path]

    def informative_of(trio):
        return informative.get(trio) or informative.get("*")

    manifest = collections.OrderedDict()
    manifest["svreview"] = __version__
    manifest["contract"] = CONTRACT_VERSION
    manifest["title"] = title or os.path.basename(os.path.dirname(os.path.abspath(candidates_path))) or "svreview"
    manifest["created"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    manifest["source"] = collections.OrderedDict(candidates=os.path.abspath(candidates_path),
                                                 bundles=os.path.abspath(bundles) if bundles else None)
    manifest["columns"] = header
    manifest["trios"] = collections.OrderedDict()
    manifest["candidates"] = []
    converted = {}          # bundle dir -> {role: (relative path, records, tagged, sample)}
    for cand in cands:
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
        m = metrics_of(cand.trio)
        if cand.trio not in manifest["trios"]:
            t = collections.OrderedDict(members=collections.OrderedDict())
            t["nhf"] = contract.nhf_block(m) if m else collections.OrderedDict(status="none", reason="no metrics given")
            t["informative_bam"] = bool(informative_of(cand.trio))
            manifest["trios"][cand.trio] = t
        regions = []
        for chrom, start, end in [(cand.chrom, cand.start, cand.end)] + cand.partners():
            r = collections.OrderedDict(chrom=chrom, start=start, end=end, own=(chrom, start, end) == (cand.chrom, cand.start, cand.end))
            mr = contract.metrics_region(m, chrom, start, end) if m else None
            if mr:
                for k, v in mr.items():
                    if k not in ("chrom", "start", "end"):
                        r[k] = v
            nr = contract.nhf_region(m, chrom, start, end) if m else None
            r["nhf_screen"] = nr
            regions.append(r)
        entry["regions"] = regions
        bd = contract.bundle_dir(bundles, cand) if bundles else None
        if bd is None:
            entry["bundle"] = None
            entry["windows"] = [collections.OrderedDict(chrom=c, start=s, end=e)
                                for c, s, e in contract.candidate_windows(cand, pad)]
            entry["members"] = collections.OrderedDict()
            entry["missing"] = "no bundle under %s" % (bundles or "(none given)")
            manifest["candidates"].append(entry)
            log("%s: no bundle" % cand.id)
            continue
        cj_path = os.path.join(bd, "candidate.json")
        cj = contract.read_candidate_json(cj_path) if os.path.isfile(cj_path) else None
        members = contract.bundle_members(bd)
        lengths = _lengths_of(members[next(iter(members))], ref) if members else None
        if cj and cj["windows0"]:
            windows = cj["windows0"]
        else:
            windows = contract.candidate_windows(cand, (cj or {}).get("pad", pad), lengths)
        entry["windows"] = [collections.OrderedDict(chrom=c, start=s, end=e) for c, s, e in windows]
        entry["bundle"] = _rel(bd, os.path.abspath(bundles))
        entry["shared_with"] = None
        if cj:
            entry["candidate_json"] = collections.OrderedDict((k, v) for k, v in cj.items() if k != "windows0")
            if cj.get("candidate") != cand.locus:
                entry["shared_with"] = cj.get("candidate")
        dest = os.path.join(out, cand.trio, os.path.basename(bd))
        if bd not in converted:
            os.makedirs(dest, exist_ok=True)
            got = collections.OrderedDict()
            inf = informative_of(cand.trio)
            dk = informative_names(inf, windows) if inf else None
            for role in ROLES:
                if role not in members:
                    got[role] = None
                    continue
                dst = os.path.join(dest, role + ".bam")
                n, tagged = convert_member(members[role], dst, ref, dk if role == "child" else None)
                got[role] = (dst, n, tagged, _sample_of(members[role], role, ref), os.path.getsize(dst))
            converted[bd] = got
            log("%s: %s" % (cand.id, "; ".join("%s %d reads%s" % (r, v[1], (" (%d informative)" % v[2]) if r == "child" and inf else "")
                                                for r, v in got.items() if v)))
        else:
            log("%s: shares %s" % (cand.id, _rel(bd, os.path.abspath(bundles))))
        entry["members"] = collections.OrderedDict()
        for role, v in converted[bd].items():
            if v is None:
                entry["members"][role] = None
                continue
            dst, n, tagged, sample, size = v
            entry["members"][role] = collections.OrderedDict(sample=sample, bam=_rel(dst, out), reads=n, bytes=size,
                                                             informative=(tagged if role == "child" else None))
            manifest["trios"][cand.trio]["members"].setdefault(role, sample)
        manifest["candidates"].append(entry)
    manifest["summary"] = collections.OrderedDict(
        candidates=len(cands), with_bundle=sum(1 for e in manifest["candidates"] if e.get("bundle")),
        tiers=collections.OrderedDict(sorted(collections.Counter(
            str(e["review_tier"]) for e in manifest["candidates"]).items())))
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    if with_viewer:
        copy_viewer(out, viewer_dir)
    log("manifest: %s (%d candidates, %d with a bundle)" % (os.path.join(out, "manifest.json"), len(cands),
                                                           manifest["summary"]["with_bundle"]))
    return manifest


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
