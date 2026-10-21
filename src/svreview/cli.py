"""cli.py - `svreview`: check a review set, build its site, serve it, read the verdicts back; make a set from one run."""

import argparse
import os
import sys

from . import __version__, contract


def _kv_paths(values):
    """``["T1=/a.json", "/b.json"]`` -> ``{"T1": "/a.json", "*": "/b.json"}``."""
    out = {}
    for v in values or []:
        if "=" in v and not os.path.exists(v):
            k, p = v.split("=", 1)
            out[k] = p
        else:
            out["*"] = v
    return out


def cmd_check(a):
    from . import review_set
    rs = review_set.ReviewSet.open(a.set)
    if a.ref:
        rs.reference = rs.reference or a.ref
    text, n_err = review_set.format_findings(rs.check(deep=a.deep))
    print(text)
    return 1 if n_err else 0


def cmd_candidates(a):
    from . import candidates
    header, rows = candidates.candidates_from_run(
        a.run, a.trio, classes=tuple(a.classes.split(",")), min_reads=a.min_reads, child_sex=a.sex,
        nhf_status=a.nhf_status, sink_multiplicity=a.sink_multiplicity, nhf_flag=a.nhf_flag,
        nhf_min_reads=a.nhf_min_reads)
    contract.write_tsv(a.out, header, rows)
    tiers = {}
    for r in rows:
        tiers[r["review_tier"]] = tiers.get(r["review_tier"], 0) + 1
    print("%d candidates (%s) -> %s" % (len(rows), ", ".join("tier %s: %d" % kv for kv in sorted(tiers.items())), a.out))


def cmd_bundle(a):
    from . import bundle
    _header, cands = contract.read_candidates(a.candidates)
    if a.tier:
        tiers = set(int(t) for t in a.tier.split(","))
        cands = [c for c in cands if c.review_tier in tiers]
    members = {"child": a.child, "father": a.father, "mother": a.mother}
    samples = dict(s.split("=", 1) for s in (a.sample or []))
    n, index = bundle.bundle_candidates(cands, members, a.out_dir, pad=a.pad, ref=a.ref, samples=samples,
                                        threads=a.threads, cram=a.cram)
    print("%d bundles for %d candidates -> %s" % (n, len(cands), index))


def cmd_site(a):
    from . import review_set, site
    if a.set:
        if a.candidates or a.bundles:
            sys.exit("svreview site: give a review set directory, or --candidates/--bundles, not both")
        rs = review_set.ReviewSet.open(a.set)
        for k, v in _kv_paths(a.metrics).items():
            rs.metrics_map[k] = v
        for k, v in _kv_paths(a.informative).items():
            rs.informative_map[k] = v
        if a.pad is not None:
            rs.pad = a.pad
    else:
        if not a.candidates:
            sys.exit("svreview site: give a review set directory, or --candidates (and --bundles)")
        rs = review_set.ReviewSet.adhoc(a.candidates, a.bundles, metrics=_kv_paths(a.metrics),
                                        informative=_kv_paths(a.informative), pad=a.pad or 10000, title=a.title,
                                        trios_path=a.trios)
    site.build_site_from_set(rs, a.out, ref=a.ref, title=a.title, with_viewer=a.with_viewer, viewer_dir=a.viewer_dir,
                             threads=a.threads)


def cmd_serve(a):
    from . import serve
    serve.serve(a.directory, port=a.port, host=a.host, quiet=a.quiet, open_browser=a.open, path=a.path)


def cmd_verdicts(a):
    verdicts = contract.read_verdicts(a.verdicts)
    print("%s: %d candidates with a verdict" % (a.verdicts, len(verdicts)))
    if a.out:
        header, rows = contract.read_tsv(a.candidates)
        out_header, out_rows, n = contract.merge_verdicts(header, rows, verdicts)
        contract.write_tsv(a.out, out_header, out_rows)
        counts = {}
        for r in out_rows:
            counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
        print("%d of %d rows got one (%s) -> %s" % (n, len(rows), ", ".join("%s %d" % kv for kv in sorted(counts.items())), a.out))


def cmd_demo(a):
    from . import demo
    demo.run_demo(a.out, simulator=a.simulator, kmer_denovo_filter=a.kmer_denovo_filter, image=a.image,
                  kmer_discovery=a.kmer_discovery, kraken2_db=a.kraken2_db, pad=a.pad, work=a.work,
                  set_dir=a.set, with_viewer=a.with_viewer, threads=a.threads)


def main(argv=None):
    p = argparse.ArgumentParser(prog="svreview", description=__doc__.split(" - ", 1)[-1])
    p.add_argument("--version", action="version", version="svreview " + __version__)
    sub = p.add_subparsers(dest="command", metavar="command")
    sub.required = True

    k = sub.add_parser("check", help="check a review set (docs/FORMAT.md): what is there, what is missing")
    k.add_argument("set", help="the review set directory")
    k.add_argument("--deep", action="store_true", help="open every member and count the child's reads over the windows")
    k.add_argument("--ref", default=None, metavar="FASTA", help="the reference for CRAMs without an embedded one")
    k.set_defaults(func=cmd_check)

    c = sub.add_parser("candidates", help="the candidate table of one kmer-discovery run, in the format's layout")
    c.add_argument("--run", required=True, metavar="PREFIX", help="the run's --out-prefix (<prefix>.bed, .sv.bedpe, .metrics.json)")
    c.add_argument("--trio", required=True, help="the trio's name (the table's trio column; the bundle directory)")
    c.add_argument("--out", required=True, metavar="TSV")
    c.add_argument("--classes", default="SV", help="region classes to list, comma-joined (default: SV)")
    c.add_argument("--min-reads", type=int, default=3, help="an AMBIGUOUS region needs this many reads (default: 3)")
    c.add_argument("--sex", default=None, help="the child's sex, 1 or 2 (chrY_female flag)")
    c.add_argument("--nhf-status", default=None, help="override the run's screen status (screened, disabled, failed)")
    c.add_argument("--sink-multiplicity", type=int, default=10)
    c.add_argument("--nhf-flag", type=float, default=0.05)
    c.add_argument("--nhf-min-reads", type=int, default=5)
    c.set_defaults(func=cmd_candidates)

    b = sub.add_parser("bundle", help="cut the review bundle of every candidate from the trio's alignments")
    b.add_argument("--candidates", required=True, metavar="TSV")
    b.add_argument("--child", required=True, metavar="BAM|CRAM")
    b.add_argument("--father", required=True, metavar="BAM|CRAM")
    b.add_argument("--mother", required=True, metavar="BAM|CRAM")
    b.add_argument("--out-dir", required=True, help="the set's bundles/ directory")
    b.add_argument("--pad", type=int, default=10000, help="bp either side of every window (default: 10000)")
    b.add_argument("--ref", default=None, metavar="FASTA", help="the CRAMs' reference (unless embedded)")
    b.add_argument("--cram", action="store_true", help="write self-contained CRAMs (an embedded reference, no file named) instead of BAMs")
    b.add_argument("--tier", default=None, help="only these review tiers, comma-joined (default: all)")
    b.add_argument("--sample", action="append", metavar="ROLE=NAME", help="name a member (default: the file's SM)")
    b.add_argument("--threads", type=int, default=1)
    b.set_defaults(func=cmd_bundle)

    s = sub.add_parser("site", help="the viewer's data from a review set: mini-BAMs, the informative reads tagged, manifest.json")
    s.add_argument("set", nargs="?", default=None, help="the review set directory (candidates.tsv, trios.tsv, bundles/)")
    s.add_argument("--out", required=True, metavar="DIR")
    s.add_argument("--candidates", default=None, metavar="TSV", help="instead of a set: the candidate table")
    s.add_argument("--bundles", default=None, metavar="DIR", help="instead of a set: the bundles (<trio>/<candidate>/)")
    s.add_argument("--trios", default=None, metavar="TSV", help="instead of a set: a trios.tsv")
    s.add_argument("--metrics", action="append", metavar="[TRIO=]JSON", help="a run's metrics.json (repeatable; without TRIO= for every trio)")
    s.add_argument("--informative", action="append", metavar="[TRIO=]BAM", help="a run's informative BAM: tags the child's reads dk:i:1")
    s.add_argument("--ref", default=None, metavar="FASTA", help="the CRAMs' reference (unless embedded or named in trios.tsv)")
    s.add_argument("--pad", type=int, default=None, help="window padding for a candidate without candidate.json (default: the set's, or 10000)")
    s.add_argument("--title", default=None)
    s.add_argument("--with-viewer", action="store_true", help="copy the viewer beside the manifest: a complete site")
    s.add_argument("--viewer-dir", default=None, help="the viewer's files (default: the package's)")
    s.add_argument("--threads", type=int, default=1)
    s.set_defaults(func=cmd_site)

    v = sub.add_parser("serve", help="serve a site directory (range requests; verdicts written to disk)")
    v.add_argument("directory")
    v.add_argument("--port", type=int, default=8765)
    v.add_argument("--host", default="127.0.0.1")
    v.add_argument("--path", default="index.html", help="the page to open (default: index.html)")
    v.add_argument("--open", action="store_true", help="open the browser")
    v.add_argument("--quiet", action="store_true")
    v.set_defaults(func=cmd_serve)

    m = sub.add_parser("verdicts", help="check a verdict TSV; merge it into the candidate table")
    m.add_argument("--verdicts", required=True, metavar="TSV", help="the viewer's export, or the served verdicts.tsv")
    m.add_argument("--candidates", default=None, metavar="TSV", help="the table to merge into")
    m.add_argument("--out", default=None, metavar="TSV", help="the table with verdict columns appended")
    m.set_defaults(func=cmd_verdicts)

    d = sub.add_parser("demo", help="the demo end to end: simulate kmer_denovo_filter's sv_demo trio, run kmer-discovery, write the set, build the site")
    d.add_argument("--out", required=True, metavar="DIR", help="the site directory (docs/demo in the repository)")
    d.add_argument("--set", default=None, metavar="DIR", help="the review set to write (default: <work>/review_set; examples/sv_demo_review in the repository)")
    d.add_argument("--kmer-denovo-filter", default=None, metavar="DIR", help="a kmer_denovo_filter checkout (its examples/sv_demo simulator)")
    d.add_argument("--simulator", default=None, metavar="PY", help="simulate_sv_trio.py (default: from the checkout)")
    d.add_argument("--kmer-discovery", default="kmer-discovery", metavar="CMD", help="the kmer-discovery to run (default: on PATH)")
    d.add_argument("--image", default=None, metavar="IMAGE", help="run the simulator and kmer-discovery in this Docker image instead")
    d.add_argument("--kraken2-db", default=None, metavar="DIR", help="a kraken2 database for the non-human screen (inside the image's /data when --image)")
    d.add_argument("--pad", type=int, default=700, help="window padding for the demo's bundles (default: 700)")
    d.add_argument("--work", default=None, metavar="DIR", help="the run's working directory (default: <out>/../work)")
    d.add_argument("--with-viewer", action="store_true")
    d.add_argument("--threads", type=int, default=2)
    d.set_defaults(func=cmd_demo)

    a = p.parse_args(argv)
    if a.command == "verdicts" and bool(a.out) != bool(a.candidates):
        p.error("--out needs --candidates (and the other way round)")
    try:
        rc = a.func(a)
    except contract.ContractError as e:
        sys.exit("svreview: %s" % e)
    return rc or 0


if __name__ == "__main__":
    sys.exit(main())
