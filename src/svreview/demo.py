"""demo.py - the demo end to end, on kmer_denovo_filter's own simulated trio (examples/sv_demo).

The trio: a random two-chromosome reference (chr1 60 kb, chr2 30 kb), 2x150 bp pairs at 30x, the parents
reference-only, the child heterozygous for a 2 kb deletion, a 1.5 kb tandem duplication, a 2 kb inversion, a
reciprocal chr1/chr2 translocation, a 55 bp and a 250 bp novel insertion and an SNV. kmer-discovery finds the
six junctions (DEL, DUP, INV x2, BND x2) and the two insertions. From its output the demo writes a complete
review set (docs/FORMAT.md) - the table, trios.tsv, CRAM bundles with the reference embedded, the run's
metrics and informative BAM under runs/ - and builds the site from it, by the same commands a real run uses:

    svreview candidates --run work/sv_demo/sv_demo --trio sv_demo --out SET/candidates.tsv
    svreview bundle --candidates SET/candidates.tsv --child work/sv_demo/child.bam ... --out-dir SET/bundles \\
        --pad 700 --cram --ref work/sv_demo/ref.fa
    svreview check SET
    svreview site SET --out docs/demo

Needs the simulator (a kmer_denovo_filter checkout, or the Docker image, which ships it) and kmer-discovery
(on PATH, or the image). With a kraken2 database the non-human screen runs too (the demo's own toy database:
the reference as human, the 250 bp insertion as E. coli; see scripts/make_demo.sh).
"""

import json
import os
import shutil
import subprocess
import sys

from . import bundle, candidates, contract, review_set, site

TRIO = "sv_demo"
RUN_FILES = (".metrics.json", ".informative.bam", ".informative.bam.bai", ".bed", ".sv.bedpe", ".summary.txt")


def _run(cmd, log, cwd=None):
    log("+ " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        raise contract.ContractError("%s failed (exit %d)" % (cmd[0], r.returncode))


def write_set(set_dir, prefix, members, ref, pad=700, title=None, threads=1, log=print):
    """A review set from one run: the table, trios.tsv, CRAM bundles (reference embedded), runs/, svreview.json."""
    set_dir = os.path.abspath(set_dir)
    os.makedirs(os.path.join(set_dir, "runs"), exist_ok=True)
    header, rows = candidates.candidates_from_run(prefix, TRIO)
    table = os.path.join(set_dir, "candidates.tsv")
    contract.write_tsv(table, header, rows)
    log("%d candidates -> %s" % (len(rows), table))
    _h, cands = contract.read_candidates(table)
    bundles = os.path.join(set_dir, "bundles")
    if os.path.isdir(bundles):
        shutil.rmtree(bundles)
    n, _index = bundle.bundle_candidates(cands, members, bundles, pad=pad, ref=ref,
                                         samples=dict((r, r) for r in contract.ROLES), threads=threads, log=log, cram=True)
    for ext in RUN_FILES:
        src = prefix + ext
        if not os.path.isfile(src):
            continue
        if ext == ".informative.bam":
            bundle.clean_copy(src, os.path.join(set_dir, "runs", TRIO + ext))      # its header names the run's paths
        elif ext != ".informative.bam.bai":
            shutil.copyfile(src, os.path.join(set_dir, "runs", TRIO + ext))
    events = os.path.join(os.path.dirname(prefix), "events.tsv")
    if os.path.isfile(events):
        shutil.copyfile(events, os.path.join(set_dir, "runs", "events.tsv"))
    review_set.write_trios(os.path.join(set_dir, "trios.tsv"), [dict(
        trio=TRIO, child="child", father="father", mother="mother", child_sex="0",
        metrics="runs/%s.metrics.json" % TRIO, informative="runs/%s.informative.bam" % TRIO)])
    with open(os.path.join(set_dir, "svreview.json"), "w") as fh:
        json.dump({"format": review_set.FORMAT_VERSION, "title": title or "the sv_demo trio (kmer_denovo_filter examples/sv_demo)",
                   "pad": pad, "reference": "the simulator's random 90 kb reference (chr1 60 kb, chr2 30 kb)"}, fh, indent=1)
    return set_dir, n


def run_demo(out, simulator=None, kmer_denovo_filter=None, image=None, kmer_discovery="kmer-discovery",
             kraken2_db=None, pad=700, work=None, set_dir=None, with_viewer=False, threads=2, log=print):
    out = os.path.abspath(out)
    work = os.path.abspath(work or os.path.join(os.path.dirname(out), "work"))
    run_dir = os.path.join(work, TRIO)
    os.makedirs(run_dir, exist_ok=True)
    prefix = os.path.join(run_dir, TRIO)
    if simulator is None and kmer_denovo_filter:
        simulator = os.path.join(kmer_denovo_filter, "examples", "sv_demo", "simulate_sv_trio.py")
    # 1. simulate and discover
    if image:
        docker = ["docker", "run", "--rm", "--platform", "linux/amd64", "-v", work + ":/data", "-w", "/data",
                  "--entrypoint", "bash", image, "-c"]
        inner = ("python3 /app/examples/sv_demo/simulate_sv_trio.py --outdir /data/%s && "
                 "kmer-discovery --child /data/%s/child.bam --mother /data/%s/mother.bam --father /data/%s/father.bam "
                 "--ref-fasta /data/%s/ref.fa --out-prefix /data/%s/%s --threads %d%s"
                 % (TRIO, TRIO, TRIO, TRIO, TRIO, TRIO, TRIO, threads,
                    (" --kraken2-db " + kraken2_db) if kraken2_db else ""))
        _run(docker + [inner], log)
    else:
        if not simulator or not os.path.isfile(simulator):
            raise contract.ContractError("the simulator was not found: give --kmer-denovo-filter DIR, --simulator PY or --image")
        _run([sys.executable, simulator, "--outdir", run_dir], log)
        if shutil.which(kmer_discovery) is None and not os.path.isfile(kmer_discovery):
            raise contract.ContractError("%s is not on PATH (a kmer_denovo_filter install), or give --image" % kmer_discovery)
        cmd = [kmer_discovery, "--child", os.path.join(run_dir, "child.bam"), "--mother", os.path.join(run_dir, "mother.bam"),
               "--father", os.path.join(run_dir, "father.bam"), "--ref-fasta", os.path.join(run_dir, "ref.fa"),
               "--out-prefix", prefix, "--threads", str(threads)]
        if kraken2_db:
            cmd += ["--kraken2-db", kraken2_db]
        _run(cmd, log)
    # 2. the review set
    set_dir = os.path.abspath(set_dir or os.path.join(work, "review_set"))
    members = dict((r, os.path.join(run_dir, r + ".bam")) for r in contract.ROLES)
    set_dir, n = write_set(set_dir, prefix, members, os.path.join(run_dir, "ref.fa"), pad=pad, threads=threads, log=log)
    # 3. the check, then the site
    rs = review_set.ReviewSet.open(set_dir)
    text, n_err = review_set.format_findings(rs.check(deep=True))
    log(text)
    if n_err:
        raise contract.ContractError("the demo's review set has errors")
    manifest = site.build_site_from_set(rs, out, with_viewer=with_viewer, threads=threads, log=log)
    log("demo: %d candidates, %d bundles; set %s; site %s" % (len(rs.candidates), n, set_dir, out))
    return manifest
