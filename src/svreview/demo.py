"""demo.py - the demo end to end, on kmer_denovo_filter's own simulated trio (examples/sv_demo).

The trio: a random two-chromosome reference (chr1 60 kb, chr2 30 kb), 2x150 bp pairs at 30x, the parents
reference-only, the child heterozygous for a 2 kb deletion, a 1.5 kb tandem duplication, a 2 kb inversion, a
reciprocal chr1/chr2 translocation, a 55 bp and a 250 bp novel insertion and an SNV. kmer-discovery finds the
six junctions (DEL, DUP, INV x2, BND x2) and the two insertions; the candidate table, the bundles and the site
follow from its output by the same commands a real run uses:

    svreview candidates --run work/sv_demo/sv_demo --trio sv_demo --out work/sv_demo/candidates.tsv
    svreview bundle --candidates ... --child work/sv_demo/child.bam ... --out-dir work/bundles --pad 700
    svreview site --candidates ... --bundles work/bundles --out docs/demo \\
        --metrics work/sv_demo/sv_demo.metrics.json --informative work/sv_demo/sv_demo.informative.bam

Needs the simulator (a kmer_denovo_filter checkout, or the Docker image, which ships it) and kmer-discovery
(on PATH, or the image). With a kraken2 database the non-human screen runs too (the demo's own toy database:
the reference as human, the 250 bp insertion as E. coli; see scripts/make_demo.sh).
"""

import os
import shutil
import subprocess
import sys

from . import bundle, candidates, contract, site

TRIO = "sv_demo"


def _run(cmd, log, cwd=None):
    log("+ " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        raise contract.ContractError("%s failed (exit %d)" % (cmd[0], r.returncode))


def run_demo(out, simulator=None, kmer_denovo_filter=None, image=None, kmer_discovery="kmer-discovery",
             kraken2_db=None, pad=700, work=None, with_viewer=False, threads=2, log=print):
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
    # 2. the candidate table
    header, rows = candidates.candidates_from_run(prefix, TRIO)
    table = os.path.join(run_dir, "candidates.tsv")
    contract.write_tsv(table, header, rows)
    log("%d candidates -> %s" % (len(rows), table))
    # 3. the bundles
    _h, cands = contract.read_candidates(table)
    members = dict((r, os.path.join(run_dir, r + ".bam")) for r in contract.ROLES)
    bundles = os.path.join(work, "bundles")
    if os.path.isdir(bundles):
        shutil.rmtree(bundles)
    n, _index = bundle.bundle_candidates(cands, members, bundles, pad=pad, samples=dict((r, r) for r in contract.ROLES),
                                         threads=threads, log=log)
    # 4. the site
    manifest = site.build_site(table, bundles, out, metrics={TRIO: prefix + ".metrics.json"},
                               informative={TRIO: prefix + ".informative.bam"}, title="svreview demo: the sv_demo trio",
                               with_viewer=with_viewer, log=log)
    log("demo: %d candidates, %d bundles -> %s" % (len(cands), n, out))
    return manifest
