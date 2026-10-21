#!/usr/bin/env python3
"""The demo's toy kraken2 library: the sv_demo reference as Homo sapiens (9606) and the child's 250 bp novel
insertion as Escherichia coli (562), with a seven-node taxonomy (the same fixture as kmer_denovo_filter's
tests/discovery/test_nhf.py). Built with the simulator's own seed, so it needs the simulator's source: a
kmer_denovo_filter checkout (--kmer-denovo-filter) or the Docker image (--simulator-image: the script is read
out of /app/examples/sv_demo). Writes <out>/library.fa and <out>/taxonomy/{nodes,names}.dmp; kraken2-build
does the rest (scripts/make_demo.sh)."""
import argparse
import importlib.util
import os
import random
import subprocess
import sys
import tempfile

import pysam

NODES = [(1, 1, "no rank"), (131567, 1, "no rank"), (2, 131567, "superkingdom"), (562, 2, "species"),
         (2759, 131567, "superkingdom"), (33208, 2759, "kingdom"), (9606, 33208, "species")]
NAMES = {1: "root", 131567: "cellular organisms", 2: "Bacteria", 562: "Escherichia coli", 2759: "Eukaryota",
         33208: "Metazoa", 9606: "Homo sapiens"}


def load_simulator(path):
    spec = importlib.util.spec_from_file_location("simulate_sv_trio", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True, help="the simulated ref.fa (to check the seed)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--kmer-denovo-filter", default=None, help="a checkout: its examples/sv_demo/simulate_sv_trio.py")
    ap.add_argument("--simulator", default=None)
    ap.add_argument("--simulator-image", default=None, help="read the simulator out of this Docker image")
    a = ap.parse_args(argv)
    sim_path = a.simulator
    if sim_path is None and a.kmer_denovo_filter:
        sim_path = os.path.join(a.kmer_denovo_filter, "examples", "sv_demo", "simulate_sv_trio.py")
    tmp = None
    if sim_path is None and a.simulator_image:
        tmp = tempfile.mkdtemp()
        sim_path = os.path.join(tmp, "simulate_sv_trio.py")
        with open(sim_path, "w") as fh:
            subprocess.run(["docker", "run", "--rm", "--platform", "linux/amd64", "--entrypoint", "cat", a.simulator_image,
                            "/app/examples/sv_demo/simulate_sv_trio.py"], stdout=fh, check=True)
    if sim_path is None or not os.path.isfile(sim_path):
        sys.exit("give --kmer-denovo-filter DIR, --simulator PY or --simulator-image IMAGE")
    sim = load_simulator(sim_path)
    rng = random.Random(1)
    ref = dict((c, sim.random_seq(rng, n)) for c, n in sim.CHROM_SIZES.items())
    with pysam.FastaFile(a.ref) as fa:
        if not all(fa.fetch(c) == ref[c] for c in ref):
            sys.exit("%s was not simulated with the default seed: the insertion cannot be reconstructed" % a.ref)
    novel = [h.seq[s:e] for h in sim.child_haplotypes(ref, rng) for (s, e, chrom, *_rest) in h.segments if chrom is None]
    (ins250,) = [s for s in novel if len(s) == 250]
    os.makedirs(os.path.join(a.out, "taxonomy"), exist_ok=True)
    with open(os.path.join(a.out, "taxonomy", "nodes.dmp"), "w") as fh:
        for t, p, r in NODES:
            fh.write("%d\t|\t%d\t|\t%s\t|\t\t|\n" % (t, p, r))
    with open(os.path.join(a.out, "taxonomy", "names.dmp"), "w") as fh:
        for t, n in NAMES.items():
            fh.write("%d\t|\t%s\t|\t\t|\tscientific name\t|\n" % (t, n))
    with open(os.path.join(a.out, "library.fa"), "w") as fh:
        for c in ref:
            fh.write(">%s|kraken:taxid|9606 human\n%s\n" % (c, ref[c]))
        fh.write(">ins250|kraken:taxid|562 Escherichia coli\n%s\n" % ins250)
    print("toy library: %d human contigs + the 250 bp insertion as E. coli -> %s" % (len(ref), a.out))


if __name__ == "__main__":
    main()
