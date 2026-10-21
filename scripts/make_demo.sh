#!/usr/bin/env bash
# The demo, end to end: kmer_denovo_filter's sv_demo trio simulated, kmer-discovery run on it (with the toy
# kraken2 database, so the non-human screen is on: the reference as human, the 250 bp insertion as E. coli),
# the candidate table, the bundles (pad 700) and the site under docs/demo. Everything runs in the pinned image.
#
#   scripts/make_demo.sh [IMAGE]          default: ghcr.io/jlanej/kmer_denovo_filter:sha-14bb559
#
# Needs docker and this repository's venv (python with pysam: `python -m venv .venv && pip install -e .`).
set -euo pipefail
cd "$(dirname "$0")/.."
IMAGE=${1:-ghcr.io/jlanej/kmer_denovo_filter:sha-14bb559}
PY=${PY:-.venv/bin/python}
WORK=$PWD/work
rm -rf "$WORK/sv_demo" "$WORK/k2" "$WORK/bundles"
mkdir -p "$WORK/sv_demo" "$WORK/k2/taxonomy"

# 1. simulate (the image ships the simulator), then the toy taxonomy and library from the simulator's seed
docker run --rm --platform linux/amd64 -v "$WORK:/data" -w /data --entrypoint bash "$IMAGE" -c \
  "python3 /app/examples/sv_demo/simulate_sv_trio.py --outdir /data/sv_demo"
$PY scripts/toy_kraken2_library.py --simulator-image "$IMAGE" --ref "$WORK/sv_demo/ref.fa" --out "$WORK/k2" || \
  { echo "the toy library needs the simulator's source: pass a checkout with --kmer-denovo-filter (see the script)"; exit 1; }

# 2. the database and the run, in the image
docker run --rm --platform linux/amd64 -v "$WORK:/data" -w /data --entrypoint bash "$IMAGE" -c "set -e
  kraken2-build --add-to-library /data/k2/library.fa --db /data/k2 --no-masking
  kraken2-build --build --db /data/k2 --threads 2
  kmer-discovery --child /data/sv_demo/child.bam --mother /data/sv_demo/mother.bam --father /data/sv_demo/father.bam \
    --ref-fasta /data/sv_demo/ref.fa --out-prefix /data/sv_demo/sv_demo --kraken2-db /data/k2 --threads 2"

# 3. the table, the bundles, the site
$PY -m svreview.cli candidates --run "$WORK/sv_demo/sv_demo" --trio sv_demo --out "$WORK/sv_demo/candidates.tsv"
$PY -m svreview.cli bundle --candidates "$WORK/sv_demo/candidates.tsv" --child "$WORK/sv_demo/child.bam" \
  --father "$WORK/sv_demo/father.bam" --mother "$WORK/sv_demo/mother.bam" --out-dir "$WORK/bundles" --pad 700 \
  --sample child=child --sample father=father --sample mother=mother
rm -rf docs/demo
$PY -m svreview.cli site --candidates "$WORK/sv_demo/candidates.tsv" --bundles "$WORK/bundles" --out docs/demo \
  --metrics "$WORK/sv_demo/sv_demo.metrics.json" --informative "$WORK/sv_demo/sv_demo.informative.bam" \
  --title "the sv_demo trio (kmer_denovo_filter examples/sv_demo)"
# the run's own files that the README points at
cp "$WORK/sv_demo/sv_demo.bed" "$WORK/sv_demo/sv_demo.sv.bedpe" "$WORK/sv_demo/sv_demo.summary.txt" "$WORK/sv_demo/events.tsv" docs/demo/
echo "demo written to docs/demo"
