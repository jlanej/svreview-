#!/usr/bin/env bash
# The demo, end to end: kmer_denovo_filter's sv_demo trio simulated, kmer-discovery run on it with the toy
# kraken2 database (the non-human screen on: the reference as human, the 250 bp insertion as E. coli), then the
# review set written to examples/sv_demo_review (docs/FORMAT.md) and the site built from it under docs/demo.
# Everything that needs the tool runs in the pinned image.
#
#   scripts/make_demo.sh [IMAGE]          default: ghcr.io/jlanej/kmer_denovo_filter:sha-14bb559
#
# Needs docker and this repository's venv (python with pysam: `python -m venv .venv && pip install -e .`).
set -euo pipefail
cd "$(dirname "$0")/.."
IMAGE=${1:-ghcr.io/jlanej/kmer_denovo_filter:sha-14bb559}
PY=${PY:-.venv/bin/python}
WORK=$PWD/work
rm -rf "$WORK/sv_demo" "$WORK/k2"
mkdir -p "$WORK/sv_demo" "$WORK/k2/taxonomy"

# 1. simulate once (the image ships the simulator; the seed is fixed), then the toy taxonomy and library
docker run --rm --platform linux/amd64 -v "$WORK:/data" -w /data --entrypoint bash "$IMAGE" -c \
  "python3 /app/examples/sv_demo/simulate_sv_trio.py --outdir /data/sv_demo"
$PY scripts/toy_kraken2_library.py --simulator-image "$IMAGE" --ref "$WORK/sv_demo/ref.fa" --out "$WORK/k2"

# 2. the database, in the image
docker run --rm --platform linux/amd64 -v "$WORK:/data" -w /data --entrypoint bash "$IMAGE" -c "set -e
  kraken2-build --add-to-library /data/k2/library.fa --db /data/k2 --no-masking
  kraken2-build --build --db /data/k2 --threads 2"

# 3. the run (simulated again, identically), the review set, the check, the site
$PY -m svreview.cli demo --image "$IMAGE" --kraken2-db /data/k2 --work "$WORK" \
  --set examples/sv_demo_review --out docs/demo --pad 700
echo "demo: examples/sv_demo_review (the set) and docs/demo (the site)"
