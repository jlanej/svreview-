# The public example: a 1000 Genomes trio (plan)

The demo under [`docs/demo`](../demo/) is a simulated trio. The example is to be a real one: a public
1000 Genomes high-coverage trio run through `kmer-discovery`, its candidates bundled and reviewed here, with the
verdicts and what they taught written down. Nothing of it is in the repository yet: this page is the plan, so
that the example is built the same way wherever the trio's files are, and is updated when it is.

## The trio

NA19139 (YRI; the 1000 Genomes "bench" trio of the k-mer screen: the child NA19139 with NA19137 and NA19138 as
the parents, from the 1000 Genomes 30x high-coverage release, CRAMs on GRCh38 with decoys and HLA, as the IGSR
mirror and ENA serve them). The whole CRAMs are being pulled to a local disk by the screen's own staging; until
they are, nothing here runs. The trio was chosen by the screen for its public status, its depth and a child
with lymphoblastoid-line rearrangements to show (IGH, IGK: every B-cell-derived line carries its own V(D)J
rearrangement, somatic and absent from the parents' clones, which the candidate table flags and the viewer
shows as "the child's reads only").

## The steps, when the CRAMs are there

```bash
# 1. the run: kmer_denovo_filter at the pin the screen uses (main 14bb559 or later), with the non-human
#    screen on (a kraken2 database: the screen's)
kmer-discovery --child NA19139.cram --father NA19137.cram --mother NA19138.cram \
  --ref-fasta GRCh38_full_analysis_set_plus_decoy_hla.fa --out-prefix bench/NA19139.discovery \
  --kraken2-db "$KRAKEN2_DB" --threads 16 --memory 64G

# 2. the candidate table, by the screen's rules (one run: no recurrence across children)
svreview candidates --run bench/NA19139.discovery --trio NA19139 --sex 2 --out bench/candidates.tsv

# 3. the bundles: every tier-1 and tier-2 candidate, 10 kb either side (the screen's pad), the three CRAMs
svreview bundle --candidates bench/candidates.tsv --tier 1,2 --pad 10000 \
  --child NA19139.cram --father NA19137.cram --mother NA19138.cram \
  --ref GRCh38_full_analysis_set_plus_decoy_hla.fa --out-dir bench/review

# 4. the site: BAMs for the browser, the informative reads tagged, one manifest
svreview site --candidates bench/candidates.tsv --bundles bench/review --out docs/example/NA19139 \
  --metrics bench/NA19139.discovery.metrics.json --informative bench/NA19139.discovery.informative.bam \
  --ref GRCh38_full_analysis_set_plus_decoy_hla.fa --title "NA19139 (1000 Genomes YRI trio)"

# 5. review: svreview serve docs/example/NA19139 --open; the verdicts land in verdicts.tsv beside the manifest
# 6. the record: svreview verdicts --verdicts docs/example/NA19139/verdicts.tsv \
#      --candidates docs/example/NA19139/candidates.tsv --out docs/example/NA19139/candidates.verdicts.tsv
```

Or, with the Analysis-side review bundles already cut (the trio screen's `bundle_candidates.py`, mini-CRAMs
with the reference embedded, `bundle.tsv`): step 4 alone, `--bundles` pointed at that directory and no `--ref`
needed.

## What goes into the repository, and what does not

* In: `docs/example/NA19139/manifest.json`, `candidates.tsv`, the per-candidate mini-BAMs of the reviewed
  candidates (a few MB each; the set is capped at what Pages serves comfortably, about 100 MB in all - the
  strongest tiers first), `verdicts.tsv` and `candidates.verdicts.tsv`, and a `README.md` of the findings: what
  each verdict was and why, which flags were right, where the viewer helped and where it did not.
* Not in: the CRAMs, the reference, the kraken2 database, the informative BAM, any local path (the commands
  above take the files from the working directory; the Analysis-side runbooks hold the site's paths).
* The reads are public (1000 Genomes, open access) and the sample IDs are public identifiers; nothing from a
  controlled-access cohort is ever written here.

## What the example is for

1. **The reader on real data.** The demo's reads are simulated (every alignment written from the read's known
   origin). Real BAMs bring novoalign/bwa tags, long CIGARs, secondary alignments, MAPQ 0 piles, supplementary
   records with hard clips: the fixtures under `tests/js/fixtures` gain a real bundle.
2. **The trio read on a lymphoblastoid line.** The IGH and IGK candidates are the child's own and absent from
   the parents: "the child's reads only" is not "de novo". The viewer must make the distinction a reviewer makes
   (flags, the locus, the verdict `artifact` with a note), not hide it.
3. **The non-human fraction on a real screen.** With the screen on, a region with a non-zero fraction and a top
   taxon; with the per-read sidecar (`<role>.taxa.tsv`, from kraken2 over the bundle's reads) the colour-by-taxon
   mode shows which reads carry it.
4. **Timing.** How long a bundle takes to open over Pages (range requests to a CDN) against a local serve.
