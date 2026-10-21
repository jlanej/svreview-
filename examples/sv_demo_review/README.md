# A review set, by example: the demo's input

This directory is a complete review set in svreview's input format ([docs/FORMAT.md](../../docs/FORMAT.md)),
the one the site under [docs/demo](../../docs/demo/) is built from. It is what an analysis hands to svreview:
whoever writes one of these for their own run - a cohort screen, a single trio - needs nothing else.

```
svreview.json        format 1, the title, the pad (700 bp here)
candidates.tsv       10 candidate regions of one trio: the seven ends of the six junctions (DEL, DUP, INV x2,
                     BND x2) and the two insertions, with evidence, junctions, flags, tier and rank
trios.tsv            the one trio: its members, its run's metrics.json and informative BAM under runs/
bundles/sv_demo/     six bundles, one per set of windows: child/father/mother as CRAMs with the reference
                     embedded (they open with no FASTA), candidate.json, the members' .json; the other end
                     of each pair is a symlink to the first end's directory (bundle.tsv lists all ten)
runs/                the run's own files: sv_demo.metrics.json, sv_demo.informative.bam (the reads carrying
                     proband-unique k-mers, dk:i:1), and for reference its BED, BEDPE and summary, and
                     events.tsv (the simulator's truth)
```

The trio is kmer_denovo_filter's `examples/sv_demo`: a random 90 kb reference, 2×150 bp pairs at 30×, the
parents reference-only, the child heterozygous for every event. The non-human screen ran with a toy kraken2
database (the reference as *Homo sapiens*, the 250 bp novel insertion as *Escherichia coli*), so the
insertion at chr2:21,867-22,078 carries `nhf` 0.5833 and the `nhf` flag (tier 4), the way a contaminant
would; the 55 bp insertion is tier 3 (clipped reads, no junction); the rest are tier 1.

```bash
svreview check examples/sv_demo_review --deep          # 0 errors
svreview site examples/sv_demo_review --out /tmp/demo   # = docs/demo
svreview serve /tmp/demo --open
```

`scripts/make_demo.sh` regenerates this directory and docs/demo in the pinned image.
