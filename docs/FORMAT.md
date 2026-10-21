# The review set: svreview's input format (version 1)

One directory holds everything a review needs. Whoever runs the discovery - the 1000 Genomes trio screen,
the GMKF arm, a single trio at home - writes one of these; svreview reads nothing else, so no arm needs an
adapter inside svreview, and `svreview check` says whether a set is complete before anyone opens it.

```
<set>/
  svreview.json           optional: {"format": 1, "title": "...", "pad": 10000, "reference": "GRCh38"}
  candidates.tsv          required: one row per candidate region (below)
  trios.tsv               required: one row per trio - who the members are, and where its files are (below)
  bundles/                per-candidate read bundles (below); or none, when trios.tsv names per-trio alignments
  runs/                   anywhere, really: the files trios.tsv points at (metrics, informative BAMs)
  verdicts.tsv            written by the review (svreview serve, or an export placed here); never by the arm
```

Paths inside the set are relative to the set's directory (forward slashes); a set is movable as a whole.
Coordinates are BED's throughout the tables: 0-based, half-open. Only `candidate.json` and the screen show
1-based samtools regions. A missing value is `.`; it is never 0.

## candidates.tsv

Tab-separated, a header line, one row per **candidate region**: a cluster of the child's informative reads
(kmer_denovo_filter's discovery BED row), with the junctions that leave it. A junction's two ends are two
rows that name each other in `junctions`; both carry `both_ends` 1. This is the unit the viewer navigates
by, and the unit a bundle is cut for.

| column | required | meaning |
|---|---|---|
| `trio` | yes | the trio's name, as in trios.tsv |
| `chrom` `start` `end` | yes | the region (0-based, half-open) |
| `class` | | `SV`, `AMBIGUOUS` or `SMALL` (the tool's) |
| `sv_type` | | `DEL` `DUP` `INV` `BND` `INS` `INTRA` or `.` |
| `reads` `unique_kmers` `split_reads` `discordant_pairs` `unmapped_mates` `breakpoint_reads` `large_indel_reads` `max_clip_len` `max_mapq` `mapq0_reads` | | the evidence as the BED has it (molecules, k-mers, bp, MAPQ) |
| `junctions` | | `partner\|strands\|type\|reads` entries joined by `;`: the partner's region as `chrom:start-end` (0-based, the partner row's own coordinates), the BEDPE line's two strands (`+-` DEL, `-+` DUP, `++`/`--` INV, `..` unknown), its type, its supporting molecules; `.` when the region has no junction |
| `locus_flags` `partner_flags` | | `;`-joined names of what the region's own locus and its partners' overlap that a reader sets aside (IGH, IGK, IGL, TRA, TRB, TRG, PAR, chrEBV, chrM, nonprimary, chrY_female); `.` when none |
| `flags` | | every reason to set the candidate aside, `;`-joined: the locus flags, the partner flags with a `partner:` prefix, `partner_recurrent`, `partner_sink`, `prior_recurrent`, `mapq0_only`, `nhf`, `nhf_unscreened`; `.` when clean |
| `review_tier` | | 1 clean with both ends, 2 clean with a junction, 3 clean without, 4 flagged |
| `rank` | | 1 = strongest, within the table |
| `both_ends` | | 1 when a junction partner is itself a row of the table |
| `nhf_reads` `nhf` | | the non-human screen's reads and fraction for the region; `.` = not screened |
| anything else | | carried unread and shown with the row (`n_trios_at_locus`, `partner_multiplicity`, `prior_recurrence`, ...) |

A table with only the four required columns is a valid set (no evidence, no junctions, every candidate
tier `.`): the viewer draws what the bundle holds. Everything the 1000 Genomes screen's table carries
(`summarize_trios.py`, then `flag_candidates.py`) is already this layout, in this column order.

**From a table of pairs** (one row per junction with `chrom1 start1 end1 chrom2 start2 end2`, as the GMKF
arm's `collate.py` writes): emit one row per end, each with `junctions` naming the other end (one entry per
orientation, the pair's strands and molecules), `both_ends` 1, the end's own evidence columns (`bp1_*` for the
first row, `bp2_*` for the second), the pair's flags on both; a region without a junction is one row with
`junctions` `.`.

## trios.tsv

Tab-separated, a header line, one row per trio.

| column | required | meaning |
|---|---|---|
| `trio` | yes | the name candidates.tsv uses |
| `child` `father` `mother` | yes | the members' sample names (what the bundles' read groups say, or any label) |
| `child_sex` | | `1` male, `2` female, `0` unknown |
| `child_alignments` `father_alignments` `mother_alignments` | | per-trio alignment files (BAM or CRAM, indexed) covering every candidate's windows of the trio - the alternative to per-candidate bundles: `svreview site` cuts each candidate's windows from them |
| `metrics` | | the run's `metrics.json` (`regions`, `nonhuman_fraction`): the region evidence and the screen |
| `informative` | | the run's informative BAM (`dk:i:1` reads): lights the k-mer evidence up |
| `reference` | | a FASTA (indexed) for CRAMs without an embedded reference |

## bundles/

`bundles/<trio>/<chrom>_<start>_<end>/`, named by the candidate's own region, holding

* `child.cram` / `father.cram` / `mother.cram` (or `.bam`), each with its index: the reads over the
  candidate's **windows** - its region and each junction partner's region, padded by the set's `pad`, merged
  where they touch - and their mates wherever they map (`samtools view -P`); a CRAM either embeds its
  reference (`--output-fmt-option embed_ref=1`) or the trio names a `reference`;
* `candidate.json`: `windows` (samtools regions, 1-based inclusive, in reference order), `members`
  (`{role: sample}`), `pad`, and `shared_by` (the other end of a pair, when its directory is a symlink to
  this one); without it the windows are computed from the row and the set's `pad`;
* optionally `<role>.taxa.tsv` (`read taxid name human`): per-read taxa for the colour-by-taxon mode;
* optionally `bundles/bundle.tsv`, the index (`candidate trio sv_type review_tier rank member sample region
  reads bytes source_url status path`); svreview rebuilds it from disk when absent.

The other end of a pair has the same windows: its directory is a symlink to the first end's (or a copy).

A candidate with neither a bundle nor per-trio alignments is listed by the viewer and marked as having
nothing to draw; it is never dropped.

## verdicts.tsv

What the review writes, and the only file an arm reads back:

```
candidate   trio   chrom   start   end   sv_type   verdict   note   reviewer   reviewed_at
chr1:9860-10000   sv_demo   chr1   9860   10000   DEL   de_novo   both ends clean   jl   2026-10-21T16:02:32
```

`verdict` is `de_novo`, `inherited`, `artifact` or `unsure`; a row with an empty verdict clears an earlier
one; the file is append-only and the last row for a candidate wins. `svreview verdicts --candidates ... --out`
appends `verdict verdict_note verdict_reviewer verdict_at` to the candidate table, keyed by
`trio chrom start end`, for the arm to read as a flag (`.` where there is none).

## Checking and using a set

```bash
svreview check <set>                     # what is there, what is missing; exit 1 on an error
svreview site <set> --out <site>         # BAMs for the browser, the informative reads tagged, one manifest
svreview serve <site> --open             # review; verdicts.tsv fills in beside the manifest
svreview verdicts --verdicts <site>/verdicts.tsv --candidates <set>/candidates.tsv --out <set>/candidates.verdicts.tsv
```

[`examples/sv_demo_review`](../examples/sv_demo_review/) is a complete set (the demo's input); the site under
[`docs/demo`](demo/) is built from it.

## Versioning

`svreview.json`'s `format` is this document's version. A column added is not a new version (columns are read
by name); a column whose meaning changes, or a file whose layout changes, is.
