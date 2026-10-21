"""The candidate table of one run, by the 1000 Genomes screen's rules."""

import os

import pytest

from svreview import candidates, contract

BED = """#filters: min_distinct_kmers=1 min_distinct_kmers_per_read=7 min_supporting_reads=1
#chrom\tstart\tend\treads\tunique_kmers\tsplit_reads\tdiscordant_pairs\tmax_clip_len\tunmapped_mates\tclass\tbreakpoint_reads\tlarge_indel_reads\tsv_type\tmax_mapq\tmapq0_reads\tnhf_reads\tnhf
chr1\t9860\t10000\t11\t28\t6\t4\t72\t0\tSV\t8\t0\tDEL\t60\t0\t11\t0.0000
chr1\t12000\t12140\t8\t28\t7\t3\t54\t0\tSV\t7\t0\tDEL\t60\t0\t8\t0.0000
chr1\t54860\t55115\t11\t60\t8\t5\t73\t0\tSV\t9\t0\tBND\t60\t0\t11\t0.0000
chr2\t4872\t5124\t10\t31\t0\t0\t0\t0\tSMALL\t0\t0\t.\t60\t0\t.\t.
chr2\t14861\t15141\t12\t60\t6\t6\t72\t0\tSV\t10\t0\tBND\t60\t0\t12\t0.0000
chr2\t21866\t22078\t10\t211\t0\t1\t114\t7\tSV\t9\t0\tINS\t60\t0\t12\t0.8333
chr2\t27000\t27100\t3\t9\t0\t1\t30\t0\tAMBIGUOUS\t0\t0\t.\t0\t3\t3\t0.0000
chr14\t105600000\t105600200\t5\t20\t3\t2\t40\t0\tSV\t3\t0\tDEL\t60\t0\t5\t0.0000
"""
BEDPE = """#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tsupporting_reads\tstrand1\tstrand2\tsv_type
chr1\t9860\t10000\tchr1\t12000\t12140\tSV_1\t15\t+\t-\tDEL
chr1\t54860\t55115\tchr2\t14861\t15141\tSV_5\t10\t+\t-\tBND
chr1\t54860\t55115\tchr2\t14861\t15141\tSV_6\t8\t-\t+\tBND
chr14\t105600000\t105600200\tchr14\t105700000\t105700100\tSV_9\t4\t+\t-\tDEL
"""


@pytest.fixture
def run(tmp_path):
    prefix = str(tmp_path / "demo")
    with open(prefix + ".bed", "w") as fh:
        fh.write(BED)
    with open(prefix + ".sv.bedpe", "w") as fh:
        fh.write(BEDPE)
    return prefix


def by_locus(rows):
    return dict(("%s:%s-%s" % (r["chrom"], r["start"], r["end"]), r) for r in rows)


def test_table_of_a_screened_run(run):
    header, rows = candidates.candidates_from_run(run, "T1", nhf_status="screened")
    assert header == list(contract.CANDIDATE_COLUMNS)
    t = by_locus(rows)
    assert set(t) == {"chr1:9860-10000", "chr1:12000-12140", "chr1:54860-55115", "chr2:14861-15141",
                      "chr2:21866-22078", "chr14:105600000-105600200"}   # SV class only; SMALL and AMBIGUOUS left out
    d = t["chr1:9860-10000"]
    assert d["junctions"] == "chr1:12000-12140|+-|DEL|15" and d["both_ends"] == 1
    assert d["review_tier"] == 1 and d["flags"] == "." and d["trio"] == "T1" and d["n_trios_at_locus"] == 1
    assert t["chr1:12000-12140"]["junctions"] == "chr1:9860-10000|+-|DEL|15"
    # two junctions of one pair, in BEDPE order, the strands as the line has them
    assert t["chr1:54860-55115"]["junctions"] == "chr2:14861-15141|+-|BND|10;chr2:14861-15141|-+|BND|8"
    assert t["chr2:14861-15141"]["junctions"] == "chr1:54860-55115|+-|BND|10;chr1:54860-55115|-+|BND|8"
    assert t["chr2:14861-15141"]["partner_multiplicity"] == 2 and t["chr2:14861-15141"]["review_tier"] == 1
    # the insertion: no junction, flagged by the screen -> tier 4 with the nhf flag
    ins = t["chr2:21866-22078"]
    assert ins["junctions"] == "." and ins["flags"] == "nhf" and ins["review_tier"] == 4 and ins["nhf"] == 0.8333
    # a junction whose partner is not a candidate: tier 2, and IGH flags both ends -> tier 4
    igh = t["chr14:105600000-105600200"]
    assert igh["both_ends"] == 0 and igh["locus_flags"] == "IGH" and igh["partner_flags"] == "IGH"
    assert igh["flags"] == "IGH;partner:IGH" and igh["review_tier"] == 4
    # rank: both ends first, then junction reads
    ranks = dict((k, r["rank"]) for k, r in t.items())
    assert ranks["chr2:14861-15141"] < ranks["chr1:9860-10000"] < ranks["chr14:105600000-105600200"]
    assert sorted(ranks.values()) == list(range(1, 7))
    # the table's order: by chromosome then start (as summarize_trios.py writes it)
    assert [(r["chrom"], r["start"]) for r in rows] == sorted((r["chrom"], r["start"]) for r in rows)
    # it reads back through the contract
    path = run + ".candidates.tsv"
    contract.write_tsv(path, header, rows)
    _h, cands = contract.read_candidates(path)
    assert [c.review_tier for c in cands].count(1) == 4
    assert cands[0].junctions[0].reads == 15


def test_unscreened_run_flags_everything(run):
    _h, rows = candidates.candidates_from_run(run, "T1", nhf_status="disabled")
    assert all("nhf_unscreened" in r["flags"] and r["review_tier"] == 4 for r in rows)
    # with no metrics file the status is "none": unscreened as well
    _h, rows = candidates.candidates_from_run(run, "T1")
    assert all("nhf_unscreened" in r["flags"] for r in rows)


def test_ambiguous_and_mapq0(run):
    _h, rows = candidates.candidates_from_run(run, "T1", nhf_status="screened", classes=("SV", "AMBIGUOUS"),
                                              min_reads=3)
    t = by_locus(rows)
    amb = t["chr2:27000-27100"]
    assert amb["class"] == "AMBIGUOUS" and amb["review_tier"] == 4 and amb["flags"] == "mapq0_only"
    _h, rows = candidates.candidates_from_run(run, "T1", nhf_status="screened", classes=("SV", "AMBIGUOUS"),
                                              min_reads=4)
    assert "chr2:27000-27100" not in by_locus(rows)


def test_bed_reader_is_by_header(tmp_path):
    """An older BED (13 columns) reads with the missing columns as '.'; a broken one is refused."""
    path = str(tmp_path / "old.bed")
    with open(path, "w") as fh:
        fh.write("#chrom\tstart\tend\treads\tunique_kmers\tsplit_reads\tdiscordant_pairs\tmax_clip_len\tunmapped_mates\tclass\tbreakpoint_reads\tlarge_indel_reads\tsv_type\n")
        fh.write("chr1\t10\t20\t5\t7\t2\t1\t30\t0\tSV\t2\t0\tDEL\n")
    rows = candidates.read_bed(path)
    assert rows[0]["start"] == 10 and rows[0]["reads"] == 5 and "nhf" not in rows[0]
    _h, out = candidates.build_candidates(rows, [], "T", "screened")
    assert out[0]["max_mapq"] == "." and out[0]["nhf"] == "." and out[0]["flags"] == "." and out[0]["review_tier"] == 3
    with open(path, "w") as fh:
        fh.write("chr1\t10\t20\n")
    with pytest.raises(contract.ContractError):
        candidates.read_bed(path)
    with open(path, "w") as fh:
        fh.write("#chrom\tstart\tend\treads\n")
        fh.write("chr1\tten\t20\t5\n")
    with pytest.raises(contract.ContractError):
        candidates.read_bed(path)


def test_locus_flags():
    assert candidates.locus_flags("chr14", 105600000, 105600200) == "IGH"
    assert candidates.locus_flags("chrX", 100000, 100100) == "PAR"
    assert candidates.locus_flags("chrY", 3000000, 3000100, child_sex="2") == "chrY_female"
    assert candidates.locus_flags("chrUn_KI270742v1", 1, 2) == "nonprimary"
    assert candidates.locus_flags("chrEBV", 1, 2) == "chrEBV" and candidates.locus_flags("chrM", 1, 2) == "chrM"
    assert candidates.locus_flags("chr1", 1, 2) == "." and candidates.locus_flags("chr22", 1, 2) == "."
