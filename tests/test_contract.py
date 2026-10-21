"""The contract: the candidate table, junctions, windows, bundles, metrics and the verdict TSV."""

import json
import os

import pytest

from svreview import contract

HERE = os.path.dirname(__file__)

CANDIDATE_HEADER = "\t".join(contract.CANDIDATE_COLUMNS)
# a tier-1 pair and a tier-3 insertion, as the 1000 Genomes screen's table has them (its exact column order)
ROWS = [
    ["chr1", "9860", "10000", "sv_demo", "SV", "DEL", "11", "28", "6", "4", "0", "8", "0", "72", "60", "0", "1",
     "chr1:12000-12140|+-|DEL|15", ".", ".", "0", "1", "1", "7", ".", ".", ".", ".", "1"],
    ["chr1", "12000", "12140", "sv_demo", "SV", "DEL", "8", "28", "7", "3", "0", "7", "0", "54", "60", "0", "1",
     "chr1:9860-10000|+-|DEL|15", ".", ".", "0", "1", "1", "8", ".", ".", ".", ".", "1"],
    ["chr2", "21866", "22078", "sv_demo", "SV", "INS", "10", "211", "0", "1", "7", "9", "0", "114", "60", "0", "1",
     ".", ".", ".", "0", "0", "0", "10", "12", "0.8333", ".", "nhf", "4"],
    ["chr1", "54860", "55115", "sv_demo", "SV", "BND", "11", "60", "8", "5", "0", "9", "0", "73", "60", "0", "1",
     "chr2:14861-15141|+-|BND|10;chr2:14861-15141|-+|BND|8", ".", ".", "0", "2", "1", "5", ".", ".", ".", ".", "1"],
]


def write_table(path, header=CANDIDATE_HEADER, rows=ROWS):
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return path


def test_junctions_parse_and_format():
    js = contract.parse_junctions("chr2:14861-15141|+-|BND|10;chr2:14861-15141|-+|BND|8")
    assert len(js) == 2
    assert js[0] == contract.Junction("chr2", 14861, 15141, "+-", "BND", 10)
    assert js[1].strands == "-+" and js[1].reads == 8
    assert contract.format_junctions(js) == "chr2:14861-15141|+-|BND|10;chr2:14861-15141|-+|BND|8"
    assert contract.parse_junctions(".") == [] and contract.parse_junctions("") == []
    assert contract.format_junctions([]) == "."
    with pytest.raises(contract.ContractError):
        contract.parse_junctions("chr2:14861-15141|+-|BND")
    with pytest.raises(contract.ContractError):
        contract.parse_junctions("chr2:14861|+-|BND|10")


def test_read_candidates(tmp_path):
    path = write_table(str(tmp_path / "c.tsv"))
    header, cands = contract.read_candidates(path)
    assert header == list(contract.CANDIDATE_COLUMNS)
    assert [c.locus for c in cands] == ["chr1:9860-10000", "chr1:12000-12140", "chr2:21866-22078", "chr1:54860-55115"]
    c = cands[0]
    assert c.id == "sv_demo/chr1_9860_10000" and c.name == "chr1_9860_10000"
    assert c.cls == "SV" and c.sv_type == "DEL" and c.review_tier == 1 and c.rank == 7 and c.both_ends == 1
    assert c.evidence["reads"] == 11 and c.evidence["breakpoint_reads"] == 8 and c.evidence["max_mapq"] == 60
    assert c.junctions[0].chrom == "chr1" and c.partners() == [("chr1", 12000, 12140)]
    assert c.flags == [] and c.locus_flags == [] and c.nhf is None and c.nhf_reads is None
    ins = cands[2]
    assert ins.junctions == [] and ins.review_tier == 4 and ins.flags == ["nhf"]
    assert ins.nhf == pytest.approx(0.8333) and ins.nhf_reads == 12
    bnd = cands[3]
    assert bnd.partners() == [("chr2", 14861, 15141)]       # two junctions, one partner
    assert bnd.row["junctions"].count(";") == 1


def test_read_candidates_minimal_layout(tmp_path):
    """Only chrom, start, end and trio are required; the rest is optional."""
    path = str(tmp_path / "m.tsv")
    with open(path, "w") as fh:
        fh.write("trio\tchrom\tstart\tend\n")
        fh.write("T1\tchr3\t100\t200\n")
    header, cands = contract.read_candidates(path)
    c = cands[0]
    assert c.trio == "T1" and c.sv_type == "." and c.junctions == [] and c.review_tier is None and c.evidence == {}
    with open(path, "w") as fh:
        fh.write("chrom\tstart\tend\n")
        fh.write("chr3\t100\t200\n")
    with pytest.raises(contract.ContractError):
        contract.read_candidates(path)
    with open(path, "w") as fh:
        fh.write("trio\tchrom\tstart\tend\n")
        fh.write("T1\tchr3\t200\t100\n")
    with pytest.raises(contract.ContractError):
        contract.read_candidates(path)


def test_candidate_windows(tmp_path):
    path = write_table(str(tmp_path / "c.tsv"))
    _h, cands = contract.read_candidates(path)
    lengths = {"chr1": 60000, "chr2": 30000}
    # own and partner, padded, kept apart when they do not touch
    assert contract.candidate_windows(cands[0], 700, lengths) == [("chr1", 9160, 10700), ("chr1", 11300, 12840)]
    # the other end of the pair gets the same windows: the bundles can be shared
    assert contract.candidate_windows(cands[1], 700, lengths) == contract.candidate_windows(cands[0], 700, lengths)
    # overlapping windows merge
    assert contract.candidate_windows(cands[0], 2000, lengths) == [("chr1", 7860, 14140)]
    # clipped to the contig; the reference's order decides the panel order
    assert contract.candidate_windows(cands[3], 20000, lengths) == [("chr1", 34860, 60000), ("chr2", 0, 30000)]
    assert contract.candidate_windows(cands[3], 20000, {"chr2": 30000, "chr1": 60000}) == \
        [("chr2", 0, 30000), ("chr1", 34860, 60000)]
    # without lengths: no clipping at the end, chromosomes in order of appearance
    assert contract.candidate_windows(cands[3], 20000) == [("chr1", 34860, 75115), ("chr2", 0, 35141)]


def test_bundle_layout(tmp_path):
    bundles = tmp_path / "review"
    d = bundles / "sv_demo" / "chr1_9860_10000"
    d.mkdir(parents=True)
    for role in contract.ROLES:
        (d / (role + ".cram")).write_bytes(b"CRAM")
        (d / (role + ".cram.crai")).write_bytes(b"")
    (d / "child.bam").write_bytes(b"BAM")     # a converted member wins over the CRAM
    (d / "child.bam.bai").write_bytes(b"")
    info = {"trio": "sv_demo", "candidate": "chr1:9860-10000", "chrom": "chr1", "start": 9860, "end": 10000,
            "sv_type": "DEL", "review_tier": 1, "rank": 7, "windows": ["chr1:9161-10700", "chr1:11301-12840"],
            "skipped_windows": [], "pad": 700, "members": {"child": "child", "father": "father", "mother": "mother"},
            "status": {"child": "ok", "father": "ok", "mother": "ok"},
            "shared_by": [{"candidate": "chr1:12000-12140", "name": "chr1_12000_12140", "sv_type": "DEL",
                           "review_tier": 1, "rank": 8}]}
    (d / "candidate.json").write_text(json.dumps(info))
    os.symlink("chr1_9860_10000", str(bundles / "sv_demo" / "chr1_12000_12140"))
    path = write_table(str(tmp_path / "c.tsv"))
    _h, cands = contract.read_candidates(path)
    assert contract.bundle_dir(str(bundles), cands[0]) == str(d)
    assert contract.bundle_dir(str(bundles), cands[1]) == str(d)          # the shared end resolves to the first
    assert contract.bundle_dir(str(bundles), cands[2]) is None
    members = contract.bundle_members(str(d))
    assert members["child"].endswith("child.bam") and members["father"].endswith("father.cram")
    assert contract.alignment_index(members["child"]).endswith("child.bam.bai")
    assert contract.alignment_index(members["father"]).endswith("father.cram.crai")
    cj = contract.read_candidate_json(str(d / "candidate.json"))
    assert cj["windows0"] == [("chr1", 9160, 10700), ("chr1", 11300, 12840)]   # samtools regions -> 0-based
    # the index
    with open(str(bundles / "bundle.tsv"), "w") as fh:
        fh.write("\t".join(contract.BUNDLE_INDEX_COLUMNS) + "\n")
        fh.write("chr1:9860-10000\tsv_demo\tDEL\t1\t7\tchild\tchild\tchr1:9161-10700,chr1:11301-12840\t100\t1000\t.\tok\tsv_demo/chr1_9860_10000/child.cram\n")
    header, rows = contract.read_bundle_index(str(bundles / "bundle.tsv"))
    assert rows[0]["status"] == "ok" and rows[0]["member"] == "child"
    with open(str(bundles / "bundle.tsv"), "w") as fh:
        fh.write("candidate\ttrio\n")
    with pytest.raises(contract.ContractError):
        contract.read_bundle_index(str(bundles / "bundle.tsv"))


def test_metrics(tmp_path):
    m = {"mode": "discovery", "candidate_regions": 1,
         "regions": [{"chrom": "chr2", "start": 21866, "end": 22078, "class": "SV", "sv_type": "INS", "reads": 10,
                      "max_mapq": 60, "mapq0_reads": 0, "nhf_reads": 12, "nhf": 0.8333}],
         "nonhuman_fraction": {"status": "screened", "reason": "", "db": "/db", "regions_screened": 9,
                               "sampled": {"reads": 143, "fractions": {"total": 143, "nonhuman": 0.07}},
                               "regions": [{"chrom": "chr2", "start": 21866, "end": 22078, "reads": 12,
                                            "placed_unmapped": 7, "fractions": {"total": 12, "nonhuman": 0.8333},
                                            "top_taxon": "Escherichia_coli(10)"}]}}
    path = tmp_path / "m.json"
    path.write_text(json.dumps(m))
    metrics = contract.read_metrics(str(path))
    assert contract.metrics_region(metrics, "chr2", 21866, 22078)["sv_type"] == "INS"
    assert contract.metrics_region(metrics, "chr2", 21866, 22079) is None
    block = contract.nhf_block(metrics)
    assert block["status"] == "screened" and block["sampled"]["nonhuman"] == 0.07 and block["regions_screened"] == 9
    reg = contract.nhf_region(metrics, "chr2", 21866, 22078)
    assert reg["top_taxon"] == "Escherichia_coli(10)" and reg["fractions"]["nonhuman"] == 0.8333 and reg["reads"] == 12
    assert contract.nhf_region(metrics, "chr1", 1, 2) is None
    path.write_text(json.dumps({"mode": "vcf"}))
    with pytest.raises(contract.ContractError):
        contract.read_metrics(str(path))
    assert contract.nhf_block({"mode": "discovery"})["status"] == "none"


def test_verdicts_round_trip(tmp_path):
    path = str(tmp_path / "verdicts.tsv")
    rows = [
        {"candidate": "chr1:9860-10000", "trio": "sv_demo", "chrom": "chr1", "start": 9860, "end": 10000,
         "sv_type": "DEL", "verdict": "de_novo", "note": "clean junction\twith a tab", "reviewer": "me",
         "reviewed_at": "2026-10-21T10:00:00"},
        {"candidate": "chr2:21866-22078", "trio": "sv_demo", "chrom": "chr2", "start": 21866, "end": 22078,
         "sv_type": "INS", "verdict": "artifact", "note": "E. coli", "reviewer": "me", "reviewed_at": "2026-10-21T10:01:00"},
    ]
    contract.write_verdicts(path, rows)
    # append: a later row for the same candidate wins, an empty verdict clears
    contract.write_verdicts(path, [dict(rows[0], verdict="unsure", note="second look")], append=True)
    v = contract.read_verdicts(path)
    assert v[("sv_demo", "chr1", 9860, 10000)]["verdict"] == "unsure"
    assert v[("sv_demo", "chr1", 9860, 10000)]["note"] == "second look"
    assert v[("sv_demo", "chr2", 21866, 22078)]["verdict"] == "artifact"
    contract.write_verdicts(path, [dict(rows[1], verdict="")], append=True)
    assert ("sv_demo", "chr2", 21866, 22078) not in contract.read_verdicts(path)
    with open(path) as fh:
        assert fh.readline().rstrip("\n") == "\t".join(contract.VERDICT_COLUMNS)
        assert "clean junction with a tab" in fh.read()
    # merged back into the candidate table as columns the analysis reads as a flag
    table = write_table(str(tmp_path / "c.tsv"))
    header, trows = contract.read_tsv(table)
    out_header, out_rows, n = contract.merge_verdicts(header, trows, contract.read_verdicts(path))
    assert n == 1 and out_header[-4:] == list(contract.MERGED_VERDICT_COLUMNS)
    assert out_rows[0]["verdict"] == "unsure" and out_rows[0]["verdict_note"] == "second look"
    assert out_rows[1]["verdict"] == "." and out_rows[2]["verdict"] == "."
    with pytest.raises(contract.ContractError):
        contract.merge_verdicts(out_header, out_rows, {})
    # a verdict outside the vocabulary is refused
    with open(path, "a") as fh:
        fh.write("x\tsv_demo\tchr1\t1\t2\tDEL\tmaybe\t\t\t\n")
    with pytest.raises(contract.ContractError):
        contract.read_verdicts(path)
