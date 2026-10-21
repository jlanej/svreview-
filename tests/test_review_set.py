"""The review set (docs/FORMAT.md): opened, checked, built into a site; per-trio alignments; the shipped example."""

import json
import os

import pysam
import pytest

from svreview import contract, review_set, site

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, os.pardir, "examples", "sv_demo_review")
REF_LEN = 6000
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": REF_LEN}, {"SN": "chr2", "LN": REF_LEN}]}
CANDIDATES = """chrom\tstart\tend\ttrio\tclass\tsv_type\tjunctions\treview_tier\trank\tboth_ends
chr1\t950\t1050\tT1\tSV\tDEL\tchr1:3000-3100|+-|DEL|6\t1\t1\t1
chr1\t3000\t3100\tT1\tSV\tDEL\tchr1:950-1050|+-|DEL|6\t1\t2\t1
chr2\t4000\t4100\tT1\tSV\tINS\t.\t3\t3\t0
"""


def write_bam(path, pairs, dk=()):
    """*pairs*: ``[(name, chrom, pos, mate_pos, proper)]``; reads named in *dk* get dk:i:1."""
    header = pysam.AlignmentHeader.from_dict(dict(HEADER, RG=[{"ID": "x", "SM": "S_" + os.path.basename(path)[:-4]}]))
    unsorted = path + ".unsorted"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as fh:
        for name, chrom, pos, mate_pos, proper in pairs:
            for i, (p, mp, rev) in enumerate(((pos, mate_pos, False), (mate_pos, pos, True))):
                r = pysam.AlignedSegment(header)
                r.query_name = name
                r.flag = 0x1 | (0x2 if proper else 0) | (0x40 if i == 0 else 0x80) | (0x10 if rev else 0x20)
                r.reference_name = chrom
                r.reference_start = p
                r.mapping_quality = 60
                r.cigarstring = "100M"
                r.query_sequence = "A" * 100
                r.query_qualities = pysam.qualitystring_to_array("I" * 100)
                r.next_reference_name = chrom
                r.next_reference_start = mp
                r.template_length = (mp + 100 - p) if i == 0 else -(mp + 100 - p)
                if name in dk:
                    r.set_tag("dk", 1)
                fh.write(r)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return path


def normal_pairs(role, chrom="chr1"):
    return [("%s:%s:n%d" % (role, chrom, i), chrom, i, i + 250, True) for i in range(0, REF_LEN - 500, 100)]


@pytest.fixture
def trio_set(tmp_path):
    """A set with per-trio alignments and no bundles/: the GMKF arm's way."""
    d = tmp_path / "set"
    (d / "runs").mkdir(parents=True)
    (d / "candidates.tsv").write_text(CANDIDATES)
    for role in contract.ROLES:
        pairs = normal_pairs(role) + normal_pairs(role, "chr2")
        if role == "child":
            pairs += [("child:d%d" % k, "chr1", 900 + k * 10, 3000 + k * 10, False) for k in range(4)]
        write_bam(str(d / "runs" / (role + ".bam")), pairs)
    write_bam(str(d / "runs" / "informative.bam"), [("child:d0", "chr1", 900, 3000, False), ("child:d1", "chr1", 910, 3010, False)],
              dk={"child:d0", "child:d1"})
    review_set.write_trios(str(d / "trios.tsv"), [dict(trio="T1", child="KID", father="DAD", mother="MOM", child_sex="2",
                                                       child_alignments="runs/child.bam", father_alignments="runs/father.bam",
                                                       mother_alignments="runs/mother.bam", informative="runs/informative.bam")])
    (d / "svreview.json").write_text(json.dumps({"format": 1, "title": "per-trio set", "pad": 400}))
    return str(d)


def test_trios_tsv_round_trip(tmp_path):
    path = str(tmp_path / "trios.tsv")
    review_set.write_trios(path, [dict(trio="A", child="a1", father="a2", mother="a3"),
                                  dict(trio="B", child="b1", father="b2", mother="b3", child_sex="1", metrics="runs/B.json")])
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    assert header == ["trio", "child", "father", "mother", "child_sex", "metrics"]
    trios = review_set.read_trios(path)
    assert list(trios) == ["A", "B"] and trios["A"].samples["child"] == "a1" and trios["A"].child_sex == "0"
    assert trios["B"].metrics == os.path.join(str(tmp_path), "runs", "B.json") and trios["B"].has_alignments is False
    with open(path, "a") as fh:
        fh.write("A\tx\ty\tz\t\t\n")
    with pytest.raises(contract.ContractError):
        review_set.read_trios(path)
    with open(path, "w") as fh:
        fh.write("trio\tchild\n")
    with pytest.raises(contract.ContractError):
        review_set.read_trios(path)


def test_per_trio_alignments_set(trio_set):
    rs = review_set.ReviewSet.open(trio_set)
    assert rs.pad == 400 and rs.title == "per-trio set" and rs.bundles_dir is None
    assert [c.id for c in rs.candidates] == ["T1/chr1_950_1050", "T1/chr1_3000_3100", "T1/chr2_4000_4100"]
    kind, src = rs.source_of(rs.candidates[0])
    assert kind == "trio" and src["child"].endswith("runs/child.bam")
    windows, cj = rs.windows_of(rs.candidates[0])
    assert windows == [("chr1", 550, 1450), ("chr1", 2600, 3500)] and cj is None
    findings = rs.check(deep=True)
    assert not [f for f in findings if f.level == "error"], findings
    assert any("3 from per-trio alignments" in f.message for f in findings)
    out = os.path.join(trio_set, "site")
    m = site.build_site_from_set(rs, out, log=lambda *a: None)
    assert m["format"] == 1 and m["summary"]["with_bundle"] == 3
    c0, c1, c2 = m["candidates"]
    assert c0["source"] == "trio" and c0["bundle"] == "trio:T1"
    assert c0["members"]["child"]["bam"] == "T1/chr1_950_1050/child.bam" and c0["members"]["child"]["sample"] == "KID"
    assert c0["members"]["child"]["informative"] == 4           # two molecules, both ends, tagged in the cut
    assert c1["shared_with"] == "chr1:950-1050" and c1["members"]["child"]["bam"] == c0["members"]["child"]["bam"]
    assert c2["windows"] == [{"chrom": "chr2", "start": 3600, "end": 4500}] and c2["shared_with"] is None
    assert m["trios"]["T1"]["members"] == {"child": "KID", "father": "DAD", "mother": "MOM"} and m["trios"]["T1"]["child_sex"] == "2"
    with pysam.AlignmentFile(os.path.join(out, "T1", "chr1_950_1050", "child.bam")) as f:
        records = list(f.fetch(until_eof=True))
    names = [r.query_name for r in records]
    assert "child:d3" in names and "child:chr1:n1000" in names
    assert "child:chr1:n4500" not in names and "child:chr2:n1000" not in names   # the windows, not the whole file
    assert sorted({r.query_name for r in records if r.has_tag("dk")}) == ["child:d0", "child:d1"]
    assert os.path.isfile(os.path.join(out, "trios.tsv"))


def test_check_finds_what_is_missing(trio_set):
    rs = review_set.ReviewSet.open(trio_set)
    os.remove(os.path.join(trio_set, "runs", "mother.bam.bai"))
    errors = [f for f in rs.check() if f.level == "error"]
    assert errors and "no index" in errors[0].message
    os.remove(os.path.join(trio_set, "trios.tsv"))
    rs = review_set.ReviewSet.open(trio_set)
    errors = [f.message for f in rs.check() if f.level == "error"]
    assert any(m.startswith("missing") for m in errors) and sum("no bundle" in m for m in errors) == 3
    with open(os.path.join(trio_set, "svreview.json"), "w") as fh:
        json.dump({"format": 2}, fh)
    with pytest.raises(contract.ContractError):
        review_set.ReviewSet.open(trio_set)
    with pytest.raises(contract.ContractError):
        review_set.ReviewSet.open(os.path.join(trio_set, "nowhere"))


def test_check_warns_on_the_table(trio_set):
    with open(os.path.join(trio_set, "candidates.tsv"), "a") as fh:
        fh.write("chr2\t4000\t4100\tT1\tSV\tINS\t.\t3\t4\t0\n")                      # twice
        fh.write("chr2\t5000\t5100\tT1\tSV\tDEL\tchr2:5500-5600|+-|DEL|3\t1\t5\t1\n")  # both_ends without the partner row
        fh.write("chr1\t100\t200\tT2\tSV\tDEL\t.\t3\t6\t0\n")                        # a trio not in trios.tsv
    rs = review_set.ReviewSet.open(trio_set)
    findings = rs.check()
    warnings = [f.message for f in findings if f.level == "warning"]
    errors = [f.message for f in findings if f.level == "error"]
    assert any("listed twice" in w for w in warnings) and any("both_ends is 1" in w for w in warnings)
    assert any("trio T2 has no row" in e for e in errors) and any("no bundle" in e and "T2" in e for e in errors)
    text, n_err = review_set.format_findings(findings)
    assert n_err == 2 and text.endswith("2 error(s), 2 warning(s)")


def test_cli_check(trio_set, capsys):
    from svreview import cli
    assert cli.main(["check", trio_set]) == 0
    assert "0 error(s)" in capsys.readouterr().out
    os.remove(os.path.join(trio_set, "runs", "child.bam"))
    assert cli.main(["check", trio_set, "--deep"]) == 1


@pytest.mark.skipif(not os.path.isdir(EXAMPLE), reason="the shipped example set is not beside the tests")
def test_shipped_example_set(tmp_path):
    """examples/sv_demo_review: the demo's input, complete, and the site it builds matches docs/demo."""
    rs = review_set.ReviewSet.open(EXAMPLE)
    findings = rs.check(deep=True)
    assert not [f for f in findings if f.level == "error"], review_set.format_findings(findings)[0]
    assert len(rs.candidates) == 10 and set(c.review_tier for c in rs.candidates) == {1, 3, 4}
    assert rs.trios["sv_demo"].metrics.endswith("runs/sv_demo.metrics.json")
    # CRAM bundles with the reference embedded: they open with no FASTA
    kind, src = rs.source_of(rs.candidates[0])
    assert kind == "bundle" and contract.bundle_members(src)["child"].endswith("child.cram")
    m = site.build_site_from_set(rs, str(tmp_path / "site"), log=lambda *a: None)
    demo = os.path.join(HERE, os.pardir, "docs", "demo", "manifest.json")
    if os.path.isfile(demo):
        with open(demo) as fh:
            shipped = json.load(fh)
        got = [(c["id"], c["review_tier"], c["windows"], c["members"]["child"]["reads"], c["members"]["child"]["informative"])
               for c in m["candidates"]]
        want = [(c["id"], c["review_tier"], c["windows"], c["members"]["child"]["reads"], c["members"]["child"]["informative"])
                for c in shipped["candidates"]]
        assert got == want
