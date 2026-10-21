"""The Python side end to end on a tiny synthetic trio: the cut, the site, the server, the verdicts."""

import json
import os
import threading
import urllib.error
import urllib.request

import pysam
import pytest

from svreview import bundle, candidates, contract, serve, site

REF_LEN = 6000
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": REF_LEN}, {"SN": "chr2", "LN": REF_LEN}]}


def _pair(header, name, chrom, pos, mate_pos, mate_chrom=None, proper=True, tags1=None, cigar="100M", read_len=100):
    """One read pair: read 1 forward at pos, read 2 reverse at mate_pos."""
    mate_chrom = mate_chrom or chrom
    out = []
    for i, (c, p, mp, mc, rev) in enumerate(((chrom, pos, mate_pos, mate_chrom, False), (mate_chrom, mate_pos, pos, chrom, True))):
        r = pysam.AlignedSegment(header)
        r.query_name = name
        r.flag = 0x1 | (0x2 if proper else 0) | (0x40 if i == 0 else 0x80) | (0x10 if rev else 0x20)
        r.reference_name = c
        r.reference_start = p
        r.mapping_quality = 60
        r.cigarstring = cigar if i == 0 else "%dM" % read_len
        r.query_sequence = "A" * read_len
        r.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        r.next_reference_name = mc
        r.next_reference_start = mp
        r.template_length = (mp + read_len - p) if (c == mc and i == 0) else (-(mp + read_len - p) if c == mc else 0)
        if i == 0:
            for k, v in (tags1 or {}).items():
                r.set_tag(k, v)
        out.append(r)
    return out


def write_trio(d):
    """child.bam / father.bam / mother.bam: normal pairs everywhere, in the child also discordant pairs and a
    split read joining chr1:1000 to chr1:3000 (a deletion), plus a pair whose mate lies far away (chr2)."""
    header = pysam.AlignmentHeader.from_dict(dict(HEADER, RG=[{"ID": "x", "SM": "SAMPLE"}]))
    paths = {}
    for role in ("child", "father", "mother"):
        path = os.path.join(d, role + ".bam")
        records = []
        for i in range(0, REF_LEN - 500, 50):
            records += _pair(header, "%s:n%d" % (role, i), "chr1", i, i + 250)
        if role == "child":
            for k in range(5):
                records += _pair(header, "child:d%d" % k, "chr1", 800 + k * 10, 3050 + k * 10, proper=False)
            records += _pair(header, "child:s0", "chr1", 950, 1150, proper=True, cigar="50M50S",
                             tags1={"SA": "chr1,3001,+,50S50M,60,0;"})
            records += _pair(header, "child:far", "chr1", 1200, 100, mate_chrom="chr2", proper=False)
        unsorted = path + ".unsorted"
        with pysam.AlignmentFile(unsorted, "wb", header=header) as fh:
            for r in records:
                fh.write(r)
        pysam.sort("-o", path, unsorted)
        os.remove(unsorted)
        pysam.index(path)
        paths[role] = path
    return paths


def write_informative(d, names):
    header = pysam.AlignmentHeader.from_dict(HEADER)
    path = os.path.join(d, "informative.bam")
    with pysam.AlignmentFile(os.path.join(d, "child.bam"), "rb") as src:
        unsorted = path + ".unsorted"
        with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
            for r in src.fetch(until_eof=True):
                if r.query_name in names:
                    r.set_tag("dk", 1)
                    out.write(r)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return path


BED = """#chrom\tstart\tend\treads\tunique_kmers\tsplit_reads\tdiscordant_pairs\tmax_clip_len\tunmapped_mates\tclass\tbreakpoint_reads\tlarge_indel_reads\tsv_type\tmax_mapq\tmapq0_reads\tnhf_reads\tnhf
chr1\t950\t1050\t6\t20\t1\t5\t50\t0\tSV\t2\t0\tDEL\t60\t0\t6\t0.0000
chr1\t3000\t3100\t6\t20\t1\t5\t50\t0\tSV\t2\t0\tDEL\t60\t0\t6\t0.0000
"""
BEDPE = """#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tsupporting_reads\tstrand1\tstrand2\tsv_type
chr1\t950\t1050\tchr1\t3000\t3100\tSV_1\t6\t+\t-\tDEL
"""


@pytest.fixture
def run(tmp_path):
    d = str(tmp_path)
    paths = write_trio(d)
    prefix = os.path.join(d, "trio.discovery")
    with open(prefix + ".bed", "w") as fh:
        fh.write(BED)
    with open(prefix + ".sv.bedpe", "w") as fh:
        fh.write(BEDPE)
    metrics = {"mode": "discovery", "candidate_regions": 2,
               "regions": [{"chrom": "chr1", "start": 950, "end": 1050, "class": "SV", "sv_type": "DEL", "reads": 6, "max_mapq": 60},
                           {"chrom": "chr1", "start": 3000, "end": 3100, "class": "SV", "sv_type": "DEL", "reads": 6, "max_mapq": 60}],
               "nonhuman_fraction": {"status": "screened", "reason": "", "db": "/db", "regions_screened": 2,
                                     "sampled": {"reads": 10, "fractions": {"total": 10, "nonhuman": 0.0}},
                                     "regions": [{"chrom": "chr1", "start": 950, "end": 1050, "reads": 6, "placed_unmapped": 0,
                                                  "fractions": {"total": 6, "nonhuman": 0.0}, "top_taxon": "."}]}}
    with open(prefix + ".metrics.json", "w") as fh:
        json.dump(metrics, fh)
    informative = write_informative(d, {"child:d0", "child:d1", "child:s0"})
    header, rows = candidates.candidates_from_run(prefix, "TRIO")
    table = os.path.join(d, "candidates.tsv")
    contract.write_tsv(table, header, rows)
    return dict(dir=d, paths=paths, prefix=prefix, table=table, informative=informative)


def test_bundle_site_and_verdicts(run):
    _h, cands = contract.read_candidates(run["table"])
    assert [c.review_tier for c in cands] == [1, 1]
    bundles = os.path.join(run["dir"], "review")
    n, index = bundle.bundle_candidates(cands, run["paths"], bundles, pad=300, log=lambda *a: None)
    assert n == 1                                            # the two ends share one bundle
    assert os.path.islink(os.path.join(bundles, "TRIO", "chr1_3000_3100"))
    header, rows = contract.read_bundle_index(index)
    assert len(rows) == 6 and {r["status"] for r in rows} == {"ok", "shared=chr1_950_1050 ok"}
    assert rows[0]["sample"] == "SAMPLE"                     # the SM of the read group
    cj = contract.read_candidate_json(os.path.join(bundles, "TRIO", "chr1_950_1050", "candidate.json"))
    assert cj["windows0"] == [("chr1", 650, 1350), ("chr1", 2700, 3400)]
    assert cj["shared_by"][0]["name"] == "chr1_3000_3100"
    # the cut holds the discordant pairs' both ends and the far mate on chr2 (samtools view -P)
    with pysam.AlignmentFile(os.path.join(bundles, "TRIO", "chr1_950_1050", "child.bam")) as f:
        names = [r.query_name for r in f.fetch(until_eof=True)]
        assert names.count("child:d0") == 2 and names.count("child:far") == 2 and names.count("child:s0") == 2
        assert sum(1 for r in f.fetch("chr2", 0, 1000)) == 1
    # the site
    out = os.path.join(run["dir"], "site")
    manifest = site.build_site(run["table"], bundles, out, metrics={"TRIO": run["prefix"] + ".metrics.json"},
                               informative={"*": run["informative"]}, title="tiny", log=lambda *a: None)
    assert manifest["contract"] == 1 and manifest["summary"]["with_bundle"] == 2
    c0, c1 = manifest["candidates"]
    assert c0["id"] == "TRIO/chr1_950_1050" and c1["shared_with"] == "chr1:950-1050"
    assert c0["members"]["child"]["bam"] == "TRIO/chr1_950_1050/child.bam"
    assert c0["members"]["child"]["informative"] == 6             # three molecules, two alignments each
    assert c0["members"]["father"]["informative"] is None
    assert c0["regions"][0]["nhf_screen"]["fractions"]["nonhuman"] == 0.0 and c0["regions"][1]["nhf_screen"] is None
    assert manifest["trios"]["TRIO"]["nhf"]["status"] == "screened"
    assert [w["start"] for w in c0["windows"]] == [650, 2700]
    with pysam.AlignmentFile(os.path.join(out, "TRIO", "chr1_950_1050", "child.bam")) as f:
        tagged = sorted({r.query_name for r in f.fetch(until_eof=True) if r.has_tag("dk")})
        assert tagged == ["child:d0", "child:d1", "child:s0"]
    with open(os.path.join(out, "manifest.json")) as fh:
        assert json.load(fh)["title"] == "tiny"
    assert os.path.isfile(os.path.join(out, "candidates.tsv"))
    # the server: a range, the manifest, a POSTed verdict
    srv = serve.Server(out, port=0, quiet=True)
    port = srv.server_address[1]
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    try:
        base = "http://127.0.0.1:%d/" % port
        bam = os.path.join(out, "TRIO", "chr1_950_1050", "child.bam")
        req = urllib.request.Request(base + "TRIO/chr1_950_1050/child.bam", headers={"Range": "bytes=10-29"})
        with urllib.request.urlopen(req) as res:
            assert res.status == 206 and res.headers["Content-Range"] == "bytes 10-29/%d" % os.path.getsize(bam)
            with open(bam, "rb") as fh:
                assert res.read() == fh.read()[10:30]
        with urllib.request.urlopen(base + "manifest.json") as res:
            assert res.status == 200 and json.load(res)["contract"] == 1
        req = urllib.request.Request(base + "TRIO/chr1_950_1050/child.bam", headers={"Range": "bytes=999999-"})
        with pytest.raises(urllib.error.HTTPError) as err:
            urllib.request.urlopen(req)
        assert err.value.code == 416
        with pytest.raises(urllib.error.HTTPError) as err:
            urllib.request.urlopen(base + "missing.json")
        assert err.value.code == 404
        with pytest.raises(urllib.error.HTTPError) as err:
            urllib.request.urlopen(base + "verdicts.tsv")
        assert err.value.code == 404                         # none yet
        body = "chr1:950-1050\tTRIO\tchr1\t950\t1050\tDEL\tde_novo\tboth ends\tme\t2026-10-21T12:00:00\n"
        req = urllib.request.Request(base + "verdicts.tsv", data=body.encode(), method="POST")
        with urllib.request.urlopen(req) as res:
            assert json.load(res)["verdicts"] == 1
        body = "\t".join(contract.VERDICT_COLUMNS) + "\nchr1:950-1050\tTRIO\tchr1\t950\t1050\tDEL\tunsure\tsecond look\tme\t2026-10-21T12:01:00\n"
        req = urllib.request.Request(base + "verdicts.tsv", data=body.encode(), method="POST")
        with urllib.request.urlopen(req) as res:
            assert json.load(res) == {"appended": 1, "verdicts": 1, "path": "verdicts.tsv"}
        req = urllib.request.Request(base + "verdicts.tsv", data=b"x\tTRIO\tchr1\t950\t1050\tDEL\tmaybe\t\t\t\n", method="POST")
        with pytest.raises(urllib.error.HTTPError) as err:
            urllib.request.urlopen(req)
        assert err.value.code == 400
        req = urllib.request.Request(base + "other.tsv", data=b"", method="POST")
        with pytest.raises(urllib.error.HTTPError) as err:
            urllib.request.urlopen(req)
        assert err.value.code == 404
        with urllib.request.urlopen(base + "verdicts.tsv") as res:
            text = res.read().decode()
        assert text.count("\n") == 3 and text.startswith("\t".join(contract.VERDICT_COLUMNS))
    finally:
        srv.shutdown()
        srv.server_close()
    # the verdicts read back and merged: the last row wins
    verdicts = contract.read_verdicts(os.path.join(out, "verdicts.tsv"))
    assert verdicts[("TRIO", "chr1", 950, 1050)]["verdict"] == "unsure"
    header, rows = contract.read_tsv(os.path.join(out, "candidates.tsv"))
    out_header, out_rows, n = contract.merge_verdicts(header, rows, verdicts)
    assert n == 1 and out_rows[0]["verdict"] == "unsure" and out_rows[0]["verdict_note"] == "second look"
    assert out_rows[1]["verdict"] == "."


def test_site_without_bundle_and_cram_needing_reference(run, tmp_path):
    """A candidate with no bundle is listed, not drawn; a CRAM without its reference is a clear error."""
    out = os.path.join(run["dir"], "site2")
    manifest = site.build_site(run["table"], os.path.join(run["dir"], "nowhere"), out, log=lambda *a: None)
    assert manifest["summary"]["with_bundle"] == 0
    assert manifest["candidates"][0]["bundle"] is None and "no bundle" in manifest["candidates"][0]["missing"]
    assert manifest["candidates"][0]["windows"][0]["start"] == 0    # padded by the default, clipped at 0
    assert manifest["trios"]["TRIO"]["nhf"]["status"] == "none"
    with pytest.raises(contract.ContractError):
        site._open_alignments(os.path.join(run["dir"], "missing.cram"))


def test_cli_round_trip(run, capsys):
    from svreview import cli
    out = os.path.join(run["dir"], "cli")
    bundles = os.path.join(run["dir"], "review")
    cli.main(["bundle", "--candidates", run["table"], "--child", run["paths"]["child"], "--father", run["paths"]["father"],
              "--mother", run["paths"]["mother"], "--out-dir", bundles, "--pad", "300", "--tier", "1"])
    cli.main(["site", "--candidates", run["table"], "--bundles", bundles, "--out", out, "--metrics",
              "TRIO=" + run["prefix"] + ".metrics.json", "--informative", run["informative"]])
    assert os.path.isfile(os.path.join(out, "manifest.json"))
    verdicts = os.path.join(run["dir"], "v.tsv")
    contract.write_verdicts(verdicts, [{"candidate": "chr1:950-1050", "trio": "TRIO", "chrom": "chr1", "start": 950, "end": 1050,
                                        "sv_type": "DEL", "verdict": "inherited", "note": "", "reviewer": "", "reviewed_at": ""}])
    merged = os.path.join(run["dir"], "merged.tsv")
    cli.main(["verdicts", "--verdicts", verdicts, "--candidates", run["table"], "--out", merged])
    header, rows = contract.read_tsv(merged)
    assert rows[0]["verdict"] == "inherited" and header[-4:] == list(contract.MERGED_VERDICT_COLUMNS)
    captured = capsys.readouterr()
    assert "1 bundles" in captured.out and "1 of 2 rows" in captured.out
    with pytest.raises(SystemExit):
        cli.main(["verdicts", "--verdicts", verdicts, "--out", merged])      # --out without --candidates
