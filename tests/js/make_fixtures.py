#!/usr/bin/env python3
"""Write the BAM-reader fixtures: a bundle member's BAM + .bai copied beside an `expected` JSON of what pysam
reads in each window (name, flag, chrom, pos, end, mapq, cigar, tlen, mate, SA/MC/MQ/NM/dk tags), so the
browser's reader can be checked against htslib. Run from the repository root after the demo is built:

    python tests/js/make_fixtures.py docs/demo/sv_demo/chr1_54860_55115 chr1:54160-55815 chr2:14161-15841
"""
import json
import os
import shutil
import sys

import pysam


def main(argv):
    bundle, windows = argv[0], argv[1:]
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(bundle.rstrip("/"))
    expected = {"bundle": name, "windows": [], "members": {}}
    for role in ("child", "father", "mother"):
        src = os.path.join(bundle, role + ".bam")
        dst = os.path.join(out, "%s.%s.bam" % (name, role))
        shutil.copyfile(src, dst)
        shutil.copyfile(src + ".bai", dst + ".bai")
        member = {"bam": os.path.basename(dst), "windows": []}
        with pysam.AlignmentFile(src, "rb") as f:
            member["refs"] = [{"name": c, "length": n} for c, n in zip(f.references, f.lengths)]
            member["total"] = sum(1 for _ in f.fetch(until_eof=True))
            for w in windows:
                chrom, span = w.rsplit(":", 1)
                start, end = [int(x) for x in span.split("-")]
                reads = []
                for r in f.fetch(chrom, start, end):
                    reads.append({
                        "name": r.query_name, "flag": r.flag, "chrom": r.reference_name, "pos": r.reference_start,
                        "end": r.reference_end if not r.is_unmapped and r.cigartuples else r.reference_start + 1,
                        "mapq": r.mapping_quality, "cigar": r.cigarstring or "*", "tlen": r.template_length,
                        "mate_chrom": r.next_reference_name, "mate_pos": r.next_reference_start, "lseq": r.query_length,
                        "seq": r.query_sequence if len(reads) < 12 else None,   # a dozen sequences a window: enough to check the decoding
                        "tags": {k: v for k, v in r.get_tags() if k in ("SA", "MC", "MQ", "NM", "dk")},
                    })
                member["windows"].append({"chrom": chrom, "start": start, "end": end, "reads": reads})
        expected["members"][role] = member
    for w in windows:
        chrom, span = w.rsplit(":", 1)
        start, end = [int(x) for x in span.split("-")]
        expected["windows"].append({"chrom": chrom, "start": start, "end": end})
    path = os.path.join(out, name + ".expected.json")
    with open(path, "w") as fh:
        json.dump(expected, fh, indent=0)
    print("%s: %s" % (path, ", ".join("%s %d reads" % (r, m["total"]) for r, m in expected["members"].items())))


if __name__ == "__main__":
    main(sys.argv[1:])
