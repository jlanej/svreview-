"""candidates.py - the candidate table of ONE kmer-discovery run, in the contract's layout.

For a cohort the analysis makes this table (kmer_denovo_filter's ``examples/1kg_trios/summarize_trios.py``
across the trios, then the review stage's ``flag_candidates.py``). For one run - the demo, a bench trio at
home, one trio of a cohort - svreview makes it by the same rules, minus what only a cohort knows
(recurrence across children: ``n_trios_at_locus`` 1, ``partner_recurrent`` 0, ``prior_recurrence`` ``.``):

* one row per region of the classes asked for (``SV`` by default; ``AMBIGUOUS`` with at least
  ``--min-reads`` reads when asked), the evidence columns as the BED has them;
* ``junctions``: every BEDPE line with an end at the region, as ``partner|strands|type|reads``;
* ``locus_flags`` / ``partner_flags``: the GRCh38 immunoglobulin and T-cell receptor loci (a lymphoblastoid
  line's own V(D)J rearrangement is its strongest "de novo" hit), the pseudoautosomal regions, ``chrEBV``,
  ``chrM``, ``nonprimary`` (any contig that is not a nuclear chromosome) and ``chrY_female``;
* ``partner_multiplicity`` (this child's junction ends in the busiest partner window: hundreds mean a mapping
  sink), ``both_ends`` (a partner that is itself a candidate), ``rank`` (1 = strongest: both ends, then a
  junction, then breakpoint reads, then junction reads, then reads);
* ``nhf_reads`` / ``nhf`` from the BED (``.`` = unscreened, not 0);
* ``flags`` and ``review_tier`` as flag_candidates.py sets them: the locus flags, the partner flags with a
  ``partner:`` prefix, ``partner_sink``, ``mapq0_only``, ``nhf`` (at or above the fraction over at least the
  reads) or ``nhf_unscreened`` (the run's screen is not ``screened``); tier 1 clean with both ends, 2 clean
  with a junction, 3 clean without, 4 flagged.

Nothing is dropped. Coordinates stay the BED's (0-based, half-open).
"""

import collections
import re

from . import contract

#: GRCh38, 0-based half-open: the loci a reader sets aside with one filter, never dropped (summarize_trios.py).
SPECIAL_LOCI = collections.OrderedDict([
    ("IGH", ("chr14", 105586437, 106879844)),
    ("IGK", ("chr2", 88857361, 90235368)),
    ("IGL", ("chr22", 22026076, 22922913)),
    ("TRA", ("chr14", 21621904, 22552132)),
    ("TRB", ("chr7", 142299011, 142813287)),
    ("TRG", ("chr7", 38240024, 38368055)),
    ("PAR1_X", ("chrX", 10000, 2781479)),
    ("PAR2_X", ("chrX", 155701382, 156030895)),
    ("PAR1_Y", ("chrY", 10000, 2781479)),
    ("PAR2_Y", ("chrY", 56887902, 57217415)),
])
PRIMARY_CONTIG_RE = re.compile(r"^chr([0-9]{1,2}|X|Y)$")
CLASSES = ("SV", "AMBIGUOUS", "SMALL")
NHF_SCREENED = "screened"
#: Columns of the table this module writes, in order (the contract's, with every column present).
HEADER = list(contract.CANDIDATE_COLUMNS)


def locus_flags(chrom, start, end, child_sex=None):
    """The flags of an interval, ``;``-joined, or ``.`` when it has none (summarize_trios.py's rule)."""
    flags = []
    for name, (c, s, e) in SPECIAL_LOCI.items():
        if chrom == c and start < e and s < end:
            flag = "PAR" if name.startswith("PAR") else name
            if flag not in flags:
                flags.append(flag)
    if chrom == "chrEBV":
        flags.append("chrEBV")
    elif chrom == "chrM":
        flags.append("chrM")
    elif not PRIMARY_CONTIG_RE.match(chrom):
        flags.append("nonprimary")
    if chrom == "chrY" and str(child_sex) == "2":
        flags.append("chrY_female")
    return ";".join(flags) if flags else contract.MISSING


def read_bed(path):
    """The discovery BED by its header names: ``[{column: value}]`` with the coordinates and counts as
    ints, ``nhf`` a float, ``.`` kept for an unscreened region or a column the BED lacks."""
    rows, header = [], None
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#filters:") or not line.strip():
                continue
            if line.startswith("#"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if header is None:
                raise contract.ContractError("%s: no '#chrom' header line before line %d" % (path, n))
            f = line.rstrip("\n").split("\t")
            if len(f) != len(header):
                raise contract.ContractError("%s line %d: %d fields, the header has %d" % (path, n, len(f), len(header)))
            r = collections.OrderedDict(zip(header, f))
            for col in ("start", "end"):
                r[col] = contract.to_int(r.get(col))
                if r[col] is None:
                    raise contract.ContractError("%s line %d: %s is not an integer" % (path, n, col))
            for col in contract.EVIDENCE_COLUMNS + ("nhf_reads",):
                if col in r and r[col] != contract.MISSING:
                    v = contract.to_int(r[col])
                    if v is None:
                        raise contract.ContractError("%s line %d: %s %r is not a count" % (path, n, col, r[col]))
                    r[col] = v
            if "nhf" in r and r["nhf"] != contract.MISSING:
                v = contract.to_float(r["nhf"])
                if v is None:
                    raise contract.ContractError("%s line %d: nhf %r is not a fraction" % (path, n, r["nhf"]))
                r["nhf"] = v
            if r.get("class", "SV") not in CLASSES:
                raise contract.ContractError("%s line %d: class %r" % (path, n, r.get("class")))
            rows.append(r)
    if header is None:
        raise contract.ContractError("%s: no '#chrom' header line" % path)
    return rows


def read_bedpe(path):
    """The discovery BEDPE: ``[{chrom1, start1, end1, chrom2, start2, end2, sv_id, reads, strand1, strand2,
    sv_type}]``; ``[]`` for a missing file."""
    import os
    links = []
    if not os.path.isfile(path):
        return links
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                raise contract.ContractError("%s line %d: %d fields, the BEDPE has 11" % (path, n, len(f)))
            try:
                links.append(collections.OrderedDict(
                    chrom1=f[0], start1=int(f[1]), end1=int(f[2]), chrom2=f[3], start2=int(f[4]), end2=int(f[5]),
                    sv_id=f[6], reads=int(f[7]), strand1=f[8], strand2=f[9], sv_type=f[10]))
            except ValueError:
                raise contract.ContractError("%s line %d: a coordinate or count is not an integer" % (path, n))
    return links


def links_for_region(region, links):
    key = (region["chrom"], region["start"], region["end"])
    out = []
    for l in links:
        if (l["chrom1"], l["start1"], l["end1"]) == key or (l["chrom2"], l["start2"], l["end2"]) == key:
            out.append(l)
    return out


def _ends_in(links, chrom, start, end):
    """Junction ends (of this child) overlapping [start, end)."""
    n = 0
    for l in links:
        for c, s, e in ((l["chrom1"], l["start1"], l["end1"]), (l["chrom2"], l["start2"], l["end2"])):
            if c == chrom and s < end and start < e:
                n += 1
    return n


def flag_and_tier(row, nhf_status, sink_multiplicity=10, nhf_flag=0.05, nhf_min_reads=5):
    """flag_candidates.py's rule on one row (strings as the table has them); sets ``flags`` and
    ``review_tier`` and returns the row."""
    flags = []
    for flag in contract.split_flags(row.get("locus_flags")):
        if flag not in flags:
            flags.append(flag)
    for flag in contract.split_flags(row.get("partner_flags")):
        if "partner:" + flag not in flags:
            flags.append("partner:" + flag)
    if contract.to_int(row.get("partner_recurrent"), 0) > 0:
        flags.append("partner_recurrent")
    if contract.to_int(row.get("partner_multiplicity"), 0) >= sink_multiplicity:
        flags.append("partner_sink")
    prior = contract.to_int(row.get("prior_recurrence"))
    if prior is not None and prior >= 3:
        flags.append("prior_recurrent")
    max_mapq = row.get("max_mapq", "")
    if max_mapq not in ("", contract.MISSING) and contract.to_int(max_mapq, -1) == 0:
        flags.append("mapq0_only")
    nhf = contract.to_float(row.get("nhf"))
    if nhf_status != NHF_SCREENED:
        flags.append("nhf_unscreened")
    elif nhf is not None and nhf >= nhf_flag and contract.to_int(row.get("nhf_reads"), 0) >= nhf_min_reads:
        flags.append("nhf")
    row["flags"] = ";".join(flags) if flags else contract.MISSING
    has_junction = row.get("junctions", contract.MISSING) not in ("", contract.MISSING)
    if flags:
        tier = 4
    elif has_junction and contract.to_int(row.get("both_ends"), 0) == 1:
        tier = 1
    elif has_junction:
        tier = 2
    else:
        tier = 3
    row["review_tier"] = tier
    return row


def build_candidates(bed, links, trio, nhf_status, classes=("SV",), min_reads=3, child_sex=None, slack=500,
                     sink_multiplicity=10, nhf_flag=0.05, nhf_min_reads=5):
    """``(HEADER, rows)``: the candidate table of one run (rows as ordered dicts of strings, in the table's
    order: by chrom, then start)."""
    picked = []
    for r in bed:
        cls = r.get("class", "SV")
        if cls not in classes:
            continue
        if cls != "SV" and contract.to_int(r.get("reads"), 0) < min_reads:
            continue
        picked.append(r)
    keys = set((r["chrom"], r["start"], r["end"]) for r in picked)
    rows = []
    for r in picked:
        key = (r["chrom"], r["start"], r["end"])
        rl = links_for_region(r, links)
        partners = []
        for l in rl:
            this_end = (l["chrom1"], l["start1"], l["end1"])
            partners.append((l["chrom2"], l["start2"], l["end2"]) if this_end == key
                            else (l["chrom1"], l["start1"], l["end1"]))
        row = collections.OrderedDict()
        row["chrom"], row["start"], row["end"] = r["chrom"], r["start"], r["end"]
        row["trio"] = trio
        row["class"], row["sv_type"] = r.get("class", contract.MISSING), r.get("sv_type", contract.MISSING)
        for col in contract.EVIDENCE_COLUMNS:
            row[col] = r.get(col, contract.MISSING)
        row["n_trios_at_locus"] = 1
        row["junctions"] = ";".join(
            "%s:%d-%d|%s%s|%s|%d" % (p[0], p[1], p[2], l["strand1"], l["strand2"], l["sv_type"], l["reads"])
            for l, p in zip(rl, partners)) or contract.MISSING
        row["locus_flags"] = locus_flags(r["chrom"], r["start"], r["end"], child_sex)
        partner_flags = []
        for p in partners:
            for flag in contract.split_flags(locus_flags(p[0], p[1], p[2], child_sex)):
                if flag not in partner_flags:
                    partner_flags.append(flag)
        row["partner_flags"] = ";".join(partner_flags) if partner_flags else contract.MISSING
        row["partner_recurrent"] = 0
        row["partner_multiplicity"] = max((_ends_in(links, p[0], p[1] - slack, p[2] + slack) for p in partners),
                                          default=0)
        row["both_ends"] = int(any(p in keys for p in partners))
        row["rank"] = ""
        for col in ("nhf_reads", "nhf"):
            row[col] = r.get(col, contract.MISSING)
        row["prior_recurrence"] = contract.MISSING
        row["_evidence"] = (row["both_ends"], int(bool(rl)), contract.to_int(r.get("breakpoint_reads"), 0),
                            max((l["reads"] for l in rl), default=0), contract.to_int(r.get("reads"), 0))
        rows.append(row)
    for rank, row in enumerate(sorted(rows, key=lambda x: x["_evidence"], reverse=True), 1):
        row["rank"] = rank
    for row in rows:
        del row["_evidence"]
        flag_and_tier(row, nhf_status, sink_multiplicity, nhf_flag, nhf_min_reads)
    rows.sort(key=lambda x: (x["chrom"], x["start"]))
    return HEADER, rows


def candidates_from_run(prefix, trio, classes=("SV",), min_reads=3, child_sex=None, nhf_status=None, **kw):
    """The table from ``<prefix>.bed``, ``<prefix>.sv.bedpe`` and ``<prefix>.metrics.json`` (the screen's
    status; *nhf_status* overrides it, e.g. ``screened`` to review a run without a database as if screened -
    say so in the note)."""
    bed = read_bed(prefix + ".bed")
    links = read_bedpe(prefix + ".sv.bedpe")
    if nhf_status is None:
        import os
        mpath = prefix + ".metrics.json"
        nhf_status = contract.nhf_block(contract.read_metrics(mpath))["status"] if os.path.isfile(mpath) else "none"
    return build_candidates(bed, links, trio, nhf_status, classes, min_reads, child_sex, **kw)
