"""contract.py - what the analysis writes and svreview reads, and what svreview writes back.

Nothing here is invented: the layouts are the ones kmer_denovo_filter's discovery mode and its 1000 Genomes
trio-screen example write (``examples/1kg_trios/summarize_trios.py``) and the review stages around them
(``flag_candidates.py``, ``bundle_candidates.py``), read by their column names so a column added later
costs nothing. Coordinates are kept as the files have them: a candidate's ``start``/``end`` are the BED's
(0-based, half-open); a bundle's ``windows`` are samtools regions (1-based, inclusive) and are converted
to 0-based half-open on the way in, so the viewer sees one convention.

The pieces, in the order a review meets them:

* **the candidate table** - one row per candidate region of a child: ``chrom start end trio class sv_type``,
  the evidence counts, ``junctions`` (``partner|strands|type|reads`` joined by ``;``), ``locus_flags``,
  ``partner_flags``, ``flags``, ``review_tier`` (1 clean with both ends, 2 clean with a junction, 3 clean
  without, 4 flagged) and ``rank`` (1 = strongest), ``nhf_reads`` / ``nhf`` (``.`` = unscreened, not 0).
  Only ``chrom start end trio`` are required; the rest is read when present.
* **a bundle** - ``<trio>/<chrom>_<start>_<end>/`` with ``child father mother`` mini-CRAMs or BAMs and their
  indexes, ``candidate.json`` (the windows cut, the members, who shares the bundle) and, in the bundle
  root, ``bundle.tsv`` (one row per member per candidate; ``status`` ``ok``, ``shared=<dir> ok``,
  ``failed: ...``, ``absent``).
* **the run's metrics** - ``<prefix>.metrics.json``: ``regions`` (the BED's rows) and ``nonhuman_fraction``
  (``status``, the sampled fractions, one entry per screened region with its fractions and top taxon).
* **the verdict TSV** - what svreview writes back: one row per reviewed candidate, keyed by
  ``trio chrom start end``, with ``verdict`` (``de_novo inherited artifact unsure``), ``note``,
  ``reviewer``, ``reviewed_at``; appended to, the last row for a candidate wins.
"""

import collections
import json
import os

ROLES = ("child", "father", "mother")

#: The 1000 Genomes screen's candidate table (summarize_trios.py, then flag_candidates.py's last three).
CANDIDATE_COLUMNS = (
    "chrom", "start", "end", "trio", "class", "sv_type",
    "reads", "unique_kmers", "split_reads", "discordant_pairs", "unmapped_mates", "breakpoint_reads",
    "large_indel_reads", "max_clip_len", "max_mapq", "mapq0_reads",
    "n_trios_at_locus", "junctions", "locus_flags", "partner_flags", "partner_recurrent",
    "partner_multiplicity", "both_ends", "rank", "nhf_reads", "nhf",
    "prior_recurrence", "flags", "review_tier",
)
CANDIDATE_REQUIRED = ("chrom", "start", "end", "trio")
EVIDENCE_COLUMNS = ("reads", "unique_kmers", "split_reads", "discordant_pairs", "unmapped_mates",
                    "breakpoint_reads", "large_indel_reads", "max_clip_len", "max_mapq", "mapq0_reads")
BUNDLE_INDEX_COLUMNS = ("candidate", "trio", "sv_type", "review_tier", "rank", "member", "sample", "region",
                        "reads", "bytes", "source_url", "status", "path")
VERDICTS = ("de_novo", "inherited", "artifact", "unsure")
VERDICT_COLUMNS = ("candidate", "trio", "chrom", "start", "end", "sv_type", "verdict", "note", "reviewer",
                   "reviewed_at")
MERGED_VERDICT_COLUMNS = ("verdict", "verdict_note", "verdict_reviewer", "verdict_at")
MISSING = "."

Junction = collections.namedtuple("Junction", "chrom start end strands sv_type reads")


class ContractError(ValueError):
    """A file is not what the contract says it is; the message names the file and the field."""


# ---- small parsers ------------------------------------------------------------------------------------------

def to_int(value, default=None):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def to_float(value, default=None):
    if value in (None, "", MISSING):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def split_flags(text):
    """``"IGH;partner:PAR"`` -> ``["IGH", "partner:PAR"]``; ``"."`` or blank -> ``[]``."""
    if text in (None, "", MISSING):
        return []
    return [f for f in text.split(";") if f and f != MISSING]


def parse_window(text):
    """``"chr2:182666976-182667119"`` -> ``("chr2", 182666976, 182667119)``, the numbers as written."""
    try:
        chrom, span = text.rsplit(":", 1)
        start, end = span.split("-")
        return chrom, int(start), int(end)
    except (AttributeError, ValueError):
        raise ContractError("not a window (chrom:start-end): %r" % (text,))


def region_str(chrom, start, end):
    return "%s:%d-%d" % (chrom, start, end)


def parse_junctions(text):
    """The ``junctions`` column: ``partner|strands|type|reads`` entries joined by ``;``; ``.`` = none.

    The partner is the other end's BED region (0-based, half-open); the strands are the BEDPE line's
    ``strand1 strand2`` in the line's order (``+-`` DEL, ``-+`` DUP, ``++`` / ``--`` INV, ``..`` unknown).
    """
    out = []
    for entry in split_flags(text):
        parts = entry.split("|")
        if len(parts) != 4:
            raise ContractError("not a junction (partner|strands|type|reads): %r" % (entry,))
        chrom, start, end = parse_window(parts[0])
        reads = to_int(parts[3])
        if reads is None:
            raise ContractError("junction reads is not a count: %r" % (entry,))
        out.append(Junction(chrom, start, end, parts[1], parts[2], reads))
    return out


def format_junctions(junctions):
    return ";".join("%s|%s|%s|%d" % (region_str(j.chrom, j.start, j.end), j.strands, j.sv_type, j.reads)
                    for j in junctions) or MISSING


def read_tsv(path):
    """``(header, rows)``; rows as ordered dicts of strings, short rows padded with ``""``."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header and header[0].startswith("#"):
            header[0] = header[0][1:]
        rows = []
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(header) - len(f))
            rows.append(collections.OrderedDict(zip(header, f)))
    return header, rows


def write_tsv(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(h, "")) for h in header) + "\n")
    os.replace(tmp, path)


# ---- the candidate table ------------------------------------------------------------------------------------

class Candidate(object):
    """One row of the candidate table, typed; ``row`` keeps every column as written."""

    def __init__(self, row, source=None, line=None):
        self.row = row
        where = "%s line %s" % (source or "candidate table", line if line is not None else "?")
        for col in CANDIDATE_REQUIRED:
            if col not in row:
                raise ContractError("%s: no column %r" % (where, col))
        self.chrom = row["chrom"]
        self.start, self.end = to_int(row["start"]), to_int(row["end"])
        if self.start is None or self.end is None or not 0 <= self.start < self.end:
            raise ContractError("%s: start/end %r/%r are not a BED interval" % (where, row["start"], row["end"]))
        self.trio = row["trio"]
        if not self.trio:
            raise ContractError("%s: empty trio" % where)
        self.cls = row.get("class", MISSING) or MISSING
        self.sv_type = row.get("sv_type", MISSING) or MISSING
        self.evidence = collections.OrderedDict(
            (c, to_int(row[c])) for c in EVIDENCE_COLUMNS if c in row)
        try:
            self.junctions = parse_junctions(row.get("junctions", MISSING))
        except ContractError as e:
            raise ContractError("%s: %s" % (where, e))
        self.locus_flags = split_flags(row.get("locus_flags"))
        self.partner_flags = split_flags(row.get("partner_flags"))
        self.flags = split_flags(row.get("flags"))
        self.review_tier = to_int(row.get("review_tier"))
        self.rank = to_int(row.get("rank"))
        self.both_ends = to_int(row.get("both_ends"))
        self.nhf_reads = to_int(row.get("nhf_reads"))
        self.nhf = to_float(row.get("nhf"))

    @property
    def name(self):
        """The bundle directory's name."""
        return "%s_%d_%d" % (self.chrom, self.start, self.end)

    @property
    def id(self):
        return "%s/%s" % (self.trio, self.name)

    @property
    def locus(self):
        return region_str(self.chrom, self.start, self.end)

    @property
    def key(self):
        return (self.trio, self.chrom, self.start, self.end)

    def partners(self):
        """The junction partners' regions, each once, in junction order."""
        seen, out = set(), []
        for j in self.junctions:
            k = (j.chrom, j.start, j.end)
            if k not in seen:
                seen.add(k)
                out.append(k)
        return out

    def __repr__(self):
        return "Candidate(%s %s %s tier %s rank %s)" % (self.trio, self.locus, self.sv_type, self.review_tier,
                                                        self.rank)


def read_candidates(path):
    """``(header, [Candidate])`` from a candidate table (any layout with the four required columns)."""
    header, rows = read_tsv(path)
    for col in CANDIDATE_REQUIRED:
        if col not in header:
            raise ContractError("%s: no column %r (the candidate table has chrom, start, end and trio)"
                                % (path, col))
    return header, [Candidate(row, path, n) for n, row in enumerate(rows, 2)]


def candidate_windows(cand, pad, lengths=None):
    """The candidate's own region and its partners', each padded by *pad* and merged where they touch:
    ``[(chrom, start, end)]`` 0-based half-open, in the order of *lengths* (a ``{chrom: length}`` in
    reference order) or of first appearance. The same rule as bundle_candidates.py's, so a bundle cut
    there and a window computed here agree."""
    regions = [(cand.chrom, cand.start, cand.end)] + cand.partners()
    order = list(lengths.keys()) if lengths else []
    for chrom, _s, _e in regions:
        if chrom not in order:
            order.append(chrom)
    rank = dict((c, i) for i, c in enumerate(order))
    padded = []
    for chrom, start, end in regions:
        limit = lengths.get(chrom) if lengths else None
        s, e = max(0, start - pad), end + pad
        if limit is not None:
            e = min(limit, e)
        padded.append((chrom, s, e))
    padded.sort(key=lambda w: (rank[w[0]], w[1], w[2]))
    merged = []
    for w in padded:
        if merged and merged[-1][0] == w[0] and w[1] <= merged[-1][2]:
            merged[-1] = (w[0], merged[-1][1], max(merged[-1][2], w[2]))
        else:
            merged.append(w)
    return merged


# ---- a bundle -----------------------------------------------------------------------------------------------

def read_bundle_index(path):
    """bundle.tsv as ``(header, rows)``, its columns checked."""
    header, rows = read_tsv(path)
    for col in BUNDLE_INDEX_COLUMNS:
        if col not in header:
            raise ContractError("%s: no column %r (bundle_candidates.py's bundle.tsv has it)" % (path, col))
    return header, rows


def read_candidate_json(path):
    """candidate.json; its windows also as ``windows0`` (0-based half-open tuples)."""
    with open(path) as fh:
        info = json.load(fh, object_pairs_hook=collections.OrderedDict)
    windows0 = []
    for w in info.get("windows", []):
        chrom, start, end = parse_window(w)
        windows0.append((chrom, max(0, start - 1), end))
    info["windows0"] = windows0
    return info


def bundle_dir(bundles, cand):
    """The candidate's bundle directory under *bundles* (a symlink resolved: a shared bundle), or None."""
    d = os.path.join(bundles, cand.trio, cand.name)
    if os.path.islink(d):
        target = os.path.join(os.path.dirname(d), os.readlink(d))
        return os.path.normpath(target) if os.path.isdir(target) else None
    return d if os.path.isdir(d) else None


def bundle_members(d):
    """``{role: path}`` of the member alignments in a bundle directory: ``<role>.bam`` when there, else
    ``<role>.cram``; a member with neither is left out."""
    out = collections.OrderedDict()
    for role in ROLES:
        for ext in (".bam", ".cram"):
            p = os.path.join(d, role + ext)
            if os.path.isfile(p):
                out[role] = p
                break
    return out


def alignment_index(path):
    """The index beside an alignment file (``.bai`` / ``.bam.bai`` / ``.crai``), or None."""
    if path.endswith(".bam"):
        names = (path + ".bai", path[:-4] + ".bai")
    elif path.endswith(".cram"):
        names = (path + ".crai",)
    else:
        names = ()
    for p in names:
        if os.path.isfile(p):
            return p
    return None


# ---- the run's metrics --------------------------------------------------------------------------------------

def read_metrics(path):
    with open(path) as fh:
        m = json.load(fh, object_pairs_hook=collections.OrderedDict)
    if m.get("mode") != "discovery":
        raise ContractError("%s: mode is %r, not discovery" % (path, m.get("mode")))
    return m


def metrics_region(metrics, chrom, start, end):
    """The BED row of a region as the metrics carry it (its coordinates exactly), or None."""
    for r in metrics.get("regions") or []:
        if r.get("chrom") == chrom and to_int(r.get("start")) == start and to_int(r.get("end")) == end:
            return r
    return None


def nhf_block(metrics):
    """The screen's status and sampled fractions: ``{status, reason, db, sampled}`` (status ``none`` when
    the run predates the screen)."""
    b = metrics.get("nonhuman_fraction") or {}
    out = collections.OrderedDict()
    out["status"] = b.get("status", "none")
    out["reason"] = b.get("reason", "") or ""
    out["db"] = b.get("db")
    out["sampled"] = (b.get("sampled") or {}).get("fractions")
    out["regions_screened"] = b.get("regions_screened")
    return out


def nhf_region(metrics, chrom, start, end):
    """The screen's entry for a region: ``{reads, placed_unmapped, fractions, top_taxon}``, or None when
    the region was not screened (which is not a measurement of 0)."""
    b = metrics.get("nonhuman_fraction") or {}
    for r in b.get("regions") or []:
        if r.get("chrom") == chrom and to_int(r.get("start")) == start and to_int(r.get("end")) == end:
            out = collections.OrderedDict()
            out["reads"] = to_int(r.get("reads"), 0)
            out["placed_unmapped"] = to_int(r.get("placed_unmapped"), 0)
            out["fractions"] = r.get("fractions")
            out["top_taxon"] = r.get("top_taxon")
            return out
    return None


# ---- the verdict TSV ----------------------------------------------------------------------------------------

def verdict_key(row):
    return (row.get("trio", ""), row.get("chrom", ""), to_int(row.get("start")), to_int(row.get("end")))


def read_verdicts(path):
    """``{(trio, chrom, start, end): row}``; the file is append-only and the last row for a candidate wins.
    A row with a verdict outside VERDICTS is a ContractError; an empty verdict clears an earlier one."""
    header, rows = read_tsv(path)
    for col in ("trio", "chrom", "start", "end", "verdict"):
        if col not in header:
            raise ContractError("%s: no column %r (svreview's verdict TSV has %s)"
                                % (path, col, ", ".join(VERDICT_COLUMNS)))
    out = collections.OrderedDict()
    for n, row in enumerate(rows, 2):
        v = row.get("verdict", "")
        if v and v not in VERDICTS:
            raise ContractError("%s line %d: verdict %r is not one of %s" % (path, n, v, ", ".join(VERDICTS)))
        key = verdict_key(row)
        if key[2] is None or key[3] is None:
            raise ContractError("%s line %d: start/end are not integers" % (path, n))
        if v:
            out[key] = row
        else:
            out.pop(key, None)
    return out


def write_verdicts(path, rows, append=False):
    """Write (or append) verdict rows; a new file gets the header."""
    new = not (append and os.path.isfile(path) and os.path.getsize(path) > 0)
    with open(path, "a" if append else "w") as fh:
        if new:
            fh.write("\t".join(VERDICT_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in VERDICT_COLUMNS) + "\n")


def merge_verdicts(header, rows, verdicts):
    """The candidate table with ``verdict verdict_note verdict_reviewer verdict_at`` appended (``.`` where
    no verdict): the columns the analysis reads back as a flag."""
    out_header = list(header)
    for c in MERGED_VERDICT_COLUMNS:
        if c in out_header:
            raise ContractError("the table already carries %r: merge into the unmerged table" % c)
        out_header.append(c)
    out_rows = []
    n = 0
    for row in rows:
        r = collections.OrderedDict(row)
        v = verdicts.get(verdict_key(row))
        r["verdict"] = v["verdict"] if v else MISSING
        r["verdict_note"] = (v.get("note") or MISSING) if v else MISSING
        r["verdict_reviewer"] = (v.get("reviewer") or MISSING) if v else MISSING
        r["verdict_at"] = (v.get("reviewed_at") or MISSING) if v else MISSING
        n += 1 if v else 0
        out_rows.append(r)
    return out_header, out_rows, n
