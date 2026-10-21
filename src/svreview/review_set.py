"""review_set.py - a review set, svreview's input (docs/FORMAT.md): candidates.tsv, trios.tsv, the bundles or the
per-trio alignments, the runs' files; opened, resolved and checked here, read by the site builder.

A set is a directory. An "ad hoc" set (``ReviewSet.adhoc``) is the same thing assembled from command-line
arguments - a candidate table, a bundle directory, metrics and informative BAMs per trio - for a run without
a trios.tsv.
"""

import collections
import json
import os

import pysam

from . import contract

FORMAT_VERSION = 1
TRIOS_REQUIRED = ("trio", "child", "father", "mother")
TRIOS_OPTIONAL = ("child_sex", "child_alignments", "father_alignments", "mother_alignments", "metrics", "informative",
                  "reference")
ROLES = contract.ROLES
DEFAULT_PAD = 10000

Finding = collections.namedtuple("Finding", "level where message")


class Trio(object):
    def __init__(self, row, base, source="trios.tsv", line=None):
        where = "%s line %s" % (source, line if line is not None else "?")
        self.row = row
        self.trio = row.get("trio", "")
        if not self.trio:
            raise contract.ContractError("%s: empty trio" % where)
        self.samples = collections.OrderedDict((r, row.get(r, "") or r) for r in ROLES)
        self.child_sex = row.get("child_sex", "") or "0"
        if self.child_sex not in ("0", "1", "2", "", "."):
            raise contract.ContractError("%s: child_sex %r is not 0, 1 or 2" % (where, self.child_sex))

        def path(col):
            v = row.get(col, "")
            return os.path.normpath(os.path.join(base, v)) if v and v != contract.MISSING else None

        self.alignments = collections.OrderedDict((r, path(r + "_alignments")) for r in ROLES)
        self.metrics = path("metrics")
        self.informative = path("informative")
        self.reference = path("reference")

    @property
    def has_alignments(self):
        return all(self.alignments[r] for r in ROLES)


def read_trios(path, base=None):
    """``{trio: Trio}`` from trios.tsv, paths resolved against *base* (default: the file's directory)."""
    base = base or os.path.dirname(os.path.abspath(path))
    header, rows = contract.read_tsv(path)
    for col in TRIOS_REQUIRED:
        if col not in header:
            raise contract.ContractError("%s: no column %r (trios.tsv has %s)" % (path, col, ", ".join(TRIOS_REQUIRED)))
    out = collections.OrderedDict()
    for n, row in enumerate(rows, 2):
        t = Trio(row, base, path, n)
        if t.trio in out:
            raise contract.ContractError("%s line %d: trio %s twice" % (path, n, t.trio))
        out[t.trio] = t
    return out


def write_trios(path, trios):
    """*trios*: ``[{column: value}]`` with at least the required columns; paths as given (relative)."""
    cols = list(TRIOS_REQUIRED) + [c for c in TRIOS_OPTIONAL if any(c in t and t[c] not in ("", None) for t in trios)]
    contract.write_tsv(path, cols, [dict((c, t.get(c, "")) for c in cols) for t in trios])


class ReviewSet(object):
    """A set opened from a directory, or assembled ad hoc."""

    def __init__(self, path=None):
        self.path = os.path.abspath(path) if path else None
        self.config = {}
        self.candidates_path = None
        self.header, self.candidates = [], []
        self.trios = collections.OrderedDict()
        self.bundles_dir = None
        self.metrics_map, self.informative_map = {}, {}
        self.pad = DEFAULT_PAD
        self.title = None
        self.reference = None
        self._metrics_cache = {}
        self._lengths_cache = {}

    # ---- construction -----------------------------------------------------------------------------------

    @classmethod
    def open(cls, path):
        rs = cls(path)
        if not os.path.isdir(rs.path):
            raise contract.ContractError("%s is not a directory" % path)
        cfg = os.path.join(rs.path, "svreview.json")
        if os.path.isfile(cfg):
            with open(cfg) as fh:
                try:
                    rs.config = json.load(fh)
                except ValueError as e:
                    raise contract.ContractError("%s: not JSON: %s" % (cfg, e))
            fmt = rs.config.get("format", FORMAT_VERSION)
            if fmt != FORMAT_VERSION:
                raise contract.ContractError("%s: format %r; this svreview reads format %d" % (cfg, fmt, FORMAT_VERSION))
        rs.pad = int(rs.config.get("pad", DEFAULT_PAD))
        rs.title = rs.config.get("title")
        rs.candidates_path = os.path.join(rs.path, "candidates.tsv")
        if not os.path.isfile(rs.candidates_path):
            raise contract.ContractError("%s: no candidates.tsv" % rs.path)
        rs.header, rs.candidates = contract.read_candidates(rs.candidates_path)
        trios = os.path.join(rs.path, "trios.tsv")
        if os.path.isfile(trios):
            rs.trios = read_trios(trios, rs.path)
        bundles = os.path.join(rs.path, "bundles")
        rs.bundles_dir = bundles if os.path.isdir(bundles) else None
        ref = rs.config.get("reference_fasta")
        rs.reference = os.path.normpath(os.path.join(rs.path, ref)) if ref else None
        return rs

    @classmethod
    def adhoc(cls, candidates_path, bundles_dir=None, metrics=None, informative=None, pad=DEFAULT_PAD, title=None,
              trios_path=None, reference=None):
        rs = cls(None)
        rs.candidates_path = os.path.abspath(candidates_path)
        rs.header, rs.candidates = contract.read_candidates(rs.candidates_path)
        rs.bundles_dir = os.path.abspath(bundles_dir) if bundles_dir else None
        rs.metrics_map = dict(metrics or {})
        rs.informative_map = dict(informative or {})
        rs.pad = pad
        rs.title = title
        rs.reference = reference
        if trios_path:
            rs.trios = read_trios(trios_path, os.path.dirname(os.path.abspath(trios_path)))
        return rs

    # ---- resolution -------------------------------------------------------------------------------------

    def trio_of(self, cand):
        return self.trios.get(cand.trio)

    def samples_of(self, cand):
        t = self.trio_of(cand)
        return t.samples if t else None

    def metrics_path_of(self, cand):
        t = self.trio_of(cand)
        if t and t.metrics:
            return t.metrics
        return self.metrics_map.get(cand.trio) or self.metrics_map.get("*")

    def metrics_of(self, cand):
        path = self.metrics_path_of(cand)
        if not path:
            return None
        if path not in self._metrics_cache:
            self._metrics_cache[path] = contract.read_metrics(path)
        return self._metrics_cache[path]

    def informative_of(self, cand):
        t = self.trio_of(cand)
        if t and t.informative:
            return t.informative
        return self.informative_map.get(cand.trio) or self.informative_map.get("*")

    def reference_of(self, cand):
        t = self.trio_of(cand)
        return (t.reference if t and t.reference else None) or self.reference

    def bundle_of(self, cand):
        return contract.bundle_dir(self.bundles_dir, cand) if self.bundles_dir else None

    def source_of(self, cand):
        """``("bundle", dir)``, ``("trio", {role: path})`` or ``(None, None)``."""
        d = self.bundle_of(cand)
        if d:
            return "bundle", d
        t = self.trio_of(cand)
        if t and t.has_alignments:
            return "trio", t.alignments
        return None, None

    def lengths_of(self, path, ref=None):
        if path not in self._lengths_cache:
            with open_alignments(path, ref) as f:
                self._lengths_cache[path] = collections.OrderedDict(zip(f.references, f.lengths))
        return self._lengths_cache[path]

    def windows_of(self, cand):
        """The candidate's windows (0-based, half-open): candidate.json's when the bundle has one, else computed
        from the row with the set's pad over the alignments' contigs."""
        kind, src = self.source_of(cand)
        if kind == "bundle":
            cj = os.path.join(src, "candidate.json")
            if os.path.isfile(cj):
                info = contract.read_candidate_json(cj)
                if info["windows0"]:
                    return info["windows0"], info
            members = contract.bundle_members(src)
            lengths = self.lengths_of(members[next(iter(members))], self.reference_of(cand)) if members else None
            return contract.candidate_windows(cand, self.pad, lengths), None
        if kind == "trio":
            lengths = self.lengths_of(src["child"], self.reference_of(cand))
            return contract.candidate_windows(cand, self.pad, lengths), None
        return contract.candidate_windows(cand, self.pad), None

    # ---- the check --------------------------------------------------------------------------------------

    def check(self, deep=False):
        """Findings (error / warning / info) about the set; an arm's exporter is done when none is an error."""
        F = []

        def err(where, msg):
            F.append(Finding("error", where, msg))

        def warn(where, msg):
            F.append(Finding("warning", where, msg))

        def info(where, msg):
            F.append(Finding("info", where, msg))

        name = os.path.basename(self.path) if self.path else "the set"
        info(name, "%d candidates in %s" % (len(self.candidates), os.path.basename(self.candidates_path)))
        missing_cols = [c for c in contract.CANDIDATE_COLUMNS if c not in self.header]
        if missing_cols:
            info("candidates.tsv", "optional columns absent: " + ", ".join(missing_cols))
        seen = {}
        for c in self.candidates:
            if c.key in seen:
                warn("candidates.tsv", "%s listed twice" % c.id)
            seen[c.key] = c
        for c in self.candidates:
            for p in c.partners():
                if c.both_ends == 1 and (c.trio,) + p not in seen:
                    warn(c.id, "both_ends is 1 but the partner %s is not a row" % contract.region_str(*p))
        if self.path and not self.trios:
            err("trios.tsv", "missing (one row per trio: trio child father mother [child_sex alignments metrics informative])")
        for t in self.trios.values():
            for role in ROLES:
                p = t.alignments[role]
                if p and not os.path.isfile(p):
                    err("trios.tsv %s" % t.trio, "%s_alignments %s is not a file" % (role, p))
                elif p and not contract.alignment_index(p):
                    err("trios.tsv %s" % t.trio, "%s has no index beside it" % p)
            if any(t.alignments.values()) and not t.has_alignments:
                err("trios.tsv %s" % t.trio, "per-trio alignments need all three members")
            for col, p in (("metrics", t.metrics), ("informative", t.informative), ("reference", t.reference)):
                if p and not os.path.isfile(p):
                    err("trios.tsv %s" % t.trio, "%s %s is not a file" % (col, p))
            if t.metrics and os.path.isfile(t.metrics):
                try:
                    m = contract.read_metrics(t.metrics)
                    block = contract.nhf_block(m)
                    info("trios.tsv %s" % t.trio, "metrics: %d regions, non-human screen %s" % (
                        len(m.get("regions") or []), block["status"]))
                except (contract.ContractError, ValueError) as e:
                    err("trios.tsv %s" % t.trio, "metrics unreadable: %s" % e)
            if t.informative and os.path.isfile(t.informative) and not contract.alignment_index(t.informative):
                err("trios.tsv %s" % t.trio, "informative BAM %s has no index" % t.informative)
        trios_named = set(c.trio for c in self.candidates)
        for trio in sorted(trios_named):
            if self.trios and trio not in self.trios:
                err("candidates.tsv", "trio %s has no row in trios.tsv" % trio)
        n_bundle = n_trio = n_none = n_shared = 0
        for c in self.candidates:
            kind, src = self.source_of(c)
            if kind is None:
                n_none += 1
                err(c.id, "no bundle under bundles/%s/%s and no per-trio alignments" % (c.trio, c.name))
                continue
            if kind == "trio":
                n_trio += 1
            else:
                n_bundle += 1
                if self.bundles_dir and os.path.islink(os.path.join(self.bundles_dir, c.trio, c.name)):
                    n_shared += 1
                members = contract.bundle_members(src)
                for role in ROLES:
                    if role not in members:
                        err(c.id, "bundle has no %s.bam/.cram" % role)
                    elif not contract.alignment_index(members[role]):
                        err(c.id, "%s has no index" % os.path.basename(members[role]))
                if not os.path.isfile(os.path.join(src, "candidate.json")):
                    warn(c.id, "no candidate.json: windows computed from the row with pad %d" % self.pad)
            if deep:
                try:
                    windows, _cj = self.windows_of(c)
                    paths = contract.bundle_members(src) if kind == "bundle" else src
                    for role, p in paths.items():
                        lengths = self.lengths_of(p, self.reference_of(c))
                        for chrom, s, e in windows:
                            if chrom not in lengths:
                                err(c.id, "%s: window contig %s is not in the file" % (role, chrom))
                        if kind == "bundle":
                            from . import bundle as _bundle
                            hits = _bundle.header_paths(p, self.reference_of(c))
                            if hits:
                                warn(c.id, "%s's header carries file paths (%s): strip them before the set is shared"
                                     % (role, ", ".join(sorted(set(h[1][:2] if ":" in h[1] else h[1] for h in hits)))))
                        if role == "child":
                            with open_alignments(p, self.reference_of(c)) as f:
                                n = sum(1 for chrom, s, e in windows if chrom in lengths for _ in f.fetch(chrom, s, e))
                            if n == 0:
                                err(c.id, "the child's file holds no read over the windows")
                except contract.ContractError as e:
                    err(c.id, str(e))
                except (OSError, ValueError) as e:
                    err(c.id, "cannot read the alignments: %s" % e)
        info(name, "sources: %d per-candidate bundles (%d shared ends), %d from per-trio alignments, %d without" % (
            n_bundle, n_shared, n_trio, n_none))
        tiers = collections.Counter(str(c.review_tier) if c.review_tier is not None else "." for c in self.candidates)
        info(name, "tiers: " + ", ".join("%s: %d" % kv for kv in sorted(tiers.items())))
        if self.path:
            v = os.path.join(self.path, "verdicts.tsv")
            if os.path.isfile(v):
                try:
                    verdicts = contract.read_verdicts(v)
                    unknown = [k for k in verdicts if k not in seen]
                    info("verdicts.tsv", "%d verdicts%s" % (len(verdicts), (", %d for candidates not in the table" % len(unknown)) if unknown else ""))
                except contract.ContractError as e:
                    err("verdicts.tsv", str(e))
        return F


def open_alignments(path, ref=None):
    """A pysam AlignmentFile for a BAM, or a CRAM with its embedded reference or *ref*."""
    if path.endswith(".cram"):
        try:
            return pysam.AlignmentFile(path, "rc", reference_filename=ref) if ref else pysam.AlignmentFile(path, "rc")
        except (OSError, ValueError) as e:
            raise contract.ContractError("%s: cannot open the CRAM%s: %s" % (
                path, "" if ref else " (no embedded reference: name one in trios.tsv or give --ref)", e))
    return pysam.AlignmentFile(path, "rb")


def format_findings(findings):
    lines = []
    for f in findings:
        lines.append("%-7s %s: %s" % (f.level.upper(), f.where, f.message))
    n_err = sum(1 for f in findings if f.level == "error")
    n_warn = sum(1 for f in findings if f.level == "warning")
    lines.append("%d error(s), %d warning(s)" % (n_err, n_warn))
    return "\n".join(lines), n_err
