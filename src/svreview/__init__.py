"""svreview - review the structural-variant candidates of kmer_denovo_filter's discovery mode.

The Python side turns what the analysis produced (a candidate table, per-candidate read bundles, a run's
metrics) into what the viewer reads (one manifest, mini-BAMs), serves it, and reads the verdicts back.
The viewer itself is the static site under docs/.
"""

__version__ = "0.1.0"

#: The manifest layout the viewer understands; bumped when a field changes meaning.
CONTRACT_VERSION = 1
