"""The `baseline_startup` suite reads exactly like `startup`, because it is that suite at
`kmer_k = 0` with the fbc9328 arm alongside.

UNTRACKED, alongside `suites/baseline_startup.toml`. A separate suite name rather than a variant of
`startup`: a session directory is keyed by suite name, and these rows are not interchangeable with
`startup`'s (that one runs at the shipped 5-mer table, this one holds it out so fbc9328 can appear).
"""

from __future__ import annotations

from .startup import analyse  # noqa: F401
