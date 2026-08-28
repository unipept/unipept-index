"""The `baseline` suite reads exactly like `defaults`, because it is that grid at `kmer_k = 0`.

A separate name rather than a variant of `defaults`: the two are not interchangeable (one holds the
k-mer table at 5, the other has none), and a session directory is keyed by suite name, so sharing
the name would let two incomparable runs land in one directory.
"""

from __future__ import annotations

from .defaults import analyse  # noqa: F401
