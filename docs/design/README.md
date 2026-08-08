# Design notes

Historical records from the performance work on the suffix-array index. They exist to answer
"why is the code shaped like this?" — the measurements that justified each optimization, and the
ones that killed the optimizations you will not find in the code.

Every load-bearing conclusion is *also* stated as a comment next to the code it explains; these
documents are the long form, not the source of truth. They are pinned to the `cleanup-baseline`
tag and to the machines and index releases named inside them, and are not re-run by CI.

| Document | What it answers |
|---|---|
| [mlp-batching-design.md](mlp-batching-design.md) | The ranked list of optimization levers, and why cross-query MLP batching was the one worth building. Start here. |
| [matrix-analysis.md](matrix-analysis.md) | The full preloaded-vs-mmap configuration matrix: which knobs moved throughput and which turned out to be noise. |
| [refactor-comparison.md](refactor-comparison.md) | Throughput across the refactor sequence, used to confirm the mmap/preloaded split cost nothing. |

## Things these documents settled

- **The k-mer bounds table pays off**, but far less than the "~60 %" figure quoted early on; the
  effect is concentrated in short peptides. See `matrix-analysis.md`.
- **Cross-query batched protein retrieval and scalar k-mer prefetch are dead** (median +1.7 % and
  +0.3 %, both inside the noise floor). They were removed rather than left as knobs — see the
  `SearchTuning` doc comment in `sa-index/src/sa_searcher/mod.rs`.
- **`MADV_WILLNEED` on k-mer SA ranges regressed the mmap backend by 16.8 %** and was reverted; the
  comment at `sa-index/src/array/mmap.rs` records why the hook is deliberately left unimplemented.
- **`validate_batch` is a cliff, not a peak**: 16 → 32 buys ~10 %, then it plateaus. Do not lower it.
