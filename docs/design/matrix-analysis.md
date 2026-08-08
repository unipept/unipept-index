> **Historical record.** Captured during the performance work that landed on
> `feature/preloaded-sa-improvements`, and preserved verbatim as the rationale behind
> decisions now encoded in the code. Measurements are pinned to the `cleanup-baseline`
> tag and to the machine and index named below; they are not re-run by CI and will drift.
> Where a conclusion is load-bearing it is also stated as a comment at the code it explains.

# Full-DB matrix: preloaded vs mmap

Full-DB (`uniprot-2025-04`) in-process matrix, median qps over 20 runs, `max_matches = 10000`.
Grid: `equate_il {T,F}` × `tryptic {T,F}` × `batch {scalar,8,16,32,64}` × `kmer {none,5,6}` ×
`{small,medium,large}` peptide files. Two index loads total (one per backend).

## Conclusions

### 1. preloaded always beats mmap (0.59–0.95×)
mmap never wins. The gap tracks what is bottlenecking:

- **Widest (~0.60×)** on *fast, search-bound* cases — large peptides + k-mer + batching + equate_il.
  That is where mmap pays the page-walk cost on the binary-search probes.
- **Narrowest (~0.90–0.95×)** on *slow, retrieval-bound* cases — small + tryptic. When
  match-volume retrieval dominates, both backends do the same random work and mmap's relative
  penalty shrinks.

Practical: on the cases that actually hurt (short peptides), mmap is only ~10% behind. The
RAM/stability tradeoff for mmap is cheapest exactly where throughput matters most.

### 2. Peptide length dominates absolute throughput
`large ~150–340k` › `medium ~20–190k` › `small ~0.6–8.6k qps`. Short peptides are 30–40× slower
(match volume). Absolute floor: small + tryptic ≈ 660 qps on both backends.

### 3. `tryptic=true` is the single biggest cost knob
6–7× slower on small/medium (`124k→23k`, `4.9k→0.7k`), ~25% on large. It validates every candidate
across huge match ranges. Unfixable by any other knob.

### 4. `equate_il=true` is *faster* (+15–22% on non-tryptic large/medium)
Takes the no-validation fast path (`equate_il && !tryptic && skip==0` → iterate+collect). Also the
common case, so the default is both correct and fast.

### 5. MLP batching: real on long peptides, nothing on short
- large **+25–30%** (`217k→273k`), medium non-tryptic **+15%**, small / tryptic **~0%**.
- **Knee at B=16** — 32 and 64 add essentially nothing.
- Safe default **B≈16**: a win for long-peptide datasets, harmless elsewhere.

### 6. k-mer table: modest, and redundant with batching on mmap
- **preloaded:** 6-mer helps search-bound cases **+15–21%** (5-mer less). Batching and 6-mer roughly
  **stack**: large F/F `186k` → +batch `223k` → +6-mer `262k`.
- **mmap:** 6-mer helps *scalar* (`149k→171k`, +14%), but once batching is on it adds **~+2%**
  (`190k→195k`), and the **5-mer actively hurts** batched mmap (−11%). Batching and the k-mer table
  hit the same page-walk-bound binary search, so on mmap they are redundant — either helps, both ≈ one.
- At full-DB scale the binary search is a small slice of total time, so k-mer's value is nowhere near
  the swissprot "+96%" (cache-resident, search-dominated).

**→ For the mmap production backend with batching on, the 6-mer table's ~3 GB buys almost nothing and
the 5-mer can hurt. Worth A/B-ing whether to carry it.**

### Best / worst
- Best: large + equate_il + non-tryptic + batch + 6-mer → preloaded **338k** / mmap **234k**.
- Worst: small + tryptic ≈ **660** on both.

## Production recommendation (mmap)
- Enable batching (**B≈16**) — win on long peptides, neutral elsewhere.
- Reconsider the k-mer table for mmap — near-redundant with batching; the 3 GB may be better dropped.
- Keep `equate_il=true` default.
- Short-peptide + tryptic is the real floor; only reducing match volume (lower `max_matches`, dedup)
  moves it.

## Measurement stability — noticed noise
The 5-mer/batched mmap cells swing more than 20-run medians should (e.g. large F/F batch16 mmap 5-mer
`169,951` sits *below* none `190,491`, but at scalar it is *above*). That is measurement noise on the
mmap page-walk path, not a real inversion.

Two fixes are now in the harness (rerun to get them):

- **Per-config spread (#1).** Matrix mode aggregates the `runs` reps into one record with
  `stats.qps_{min,p10,p50,p90,max}`; `matrix_bench.sh` shows the median plus a `noise` column
  (`±(p90−p10)/2 / p50`). Treat sub-band deltas as ties — the 5-mer wobble sits inside its band; the
  real stories (batching +30%, tryptic 6×) clear it easily.
- **Reproducible query stream (#5).** `--seed <u64>` fixes random-mode peptide generation across runs
  and invocations. (The matrix already fixes its stream by reusing the same file lines for every
  config, so cross-cell input variance there was already zero.)

Still worth doing on the server: pin the bench to a core + `performance` governor + turbo off,
interleave configs instead of running each backend's reps in a block, and pre-fault pages
(`MAP_POPULATE`/`madvise(WILLNEED)`) so mmap cells don't differ on cold-page luck.
