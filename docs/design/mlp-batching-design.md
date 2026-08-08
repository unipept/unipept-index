> **Historical record.** Captured during the performance work that landed on
> `feature/preloaded-sa-improvements`, and preserved verbatim as the rationale behind
> decisions now encoded in the code. Measurements are pinned to the `cleanup-baseline`
> tag and to the machine and index named below; they are not re-run by CI and will drift.
> Where a conclusion is load-bearing it is also stated as a comment at the code it explains.

# Optimization levers — summary (2026-07-31)

Ranked by impact for the **mmap** backend (production: lower RAM pressure, less crash-prone).
Swissprot medians unless noted; "full DB" = `uniprot-2025-04` on the server (300 GB RAM, 12 cores).

| lever | mmap | preloaded | notes |
|---|---|---|---|
| **k-mer bounds cache** (k=5) | **+96%** | +50% | *already implemented*; fewer probes → fewer page walks. **Biggest mmap win.** Enable via `sa-builder --output-kmer-table` + `sa-server --kmer-table-file`. ~5 MB (k=4) / ~127 MB (k=5) RAM. |
| **huge pages** (`SA_MADV_HUGEPAGE`) | ✗ no-op | +15.5% (full DB) | file-backed mmap can't get THP; preloaded's anonymous `Vec`s can. |
| **MLP batching** (`SA_MLP_BATCH`) | +57% swissprot / +4% full DB | +51% swissprot / +6% full DB | full DB is TLB/page-walk-bound, so batching (data-load MLP) barely helps there; shines when DRAM-data-latency bound. |

**Status:** the k-mer table file already exists on the server and production `sa-server`
loads it (`--kmer-table-file`), so production mmap **already realizes the ~2×** — it's the
current baseline, not a new win. The earlier benchmark sweeps ran *without* the table
(mlp_sweep didn't pass one), so those ~31–32k mmap numbers understate production by ~2×.
To benchmark like production, pass `MLP_KMER_FILE=<path>` to `mlp_sweep.sh`.

**Remaining mmap headroom is small:** with k-mer already on, batching adds ~little (short
binary search), and huge pages don't apply to file-backed mmap. Absent a heavy hugetlbfs
change, mmap is near its practical ceiling for these constraints. Bigger absolute
throughput lives on **preloaded** (k-mer + `SA_MADV_HUGEPAGE` + modest `SA_MLP_BATCH`) — at
the RAM/crash-risk cost you're deliberately avoiding in production.

Caveat: k-mer and batching partly overlap — with the k-mer table the binary search is short,
so batching adds little on top (measured +1% on swissprot preloaded with k=5). k-mer first.

---

# Cross-query batching for memory-level parallelism (MLP)

Goal: hide the random-access DRAM latency that caps the **mmap** backend, by keeping
many independent memory misses in flight per core. (Huge pages were the wrong lever —
they don't apply to a file-backed mmap; see `refactor-comparison.md`.)

## Where the latency is (current code)

- `search_all_peptides` → `peptides.par_iter().filter_map(search_peptide)`.
  Rayon gives **thread-level** parallelism (one peptide per core at a time). Within a
  core, a peptide is processed to completion before the next.
- Per peptide: `search_matching_suffixes` loops `skip = 0..sample_rate` →
  `search_bounds(&s[skip..])` → **two** `binary_search_bound_in_range` calls (Min, Max).
- `binary_search_bound_in_range` (sa_searcher.rs:197) is the hot **dependent chain**:
  ```rust
  while hi - lo > 1 {
      let center = (lo + hi) / 2;
      let (retval, lcp) = self.compare(search_string, self.sa.get(center), skip, bound);
      if ... { hi = center } else { lo = center }   // next center depends on this result
  }
  ```
  Each step has two random misses: `sa.get(center)` (SA) and the first `text.get(suffix)`
  inside `compare()`. The next `center` can't be computed until `compare` returns, so a
  core that can sustain ~10 outstanding misses runs at **1**. That underused MLP is the
  entire opportunity.

The intra-query `prefetch_binary_search_pivots` (prefetch both children before reading
`center`) only partially helps a single stream; it can't overlap the `sa→text`
dependency. Real overlap needs *independent* work — other queries.

## The idea — AMAC-style multi-stream interleaving

Process **B** independent peptides per core and advance their binary searches in
lockstep, so ~B misses are in flight at once. Three stages per level:

| stage | for each active stream k | overlaps |
|---|---|---|
| A | `center[k] = (lo[k]+hi[k])/2;` `sa.prefetch_sa_index(center[k])` | issues B SA prefetches |
| B | `suffix[k] = sa.get(center[k]);` `text.prefetch_at(suffix[k]+skip[k])` | B SA reads land; issues B text prefetches |
| C | `compare(...)` (text warm) → update `lo/hi/lcp/found`; retire if `hi-lo<=1` | B text reads land |

Stream k+1's prefetch in each stage overlaps stream k's miss → throughput scales with
B until memory bandwidth (not latency) binds. `compare()`'s inner loop stays scalar:
only its *first* `text.get` is the random miss (prefetched in B); the rest are
sequential and covered by the HW prefetcher.

## Integration points

1. **`binary_search_bounds_batched(bound, &mut [StreamState])`** — new fn mirroring
   `binary_search_bound_in_range` *exactly*, including the `hi == left+1 && lo == left`
   edge case, `lcp_left/lcp_right` tracking, and the `found` flag — but over an array of
   `StreamState { lo, hi, lcp_left, lcp_right, found, active, search_string, skip }`.
2. **Batched `search_bounds`** — for a batch: (a) batch the k-mer `table.lookup`s (also a
   random access — prefetch/pipeline them), (b) run Min for all streams, (c) run Max for
   the streams whose Min was found.
3. **`search_all_peptides`** — swap `par_iter()` for `par_chunks(B)`; each chunk runs the
   batched bound search, then per-peptide retrieval (retrieval already has its own
   two-pass prefetch batching — leave it, or batch across the chunk too). Handle the odd
   remainder with the existing scalar path.
4. **Tunable B** (start 8–16). Keep the scalar path as the B=1 fallback and for the
   remainder.

## Correctness

The batched state machine must return **identical** `(min,max)` bounds as the scalar
version. Guard with:
- the existing `sa_searcher` tests (`test_search_simple/sparse/dense`, `test_il_*`), and
- a **differential test**: for a corpus of peptides, assert
  `batched_search(p) == scalar_search(p)` for every p (bounds and final suffix sets).

## Expected benefit & caveats

- Targets `search_bounds`, a large share of mmap query time (full-DB timing breakdown).
  If per-core latency-bound — likely for random 8-byte reads over a 242 GB mapping —
  B-way interleaving can cut the bound-search phase by a meaningful factor.
- Helps **mmap more than preloaded** (more latency to hide); preloaded may see a smaller
  gain or none if it's already closer to bandwidth-bound.
- **Won't show on swissprot** (cache-resident, ~730k qps) — correctness can be verified
  locally, but the latency win must be measured on the **full-DB server**.
- If the machine is already memory-bandwidth-saturated across all cores under rayon, the
  ceiling is lower — worth measuring B∈{1,4,8,16} to find the knee.

## Prototype plan

1. Implement in a worktree behind a tuned constant `B` (or env/CLI toggle for A/B).
2. `cargo test -p sa-index` + the new differential test (correctness gate).
3. Build mmap + preloaded; A/B batched-vs-scalar and sweep B on the full DB via the harness.

## Prototype results (built; correctness-gated)

Implemented as `Searcher::search_matching_suffixes_batched` (3-stage interleaved
binary search + batched bounds), wired into the benchmark behind `SA_MLP_BATCH=B`
(order-preserving `par_chunks`). Patch: `benchmarks/mlp-batching.patch`; sweep harness:
`benchmarks/ab_mlp_batch.sh`.

- **Correctness:** `cargo test -p sa-index` 51/51 pass, incl. `test_batched_matches_scalar`
  (batched == scalar for both sample_rate 1 and 3, across equate_il and max_matches).
  End-to-end on swissprot, scalar vs B=8 gave **identical** `suffix_hit_count`
  (1,779,948) for both backends.
- **Throughput — swissprot (cache-resident, ~560 MB; only L3/DRAM latency to hide):**

  | backend | B=1 | B=4 | B=8 | B=16 | B=32 | B=64 |
  |---|---|---|---|---|---|---|
  | preloaded | 751k | +12% | +28% | +30% | +46% | **+51%** |
  | mmap | 679k | +20% | +23% | +44% | +52% | **+57%** |

  Still climbing at B=64 (knee not reached locally). On the **full-DB server** the
  latency being hidden is page-fault-scale, not L3 — so the win should be at least as
  large; sweep B there with `ab_mlp_batch.sh`.

*Caveats:* quick 10-run smokes, not rigorous ABBA (effect is large + monotonic, clearly
real). The prototype allocates a few small Vecs per skip-step per chunk (active/sub/
bounds) — reusable buffers would trim overhead but weren't needed to show the win.

## Full-DB result (`uniprot-2025-04`, mmap) — batching does NOT carry over

| B | 1 | 4 | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|---|
| mmap qps | 31,610 | +1% | +3% | **+4%** | +3% | +3% | +1% |

**Only +4% (peaks at B=16, then declines)** — vs +57% on swissprot. Correctness is
fine; the full-DB bottleneck is simply not what batching hides.

Leading explanation: swissprot (560 MB) is cache-resident, so random accesses are pure
DRAM-data-latency — batching's sweet spot. Full-DB random access over 242 GB is bound
by **4 KB-page TLB misses / page walks** (each access misses the TLB → multi-level walk;
only ~2–4 hardware page-walkers → per-core MLP capped there, independent of how many
*data* prefetches we issue). The B=1 qps is ~21× slower per query than swissprot,
consistent with every access paying a TLB-miss + DRAM hit. The decline past B=16 is the
prototype's per-step allocation overhead exceeding the tiny benefit. (If RAM < ~250 GB
the index also isn't fully resident → a disk-I/O component `_mm_prefetch` can't touch.)

This ties back to the huge-page null result: file-backed mmap couldn't get huge pages,
so its TLB pressure stayed maximal — and that TLB wall is now what throttles batching.

**Diagnostic confirmed.** Preloaded full-DB sweep (300+ GB RAM, 12 cores) also caps at
**~+6% (B=8), no scaling** — same as mmap. So: fully RAM-resident (not disk I/O), not
mmap-specific, and ~140 MB/s ⇒ not bandwidth-bound. Both backends are **TLB/page-walk
bound**. THP is in `madvise` mode, so the preloaded anonymous `Vec`s sit on 4 KB pages.

**The lever: huge pages on the preloaded `Vec`s** (anonymous ⇒ THP-eligible, unlike the
file-backed mmap where `MADV_HUGEPAGE` was a no-op). Two tests:
- zero-code: `echo always > /sys/kernel/mm/transparent_hugepage/enabled`, re-sweep;
- shippable: env-gated `madvise(MADV_HUGEPAGE)` on the preloaded SA + text `Vec`s.

If page walks are the wall, huge pages should lift the scalar B=1 baseline *and* let
batching resume scaling (data-load MLP is no longer throttled by walker scarcity).
Batching stays valuable wherever the workload is DRAM-data-latency bound (smaller DBs,
or once TLB is addressed) — keep it (env-gated, default off).
