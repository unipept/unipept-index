# sa-index benchmark record

Measurements behind the design decisions in `sa-index`. The crate docs carry the conclusions and
the ratios, which hold generally; this file carries the runs they came from, which do not.

**Read the absolute figures as one machine's, not as targets.** Throughput, resident sizes and
memory ceilings below are specific to the hardware and index named in each session. What transfers
is the shape: which arm wins, by roughly how much, and where the crossovers are.

Nothing here is reproducible by rebuilding the crate. `sa-index` declares no features and carries
no instrumentation; figures attributed to instrumentation come from builds that no longer exist.
Whole-run measurement lives in `sa-benchmarks`, which is excluded from the workspace's
`default-members`.

## Session: memory-ceiling sweep

Index and host: the 223 GB UniProt index on a 295 GB / 12-core Xeon Silver 4410Y. Ceilings imposed
with cgroup v2 `MemoryMax`, swap off, page cache dropped before every cell, 100 reps x 10,000
mixed-length peptides per cell, a 6-mer table attached throughout. Every cell runs
`tryptic=false`, which is the regime that flatters preloading — see the caveat below.

For scale, the index divides as SA 160.2 GB, text 43.3 GB, mapping 10.8 GB and protein metadata
8.3 GB, the last becoming ~24 GB once preloaded (the RSS delta when `preloaded-proteins` moves it
to the heap). The ratio is the durable part: the SA is 72% of the index, so it dominates
residency, and nothing can preload it — there is no `preloaded-sa`.

## When the index does not fit in RAM

**The degradation is a concurrency limit, not a bandwidth one.** A major fault blocks its
thread, and `prefetch_read` cannot help: a CPU hint instruction cannot fault, so every prefetch
in this crate is inert against an absent page. With rayon at the core count, each faulting
thread idles a core. Raising `RAYON_NUM_THREADS` leaves the fault *count* unchanged to within
0.4% and still buys, on the `mmap` arm:

| ceiling | major faults/rep | default threads | tuned | change |
|---|---|---|---|---|
| none | 0 | **38,081** | 34,264 @ 48 | **-10.0%** |
| 167 GB (75%) | 23,828 | 15,518 | **25,734** @ 96 | **+65.8%** |
| 112 GB (50%) | 45,572 | 10,503 | **20,473** @ 96 | **+94.9%** |

Unconstrained the default is already the peak, and every tuned value is a loss; under a ceiling
96 threads is the best of the three everywhere, and the peak is at or below it (an earlier
sweep found the curve monotonically downward from 96 through 128 and 192). Every arm has the
same shape. At 96 threads under a 167 GB ceiling the gains are +65.8% (`mmap`), +86.2%
(`preloaded-proteins`), +102.7% (`preloaded-text,preloaded-proteins`) and +94.9% (that plus
`preloaded-mapping`); at 112 GB, +94.9%, +116.4%, +112.8% and +116.4%. The fully preloaded
build did not fit under either ceiling. Unconstrained every arm loses — between -5.4%
(preloaded) and -18.1% (`preloaded-proteins`) at 96 threads — so this is a deployment knob, not
a default. It is also the largest single effect anywhere in this investigation, larger than any
storage-backend choice.

Reproduced across three months and a rewritten harness: at 3259427 (2026-08-10/11) the `mmap`
gains were +65.6% and +86.1% at the same two ceilings, at 2dfa6517b7 (2026-08-16) +62.6% and
+92.2%, each time with the cost unconstrained and major faults flat across thread counts. Every
sign and rough magnitude has held.

**Every `preloaded-*` feature is a bet on full residency.** With the index resident they all
pay: against plain `mmap`, `preloaded-proteins` is +22.1%, `preloaded-text,preloaded-proteins`
is +46.1% and adding `preloaded-mapping` is +52.7%, with the fully preloaded build +57.6%. But
preloaded memory is non-evictable anonymous memory, so under pressure it cannot be reclaimed
and instead displaces the file-backed page cache the mapped structures live in — both bigger
and hotter per query. Sweeping the ceiling down, with each column adding one structure to the
one on its left:

| ceiling | mmap | +proteins | +text | +mapping | preloaded |
|---|---|---|---|---|---|
| none | 35,725 | 43,603 | 52,198 | 54,553 | **56,297** |
| 223 GB | 27,740 | 19,045 | 29,187 | 30,202 | did not fit |
| 167 GB | 14,748 | 13,832 | 15,161 | 15,617 | did not fit |
| 140 GB | 12,724 | 11,683 | 12,497 | 12,658 | did not fit |
| 112 GB | 10,679 | 9,582 | 10,212 | 10,425 | did not fit |
| 78 GB | **7,411** | 292 | 169 | did not fit | did not fit |

The unconstrained `mmap` figure here (35,725) and the one in the thread table above (38,081)
are the same configuration measured by two different suites in the same session. They differ by
6.6%, which is inside both cells' own resolution floors; neither suite can tell them apart, and
nothing should be read into the gap. Compare within a table, never across.

There is no crossover to find above 78 GB: from the first ceiling that binds, no preloading arm
is ahead of plain `mmap` by more than the floor, and `preloaded-proteins` is behind it at every
one. At 78 GB — roughly a third of the index — the preloaded arms do not degrade, they
collapse: 31x and 55x `mmap`'s fault rate (2,310,152 and 4,173,435 major faults per rep against
75,383) for a 96% and 98% loss. The fully preloaded build was OOM-killed at every ceiling in
the sweep. So preloading is worth having exactly when the whole index is guaranteed resident,
and is the wrong default anywhere the ceiling might move.

Note the run's own caveats: the `mmap` and `preloaded-proteins` cells at 223 GB and the `mmap`
cell at 78 GB had not reached steady state (drift +25.2%, +13.0% and -62.5% from the first
quarter of reps to the last), so those three are softer than the rows between them.

**Even when residency is guaranteed, `preloaded-proteins` is not the arm to reach for.** This
sweep is its best case: at `tryptic=false` the search accepts 9-13% of the candidates it
examines, so retrieval is a third of the work and the structure it preloads is hot. Under
`tryptic` acceptance falls to ~0.5%, retrieval drops to 1-7% of the work, and the same
session's 16-cell throughput sweep put `mmap,preloaded-text,preloaded-proteins` ahead of it in
all sixteen, at a lower resident footprint (242 GB against 250 GB — `preloaded-proteins` has
the highest of any arm). Preloading the metadata alone closes the retrieval gap and none of the
search one, which is the larger half.

**A 6-mer k-mer table is worth its 3.06 GB here**, against the resident-case measurement that
cannot separate it from a 5-mer. Established by the ceiling sweep at 3259427 and carried as
background since: at a 167 GB ceiling the 6-mer is +18.4% and -27.9% faults versus no table,
where a 5-mer is +3.2% and -6.2%, i.e. barely distinguishable from nothing. The difference is
working-set size, not probe count — a 5-mer narrows the search to ~7 SA pages per query, a
6-mer to ~1 — and that only matters once pages can be evicted. Every cell of the sweeps above
therefore has a 6-mer attached, since the table removes exactly the phase that degrades under a
cap and a sweep without it would answer a different question. This is why `sa-builder
--kmer-size` documents 6 as the tuning step for constrained deployments; the default stays 5,
which is the size the resident measurements support.

**All of the loss is in the search phase**, measured at 3259427 by instrumentation since
removed: retrieval was flat at ~147 ms per rep across every ceiling while search went
135 ms -> 1127 ms, so the two-pass prefetch pipeline in `sa_searcher::retrieval` keeps working
under paging and needs nothing. Within search, the split is roughly even between the dependent
binary-search chain and the contiguous SA range scan (52% / 48% of thread-time at a 167 GB
ceiling).

Two further ideas were measured and rejected. **`MADV_WILLNEED` over the SA range about to be
scanned does not pay**: the advice lands (major faults -23-25% under a ceiling) but the
throughput decays from +12.0% at the core count to ~0% at 96 threads, since oversubscription
already overlaps those faults, and it costs -3.7% resident. The comment in `array::mmap`
carries the reasoning in full. And **sorting queries by k-mer prefix to create page
locality does not work either**: it changed the fault count by -0.1% and cost 4.4% throughput.
With 10,000 queries per rep drawn against 20^6 possible
6-mers, the expected number of queries sharing a prefix is under one, so there is no page reuse
for sorting to expose. Locality needs reuse, and this workload has none.
