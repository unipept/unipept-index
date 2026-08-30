//! Best-effort hints to the memory subsystem.
//!
//! The index is far larger than any cache and larger than the TLB can map at 4 KB granularity — a
//! full UniProt suffix array runs to ~160 GB out of 223 GB total; see the crate docs of `sa-index`
//! for the breakdown — so both backends spend most of their time waiting on memory. The two hints
//! here attack that from opposite ends:
//!
//! * [`prefetch::prefetch_read`] starts a load early, so its ~80-100 ns DRAM latency overlaps
//!   with useful work. See [`prefetch`].
//! * [`hugepages::advise_capacity`] asks for 2 MB pages, so walking a multi-gigabyte buffer costs
//!   far fewer TLB misses and page-walks. See [`hugepages`].
//! * [`warmup::touch_all_pages`] faults a mapped section in up front, so the reads that follow
//!   find it resident. See [`warmup`].
//!
//! None of them affects correctness. The first two are advisory — the CPU and the kernel are free
//! to ignore them — and are silent no-ops off their target platform, so callers need no `cfg`
//! guards. [`warmup`] is the exception: it does the faulting itself and costs real time, which is
//! why it runs once at startup rather than per query.
//!
//! # Two hints, opposite disciplines
//!
//! They are grouped here because they share a rationale, not a usage pattern, and getting either
//! one wrong is silent. The mistakes are not the same mistake:
//!
//! * [`prefetch::prefetch_read`] is issued **per access, on the hot path** — thousands of times
//!   per query, from the innermost loops of binary search and protein retrieval. What matters is
//!   that it costs nothing, which is why it is `#[inline(always)]`: a `call`/`ret` pair around one
//!   hint instruction costs more than the hint saves.
//! * [`hugepages::advise_capacity`] is issued **once per allocation, at load time**. What matters
//!   is *when* — between reserving the allocation and first writing to it. Issued one line too
//!   late it is not merely weaker, it is worthless, and it still looks like it is doing something.
//! * [`warmup::touch_all_pages`] is issued **once per mapped section, before serving or timing**.
//!   What matters is that its read is laundered through `black_box`: without that the optimizer
//!   deletes the loop and the warmup silently does nothing.
//!
//! So the thing to read before using `prefetch_read` is its `inline(always)` note; the thing to
//! read before using `advise_capacity` is the ordering argument in [`hugepages`].
//!
//! # Why a separate crate
//!
//! `protein-text`, `protein-metadata`, `sa-index` and `bitarray` all need these hints without
//! depending on one another. Both modules were previously separate crates for exactly that reason
//! — `prefetch`, and `hugepages` inside `bitarray` — and were merged once it was clear the second
//! had nothing to do with bit packing: four of its six call sites advise plain `Vec`s in
//! `sa-index` with no bit array in sight.
#![warn(missing_docs)]

pub mod hugepages;
pub mod prefetch;
pub mod warmup;
