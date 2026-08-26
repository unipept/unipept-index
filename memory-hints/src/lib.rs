//! Best-effort hints to the memory subsystem.
//!
//! The index is far larger than any cache and larger than the TLB can map at 4 KB granularity — a
//! full UniProt suffix array runs to ~149 GB out of 223 GB total; see the crate docs of `sa-index`
//! for the breakdown — so both backends spend most of their time waiting on memory. The two hints
//! here attack that from opposite ends:
//!
//! * [`prefetch::prefetch_read`] starts a load early, so its ~80-100 ns DRAM latency overlaps
//!   with useful work. See [`prefetch`].
//! * [`hugepages::advise_capacity`] asks for 2 MB pages, so walking a multi-gigabyte buffer costs
//!   far fewer TLB misses and page-walks. See [`hugepages`].
//!
//! Neither affects correctness. Both are advisory — the CPU and the kernel are free to ignore them
//! — and both are silent no-ops off their target platform, so callers need no `cfg` guards.
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
//!
//! So the thing to read before using `prefetch_read` is its `inline(always)` note; the thing to
//! read before using `advise_capacity` is the ordering argument in [`hugepages`].
//!
//! # Why a separate crate
//!
//! `text-compression`, `sa-mappings`, `sa-index` and `bitarray` all need these hints without
//! depending on one another. Both modules were previously separate crates for exactly that reason
//! — `prefetch`, and `hugepages` inside `bitarray` — and were merged once it was clear the second
//! had nothing to do with bit packing: four of its six call sites advise plain `Vec`s in
//! `sa-index` with no bit array in sight.
#![warn(missing_docs)]

pub mod hugepages;
pub mod prefetch;
