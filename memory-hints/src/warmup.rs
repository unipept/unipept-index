//! Faulting a memory mapping into the page cache, ahead of the reads that need it.
//!
//! The third hint in this crate, and the one with the coarsest grain: [`prefetch`](crate::prefetch)
//! starts one load, [`hugepages`](crate::hugepages) changes how a buffer is paged, and this walks a
//! whole mapped section so the reads that follow find it resident.
//!
//! Unlike the other two this is **not** advisory — it does the faulting itself, and it is not free:
//! it is warmup, run once before timing or serving, never per query.

use memmap2::Mmap;

/// Page size assumed when warming a mapping. Touching one byte per this many bytes is enough to
/// fault in every page; a larger real page size just means some touches are redundant.
const ASSUMED_PAGE_SIZE: usize = 4096;

/// Reads every page of `mmap[range]` into the page cache.
///
/// Warmup only — never called per query. Serving from a cold mapping means the first requests pay
/// the page faults; `sa-benchmarks` sweeps every structure before it times anything, which is what
/// makes a throughput figure a steady-state one. `sa-server` does **not** call this today, so a
/// freshly started server warms itself on live traffic.
///
/// All three steps matter:
///
/// 1. [`memmap2::Advice::Sequential`] tells the kernel to read far ahead, so the sweep faults in long runs
///    instead of one page at a time.
/// 2. Touching one byte per page forces the fault. The read **must** be laundered through
///    [`std::hint::black_box`]: without it the optimizer deletes a loop whose result is unused,
///    and the warmup silently does nothing.
/// 3. [`memmap2::Advice::Random`] restores the steady-state pattern. The index is probed in an order the
///    kernel cannot predict, so leaving readahead enabled would make every later miss drag in
///    neighbouring pages that will not be used.
///
/// `range` is a byte range into `mmap`; callers pass only their own section, so a structure
/// sharing a file with others does not warm its neighbours.
///
/// Returns the number of bytes swept. Every caller passes it back up so the benchmark harness can
/// divide it by the elapsed time: a sweep running at disk bandwidth and one running at memcpy
/// bandwidth do the same work and take an order of magnitude apart, and without the byte count the
/// two are indistinguishable in a report.
///
/// It previously existed as five near-identical copies, then as one copy in the protein-text
/// crate — which every caller happened to depend on, but which has nothing to do with warming a
/// mapping. It sits here with the other memory hints instead.
pub fn touch_all_pages(mmap: &Mmap, range: std::ops::Range<usize>) -> u64 {
    #[cfg(unix)]
    let _ = mmap.advise(memmap2::Advice::Sequential);

    let swept = mmap[range.clone()].len() as u64;
    for chunk in mmap[range].chunks(ASSUMED_PAGE_SIZE) {
        std::hint::black_box(chunk[0]);
    }

    #[cfg(unix)]
    let _ = mmap.advise(memmap2::Advice::Random);

    swept
}
