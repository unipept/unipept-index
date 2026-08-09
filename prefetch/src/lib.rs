#![warn(missing_docs)]
//! A single portable software-prefetch hint.
//!
//! The index is far larger than any cache — a full UniProt suffix array runs to hundreds of
//! megabytes even compressed, and the suffix-to-protein mapping to over a gigabyte — so both
//! backends spend most of their time waiting on DRAM. Binary search and protein retrieval walk
//! these structures in an order the hardware prefetcher cannot predict, which is exactly the
//! case software prefetching exists for: the address of the *next* access is known several
//! iterations before the value is needed, so the load can be started early and its ~80-100 ns
//! latency overlapped with useful work.
//!
//! Callers pair this with two-pass batching (fill a batch and issue hints, then process the
//! batch) rather than prefetching one element ahead; see `sa_searcher::iterate_sa_range` and
//! `sa_searcher::retrieval` in `sa-index`.
//!
//! This lives in its own crate so that `bitarray`, `text-compression`, `sa-mappings` and
//! `sa-index` can all issue hints without depending on one another.

/// Issues a non-blocking hardware prefetch hint for the cache line at `ptr`.
///
/// `ptr` is never dereferenced and need not be valid, aligned, or even mapped: on both supported
/// architectures this compiles to a hint instruction that the CPU is free to ignore and that
/// cannot fault. That is what makes the function safe despite taking a raw pointer, and it is why
/// callers may prefetch one element past the end of a slice without a bounds check.
///
/// On architectures other than x86-64 and aarch64 this is a no-op. The function is always defined
/// so that callers need no `cfg` guards.
///
/// # Why `inline(always)` rather than `#[inline]`
///
/// The body is a single instruction, and every caller is in another crate. The workspace sets no
/// `[profile.release]`, so there is no LTO and a cross-crate call is a real call unless the
/// function is inlined. A `call`/`ret` pair around one hint instruction costs more than the hint
/// saves, so an un-inlined `prefetch_read` is strictly worse than not prefetching at all — this
/// attribute is load-bearing, not a micro-optimization.
#[inline(always)]
pub fn prefetch_read<T>(ptr: *const T) {
    // The two architectures use different hint strengths on purpose.
    //
    // x86-64 `_MM_HINT_T0` requests the line into *all* cache levels including L1. The x86
    // prefetcher drops hints that would cause pressure, so asking for L1 is safe and gives the
    // shortest latency when the value is used soon after — which is the case here, since the
    // two-pass loops consume a batch immediately after issuing its hints.
    //
    // aarch64 `pldl1keep` is the closest equivalent: preload into L1, "keep" (temporal) locality,
    // meaning the line should be retained rather than streamed past. The alternative `pldl1strm`
    // would mark it as streaming and evict it early, which is wrong for data that is about to be
    // read.
    #[cfg(target_arch = "x86_64")]
    // SAFETY: `_mm_prefetch` is a pure hint — it never faults and never reads.
    unsafe {
        std::arch::x86_64::_mm_prefetch(ptr as *const i8, std::arch::x86_64::_MM_HINT_T0)
    }
    #[cfg(target_arch = "aarch64")]
    // SAFETY: `prfm` is a pure hint — it never faults and never reads. `readonly` and `nostack`
    // hold because the instruction touches neither memory nor the stack.
    unsafe {
        std::arch::asm!(
            "prfm pldl1keep, [{p}]",
            p = in(reg) ptr,
            options(nostack, preserves_flags, readonly)
        )
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = ptr;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The crate's whole safety claim is that a hint never faults, whatever the pointer. Callers
    /// depend on this: the two-pass loops prefetch a lookahead index that may sit past the end of
    /// the data, and skipping the bounds check is the point.
    #[test]
    fn prefetching_never_faults() {
        let data: Vec<u64> = (0..64).collect();

        prefetch_read(data.as_ptr());
        prefetch_read(unsafe { data.as_ptr().add(data.len() - 1) });

        // One past the end: a valid address to *form*, never to read.
        prefetch_read(unsafe { data.as_ptr().add(data.len()) });

        // Far past the end, and an unmapped low address. Both are hints the CPU will drop.
        prefetch_read((data.as_ptr() as usize + (1 << 20)) as *const u64);
        prefetch_read(0x1000 as *const u64);
    }

    /// Works for any `T`, including zero-sized types, since the pointer is only ever an address.
    #[test]
    fn prefetching_is_generic_over_the_pointee() {
        let bytes = [0u8; 8];
        let unit = ();

        prefetch_read(bytes.as_ptr());
        prefetch_read(&unit as *const ());
        prefetch_read(&bytes as *const [u8; 8]);
    }
}
