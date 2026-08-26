//! A single portable software-prefetch hint.
//!
//! Binary search and protein retrieval walk the index in an order the hardware prefetcher cannot
//! predict, which is exactly the case software prefetching exists for: the address of the *next*
//! access is known several iterations before the value is needed, so the load can be started
//! early and its ~80-100 ns latency overlapped with useful work. See the [crate docs](crate) for
//! how large the structures being walked are, and why the sibling [`hugepages`](crate::hugepages)
//! hint attacks the same problem from the other end.
//!
//! Callers pair this with two-pass batching (fill a batch and issue hints, then process the
//! batch) rather than prefetching one element ahead. Both loops are methods on `sa_searcher`'s
//! `Searcher`: the private `iterate_sa_range`, and `retrieve_proteins`, which is public but is
//! written in the private `sa_searcher::retrieval` module.

/// Issues a non-blocking hardware prefetch hint for the cache line at `ptr`.
///
/// `ptr` is never dereferenced and need not be valid, aligned, or even mapped: on both supported
/// architectures this compiles to a hint instruction that the CPU is free to ignore and that
/// cannot fault. That is what makes the function safe despite taking a raw pointer.
///
/// It does not remove the caller's bounds check. Every call site in the workspace forms its
/// pointer by indexing — `&self.blocks[word]`, `&self.mmap[off]` — and that panics on an
/// out-of-range index in safe Rust, before the hint instruction the CPU would have dropped is
/// ever reached. What fault-freedom buys is that a *useless* hint costs nothing, so each guard
/// can be a cheap `<` compare that skips the hint rather than a check that has to be exact.
///
/// The lookahead position that genuinely may sit past the end is absorbed one level up: the
/// `prefetch_at`, `prefetch_sa_index` and `prefetch_for_suffix` trait methods each bounds-check
/// in their own address domain and return silently, so the two-pass loops calling them never
/// check. Those guards belong at the backend because they are not the same check — an element
/// index, a byte offset, a 16-byte span and a residue count that indexes a *word* array are four
/// different domains, and only the backend knows its own address arithmetic. See
/// `text_compression::ProteinTextBackend::prefetch_at` for that contract.
///
/// On architectures other than x86-64 and aarch64 this is a no-op. The function is always defined
/// so that callers need no `cfg` guards.
///
/// # Why `inline(always)` rather than `#[inline]`
///
/// The body is a single instruction, and every caller is in another crate. The workspace sets no
/// `[profile.release]`, so there is no cross-crate LTO: a call into this crate is a real call
/// unless the callee's body reaches the caller's codegen unit. Being generic already achieves
/// that much — a generic function's MIR is exported and monomorphised in the calling crate, so
/// LLVM sees this body and would almost certainly inline it unprompted. What `inline(always)`
/// adds is the *guarantee*, and at every opt-level rather than only where LLVM's cost model
/// happens to agree. That guarantee is worth having: a `call`/`ret` pair around one hint
/// instruction costs more than the hint saves, so an un-inlined `prefetch_read` would be
/// strictly worse than not prefetching at all.
///
/// Note that LTO would not change this reasoning — it was measured on the full index and
/// rejected; see the crate docs of `sa-index` for the numbers.
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

    /// The crate's whole safety claim is that a hint never faults, whatever the pointer. That is
    /// what lets a backend guard its hints with a cheap `<` compare and drop the ones that fall
    /// out of range, rather than having to prove every address is one it could legally read.
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
