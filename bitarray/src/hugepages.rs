//! Best-effort transparent huge page advice for large anonymous buffers.
//!
//! WHEN the advice is issued decides whether it does anything at all, so both `with_capacity`
//! constructors call [`advise_capacity`] on the reserved allocation *before* zeroing it, and no
//! caller should need to. `MADV_HUGEPAGE` on an untouched region makes the *page faults that
//! populate it* allocate 2 MB pages directly. On a region that is already populated it does
//! nothing of the sort: the 4 KB pages are already there, and all the advice buys is eligibility
//! for khugepaged to collapse them in the background — at a default 4096 pages per 10 s scan,
//! which for a 160 GB suffix array is on the order of a day.
//!
//! Two versions of that mistake have already been made here. The advice was first issued by the
//! loaders, after `read_binary` had filled the buffer. Moving it into the constructors did not fix
//! it: they reserved the allocation and then `resize`d it to zero, and that `resize` memsets the
//! whole buffer — 1 GiB of `Vec<u64>` goes from 1.5 MB resident after `try_reserve_exact` to
//! 1025 MB after `resize`, so every page was faulted in before the advice was issued. Anything
//! that writes to the buffer, zeroes included, counts as populating it.
//!
//! Consequence worth knowing when reading a benchmark: on a box with
//! `/sys/kernel/mm/transparent_hugepage/enabled` at `[always]` this is all moot, because every
//! anonymous mapping gets huge pages regardless. It only matters at `[madvise]`.

/// Requests transparent huge pages over a `Vec`'s whole allocation, reserved capacity included.
///
/// Deliberately the capacity and not the initialised prefix: a `Vec` that has only been reserved
/// has no initialised prefix at all, and "reserved but not yet written" is precisely when the
/// advice is worth issuing — see the [module docs](self). So the call goes between the reservation
/// and the loop or `resize` that fills it.
///
/// The caller must not let the `Vec` reallocate afterwards, or the advice is left behind on the old
/// allocation. Reserving the final size up front, as every caller here does, is what guarantees
/// that.
pub fn advise_capacity<T>(vec: &Vec<T>) {
    advise_raw(vec.as_ptr() as usize, vec.capacity() * std::mem::size_of::<T>());
}

/// Cutting TLB/page-walk cost on the large random-access buffers the index preloads. No-op off
/// Linux; errors are ignored (the kernel may lack THP). This anonymous memory is THP-eligible,
/// unlike a file-backed mmap.
#[cfg(target_os = "linux")]
fn advise_raw(start: usize, len: usize) {
    // SAFETY: sysconf is always safe to call; a non-positive result is guarded below.
    let page = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
    if page <= 0 {
        return;
    }
    let page = page as usize;
    let end = start + len;
    // Advise only the page-aligned interior (madvise requires page-aligned bounds).
    let aligned_start = (start + page - 1) & !(page - 1);
    let aligned_end = end & !(page - 1);
    if aligned_end > aligned_start {
        // SAFETY: advises a sub-range of a live, page-aligned allocation; MADV_HUGEPAGE is
        // a hint that never reads, frees, or moves the memory.
        unsafe {
            libc::madvise(aligned_start as *mut libc::c_void, aligned_end - aligned_start, libc::MADV_HUGEPAGE);
        }
    }
}

#[cfg(not(target_os = "linux"))]
fn advise_raw(_start: usize, _len: usize) {}
