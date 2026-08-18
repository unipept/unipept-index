//! Best-effort transparent huge page advice for large anonymous buffers.
//!
//! WHEN the advice is issued decides whether it does anything at all, so both `with_capacity`
//! constructors call it on the fresh allocation and no caller should need to. `MADV_HUGEPAGE` on
//! an untouched region makes the *page faults that populate it* allocate 2 MB pages directly. On a
//! region that is already populated it does nothing of the sort: the 4 KB pages are already there,
//! and all the advice buys is eligibility for khugepaged to collapse them in the background — at a
//! default 4096 pages per 10 s scan, which for a 160 GB suffix array is on the order of a day. The
//! advice used to be issued by the loaders after `read_binary` had filled the buffer, which is
//! exactly that useless case.
//!
//! Consequence worth knowing when reading a benchmark: on a box with
//! `/sys/kernel/mm/transparent_hugepage/enabled` at `[always]` this is all moot, because every
//! anonymous mapping gets huge pages regardless. It only matters at `[madvise]`.

/// Requests transparent huge pages over a bit array's backing words.
///
/// Call this before the region is written to; see the [module docs](self).
pub fn advise(data: &[u64]) {
    advise_raw(data.as_ptr() as usize, std::mem::size_of_val(data));
}

/// Requests transparent huge pages over a `Vec`'s whole allocation, reserved capacity included.
///
/// [`advise`] only covers the initialised prefix, which for a `Vec::with_capacity` that has not
/// been filled yet is nothing at all — and "not filled yet" is precisely when the advice is worth
/// issuing. The caller must not let the `Vec` reallocate afterwards, or the advice is left behind
/// on the old allocation.
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
