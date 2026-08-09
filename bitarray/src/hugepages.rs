//! Best-effort transparent huge page advice for large anonymous bit-array buffers.

/// Requests transparent huge pages (2 MB) over `data`, cutting TLB/page-walk cost on the
/// large random-access buffers this crate backs (the preloaded suffix array + protein
/// text). No-op off Linux; errors are ignored (the kernel may lack THP). This anonymous
/// memory is THP-eligible, unlike a file-backed mmap.
#[cfg(target_os = "linux")]
pub(crate) fn advise(data: &[u64]) {
    // SAFETY: sysconf is always safe to call; a non-positive result is guarded below.
    let page = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
    if page <= 0 {
        return;
    }
    let page = page as usize;
    let start = data.as_ptr() as usize;
    let end = start + std::mem::size_of_val(data);
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
pub(crate) fn advise(_data: &[u64]) {}
