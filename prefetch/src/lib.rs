/// Issues a non-blocking hardware prefetch hint for the cache line at `ptr`.
/// On unsupported platforms this is a no-op. The function itself is always
/// defined so callers don't need `cfg` guards.
#[inline(always)]
pub fn prefetch_read<T>(ptr: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: `_mm_prefetch` is a pure hint — it never faults and never reads.
    unsafe { std::arch::x86_64::_mm_prefetch(ptr as *const i8, std::arch::x86_64::_MM_HINT_T0) }
    #[cfg(target_arch = "aarch64")]
    // SAFETY: `prfm` is a pure hint — it never faults and never reads.
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
