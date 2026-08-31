// ignore errors because of different style in c code and import the c bindings
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use std::ptr::null_mut;

use crate::bitpacking::{BITS_PER_CHAR, bitpack_text_8, bitpack_text_16, bitpack_text_32};
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub mod bitpacking;

/// Builds the suffix array over the `text` using the libsais algorithm
///
/// # Arguments
/// * `text` - The text used for suffix array construction
///
/// # Returns
///
/// Returns Some with the suffix array build over the text if construction succeeds
/// Returns None if construction of the suffix array failed
pub fn sais64(text: Vec<u8>, libsais_sparseness: usize) -> Result<Vec<i64>, &'static str> {
    let mut sa;

    let required_bits = libsais_sparseness * BITS_PER_CHAR;
    let exit_code = if required_bits <= 8 {
        // bitpacked values fit in uint8_t
        let packed_text = if libsais_sparseness == 1 { text } else { bitpack_text_8(text, libsais_sparseness) };

        sa = vec![0; packed_text.len()];
        // SAFETY: `T` points at `packed_text`'s `n` initialised `u8`s and libsais only reads them;
        // `SA` points at `sa`, which was just allocated with exactly `n` `i64`s, and libsais writes
        // `n + fs` of them with `fs = 0`. Both vectors outlive the call, and `freq` is documented as
        // optional, so a null pointer means "do not collect frequencies". `n` cannot overflow `i64`
        // because a `Vec`'s length is at most `isize::MAX`.
        unsafe { libsais64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, 0, null_mut()) }
    } else if required_bits <= 16 {
        // bitpacked values fit in uint16_t
        let packed_text = bitpack_text_16(text, libsais_sparseness);
        sa = vec![0; packed_text.len()];
        // SAFETY: as above, with `T` a `*const u16` matching `bitpack_text_16`'s `Vec<u16>`.
        unsafe { libsais16x64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, 0, null_mut()) }
    } else {
        let packed_text = bitpack_text_32(text, libsais_sparseness);
        sa = vec![0; packed_text.len()];
        let k = 1 << (libsais_sparseness * BITS_PER_CHAR);
        // SAFETY: as above, with `T` a `*const u32` matching `bitpack_text_32`'s `Vec<u32>`. This
        // entry point also takes the alphabet size `k`, and libsais indexes a `k`-element bucket
        // array with the values it reads: `bitpack_text_32` packs `libsais_sparseness` characters of
        // `BITS_PER_CHAR` bits each, so every value is below the `k` computed on the line above.
        unsafe { libsais32x64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, k, 0, null_mut()) }
    };

    if exit_code == 0 {
        for elem in sa.iter_mut() {
            let libsais_sparseness = libsais_sparseness as i64;
            *elem *= libsais_sparseness;
        }
        Ok(sa)
    } else {
        Err("Failed building suffix array")
    }
}

#[cfg(test)]
mod tests {
    use crate::sais64;

    #[test]
    fn check_build_sa_with_libsais64() {
        let sparseness_factor = 4;
        let text = "BANANA-BANANA$".as_bytes().to_vec();
        let sa = sais64(text, sparseness_factor);
        let correct_sa: Vec<i64> = vec![12, 8, 0, 4];
        assert_eq!(sa, Ok(correct_sa));
    }
}
