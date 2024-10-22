// ignore errors because of different style in c code and import the c bindings
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use std::ptr::null_mut;
use crate::bitpacking::{bitpack_text_16, bitpack_text_32, BITS_PER_CHAR};
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub mod bitpacking;

/// Builds the suffix array over the `text` using the libsais64 algorithm
///
/// # Arguments
/// * `text` - The text used for suffix array construction
///
/// # Returns
///
/// Returns Some with the suffix array build over the text if construction succeeds
/// Returns None if construction of the suffix array failed
pub fn sais64(text: &Vec<u8>, sparseness_factor: u8) -> Result<Vec<i64>, &str> {
    let sparseness_factor = sparseness_factor as usize;
    let mut libsais_sparseness = sparseness_factor;
    let mut sa;
    let exit_code;

    if sparseness_factor * BITS_PER_CHAR <= 16 {
        // bitpacked values fit in uint16_t
        let packed_text = bitpack_text_16(text, libsais_sparseness);
        sa = vec![0; packed_text.len()];
        exit_code = unsafe { libsais16x64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, 0, null_mut()) };
    } else {
        // bitpacked values do not fit in uint16_t, use 32-bit text
        // set libsais_sparseness to highest sparseness factor fitting in 32-bit value and sparseness factor divisible by libsais sparseness
        // max 28 out of 32 bits used, because a bucket is created for every element of the alfabet 8 * 2^28).
        libsais_sparseness = 28 / BITS_PER_CHAR;
        while sparseness_factor % libsais_sparseness != 0 && libsais_sparseness * BITS_PER_CHAR > 16 {
            libsais_sparseness -= 1;
        }

        if libsais_sparseness * BITS_PER_CHAR <= 16 {
            return Err("invalid sparseness factor");
        }

        let packed_text = bitpack_text_32(text, libsais_sparseness);
        sa = vec![0; packed_text.len()];
        let k = 1 << (libsais_sparseness * BITS_PER_CHAR);
        exit_code = unsafe { libsais32x64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, k, 0, null_mut()) };
    }
    
    if exit_code == 0 {
        for elem in sa.iter_mut() {
            let libsais_sparseness = libsais_sparseness as i64;
            *elem *= libsais_sparseness;
        }
        Ok(sa) 
    } else { Err("Failed building suffix array") }
}

#[cfg(test)]
mod tests {
    use crate::sais64_long;

    #[test]
    fn check_build_sa_with_libsais64() {
        let bits_per_char = 5;
        let sparseness_factor = 4;
        let mut text = [100834,   // BANA
            493603,                     // NA-B
            80975,                      // ANAN
            65536                       // A$
        ].to_vec();
        let sa = sais64_long(&mut text, 1 << (bits_per_char * sparseness_factor), sparseness_factor);
        assert_eq!(sa, Some(vec![12, 8, 0, 4]));
    }
}
