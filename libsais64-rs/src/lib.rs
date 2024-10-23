// ignore errors because of different style in c code and import the c bindings
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use std::ptr::null_mut;
use crate::bitpacking::{bitpack_text_16, bitpack_text_32, bitpack_text_8, BITS_PER_CHAR};
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
pub fn sais64(text: &Vec<u8>, libsais_sparseness: usize) -> Result<Vec<i64>, &str> {
    let exit_code;
    let mut sa;

    let required_bits = libsais_sparseness * BITS_PER_CHAR;
    if required_bits <= 8 {
        // bitpacked values fit in uint8_t
        let packed_text = bitpack_text_8(text, libsais_sparseness);
        sa = vec![0; packed_text.len()];
        exit_code = unsafe { libsais64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, 0, null_mut()) };
    } else if required_bits <= 16 {
        // bitpacked values fit in uint16_t
        let packed_text = bitpack_text_16(text, libsais_sparseness);
        sa = vec![0; packed_text.len()];
        exit_code = unsafe { libsais16x64(packed_text.as_ptr(), sa.as_mut_ptr(), packed_text.len() as i64, 0, null_mut()) };
    } else {
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
    use crate::sais64;

    #[test]
    fn check_build_sa_with_libsais64() {
        let bits_per_char = 5;
        let sparseness_factor = 4;
        let mut text = [100834,   // BANA
            493603,                     // NA-B
            80975,                      // ANAN
            65536                       // A$
        ].to_vec();
        let sa = sais64(&mut text, sparseness_factor);
        assert_eq!(sa, Some(vec![12, 8, 0, 4]));
    }
}
