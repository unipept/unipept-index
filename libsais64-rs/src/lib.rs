// ignore errors because of different style in c code and import the c bindings
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use std::ptr::null_mut;
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

/// Builds the suffix array over the `text` using the libsais64 algorithm
///
/// # Arguments
/// * `text` - The text used for suffix array construction
///
/// # Returns
///
/// Returns Some with the suffix array build over the text if construction succeeds
/// Returns None if construction of the suffix array failed
pub fn sais64(text: &Vec<u16>, sparseness_factor: u8) -> Option<Vec<i64>> {
    let mut sa = vec![0; text.len()];
    let exit_code = unsafe { libsais16x64(text.as_ptr(), sa.as_mut_ptr(), text.len() as i64, 0, null_mut()) };
    if exit_code == 0 {
        let sparseness_factor = sparseness_factor as i64;
        for elem in sa.iter_mut() {
            *elem *= sparseness_factor;
        }
        Some(sa) 
    } else { None }
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
