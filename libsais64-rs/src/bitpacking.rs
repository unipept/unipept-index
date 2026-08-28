/// Maps a text byte to its 5-bit rank, or `None` if the byte is not one the text can hold.
///
/// The alphabet is `$`, `-`, and `A`–`Z`, giving ranks 0..=27, which is what fits in
/// [`BITS_PER_CHAR`]. Everything else is rejected rather than computed:
///
/// * a byte below `b'A'` (every digit, space and most punctuation) underflows `c - b'A'` — a panic
///   in debug, a wrapped rank in release;
/// * a byte at or above 95, which includes every non-ASCII UTF-8 byte, produces a rank above 31.
///   The packing loops OR that in, so the excess bits land in the *neighbouring* residue's field
///   and silently corrupt a character the caller never supplied.
///
/// In practice `sa-builder` cannot reach either case — it packs the text *after* it has been
/// through the 5-bit protein encoding, so every byte is already in the alphabet. This is a `pub`
/// module though, and a caller packing raw FASTA would hit both.
fn get_rank(c: u8) -> Option<u8> {
    match c {
        b'$' => Some(0),
        b'-' => Some(1),
        b'A'..=b'Z' => Some(2 + (c - b'A')),
        _ => None
    }
}

/// The error a packing function returns when the text holds a byte outside the alphabet.
#[derive(Debug, PartialEq, Eq)]
pub struct UnsupportedCharacter {
    /// The offending byte.
    pub byte: u8,
    /// Where it sits in the text.
    pub index: usize
}

impl std::fmt::Display for UnsupportedCharacter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "byte {:#04x} at text position {} is not in the protein alphabet ($, -, A-Z)", self.byte, self.index)
    }
}

impl std::error::Error for UnsupportedCharacter {}

/// Looks up one rank, reporting where an unsupported byte was found.
#[inline]
fn rank_at(text: &[u8], index: usize) -> Result<u8, UnsupportedCharacter> {
    get_rank(text[index]).ok_or(UnsupportedCharacter { byte: text[index], index })
}

// Amount of bits necessary to represent one character in the protein text.
pub const BITS_PER_CHAR: usize = 5;

// Bitpack text in a vector of u16 elements. BITS_PER_CHAR * sparseness_factor <= 16.
pub fn bitpack_text_16(text: Vec<u8>, sparseness_factor: usize) -> Result<Vec<u16>, UnsupportedCharacter> {
    assert!(BITS_PER_CHAR * sparseness_factor <= 16);

    let num_ints = text.len().div_ceil(sparseness_factor);
    let mut text_packed = vec![0; num_ints];

    if text.is_empty() {
        return Ok(text_packed);
    }

    for (i, element) in text_packed.iter_mut().enumerate().take(num_ints - 1) {
        let ti = i * sparseness_factor;
        *element = 0u16;
        for j in 0..sparseness_factor {
            let rank_c = rank_at(&text, ti + j)? as u16;
            *element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - j));
        }
    }

    // Handle the last element
    let mut last_element = 0u16;
    let last_el_start = sparseness_factor * (num_ints - 1);
    for i in 0..((text.len() - 1) % sparseness_factor + 1) {
        let rank_c = rank_at(&text, last_el_start + i)? as u16;
        last_element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - i));
    }
    text_packed[num_ints - 1] = last_element;

    Ok(text_packed)
}

// Bitpack text in a vector of u32 elements. BITS_PER_CHAR * sparseness_factor <= 32.
pub fn bitpack_text_32(text: Vec<u8>, sparseness_factor: usize) -> Result<Vec<u32>, UnsupportedCharacter> {
    assert!(BITS_PER_CHAR * sparseness_factor <= 32);

    let num_ints = text.len().div_ceil(sparseness_factor);
    let mut text_packed = vec![0; num_ints];

    if text.is_empty() {
        return Ok(text_packed);
    }

    for (i, element) in text_packed.iter_mut().enumerate().take(num_ints - 1) {
        let ti = i * sparseness_factor;
        *element = 0u32;
        for j in 0..sparseness_factor {
            let rank_c = rank_at(&text, ti + j)? as u32;
            *element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - j));
        }
    }

    // Handle the last element
    let mut last_element = 0u32;
    let last_el_start = sparseness_factor * (num_ints - 1);
    for i in 0..((text.len() - 1) % sparseness_factor + 1) {
        let rank_c = rank_at(&text, last_el_start + i)? as u32;
        last_element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - i));
    }
    text_packed[num_ints - 1] = last_element;

    Ok(text_packed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every byte the protein text can actually hold maps to a rank that fits in `BITS_PER_CHAR`.
    #[test]
    fn every_alphabet_byte_fits_the_field() {
        let max = (1u8 << BITS_PER_CHAR) - 1;
        for c in (*b"$-").into_iter().chain(b'A'..=b'Z') {
            let rank = get_rank(c).unwrap_or_else(|| panic!("{} should be in the alphabet", c as char));
            assert!(rank <= max, "{} -> rank {rank} does not fit in {BITS_PER_CHAR} bits", c as char);
        }
    }

    /// Everything else is rejected rather than computed. Unchecked, the two halves of this range
    /// failed differently: below `A` the subtraction underflowed, and from 95 up the rank exceeded
    /// the field and spilled into the neighbouring residue.
    #[test]
    fn every_other_byte_is_rejected() {
        let alphabet: Vec<u8> = (*b"$-").into_iter().chain(b'A'..=b'Z').collect();
        for c in 0..=255u8 {
            if alphabet.contains(&c) {
                continue;
            }
            assert_eq!(get_rank(c), None, "byte {c:#04x} should not be in the alphabet");
        }
    }

    /// The packers report the offending byte and its position rather than corrupting the output.
    #[test]
    fn packing_reports_an_unsupported_byte() {
        assert!(bitpack_text_16(b"ACGT-ACGT$".to_vec(), 3).is_ok());

        let err = bitpack_text_16(b"ACG*T$".to_vec(), 3).expect_err("'*' is not in the alphabet");
        assert_eq!(err.byte, b'*');
        assert_eq!(err.index, 3);

        // Non-ASCII: every byte of a multi-byte character is >= 128, the class that used to
        // overflow the 5-bit field.
        let err = bitpack_text_32("ACGé$".as_bytes().to_vec(), 5).expect_err("non-ASCII is not in the alphabet");
        assert!(err.byte >= 128, "expected a non-ASCII byte, got {:#04x}", err.byte);
    }

    /// Empty input is not an error, and neither packer indexes anything.
    #[test]
    fn empty_text_packs_to_nothing() {
        assert_eq!(bitpack_text_16(Vec::new(), 3), Ok(Vec::new()));
        assert_eq!(bitpack_text_32(Vec::new(), 6), Ok(Vec::new()));
    }
}
