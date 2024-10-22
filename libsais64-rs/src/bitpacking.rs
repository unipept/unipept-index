

// Function to get the rank of a character
fn get_rank(c: u8) -> u8 {
    match c {
        b'$' => 0,
        b'-' => 1,
        _ => 2 + (c - b'A'),
    }
}

pub const BITS_PER_CHAR: usize = 5;
pub fn bitpack_text_16(text: &Vec<u8>, sparseness_factor: usize) -> Vec<u16> {

    let num_ints = (text.len() + (sparseness_factor-1)) / sparseness_factor;
    let mut text_packed = vec![0; num_ints];

    if text.len() == 0 {
        return text_packed;
    }

    for i in 0..(num_ints-1) {
        let ti = i * sparseness_factor;
        let mut element = 0u16;
        for j in 0..sparseness_factor {
            let rank_c = get_rank(text[ti + j]) as u16;
            element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - j));
        }
        text_packed[i] = element;
    }

    // Handle the last element
    let mut last_element = 0u16;
    let last_el_start = sparseness_factor * (num_ints - 1);
    for i in 0..((text.len() - 1) % sparseness_factor + 1) {
        let rank_c = get_rank(text[last_el_start + i]) as u16;
        last_element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - i));
    }
    text_packed[num_ints - 1] = last_element;

    text_packed

}

pub fn bitpack_text_32(text: &Vec<u8>, sparseness_factor: usize) -> Vec<u32> {

    let num_ints = (text.len() + (sparseness_factor-1)) / sparseness_factor;
    let mut text_packed = vec![0; num_ints];

    if text.len() == 0 {
        return text_packed;
    }

    for i in 0..(num_ints-1) {
        let ti = i * sparseness_factor;
        let mut element = 0u32;
        for j in 0..sparseness_factor {
            let rank_c = get_rank(text[ti + j]) as u32;
            element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - j));
        }
        text_packed[i] = element;
    }

    // Handle the last element
    let mut last_element = 0u32;
    let last_el_start = sparseness_factor * (num_ints - 1);
    for i in 0..((text.len() - 1) % sparseness_factor + 1) {
        let rank_c = get_rank(text[last_el_start + i]) as u32;
        last_element |= rank_c << (BITS_PER_CHAR * (sparseness_factor - 1 - i));
    }
    text_packed[num_ints - 1] = last_element;

    text_packed

}