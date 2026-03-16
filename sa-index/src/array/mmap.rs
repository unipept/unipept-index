use memmap2::Mmap;

/// Reads a u64 value in little-endian byte order from the given mmap at the given byte offset.
pub(super) fn read_u64_le(mmap: &Mmap, byte_offset: usize) -> u64 {
    let bytes: [u8; 8] = mmap[byte_offset..byte_offset + 8].try_into().unwrap();
    u64::from_le_bytes(bytes)
}

/// Returns the suffix array value at the given index from a memory-mapped file.
pub(super) fn get_mmap(mmap: &Mmap, data_offset: usize, bits_per_value: usize, index: usize) -> i64 {
    if bits_per_value == 64 {
        let offset = data_offset + index * 8;
        let bytes: [u8; 8] = mmap[offset..offset + 8].try_into().unwrap();
        i64::from_le_bytes(bytes)
    } else {
        let mask: u64 = (1u64 << bits_per_value) - 1;
        let bit_offset = index * bits_per_value;
        let start_block = bit_offset / 64;
        let start_block_offset = bit_offset % 64;
        let block_byte_offset = data_offset + start_block * 8;
        let start_val = read_u64_le(mmap, block_byte_offset);
        if start_block_offset + bits_per_value <= 64 {
            ((start_val >> (64 - start_block_offset - bits_per_value)) & mask) as i64
        } else {
            let end_block_offset = (index + 1) * bits_per_value % 64;
            let end_val = read_u64_le(mmap, block_byte_offset + 8);
            let a = start_val << end_block_offset;
            let b = end_val >> (64 - end_block_offset);
            ((a | b) & mask) as i64
        }
    }
}
