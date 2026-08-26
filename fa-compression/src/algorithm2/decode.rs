//! This module provides a function to decode a byte array into a string representation of
//! annotations.

use super::CompressionTable;

/// Decodes a byte slice using a compression table and returns the corresponding string.
///
/// Each group of 3 bytes is read as a little-endian index into `compression_table`, and the
/// annotations it names are joined with `;`.
///
/// # Arguments
///
/// * `input` - The byte slice to decode. Must be the output of [`encode()`](super::encode()) against
///   **this same table**; the bytes carry no way to check that.
/// * `compression_table` - The compression table used for decoding. **Consumed**: the table is
///   taken by value, so decoding twice against the same table needs two tables.
///
/// # Returns
///
/// The decoded string.
///
/// # Panics
///
/// If any 3-byte group names an index the table does not have. There is no validation step and no
/// error return: decoding with a table that is not the one used to encode either panics or, worse,
/// silently yields the wrong annotations.
///
/// # A trailing partial group is dropped
///
/// Input is read in exact 3-byte chunks, so a length that is not a multiple of 3 has its final 1 or
/// 2 bytes **silently ignored** rather than reported.
///
/// ```
/// use fa_compression::algorithm2::{CompressionTable, decode};
///
/// let mut table = CompressionTable::new();
/// table.add_entry("IPR:IPR000001".to_string());
///
/// // The stray 4th byte is discarded.
/// assert_eq!(decode(&[0, 0, 0, 0], table), "IPR:IPR000001");
/// ```
///
/// # Examples
///
/// ```
/// use fa_compression::algorithm2::decode;
/// use fa_compression::algorithm2::CompressionTable;
///
/// let input = &[0, 0, 0, 1, 0, 0];
/// let mut compression_table = CompressionTable::new();
/// compression_table.add_entry("IPR:IPR000001".to_string());
/// compression_table.add_entry("IPR:IPR000002".to_string());
///
/// let decoded_string = decode(input, compression_table);
/// assert_eq!(decoded_string, "IPR:IPR000001;IPR:IPR000002");
/// ```
pub fn decode(input: &[u8], compression_table: CompressionTable) -> String {
    if input.is_empty() {
        return String::new();
    }

    let mut result = String::with_capacity(input.len() / 3 * 15);
    for bytes in input.chunks_exact(3) {
        // Convert the first 3 bytes to a u32 and use it as an index in the compression table
        let index = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], 0]) as usize;
        result.push_str(&compression_table[index].annotation);
        result.push(';');
    }

    // Remove the trailing semicolon
    result.pop();

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_compresion_table() -> CompressionTable {
        let mut table = CompressionTable::new();

        table.add_entry("IPR:IPR000001".to_string());
        table.add_entry("IPR:IPR000002".to_string());
        table.add_entry("IPR:IPR000003".to_string());
        table.add_entry("IPR:IPR000004".to_string());
        table.add_entry("GO:0000001".to_string());
        table.add_entry("GO:0000002".to_string());
        table.add_entry("GO:0000003".to_string());
        table.add_entry("EC:1.1.1.-".to_string());
        table.add_entry("EC:2.12.3.7".to_string());
        table.add_entry("EC:2.2.-.-".to_string());

        table
    }

    #[test]
    fn test_decode_empty() {
        let table = create_compresion_table();
        assert_eq!(decode(&[], table), "")
    }

    #[test]
    fn test_decode_single_ec() {
        let table = create_compresion_table();
        assert_eq!(decode(&[8, 0, 0], table), "EC:2.12.3.7");
    }

    #[test]
    fn test_decode_single_go() {
        let table = create_compresion_table();
        assert_eq!(decode(&[6, 0, 0], table), "GO:0000003");
    }

    #[test]
    fn test_decode_single_ipr() {
        let table = create_compresion_table();
        assert_eq!(decode(&[0, 0, 0], table), "IPR:IPR000001");
    }

    #[test]
    fn test_decode_all() {
        let table = create_compresion_table();
        assert_eq!(
            decode(&[0, 0, 0, 7, 0, 0, 2, 0, 0, 5, 0, 0], table),
            "IPR:IPR000001;EC:1.1.1.-;IPR:IPR000003;GO:0000002"
        )
    }
}
