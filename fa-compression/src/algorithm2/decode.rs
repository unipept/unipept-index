use super::CompressionTable;

pub fn decode(input: &[u8], compression_table: CompressionTable) -> String {
    if input.is_empty() {
        return String::new();
    }

    let mut result = String::with_capacity(input.len() / 3 * 15);
    for bytes in input.chunks_exact(3) {
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
