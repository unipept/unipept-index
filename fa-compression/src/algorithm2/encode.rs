use super::CompressionTable;

pub fn encode(input: &str, compression_table: CompressionTable) -> Vec<u8> {
    if input.is_empty() {
        return Vec::new();
    }

    let mut encoded: Vec<u8> = Vec::with_capacity(input.len() / 3);

    for annotation in input.split(';') {
        if let Some(index) = compression_table.index_of(annotation) {
            encoded.extend_from_slice(&index.to_le_bytes()[0..3])
        }
    }

    encoded
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
    fn test_encode_empty() {
        let table = create_compresion_table();
        assert_eq!(encode("", table), vec![])
    }

    #[test]
    fn test_encode_single_ec() {
        let table = create_compresion_table();
        assert_eq!(encode("EC:2.12.3.7", table), vec![8, 0, 0])
    }

    #[test]
    fn test_encode_single_go() {
        let table = create_compresion_table();
        assert_eq!(encode("GO:0000003", table), vec![6, 0, 0])
    }

    #[test]
    fn test_encode_single_ipr() {
        let table = create_compresion_table();
        assert_eq!(encode("IPR:IPR000002", table), vec![1, 0, 0])
    }

    #[test]
    fn test_encode_all() {
        let table = create_compresion_table();
        assert_eq!(encode("IPR:IPR000001;EC:1.1.1.-;IPR:IPR000003;GO:0000002", table), vec![ 0, 0, 0, 7, 0, 0, 2, 0, 0, 5, 0, 0 ])
    }
}
