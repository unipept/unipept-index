mod decode;
mod encode;

use std::ops::Index;

pub use decode::decode;
pub use encode::encode;

pub struct CompressionTableEntry {
    pub annotation: String
}

pub struct CompressionTable {
    pub entries: Vec<CompressionTableEntry>
}

impl CompressionTable {
    pub fn new() -> CompressionTable {
        CompressionTable {
            entries: Vec::new()
        }
    }

    pub fn add_entry(&mut self, annotation: String) {
        self.entries.push(CompressionTableEntry {
            annotation
        });
    }

    pub fn index_of(&self, annotation: &str) -> Option<usize> {
        self.entries
            .iter()
            .position(|entry| entry.annotation == annotation)
    }
}

impl Default for CompressionTable {
    fn default() -> Self {
        Self::new()
    }
}

impl Index<usize> for CompressionTable {
    type Output = CompressionTableEntry;

    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_compresion_table() -> CompressionTable {
        let mut table = CompressionTable::new();

        table.add_entry("IPR:IPR000001".to_string());
        table.add_entry("IPR:IPR000002".to_string());
        table.add_entry("GO:0000001".to_string());
        table.add_entry("GO:0000002".to_string());
        table.add_entry("EC:1.1.1.-".to_string());

        table
    }

    #[test]
    fn test_add_entry() {
        assert_eq!(create_compresion_table().entries.len(), 5);
    }

    #[test]
    fn test_index_of() {
        let table = create_compresion_table();

        assert_eq!(table.index_of("IPR:IPR000001"), Some(0));
        assert_eq!(table.index_of("IPR:IPR000002"), Some(1));
        assert_eq!(table.index_of("GO:0000001"), Some(2));
        assert_eq!(table.index_of("GO:0000002"), Some(3));
        assert_eq!(table.index_of("EC:1.1.1.-"), Some(4));
    }

    #[test]
    fn test_index_of_not_found() {
        let table = create_compresion_table();

        assert_eq!(table.index_of("IPR:IPR000003"), None);
        assert_eq!(table.index_of("GO:0000003"), None);
        assert_eq!(table.index_of("EC:2.2.2.-"), None);
    }

    #[test]
    fn test_index() {
        let table = create_compresion_table();

        assert_eq!(table[0].annotation, "IPR:IPR000001");
        assert_eq!(table[1].annotation, "IPR:IPR000002");
        assert_eq!(table[2].annotation, "GO:0000001");
        assert_eq!(table[3].annotation, "GO:0000002");
        assert_eq!(table[4].annotation, "EC:1.1.1.-");
    }
}
