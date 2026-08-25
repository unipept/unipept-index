//! Table-driven compression against a global annotation dictionary.
//!
//! Every distinct annotation in the database gets one [`CompressionTable`] entry, and an annotation
//! is then stored as its index in that table: 3 bytes each, for a compression ratio around **76%**
//! — better than [`algorithm1`](super::algorithm1)'s 68-71%. The costs are what decide between
//! them:
//!
//! * **The table is not optional.** Encoded bytes are meaningless without the exact table that
//!   produced them, so it has to be built over the whole database, kept alongside the data, and
//!   passed to every call. A single entry cannot be decoded in isolation.
//! * **Encoding is a linear scan.** Looking an annotation up walks every entry, once per
//!   annotation, so encoding costs *table size × annotations encoded*. On a database-sized table
//!   that is a different order of magnitude, not a constant factor.
//! * **The index space is 24 bits.** Only three of the four index bytes are stored, so a table
//!   beyond `2^24` entries cannot be addressed.
//!
//! [`algorithm1`](super::algorithm1) is what the index and the server use; this module is the
//! higher-ratio alternative for cases where a whole-database table is acceptable.

mod decode;
mod encode;

use std::ops::Index;

pub use decode::decode;
pub use encode::encode;

/// A single annotation stored in a [`CompressionTable`].
///
/// Obtained by indexing the table with an annotation's index, which is exactly what an encoded
/// 3-byte value is.
pub struct CompressionTableEntry {
    annotation: String
}

impl CompressionTableEntry {
    /// The annotation this entry holds, prefix included (`"IPR:IPR000001"`, `"EC:1.1.1.-"`, ...).
    ///
    /// # Examples
    ///
    /// ```
    /// use fa_compression::algorithm2::CompressionTable;
    ///
    /// let mut table = CompressionTable::new();
    /// table.add_entry("IPR:IPR000001".to_string());
    ///
    /// assert_eq!(table[0].annotation(), "IPR:IPR000001");
    /// ```
    pub fn annotation(&self) -> &str {
        &self.annotation
    }
}

/// Represents a compression table.
pub struct CompressionTable {
    /// List of annotations in the compression table.
    entries: Vec<CompressionTableEntry>
}

impl CompressionTable {
    /// Creates a new compression table.
    ///
    /// # Returns
    ///
    /// An empty compression table.
    ///
    /// # Examples
    ///
    /// ```
    /// use fa_compression::algorithm2::CompressionTable;
    ///
    /// let table = CompressionTable::new();
    /// ```
    pub fn new() -> CompressionTable {
        CompressionTable { entries: Vec::new() }
    }

    /// Adds a new entry to the compression table.
    ///
    /// # Arguments
    ///
    /// * `annotation` - The annotation to add to the compression table.
    ///
    /// # Examples
    ///
    /// ```
    /// use fa_compression::algorithm2::CompressionTable;
    ///
    /// let mut table = CompressionTable::new();
    /// table.add_entry("IPR:IPR000001".to_string());
    /// table.add_entry("IPR:IPR000002".to_string());
    /// ```
    pub fn add_entry(&mut self, annotation: String) {
        self.entries.push(CompressionTableEntry { annotation });
    }

    /// Returns the index of the given annotation in the compression table, if it exists.
    ///
    /// A linear scan over every entry. This is what makes [`encode`] cost *table size × annotations
    /// encoded*; see the module docs.
    fn index_of(&self, annotation: &str) -> Option<usize> {
        self.entries.iter().position(|entry| entry.annotation == annotation)
    }
}

impl Default for CompressionTable {
    /// Creates a default compression table.
    fn default() -> Self {
        Self::new()
    }
}

impl Index<usize> for CompressionTable {
    type Output = CompressionTableEntry;

    /// Returns a reference to the compression table entry at the given index.
    ///
    /// # Panics
    ///
    /// If `index` is out of bounds, like any other slice index. [`decode`] indexes the table with
    /// values taken straight from its input, so decoding with the wrong table panics here.
    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Creates a compression table with some predefined entries for testing.
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
    fn test_default() {
        assert_eq!(CompressionTable::default().entries.len(), 0);
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
