//! Module for managing partitioned search patterns.
//!
//! `SearchPattern` is used to divide a pattern into evenly (or almost evenly)
//! distributed parts for multi-phase or parallel search. This is particularly
//! useful when applying bidirectional or block-based pattern matching algorithms.

use std::error::Error;

/// A pattern split into multiple parts for segmented searching.
///
/// `SearchPattern` splits a given byte pattern into approximately equal-length
/// chunks. This is helpful in FM-index search scenarios where patterns are
/// extended part by part.
///
/// # Fields
/// - `parts`: Vector of pattern parts, in original order.
/// - `length`: Total length of the original pattern.
pub struct SearchPattern {
    parts: Vec<Vec<u8>>,
    length: usize
}

impl SearchPattern {
    
    /// Creates a new `SearchPattern` by splitting a pattern into `parts_amount` segments.
    ///
    /// # Arguments
    /// * `pattern` - A byte vector representing the full pattern to search.
    /// * `parts_amount` - Number of parts to split the pattern into.
    ///
    /// # Errors
    /// Returns an error if:
    /// - `parts_amount` is zero.
    /// - `parts_amount` is greater than the length of the pattern.
    pub fn new(pattern: Vec<u8>, parts_amount: usize) -> Result<Self, Box<dyn Error + Send + Sync>> {

        if parts_amount == 0 {
            return Err("parts_amount must be greater than 0".into());
        }

        let total_len = pattern.len();

        if total_len < parts_amount {
            return Err("parts_amount is greater than pattern length".into());
        }

        let base_size = total_len / parts_amount;
        let remainder = total_len % parts_amount;

        let mut parts = Vec::with_capacity(parts_amount);
        let mut start = 0;

        for i in 0..parts_amount {
            let extra = if i < remainder { 1 } else { 0 };
            let end = start + base_size + extra;
            parts.push(pattern[start..end].to_vec());
            start = end;
        }

        let length = total_len;

        Ok(SearchPattern { parts, length })


    }

    /// Returns the part at the given index, optionally reversed.
    ///
    /// # Arguments
    /// * `index` - The index of the part to retrieve.
    /// * `direction_left` - Whether the part should be reversed for right-to-left search.
    ///
    /// # Panics
    /// Panics if `index` is out of bounds.
    pub fn get_part(&self, index: u8, direction_left: bool) -> Vec<u8> {
        let mut part = self.parts.get(index as usize).unwrap().clone();
        if direction_left { part.reverse() };
        part
    }

    /// Returns the length of the part at the specified index.
    ///
    /// # Panics
    /// Panics if `index` is out of bounds.
    pub fn get_part_len(&self, index: u8) -> usize {
        self.parts.get(index as usize).unwrap().len()
    }

    /// Returns the total length of the original pattern.
    pub fn len(&self) -> usize {
        self.length
    }
}