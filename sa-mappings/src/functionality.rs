//! This module contains the FunctionAggregator struct that is responsible for aggregating the functional annotations of proteins.

use crate::proteins::Protein;

/// A struct that represents a function aggregator
pub struct FunctionAggregator {}

impl FunctionAggregator {
    /// Aggregates the functional annotations of proteins
    /// 
    /// # Arguments
    /// * `proteins` - A vector of proteins
    /// 
    /// # Returns
    /// 
    /// Returns a string containing the aggregated functional annotations
    pub fn aggregate(&self, proteins: Vec<Protein>) -> String {
        proteins
            .iter()
            .map(|protein| protein.get_functional_annotations())
            .collect::<Vec<String>>()
            .join(";")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aggregate() {
        let proteins = vec![
            Protein {
                uniprot_id: "uniprot1".to_string(),
                sequence: (0, 3),
                taxon_id: 1,
                functional_annotations: vec![0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27],
            },
            Protein {
                uniprot_id: "uniprot2".to_string(),
                sequence: (4, 3),
                taxon_id: 2,
                functional_annotations: vec![0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27],
            },
        ];

        let function_aggregator = FunctionAggregator {};

        assert_eq!(function_aggregator.aggregate(proteins), "GO:0009279;IPR:IPR016364;IPR:IPR008816;GO:0009279;IPR:IPR016364;IPR:IPR008816");
    }
}
