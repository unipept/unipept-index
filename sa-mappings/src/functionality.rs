//! This module contains the FunctionAggregator struct that is responsible for aggregating the
//! functional annotations of proteins.

use std::collections::{
    HashMap,
    HashSet
};

use serde::Serialize;

use crate::proteins::Protein;

/// A struct that represents the functional annotations once aggregated
#[derive(Debug, Serialize)]
pub struct FunctionalAggregation {
    /// A HashMap representing how many GO, EC and IPR terms were found
    pub counts: HashMap<String, usize>,
    /// A HashMap representing how often a certain functional annotation was found
    pub data:   HashMap<String, u32>
}

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
    /// Returns a JSON string containing the aggregated functional annotations
    pub fn aggregate(&self, proteins: Vec<&Protein>) -> FunctionalAggregation {
        // Keep track of the proteins that have a certain annotation
        let mut proteins_with_ec: HashSet<String> = HashSet::new();
        let mut proteins_with_go: HashSet<String> = HashSet::new();
        let mut proteins_with_ipr: HashSet<String> = HashSet::new();

        // Keep track of the counts of the different annotations
        let mut data: HashMap<String, u32> = HashMap::new();

        for protein in proteins.iter() {
            for annotation in protein.get_functional_annotations().split(';') {
                match annotation.chars().next() {
                    Some('E') => proteins_with_ec.insert(protein.uniprot_id.clone()),
                    Some('G') => proteins_with_go.insert(protein.uniprot_id.clone()),
                    Some('I') => proteins_with_ipr.insert(protein.uniprot_id.clone()),
                    _ => false
                };

                data.entry(annotation.to_string())
                    .and_modify(|c| *c += 1)
                    .or_insert(1);
            }
        }

        let mut counts: HashMap<String, usize> = HashMap::new();
        counts.insert("all".to_string(), proteins.len());
        counts.insert("EC".to_string(), proteins_with_ec.len());
        counts.insert("GO".to_string(), proteins_with_go.len());
        counts.insert("IPR".to_string(), proteins_with_ipr.len());

        data.remove("");

        FunctionalAggregation {
            counts,
            data
        }
    }

    /// Aggregates the functional annotations of proteins
    ///
    /// # Arguments
    /// * `proteins` - A vector of proteins
    ///
    /// # Returns
    ///
    /// Returns a list of lists with all the functional annotations per protein
    pub fn get_all_functional_annotations(&self, proteins: &[&Protein]) -> Vec<Vec<String>> {
        proteins
            .iter()
            .map(|&prot| {
                prot.get_functional_annotations()
                    .split(';')
                    .map(|ann| ann.to_string())
                    .filter(|s| !s.is_empty())
                    .collect()
            })
            .collect::<Vec<Vec<String>>>()
    }
}

#[cfg(test)]
mod tests {
    use fa_compression::algorithm1::encode;

    use super::*;

    #[test]
    fn test_aggregate() {
        let mut proteins: Vec<Protein> = Vec::new();
        proteins.push(Protein {
            uniprot_id:             "P12345".to_string(),
            taxon_id:               9606,
            functional_annotations: encode("GO:0001234;GO:0005678")
        });
        proteins.push(Protein {
            uniprot_id:             "P23456".to_string(),
            taxon_id:               9606,
            functional_annotations: encode("EC:1.1.1.-")
        });

        let function_aggregator = FunctionAggregator {};

        let result = function_aggregator.aggregate(proteins.iter().collect());

        assert_eq!(result.counts.get("all"), Some(&2));
        assert_eq!(result.counts.get("EC"), Some(&1));
        assert_eq!(result.counts.get("GO"), Some(&1));
        assert_eq!(result.counts.get("IPR"), Some(&0));
        assert_eq!(result.counts.get("NOTHING"), None);

        assert_eq!(result.data.get("GO:0001234"), Some(&1));
        assert_eq!(result.data.get("GO:0005678"), Some(&1));
        assert_eq!(result.data.get("EC:1.1.1.-"), Some(&1));
        assert_eq!(result.data.get("EC:1.1.2.-"), None);
    }

    #[test]
    fn test_get_all_functional_annotations() {
        let mut proteins: Vec<&Protein> = Vec::new();

        let protein1 = Protein {
            uniprot_id:             "P12345".to_string(),
            taxon_id:               9606,
            functional_annotations: encode("GO:0001234;GO:0005678")
        };
        let protein2 = Protein {
            uniprot_id:             "P23456".to_string(),
            taxon_id:               9606,
            functional_annotations: encode("EC:1.1.1.-")
        };

        proteins.push(&protein1);
        proteins.push(&protein2);

        let function_aggregator = FunctionAggregator {};

        let result = function_aggregator.get_all_functional_annotations(proteins.as_slice());

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].len(), 2);
        assert_eq!(result[1].len(), 1);
    }

    #[test]
    fn test_serialize_functional_aggregation() {
        let mut proteins: Vec<Protein> = Vec::new();
        proteins.push(Protein {
            uniprot_id:             "P12345".to_string(),
            taxon_id:               9606,
            functional_annotations: encode("GO:0001234;GO:0005678")
        });
        proteins.push(Protein {
            uniprot_id:             "P23456".to_string(),
            taxon_id:               9606,
            functional_annotations: encode("EC:1.1.1.-")
        });

        let function_aggregator = FunctionAggregator {};

        let result = function_aggregator.aggregate(proteins.iter().collect());

        let generated_json = serde_json::to_string(&result).unwrap();
        let expected_json = "{\"counts\":{\"all\":2,\"GO\":1,\"EC\":1,\"IPR\":0},\"data\":{\"GO:0001234\":1,\"GO:0005678\":1,\"EC:1.1.1.-\":1}}";

        assert_eq!(
            generated_json.parse::<serde_json::Value>().unwrap(),
            expected_json.parse::<serde_json::Value>().unwrap(),
        );
    }
}
