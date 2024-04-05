use std::error::Error;

mod proteins;
mod taxonomy;

#[derive(Debug)]
struct DatabaseFormatError {
    error: Vec<String>
}

impl DatabaseFormatError {
    fn new(error: Vec<String>) -> Self {
        Self { error }
    }
}

impl std::fmt::Display for DatabaseFormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Expected the protein database file to have the following fields separated by a tab: <Uniprot_accession> <protein id> <sequence>\nBut tried to unpack following vector in 3 variables: {:?}", self.error)
    }
}

impl Error for DatabaseFormatError {}
