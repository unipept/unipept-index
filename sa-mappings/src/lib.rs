#![warn(missing_docs)]

//! Protein metadata: accessions, taxon ids and functional annotations, addressed by index.
//!
//! The suffix array resolves a peptide to positions in the concatenated protein text; the
//! suffix-to-protein mapping turns those positions into protein indices; and this crate turns an
//! index into the metadata a search result actually reports.
//!
//! Like the rest of the index it has two backends, selected by the `mmap` feature and resolved
//! through the [`proteins::Proteins`] alias. Both hand out [`proteins::ProteinRef`], which
//! borrows rather than copies, so the two differ only in where the bytes live.

pub mod proteins;
