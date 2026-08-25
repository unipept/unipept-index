//! Protein metadata: accessions, taxon ids and functional annotations, addressed by index.
//!
//! The suffix array resolves a peptide to positions in the concatenated protein text; the
//! suffix-to-protein mapping turns those positions into protein indices; and this crate turns an
//! index into the metadata a search result actually reports.
//!
//! Like the rest of the index it has two backends — one holding owned memory, one borrowing a
//! memory mapping. Both are always compiled; a caller chooses by naming one of the two structs.
//! Both hand out [`proteins::ProteinRef`], which borrows rather than copies, so the two differ
//! only in where the bytes live.
//!
//! Each is in turn generic over its text backend, so metadata storage and text storage are two
//! independent axes rather than one switch; see [`proteins`].
#![warn(missing_docs)]

pub mod proteins;
