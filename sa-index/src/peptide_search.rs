//! The public entry point: peptides in, proteins out.
//!
//! Wraps the searcher with the parts a caller actually wants — normalising the query, dropping
//! peptides too short to be searchable, resolving suffixes to proteins, and decoding annotations.
//! `sa-server` calls [`search_all_peptides_json`], which hands back the response body already
//! serialised.
//!
//! # Why the results borrow, and why serialisation happens here
//!
//! A result set is large — a 10,000-peptide request at the default cutoff can reach millions of
//! protein hits and hundreds of megabytes of JSON — so the two obvious conveniences are both
//! expensive at that scale:
//!
//! * Owning the accession and the decoded annotations costs **two allocations per hit**, purely so
//!   `serde` has a `&str` to copy out of a moment later. [`ProteinInfo`] therefore borrows the
//!   accession from the index and keeps the annotations *encoded*, decoding them straight into the
//!   serialiser's output buffer via [`fa_compression::algorithm1::decoded`].
//! * Serialising the whole `Vec<SearchResult>` in one go is single-threaded, and it was the largest
//!   serial stretch of a request. [`search_all_peptides_json`] instead serialises each peptide's
//!   result on the rayon worker that already retrieved it, and the caller only has to frame the
//!   chunks.
//!
//! Borrowing has a consequence worth stating: [`SearchResult`] holds a reference into the
//! [`Searcher`], so it can never satisfy the `'static` bound an HTTP body needs. Serialising inside
//! the request handler is not merely an optimisation here — it is what the lifetimes require.

use rayon::prelude::*;
use sa_mappings::proteins::{ProteinRef, ProteinsBackend};
use serde::{Serialize, Serializer};

use crate::{
    array::SuffixArrayBackend,
    sa_searcher::{SearchAllSuffixesResult, Searcher},
    suffix_to_protein_index::SuffixToProteinMappingBackend
};

/// Everything found for one peptide. Serialised straight to JSON by the server.
#[derive(Debug, Serialize)]
pub struct SearchResult<'a> {
    /// The peptide as the caller wrote it, before normalisation.
    pub sequence: &'a str,
    /// Every protein containing the peptide.
    pub proteins: Vec<ProteinInfo<'a>>,
    /// Whether the match cutoff was hit, meaning `proteins` is a truncated sample rather than the
    /// complete set.
    pub cutoff_used: bool
}

/// One matching protein, borrowed from the index.
///
/// Nothing here is owned and nothing is decoded until serialisation: see the module docs for why.
#[derive(Debug, Serialize)]
pub struct ProteinInfo<'a> {
    /// NCBI taxon id.
    pub taxon: u32,
    /// UniProt accession, e.g. `P12345`.
    pub uniprot_accession: &'a str,
    /// Functional annotations, **still `fa-compression`-encoded**.
    ///
    /// Serialised as the decoded `GO:`/`EC:`/`IPR:` text under the key `functional_annotations`,
    /// which is why the field is not called that: it holds bytes that are not yet that string.
    #[serde(rename = "functional_annotations", serialize_with = "serialize_annotations")]
    pub annotations: &'a [u8]
}

/// Writes the annotations as their decoded text, without building that text first.
///
/// `serde_json` overrides `collect_str` to stream a `Display` value through its string-escaping
/// writer, so the decode lands directly in the output buffer. Other serialisers fall back to
/// serde's default, which materialises a `String` — correct, just not free.
fn serialize_annotations<S: Serializer>(annotations: &&[u8], serializer: S) -> Result<S::Ok, S::Error> {
    serializer.collect_str(&fa_compression::algorithm1::decoded(annotations))
}

impl<'a> From<ProteinRef<'a>> for ProteinInfo<'a> {
    fn from(protein: ProteinRef<'a>) -> Self {
        ProteinInfo {
            taxon: protein.taxon_id,
            uniprot_accession: protein.uniprot_id,
            annotations: protein.functional_annotations
        }
    }
}

/// Searches the `peptide` in the index multithreaded and retrieves the matching proteins
///
/// # Arguments
/// * `searcher` - The Searcher which contains the protein database
/// * `peptide` - The peptide that is being searched in the index
/// * `cutoff` - The maximum amount of matches we want to process from the index
/// * `equate_il` - Boolean indicating if we want to equate I and L during search
/// * `tryptic` - Boolean indicating if we only want tryptic matches.
///
/// # Returns
///
/// Returns Some if matches are found.
/// The first argument is true if the cutoff is used, otherwise false
/// The second argument is a list of all matching proteins for the peptide
/// Returns None if the peptides does not have any matches, or if the peptide is shorter than the
/// sparseness factor k used in the index
pub fn search_proteins_for_peptide<
    'a,
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend
>(
    searcher: &'a Searcher<SA, P, STPM>,
    peptide: &str,
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Option<(bool, Vec<ProteinRef<'a>>)> {
    let peptide = peptide.trim_end().to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sa.sample_rate() as usize {
        return None;
    }

    let suffix_search = searcher.search_matching_suffixes_scalar(peptide.as_bytes(), cutoff, equate_il, tryptic);
    let (suffixes, cutoff_used) = match suffix_search {
        SearchAllSuffixesResult::MaxMatches(matched_suffixes) => Some((matched_suffixes, true)),
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => Some((matched_suffixes, false)),
        SearchAllSuffixesResult::NoMatches => None
    }?;

    let proteins = searcher.retrieve_proteins(&suffixes);

    Some((cutoff_used, proteins))
}

/// Searches one peptide and packages the result.
///
/// The single-peptide path. `sa-server` uses [`search_all_peptides`] instead, which batches the
/// suffix searches across the whole request; this one is kept for callers with a single query.
///
/// Returns `None` when the peptide has no matches or is shorter than the sparseness factor.
pub fn search_peptide<'a, SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &'a Searcher<SA, P, STPM>,
    peptide: &'a str,
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Option<SearchResult<'a>> {
    let (cutoff_used, proteins) = search_proteins_for_peptide(searcher, peptide, cutoff, equate_il, tryptic)?;

    Some(SearchResult {
        sequence: peptide,
        proteins: proteins.into_iter().map(|protein| protein.into()).collect(),
        cutoff_used
    })
}

/// Searches the list of `peptides` in the index and retrieves all related information about the
/// found proteins This does NOT perform any of the analyses
///
/// # Arguments
/// * `searcher` - The Searcher which contains the protein database
/// * `peptides` - List of peptides we want to search in the index
/// * `cutoff` - The maximum amount of matches we want to process from the index
/// * `equate_il` - Boolean indicating if we want to equate I and L during search
/// * `tryptic` - Boolean indicating if we only want tryptic matches.
///
/// # Returns
///
/// Returns an `OutputData<SearchOnlyResult>` object with the search results for the peptides
pub fn search_all_peptides<'a, SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &'a Searcher<SA, P, STPM>,
    peptides: &'a [String],
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Vec<SearchResult<'a>> {
    search_all_peptides_with(searcher, peptides, cutoff, equate_il, tryptic, |result| result)
}

/// The shared pipeline behind [`search_all_peptides`] and [`search_all_peptides_json`]: search,
/// retrieve, and hand each peptide's [`SearchResult`] to `f` on the rayon worker that built it.
///
/// Existing so the two entry points cannot drift — whatever the server does to a result, the
/// benchmark harness can do to the same one.
fn search_all_peptides_with<'a, SA, P, STPM, T, F>(
    searcher: &'a Searcher<SA, P, STPM>,
    peptides: &'a [String],
    cutoff: usize,
    equate_il: bool,
    tryptic: bool,
    f: F
) -> Vec<T>
where
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend,
    T: Send,
    F: Fn(SearchResult<'a>) -> T + Sync + Send
{
    let sample_rate = searcher.sa.sample_rate() as usize;

    // Normalise once and keep only searchable peptides (anything shorter than the sample rate
    // cannot be searched), remembering each peptide's original index so the returned `sequence`
    // stays verbatim — matching the previous per-peptide behaviour.
    let prepared: Vec<(usize, String)> = peptides
        .iter()
        .enumerate()
        .map(|(index, peptide)| (index, peptide.trim_end().to_uppercase()))
        .filter(|(_, normalised)| normalised.len() >= sample_rate)
        .collect();

    let byte_peptides: Vec<&[u8]> = prepared.iter().map(|(_, peptide)| peptide.as_bytes()).collect();

    // One batched (MLP) search over the whole list — the same code path the benchmark measures
    // and the only place the batch size is chosen.
    let suffix_results = searcher.search_all_matching_suffixes(&byte_peptides, cutoff, equate_il, tryptic);

    // Retrieve the proteins for each hit and build the result, dropping peptides with no matches
    // (preserves the previous filter_map semantics and result ordering).
    //
    // Retrieval is per peptide by design: a cross-query batched variant was built and measured
    // (run3) and moved throughput by a median of +1.7%, never clearing the noise floor. See the
    // module comment on `sa_searcher::batched` for the full result.
    prepared
        .par_iter()
        .zip(suffix_results.par_iter())
        .filter_map(|((original_index, _), suffix_result)| {
            let (suffixes, cutoff_used) = match suffix_result {
                SearchAllSuffixesResult::MaxMatches(matched) => (matched, true),
                SearchAllSuffixesResult::SearchResult(matched) => (matched, false),
                SearchAllSuffixesResult::NoMatches => return None
            };

            let proteins = searcher.retrieve_proteins(suffixes);

            Some(f(SearchResult {
                sequence: &peptides[*original_index],
                proteins: proteins.into_iter().map(|protein| protein.into()).collect(),
                cutoff_used
            }))
        })
        .collect()
}

/// Searches the list of `peptides` and returns the response body already serialised, one chunk per
/// peptide with matches.
///
/// The server path. Each chunk is one `SearchResult` as JSON, prefixed with a `,` so
/// [`frame_chunks`] can turn the collection into a JSON array by overwriting the first byte —
/// see there. Serialising per peptide is what keeps this phase parallel; a single
/// `serde_json::to_vec` over the whole answer would put hundreds of megabytes of JSON writing on
/// one thread.
///
/// `chunks.concat()` is byte-for-byte what `serde_json::to_vec(&search_all_peptides(..))` produces.
pub fn search_all_peptides_json<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    peptides: &[String],
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Vec<Vec<u8>> {
    let mut chunks =
        search_all_peptides_with(searcher, peptides, cutoff, equate_il, tryptic, |result| json_chunk(&result));
    frame_chunks(&mut chunks);
    chunks
}

/// Serialises one result, prefixed with the `,` that separates it from the previous element.
///
/// The buffer is sized up front because it can reach megabytes for a single peptide at the default
/// cutoff, and `serde_json`'s own writer starts at 128 bytes and doubles.
pub fn json_chunk(result: &SearchResult<'_>) -> Vec<u8> {
    let estimate = 64
        + result.sequence.len()
        + result
            .proteins
            .iter()
            .map(|protein| 64 + protein.uniprot_accession.len() + protein.annotations.len() * 3)
            .sum::<usize>()
        // Room for the closing `]` that `frame_chunks` pushes onto the last chunk.
        + 1;
    let mut chunk = Vec::with_capacity(estimate);
    chunk.push(b',');
    // Infallible: the writer is a `Vec` and every value in the tree serialises.
    serde_json::to_writer(&mut chunk, result).expect("serialising a SearchResult into a Vec cannot fail");
    chunk
}

/// Turns [`json_chunk`] output into a JSON array, in place.
///
/// Each chunk already starts with the `,` that precedes it, so the array is finished by overwriting
/// the very first one with `[` and appending `]`. That is exactly `serde_json`'s own compact
/// framing — `[` , element , `,` , element , `]`, no whitespace — and an empty collection becomes
/// `[]`.
pub fn frame_chunks(chunks: &mut Vec<Vec<u8>>) {
    match chunks.first_mut() {
        None => chunks.push(b"[]".to_vec()),
        Some(first) => {
            first[0] = b'[';
            chunks.last_mut().expect("non-empty").push(b']');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // `ProteinsBackend` (for `proteins.text()`) comes in through `use super::*`.
    use sa_mappings::proteins::{InMemoryProteins, Protein};
    use text_compression::InMemoryProteinText;

    use crate::{
        array::{InMemorySA, OriginalSA},
        suffix_to_protein_index::{BitVecSuffixToProtein, InMemorySuffixToProteinMapping}
    };

    type TestSearcher = Searcher<InMemorySA, InMemoryProteins<InMemoryProteinText>, InMemorySuffixToProteinMapping>;

    /// Built from the owned types directly rather than through `sa_searcher::test_utils`: these
    /// test the peptide-level API, not the storage backends, and the backends are already proved
    /// interchangeable by `sa_searcher::backend_agreement`.
    ///
    /// `annotations` decides whether the proteins carry real `fa-compression` payloads. Empty ones
    /// keep the fixture readable, but they never exercise the decode path — which is most of what
    /// the response phase does — so the JSON tests run against both.
    fn test_searcher(annotations: bool) -> TestSearcher {
        // Example DB "AI-CLACVAA-AC-KCRLY$", sample rate 1 (so single characters are searchable).
        let text = InMemoryProteinText::from_string("AI-CLACVAA-AC-KCRLY$");
        let encoded = |text: &str| {
            if annotations {
                fa_compression::algorithm1::encode(text)
            } else {
                vec![]
            }
        };
        let proteins = InMemoryProteins::new(text, vec![
            Protein {
                uniprot_id: "P0".to_string(),
                taxon_id: 10,
                functional_annotations: encoded("EC:1.1.1.-;GO:0009279;IPR:IPR016364")
            },
            Protein {
                uniprot_id: "P1".to_string(),
                taxon_id: 20,
                functional_annotations: encoded("GO:0046782")
            },
            Protein {
                uniprot_id: "P2".to_string(),
                taxon_id: 30,
                functional_annotations: encoded("")
            },
            Protein {
                uniprot_id: "P3".to_string(),
                taxon_id: 40,
                functional_annotations: encoded("IPR:IPR008816;IPR:IPR032635")
            },
        ]);
        let sa = InMemorySA::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18],
            1
        ));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        Searcher::new(sa, proteins, InMemorySuffixToProteinMapping::BitVec(stp))
    }

    /// Mixed case (normalisation), whitespace + empty (trim/length filter), a no-match ("ZZZ").
    fn test_peptides() -> Vec<String> {
        ["A", "ai", "CLA", "kcrly", "VAA", "ZZZ", "", "  ", "AC"].iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_search_all_peptides_matches_scalar_reference() {
        let searcher = test_searcher(true);
        let peptides = test_peptides();

        let as_json = |v: &[SearchResult]| serde_json::to_string(v).unwrap();
        for equate_il in [true, false] {
            for tryptic in [true, false] {
                // Batched production path vs the scalar per-peptide reference (the old behaviour).
                let got = search_all_peptides(&searcher, &peptides, 1000, equate_il, tryptic);
                let reference: Vec<SearchResult> =
                    peptides.iter().filter_map(|p| search_peptide(&searcher, p, 1000, equate_il, tryptic)).collect();
                assert_eq!(as_json(&got), as_json(&reference), "il={} tryptic={}", equate_il, tryptic);
            }
        }
    }

    /// The claim the server depends on: the parallel per-peptide chunks, concatenated, are
    /// byte-for-byte what one `serde_json::to_vec` over the whole answer produces.
    #[test]
    fn json_chunks_concatenate_to_the_same_bytes_as_one_pass() {
        for annotations in [true, false] {
            let searcher = test_searcher(annotations);
            let peptides = test_peptides();

            for equate_il in [true, false] {
                for tryptic in [true, false] {
                    let chunks = search_all_peptides_json(&searcher, &peptides, 1000, equate_il, tryptic);
                    let reference =
                        serde_json::to_vec(&search_all_peptides(&searcher, &peptides, 1000, equate_il, tryptic))
                            .unwrap();
                    assert_eq!(
                        chunks.concat(),
                        reference,
                        "annotations={annotations} il={equate_il} tryptic={tryptic}"
                    );
                }
            }
        }
    }

    /// Nothing matches, so there are no chunks to frame — the array still has to be `[]` rather
    /// than empty, which is the one case the leading-comma trick has to special-case.
    #[test]
    fn json_chunks_of_an_empty_result_are_an_empty_array() {
        let searcher = test_searcher(true);

        let no_matches = vec!["ZZZ".to_string()];
        assert_eq!(search_all_peptides_json(&searcher, &no_matches, 1000, false, false).concat(), b"[]");

        let no_peptides: Vec<String> = vec![];
        assert_eq!(search_all_peptides_json(&searcher, &no_peptides, 1000, false, false).concat(), b"[]");
    }

    /// One match is the other edge of the framing: a single chunk is both the first (its comma
    /// becomes `[`) and the last (it gets the `]`).
    #[test]
    fn json_chunks_of_a_single_result_are_framed_at_both_ends() {
        let searcher = test_searcher(true);

        let one = vec!["VAA".to_string()];
        let chunks = search_all_peptides_json(&searcher, &one, 1000, false, false);
        let body = chunks.concat();
        assert_eq!(body.first(), Some(&b'['));
        assert_eq!(body.last(), Some(&b']'));
        assert_eq!(body, serde_json::to_vec(&search_all_peptides(&searcher, &one, 1000, false, false)).unwrap());
    }

    fn assert_json_eq(generated_json: &str, expected_json: &str) {
        assert_eq!(
            generated_json.parse::<serde_json::Value>().unwrap(),
            expected_json.parse::<serde_json::Value>().unwrap(),
        );
    }

    #[test]
    fn test_serialize_protein_info() {
        let annotations = fa_compression::algorithm1::encode("GO:0001234;GO:0005678");
        let protein_info = ProteinInfo {
            taxon: 1,
            uniprot_accession: "P12345",
            annotations: &annotations
        };

        let generated_json = serde_json::to_string(&protein_info).unwrap();
        let expected_json =
            "{\"taxon\":1,\"uniprot_accession\":\"P12345\",\"functional_annotations\":\"GO:0001234;GO:0005678\"}";

        assert_json_eq(&generated_json, expected_json);
    }

    #[test]
    fn test_serialize_search_result() {
        let search_result = SearchResult {
            sequence: "MSKIAALLPSV",
            proteins: vec![],
            cutoff_used: true
        };

        let generated_json = serde_json::to_string(&search_result).unwrap();
        let expected_json = "{\"sequence\":\"MSKIAALLPSV\",\"proteins\":[],\"cutoff_used\":true}";

        assert_json_eq(&generated_json, expected_json);
    }

    /// `serialize_annotations` streams through `collect_str` rather than handing serde a `&str`,
    /// so the escaping is `serde_json`'s `collect_str` path rather than its `serialize_str` one.
    /// The two are supposed to be identical; this pins that, and pins that the borrowed accession
    /// and sequence still go through ordinary string escaping.
    #[test]
    fn borrowed_fields_escape_exactly_like_owned_ones() {
        /// What the struct looked like before it borrowed: plain derived `Serialize` over owned
        /// `String`s, which is the behaviour that must not have changed.
        #[derive(Serialize)]
        struct OwnedProteinInfo {
            taxon: u32,
            uniprot_accession: String,
            functional_annotations: String
        }

        #[derive(Serialize)]
        struct OwnedSearchResult {
            sequence: String,
            proteins: Vec<OwnedProteinInfo>,
            cutoff_used: bool
        }

        // Nothing `fa-compression` emits needs escaping, but the accession and the sequence are
        // arbitrary text from the index and the caller.
        let awkward = ["", "P12345", "quote\"inside", "back\\slash", "new\nline", "nul\u{0}byte", "ünïcøde", "🧬"];

        for accession in awkward {
            for sequence in awkward {
                let annotations = fa_compression::algorithm1::encode("EC:1.1.1.-;GO:0009279;IPR:IPR016364");
                let decoded = fa_compression::algorithm1::decode(&annotations);

                let borrowed = SearchResult {
                    sequence,
                    proteins: vec![ProteinInfo {
                        taxon: 7,
                        uniprot_accession: accession,
                        annotations: &annotations
                    }],
                    cutoff_used: false
                };
                let owned = OwnedSearchResult {
                    sequence: sequence.to_string(),
                    proteins: vec![OwnedProteinInfo {
                        taxon: 7,
                        uniprot_accession: accession.to_string(),
                        functional_annotations: decoded
                    }],
                    cutoff_used: false
                };

                assert_eq!(
                    serde_json::to_string(&borrowed).unwrap(),
                    serde_json::to_string(&owned).unwrap(),
                    "accession={accession:?} sequence={sequence:?}"
                );
            }
        }
    }

    /// Empty annotations are the common case in a real index and the one where an off-by-one in the
    /// framing would be invisible in the tests above.
    #[test]
    fn empty_annotations_serialise_as_an_empty_string() {
        let protein_info = ProteinInfo {
            taxon: 1,
            uniprot_accession: "P12345",
            annotations: &[]
        };
        assert_eq!(
            serde_json::to_string(&protein_info).unwrap(),
            "{\"taxon\":1,\"uniprot_accession\":\"P12345\",\"functional_annotations\":\"\"}"
        );
    }
}
