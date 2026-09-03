//! The public entry point: peptides in, proteins out.
//!
//! Wraps the searcher with the parts a caller actually wants — normalising the query, dropping
//! peptides too short to be searchable, resolving suffixes to proteins, and decoding annotations.
//! `sa-server` calls [`search_all_peptides_json`], which hands back the response body already
//! serialised — and so does anything else building an endpoint on this crate, which is how the
//! search actually reaches production. `sa-server` is a testing and direct-serving tool, not the
//! path a deployed request takes.
//!
//! # The taxa-only path
//!
//! [`search_all_peptides_taxa`] and [`search_peptide_taxa`] answer the same query with taxon ids
//! alone. `pept2taxa` and `pept2lca` want nothing else, and the protein path makes them pay for an
//! accession and an encoded annotation blob per hit — plus the decode below — before they throw
//! all of it away. Both paths run the same search over the same suffixes and drop the same
//! peptides; they differ only in what is retrieved once the suffixes are known.
//!
//! # Resource limits are the caller's, and there are none here
//!
//! Nothing in this module bounds the work a request can ask for. A caller exposing these functions
//! to untrusted input owns that, and two things multiply:
//!
//! * **The peptide list.** Every entry costs an independent search plus per-hit protein retrieval,
//!   and short peptides are the expensive ones — a single residue matches a large fraction of the
//!   index.
//! * **`cutoff`.** It is *not* an upper bound on work in the direction that matters. It caps a
//!   result set only while the match range is larger than it; once the range is smaller, the whole
//!   range is collected regardless (see `sa_searcher::scalar`). A very large `cutoff` therefore
//!   means "return everything", not "return at most this many".
//!
//! So the cost of one request is roughly the sum over peptides of that peptide's match count, each
//! match carrying a random access into the protein store and an annotation decode. Bound the
//! peptide count and clamp `cutoff` at the boundary that accepts the request; neither is done here,
//! and the alphabet and length filters below are correctness filters, not limits.
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

use protein_metadata::{ProteinRef, ProteinsBackend};
use rayon::prelude::*;
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

/// Every taxon found for one peptide, sorted and deduplicated.
///
/// The lightweight answer: `pept2taxa` and `pept2lca` want taxon ids and nothing else, and
/// producing a [`SearchResult`] for them means retrieving an accession and an annotation blob per
/// hit that are then discarded. Sorting and deduplicating here rather than in the caller is what
/// makes the difference worth having — a peptide matching thousands of suffixes usually spans far
/// fewer taxa.
///
/// Borrows `sequence` for the same reason [`SearchResult`] does; see the module docs.
#[derive(Debug, Serialize)]
pub struct TaxaSearchResult<'a> {
    /// The peptide as the caller wrote it, before normalisation.
    pub sequence: &'a str,
    /// Every distinct taxon id among the matching proteins, ascending.
    pub taxa: Vec<u32>,
    /// Whether the match cutoff was hit, meaning `taxa` is drawn from a truncated sample of the
    /// matches rather than all of them.
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

/// Normalises one peptide for search, or rejects it.
///
/// Uppercases ASCII in the same pass that validates, and returns `None` for any byte the index
/// cannot contain. Three groups are rejected:
///
/// * `-` and `$` — structural characters in the text (the protein separator and the terminator),
///   never part of a sequence. A query containing either could otherwise match across a protein
///   boundary, and a trailing `$` is what drives a match to `text.len()`.
/// * `J` — absent from the index alphabet (`BIT5_TO_CHAR`), so it can never match.
/// * everything else — digits, punctuation, and every non-ASCII byte. `to_uppercase` is
///   Unicode-aware, so a character like `é` survives an uppercase-only normalisation as multi-byte
///   UTF-8 whose every byte is >= 128.
///
/// Used by both entry points in place of `to_uppercase()`: it walks the same bytes and makes the
/// same single allocation, without the Unicode table lookups.
///
/// Note this does *not* make the checks in `KmerTable::lookup` or `check_tryptic_c_term`
/// redundant. `search_matching_suffixes_scalar` and `search_all_matching_suffixes_batched` are
/// public and are called directly (by `sa-benchmarks`, among others), so those paths never pass
/// through here.
fn normalise_peptide(peptide: &str) -> Option<String> {
    let trimmed = peptide.trim_end();
    let mut out = String::with_capacity(trimmed.len());
    for &b in trimmed.as_bytes() {
        let upper = b.to_ascii_uppercase();
        if !upper.is_ascii_uppercase() || upper == b'J' {
            return None;
        }
        out.push(upper as char);
    }
    Some(out)
}

/// Searches one `peptide` in the index and retrieves the matching proteins.
///
/// Single-threaded: this is the scalar path, and it is what [`search_peptide`] is built on.
/// [`search_all_peptides`] is the batched, rayon-parallel path for a whole list.
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
/// Returns None if the peptide does not have any matches, if it is shorter than the sparseness
/// factor k used in the index, or if it contains a character the index cannot hold: the
/// structural `-` and `$`, `J` (absent from the index alphabet), or anything else outside the
/// ASCII amino-acid letters. ASCII case is normalised first, so lowercase input is fine.
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
    let (cutoff_used, suffixes) = search_suffixes_for_peptide(searcher, peptide, cutoff, equate_il, tryptic)?;

    let proteins = searcher.retrieve_proteins(&suffixes);

    Some((cutoff_used, proteins))
}

/// Normalises one peptide and runs the scalar suffix search, returning the matched suffixes and
/// whether the cutoff truncated them.
///
/// The half of the single-peptide path that does not depend on what is retrieved afterwards, so
/// the protein and taxa entry points share it and cannot disagree about which peptides are
/// searchable. The batched path has its own equivalent inside [`search_all_suffixes_with`].
///
/// `None` means the peptide is unsearchable — rejected by [`normalise_peptide`], or shorter than
/// the sparseness factor — or that it simply matched nothing. Callers cannot tell those apart, and
/// none of them need to.
fn search_suffixes_for_peptide<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    peptide: &str,
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Option<(bool, Vec<i64>)> {
    // Rejects anything outside the index alphabet before it reaches the searcher.
    let peptide = normalise_peptide(peptide)?;

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sa.sample_rate() as usize {
        return None;
    }

    let suffix_search = searcher.search_matching_suffixes_scalar(peptide.as_bytes(), cutoff, equate_il, tryptic);
    match suffix_search {
        SearchAllSuffixesResult::MaxMatches(matched_suffixes) => Some((true, matched_suffixes)),
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => Some((false, matched_suffixes)),
        SearchAllSuffixesResult::NoMatches => None
    }
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
/// Returns one [`SearchResult`] per peptide that matched, in input order. Peptides that match
/// nothing, are shorter than the sparseness factor, or fall outside the index alphabet are dropped
/// rather than represented by an empty entry.
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
    search_all_suffixes_with(searcher, peptides, cutoff, equate_il, tryptic, |sequence, suffixes, cutoff_used| {
        let proteins = searcher.retrieve_proteins(suffixes);

        f(SearchResult {
            sequence,
            proteins: proteins.into_iter().map(|protein| protein.into()).collect(),
            cutoff_used
        })
    })
}

/// The batched search itself, with retrieval left to the caller: normalise the list, run one
/// batched suffix search over all of it, and hand each matching peptide's suffixes to `f` on the
/// rayon worker that produced them.
///
/// Sits under [`search_all_peptides_with`] and [`search_all_peptides_taxa`] so that batching is
/// decided in exactly one place. What differs between the protein and taxa paths is only what `f`
/// does with the suffixes; nothing above this function chooses a batch size, normalises a peptide,
/// or decides which peptides are searchable.
fn search_all_suffixes_with<'a, SA, P, STPM, T, F>(
    searcher: &Searcher<SA, P, STPM>,
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
    F: Fn(&'a str, &[i64], bool) -> T + Sync + Send
{
    let sample_rate = searcher.sa.sample_rate() as usize;

    // Normalise once and keep only searchable peptides — anything shorter than the sample rate
    // cannot be searched, and anything outside the index alphabet cannot match — remembering each
    // peptide's original index so the returned `sequence` stays verbatim.
    let prepared: Vec<(usize, String)> = peptides
        .iter()
        .enumerate()
        .filter_map(|(index, peptide)| Some((index, normalise_peptide(peptide)?)))
        .filter(|(_, normalised)| normalised.len() >= sample_rate)
        .collect();

    let byte_peptides: Vec<&[u8]> = prepared.iter().map(|(_, peptide)| peptide.as_bytes()).collect();

    // One batched (MLP) search over the whole list — the same code path the benchmark measures
    // and the only place the batch size is chosen.
    let suffix_results = searcher.search_all_matching_suffixes_batched(&byte_peptides, cutoff, equate_il, tryptic);

    // Retrieve for each hit and build the result, dropping peptides with no matches (preserves the
    // previous filter_map semantics and result ordering).
    //
    // Retrieval is per peptide by design: a cross-query batched variant was built and measured,
    // and moved throughput by a median well inside the noise floor. See the module comment on
    // `sa_searcher::batched` for the full result.
    prepared
        .par_iter()
        .zip(suffix_results.par_iter())
        .filter_map(|((original_index, _), suffix_result)| {
            let (suffixes, cutoff_used) = match suffix_result {
                SearchAllSuffixesResult::MaxMatches(matched) => (matched, true),
                SearchAllSuffixesResult::SearchResult(matched) => (matched, false),
                SearchAllSuffixesResult::NoMatches => return None
            };

            Some(f(&peptides[*original_index], suffixes, cutoff_used))
        })
        .collect()
}

/// Searches one peptide and returns only its taxa.
///
/// The taxa counterpart of [`search_peptide`], and the single-peptide sibling of
/// [`search_all_peptides_taxa`]. Nothing is retrieved that a taxon id does not need: no accession,
/// no annotation blob, and so no decode at serialisation time either.
///
/// Returns `None` on the same terms as [`search_peptide`] — no matches, or a peptide that is
/// unsearchable because it is too short or outside the index alphabet.
pub fn search_peptide_taxa<'a, SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    peptide: &'a str,
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Option<TaxaSearchResult<'a>> {
    let (cutoff_used, suffixes) = search_suffixes_for_peptide(searcher, peptide, cutoff, equate_il, tryptic)?;

    Some(TaxaSearchResult {
        sequence: peptide,
        taxa: sorted_taxa(searcher, &suffixes),
        cutoff_used
    })
}

/// Searches the list of `peptides` and returns only their taxa.
///
/// The taxa counterpart of [`search_all_peptides`], down the same batched path — both are built
/// on the private `search_all_suffixes_with`. For `pept2taxa` and `pept2lca` this is
/// the entry point to reach for: a [`SearchResult`] would carry an accession and an encoded
/// annotation blob per hit that those consumers immediately discard.
///
/// Returns one [`TaxaSearchResult`] per peptide that matched, in input order; peptides that match
/// nothing are omitted, exactly as in [`search_all_peptides`].
pub fn search_all_peptides_taxa<'a, SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    peptides: &'a [String],
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Vec<TaxaSearchResult<'a>> {
    search_all_suffixes_with(searcher, peptides, cutoff, equate_il, tryptic, |sequence, suffixes, cutoff_used| {
        TaxaSearchResult { sequence, taxa: sorted_taxa(searcher, suffixes), cutoff_used }
    })
}

/// Retrieves the taxa for `suffixes` and reduces them to the ascending distinct set.
///
/// `sort_unstable` because the values are `u32` with no meaning attached to equal elements, and
/// because the input is far longer than the output: a peptide matching thousands of suffixes
/// typically spans orders of magnitude fewer taxa.
fn sorted_taxa<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    suffixes: &[i64]
) -> Vec<u32> {
    let mut taxa = searcher.retrieve_taxa(suffixes);
    taxa.sort_unstable();
    taxa.dedup();
    taxa
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
///
/// # Panics
///
/// Every chunk must be [`json_chunk`] output — non-empty and comma-prefixed. Both functions are
/// `pub`, so the precondition is stated and checked rather than assumed: unchecked, a violation is
/// silent rather than loud — a chunk of `b"x"` comes back as `"[]"`, the caller's byte overwritten
/// and the result quietly wrong, and an empty first chunk panics with an index-out-of-bounds and no
/// explanation. Checking turns both into one clear message; this is a programming error on the
/// caller's side, not untrusted input, so it is an assertion rather than a `Result`.
pub fn frame_chunks(chunks: &mut Vec<Vec<u8>>) {
    match chunks.first_mut() {
        None => chunks.push(b"[]".to_vec()),
        Some(first) => {
            assert_eq!(
                first.first(),
                Some(&b','),
                "frame_chunks expects json_chunk output: every chunk must be non-empty and start with ','"
            );
            first[0] = b'[';
            chunks.last_mut().expect("non-empty").push(b']');
        }
    }
}

#[cfg(test)]
mod tests {
    // `ProteinsBackend` (for `proteins.text()`) comes in through `use super::*`.
    use protein_metadata::{InMemoryProteins, Protein};
    use protein_text::InMemoryProteinText;

    use super::*;
    use crate::{
        array::{InMemorySA, OriginalSA},
        suffix_to_protein_index::{BitVecSuffixToProtein, InMemorySuffixToProteinMapping}
    };

    type TestSearcher = Searcher<InMemorySA, InMemoryProteins<InMemoryProteinText>, InMemorySuffixToProteinMapping>;

    /// Built from the owned types directly rather than through `sa_searcher::test_utils`: these
    /// test the peptide-level API, not the storage backends, and the backends are already proved
    /// interchangeable by `sa_searcher::tests::every_backend_combination_returns_identical_results`.
    ///
    /// `annotations` decides whether the proteins carry real `fa-compression` payloads. Empty ones
    /// keep the fixture readable, but they never exercise the decode path — which is most of what
    /// the response phase does — so the JSON tests run against both.
    fn test_searcher(annotations: bool) -> TestSearcher {
        // Example DB "AI-CLACVAA-AC-KCRLY$", sample rate 1 (so single characters are searchable).
        let text = InMemoryProteinText::from_string("AI-CLACVAA-AC-KCRLY$");
        let encoded = |text: &str| {
            if annotations { fa_compression::algorithm1::encode(text) } else { vec![] }
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

    /// Every byte the index cannot hold is rejected before it reaches the searcher — including the
    /// two that would otherwise reach an index guarded nowhere else: a non-ASCII character, whose
    /// UTF-8 bytes are all at least 128 (the `KmerTable::lookup` case), and a trailing `$`, which
    /// drives a match to `text.len()` (the `check_tryptic_c_term` case).
    #[test]
    fn test_normalise_peptide_rejects_everything_outside_the_alphabet() {
        // Accepted, and uppercased in the same pass.
        assert_eq!(normalise_peptide("aik").as_deref(), Some("AIK"));
        assert_eq!(normalise_peptide("CLA  ").as_deref(), Some("CLA"));
        assert_eq!(normalise_peptide("").as_deref(), Some(""));

        // Structural characters: never part of a sequence.
        assert!(normalise_peptide("AC$").is_none(), "terminator");
        assert!(normalise_peptide("AC-AC").is_none(), "separator");

        // Absent from BIT5_TO_CHAR, so it can never match.
        assert!(normalise_peptide("AJK").is_none(), "J is not in the index alphabet");

        // Non-ASCII: 'é' is 0xC3 0xA9, both >= 128.
        assert!(normalise_peptide("ACé").is_none(), "non-ASCII");
        assert!(normalise_peptide("Ω").is_none());

        // Everything else.
        assert!(normalise_peptide("AC1").is_none(), "digit");
        assert!(normalise_peptide("AC*").is_none(), "punctuation");
        assert!(normalise_peptide("A C").is_none(), "interior whitespace");
    }

    /// A rejected peptide is dropped, not searched — the same treatment the length filter already
    /// gives a too-short one, so it simply does not appear in the results.
    ///
    /// The cases here are chosen so the test can actually fail: `A-A` sits at index 9 of
    /// "AI-CLACVAA-AC-KCRLY$" and `Y$` at index 18, so both *would* be found if the structural
    /// characters were searchable. A peptide that merely does not occur would be dropped by the
    /// no-match path whether or not it was filtered, and would prove nothing.
    #[test]
    fn test_rejected_peptides_are_dropped_from_both_paths() {
        let searcher = test_searcher(true);

        // Control: strip the structural character and each one is findable, so the assertions
        // below fail for the right reason — rejected, not simply absent.
        for present in ["AA", "Y", "AC"] {
            assert!(
                search_peptide(&searcher, present, 1000, false, false).is_some(),
                "{present} should be found in the example text"
            );
        }

        // Each of these would match the raw text if it reached the searcher.
        for rejected in ["A-A", "Y$", "AA-AC"] {
            assert!(
                search_peptide(&searcher, rejected, 1000, false, false).is_none(),
                "{rejected} contains a structural character and must be rejected"
            );
        }

        // Outside the alphabet for other reasons: unmatchable, but they must not reach the
        // searcher's unchecked table indices either.
        for rejected in ["ACé", "AJ", "AC1"] {
            assert!(search_peptide(&searcher, rejected, 1000, false, false).is_none(), "{rejected}");
        }

        // The batch path drops exactly the same peptides, and keeps input order for the rest.
        let peptides: Vec<String> = ["AA", "A-A", "Y", "Y$", "AC", "ACé"].iter().map(|s| s.to_string()).collect();
        let got = search_all_peptides(&searcher, &peptides, 1000, false, false);
        let found: Vec<&str> = got.iter().map(|r| r.sequence).collect();
        assert_eq!(found, vec!["AA", "Y", "AC"], "only the alphabet-clean peptides survive");
    }

    /// `cutoff_used` must mean "matches were dropped", not "the cutoff was reached".
    ///
    /// Every cutoff test used `>=`, so a complete set of exactly `cutoff` matches came back flagged
    /// as a truncated sample. That flag is what an API consumer reads to decide whether a protein
    /// list is exhaustive, so it was a user-visible wrong answer on a healthy index — and the
    /// smaller the client's `cutoff`, the more often it landed on the boundary.
    ///
    /// The taxa a peptide reports must be exactly the distinct taxa of the proteins the protein
    /// path reports for it — that is the whole contract of the lightweight path, and the only way
    /// it can be wrong without any test noticing is if the two disagree.
    ///
    /// Checked for both entry points, since they reach the suffixes by different routes: the
    /// batched pipeline and the scalar one.
    #[test]
    fn the_taxa_path_reports_exactly_the_taxa_of_the_protein_path() {
        let searcher = test_searcher(true);
        let peptides: Vec<String> =
            ["A", "C", "AC", "CLACVAA", "KCRLY", "WWW", "ai"].iter().map(|p| p.to_string()).collect();

        let by_protein = search_all_peptides(&searcher, &peptides, 1_000_000, false, false);
        let by_taxa = search_all_peptides_taxa(&searcher, &peptides, 1_000_000, false, false);

        assert!(!by_taxa.is_empty(), "fixture must produce matches");
        assert_eq!(by_taxa.len(), by_protein.len(), "both paths must drop the same peptides");

        for (proteins, taxa) in by_protein.iter().zip(by_taxa.iter()) {
            assert_eq!(proteins.sequence, taxa.sequence, "results must stay in the same order");
            assert_eq!(proteins.cutoff_used, taxa.cutoff_used, "{}: cutoff flag differs", taxa.sequence);

            let mut expected: Vec<u32> = proteins.proteins.iter().map(|protein| protein.taxon).collect();
            expected.sort_unstable();
            expected.dedup();
            assert_eq!(taxa.taxa, expected, "{}: taxa differ from the protein path", taxa.sequence);

            // The single-peptide path must agree with the batched one it does not share.
            let single = search_peptide_taxa(&searcher, taxa.sequence, 1_000_000, false, false)
                .unwrap_or_else(|| panic!("{}: batched path matched, scalar path did not", taxa.sequence));
            assert_eq!(single.taxa, taxa.taxa, "{}: scalar and batched taxa differ", taxa.sequence);
            assert_eq!(single.cutoff_used, taxa.cutoff_used, "{}: scalar and batched cutoff differ", taxa.sequence);
        }
    }

    /// Taxa are deduplicated, which is the point: `A` matches five suffixes across three proteins
    /// carrying three distinct taxa, and the protein path reports all five.
    #[test]
    fn taxa_are_deduplicated_where_the_protein_path_repeats_them() {
        let searcher = test_searcher(false);

        let proteins = search_peptide(&searcher, "A", 1_000_000, false, false).unwrap();
        let taxa = search_peptide_taxa(&searcher, "A", 1_000_000, false, false).unwrap();

        assert!(proteins.proteins.len() > taxa.taxa.len(), "fixture must actually contain duplicate taxa");
        assert_eq!(taxa.taxa, vec![10, 20, 30]);
    }

    /// An unsearchable peptide is rejected by the taxa path on the same terms as the protein path,
    /// including the alphabet filter — a peptide carrying the structural `-` must not reach the
    /// searcher, where it could match across a protein boundary.
    #[test]
    fn the_taxa_path_rejects_what_the_protein_path_rejects() {
        let searcher = test_searcher(false);

        for peptide in ["WWW", "A-C", "A$", "J", "AJ", "a c", "\u{e9}"] {
            assert_eq!(
                search_peptide_taxa(&searcher, peptide, 1_000_000, false, false).is_some(),
                search_peptide(&searcher, peptide, 1_000_000, false, false).is_some(),
                "{peptide}: the two paths disagree about whether this is searchable"
            );
        }
    }

    #[test]
    fn a_taxa_result_serialises_to_the_documented_shape() {
        let result = TaxaSearchResult {
            sequence: "MSKIAALLPSV",
            taxa: vec![1, 9606],
            cutoff_used: false
        };

        assert_json_eq(
            &serde_json::to_string(&result).unwrap(),
            "{\"sequence\":\"MSKIAALLPSV\",\"taxa\":[1,9606],\"cutoff_used\":false}"
        );
    }

    /// The counts are discovered rather than hard-coded, so the test states the property instead of
    /// restating the fixture.
    #[test]
    fn cutoff_used_marks_only_genuinely_truncated_results() {
        let searcher = test_searcher(true);

        for peptide in ["A", "C", "V", "AC", "KCRLY"] {
            // Ground truth: the complete result, found with a cutoff nothing can reach.
            let Some(full) = search_peptide(&searcher, peptide, 1_000_000, false, false) else {
                continue;
            };
            let total = full.proteins.len();
            assert!(total > 0, "{peptide}: fixture peptide should match");
            assert!(!full.cutoff_used, "{peptide}: a cutoff far above the match count is not truncation");

            // Exactly the match count: everything is returned, so nothing was dropped.
            let at = search_peptide(&searcher, peptide, total, false, false).unwrap();
            assert_eq!(at.proteins.len(), total, "{peptide}: cutoff == total must still return everything");
            assert!(!at.cutoff_used, "{peptide}: a complete set of exactly `cutoff` is not truncated");

            // One below: a match really was dropped.
            if total > 1 {
                let below = search_peptide(&searcher, peptide, total - 1, false, false).unwrap();
                assert_eq!(below.proteins.len(), total - 1, "{peptide}: cutoff must cap the result");
                assert!(below.cutoff_used, "{peptide}: dropping a match must be reported");
            }

            // The batched path drives the same flag and must agree at every boundary — this is the
            // kind of one-element difference the two searchers drift apart on.
            let one = [peptide.to_string()];
            for cutoff in [total.saturating_sub(1).max(1), total, total + 1] {
                let batched = search_all_peptides(&searcher, &one, cutoff, false, false);
                let scalar = search_peptide(&searcher, peptide, cutoff, false, false).unwrap();
                assert_eq!(batched.len(), 1, "{peptide}: expected one result at cutoff {cutoff}");
                assert_eq!(
                    batched[0].cutoff_used, scalar.cutoff_used,
                    "{peptide}: batched and scalar disagree on cutoff_used at cutoff {cutoff}"
                );
                assert_eq!(
                    batched[0].proteins.len(),
                    scalar.proteins.len(),
                    "{peptide}: batched and scalar returned different counts at cutoff {cutoff}"
                );
            }
        }
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
        let search_result = SearchResult { sequence: "MSKIAALLPSV", proteins: vec![], cutoff_used: true };

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
        let protein_info = ProteinInfo { taxon: 1, uniprot_accession: "P12345", annotations: &[] };
        assert_eq!(
            serde_json::to_string(&protein_info).unwrap(),
            "{\"taxon\":1,\"uniprot_accession\":\"P12345\",\"functional_annotations\":\"\"}"
        );
    }
}
