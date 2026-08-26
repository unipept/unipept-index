//! The tryptic filter: deciding whether a match sits at a real trypsin cut.
//!
//! Trypsin cuts after K or R, except when proline follows. A tryptic search therefore accepts a
//! match only when both of its boundaries are either a cut site or a protein edge — which is what
//! [`Searcher::check_tryptic_boundaries`] decides.
//!
//! Two things here exist to make that cheap enough for the hot loop:
//!
//! * [`TrypticQuery`] resolves the half of the predicate that depends only on the peptide, once
//!   per query rather than once per candidate, so the per-candidate check drops from four text
//!   reads to two.
//! * [`tryptic_extension_chars`] lets the *search* find left-extended matches directly, instead of
//!   finding every match and filtering afterwards. `scalar.rs` and `batched.rs` drive that.

use protein_metadata::{ProteinsBackend, SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinTextBackend;

use super::Searcher;
use crate::{array::SuffixArrayBackend, suffix_to_protein_index::SuffixToProteinMappingBackend};

/// The query-invariant half of the tryptic-cut predicate, resolved once per peptide instead of
/// once per candidate match.
///
/// A candidate that reaches the tryptic check has already been confirmed to equal the search
/// string over `[match_start, match_end)` — either by the suffix array range itself or by the
/// `check_prefix` / `check_suffix` conjuncts that run before it. So the two text positions the
/// naive predicate reads *inside* the match are redundant: the peptide already knows them.
///
/// * `text[match_start]` is `search_string[0]`
/// * `text[match_end - 1]` is `search_string.last()`
///
/// The substitution holds under both `equate_il` settings: K, R and P are untouched by the
/// L → I normalization the index is built with, so a text position can only hold P/K/R when the
/// corresponding peptide character does.
#[derive(Clone, Copy)]
pub(super) enum TrypticQuery {
    /// Not filtering for tryptic matches — every candidate passes the boundary check.
    Off,
    /// Filtering, with both query-invariant terms already resolved.
    On {
        /// `search_string[0] != b'P'`, standing in for `text[match_start] != b'P'`.
        /// Proline directly after a K/R blocks the trypsin cut at the N-terminus.
        first_not_proline: bool,
        /// `search_string.last()` is K or R, standing in for `text[match_end - 1]` being a
        /// trypsin cut site.
        last_is_kr: bool
    },
    /// Degenerate zero-length query: `match_start == match_end`, so neither substitution above
    /// exists and we fall back to the original four-read formulation. Unreachable in production
    /// (`search_proteins_for_peptide` drops peptides shorter than the sparseness factor), but
    /// `search_matching_suffixes_scalar` is public, so stay bit-exact for it anyway.
    ZeroLength
}

/// Characters that can legally precede a tryptic match: the two trypsin cut residues plus the
/// protein separator (a match at a protein start needs no cut).
///
/// Searching `X + peptide` for each of these — instead of truncating the peptide — is what makes
/// the tryptic path cheap; see `search_matching_suffixes_scalar` for the derivation. Ordered K, R,
/// separator so the two common cases run first.
pub(super) const TRYPTIC_EXTENSION_CHARS: [u8; 3] = [b'K', b'R', b'-'];

/// The separator-only subset, used when the peptide starts with proline: `K|R` followed by P is
/// not a trypsin cut site, so those two searches are guaranteed empty and are skipped entirely.
const TRYPTIC_EXTENSION_CHARS_PROLINE: [u8; 1] = [b'-'];

/// Left-extension characters to search for this peptide.
#[inline]
pub(super) fn tryptic_extension_chars(search_string: &[u8]) -> &'static [u8] {
    debug_assert_eq!(TRYPTIC_EXTENSION_CHARS[2], SEPARATION_CHARACTER);
    if search_string.first() == Some(&b'P') { &TRYPTIC_EXTENSION_CHARS_PROLINE } else { &TRYPTIC_EXTENSION_CHARS }
}

impl TrypticQuery {
    #[inline]
    pub(super) fn new(tryptic: bool, search_string: &[u8]) -> Self {
        if !tryptic {
            return Self::Off;
        }
        match (search_string.first(), search_string.last()) {
            (Some(&first), Some(&last)) => Self::On {
                first_not_proline: first != b'P',
                last_is_kr: last == b'K' || last == b'R'
            },
            _ => Self::ZeroLength
        }
    }
}

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Check if a cut is the start of a protein.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is at the start of a protein.
    #[inline]
    fn check_start_of_protein(&self, cut_index: usize) -> bool {
        cut_index == 0 || self.proteins.text().get(cut_index - 1) == SEPARATION_CHARACTER
    }

    /// Check if a cut is the end of a protein.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is at the end of a protein.
    #[inline]
    fn check_end_of_protein(&self, cut_index: usize) -> bool {
        self.proteins.text().get(cut_index) == TERMINATION_CHARACTER
            || self.proteins.text().get(cut_index) == SEPARATION_CHARACTER
    }

    /// Check if a cut is a tryptic cut, so check if the amino acid preceding the cut is K or R and the amino acid at the cut is not P.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is a tryptic cut.
    #[inline]
    fn check_tryptic_cut(&self, cut_index: usize) -> bool {
        (self.proteins.text().get(cut_index - 1) == b'K' || self.proteins.text().get(cut_index - 1) == b'R')
            && self.proteins.text().get(cut_index) != b'P'
    }

    /// Tryptic-boundary check for a candidate match spanning `[match_start, match_end)`.
    ///
    /// Equivalent to
    /// `(check_start_of_protein(ms) || check_tryptic_cut(ms)) && (check_end_of_protein(me) || check_tryptic_cut(me))`
    /// but reads **two** text positions instead of four: only `ms - 1` and `me`. The other two
    /// (`ms` and `me - 1`) sit inside the match, where the text equals the peptide, so `query`
    /// supplies them without touching memory (see `TrypticQuery`).
    ///
    /// This is the hottest loop in the index: `tryptic` disables the fast path in `scalar.rs`,
    /// so the candidate scan has to run until `max_matches` *accepted* matches accumulate —
    /// work that scales as `max_matches / acceptance_rate`. Halving the reads also halves the
    /// prefetch hints `prefetch_match_positions` has to keep in flight, which matters more than
    /// the ALU saving (see the comment there).
    #[inline]
    pub(super) fn check_tryptic_boundaries(
        &self,
        text: &P::Text,
        match_start: usize,
        match_end: usize,
        query: TrypticQuery
    ) -> bool {
        match query {
            TrypticQuery::Off => true,
            TrypticQuery::On { first_not_proline, last_is_kr } => {
                // N-terminus. The `match_start == 0` guard has to come first: without it the
                // `match_start - 1` read underflows. (The original formulation was only safe
                // because `check_start_of_protein` short-circuits the `||` on that same test.)
                let n_term_ok = match_start == 0 || {
                    let before = text.get(match_start - 1);
                    before == SEPARATION_CHARACTER || ((before == b'K' || before == b'R') && first_not_proline)
                };

                // C-terminus. Delegated rather than repeated: this used to be an inline copy of
                // `check_tryptic_c_term`'s body *without* its `match_end >= text.len()` guard, so
                // the two halves of the same predicate disagreed about the end of the text and
                // this one panicked there. `#[inline]` keeps it a single read either way.
                n_term_ok && self.check_tryptic_c_term(text, match_end, last_is_kr)
            }
            // Zero-length query: `match_end == match_start`, so the peptide has no character to
            // stand in for either read and the original formulation is the only correct one.
            TrypticQuery::ZeroLength => {
                (self.check_start_of_protein(match_start) || self.check_tryptic_cut(match_start))
                    && (self.check_end_of_protein(match_end) || self.check_tryptic_cut(match_end))
            }
        }
    }

    /// Checks only the C-terminal half of the tryptic predicate.
    ///
    /// Called by [`Self::check_tryptic_boundaries`] for its second conjunct, and directly by the
    /// left-extended search path, where the N-terminal half holds by construction: the search
    /// string was `X + peptide` for `X` in {K, R, separator}, so the character before the match
    /// *is* `X`. Nothing to read and nothing to check.
    ///
    /// `last_is_kr` is the query-invariant stand-in for `text[match_end - 1]` — see
    /// [`TrypticQuery`] for why the substitution is sound.
    ///
    /// # On reading `text[match_end]`
    ///
    /// `match_end` is one past the match, so this reads the character *after* it. For a real
    /// peptide that is always in bounds: the text ends with [`TERMINATION_CHARACTER`], and a
    /// match reaching `text.len()` would have to have consumed the terminator — that is, the
    /// query's own last character would be `$`, which no protein sequence contains.
    ///
    /// A caller *can* construct such a query, since `$` is in the alphabet, so the bound is
    /// checked rather than assumed. It used to be neither, and the two backends then failed
    /// differently on the same input: the preloaded one panics whenever index `len` lands in a
    /// word the bit array never allocated (`5 * len % 64 >= 60`, i.e. `len % 64` in
    /// {12, 25, 38, 51}, plus `len % 64 == 0`), while the mmap one reads whatever sits inside the
    /// page-rounded mapping and can return *either* verdict — the first two clauses test `after`
    /// directly, so a stray separator or terminator byte makes the predicate accept.
    ///
    /// The guard has to live *here* rather than at each call site, because it is exactly what the
    /// inline copy in [`Self::check_tryptic_boundaries`] was missing: the extended path was fixed
    /// and the normal path was not, and nothing tested the pair together.
    ///
    /// One past the last residue is the end of the final protein, so it is reported as a valid
    /// C-terminal cut, symmetric with [`Searcher::check_start_of_protein`] treating index 0 as a
    /// protein start.
    #[inline]
    pub(super) fn check_tryptic_c_term(&self, text: &P::Text, match_end: usize, last_is_kr: bool) -> bool {
        if match_end >= text.len() {
            return true;
        }
        let after = text.get(match_end);
        after == TERMINATION_CHARACTER || after == SEPARATION_CHARACTER || (last_is_kr && after != b'P')
    }
}

#[cfg(test)]
mod tests {
    use protein_metadata::ProteinsBackend as _;
    use text_compression::ProteinTextBackend as _;

    use super::TrypticQuery;
    use crate::sa_searcher::test_utils::{PreloadedSearcher, example_searcher, searcher_over_text};

    /// `match_end == text.len()` is reachable from the server: `$` is in the query alphabet, so a
    /// query ending in `$` can match through the terminator. The read must be bounded.
    ///
    /// The text length is what makes this bite. The preloaded bit array allocates only
    /// `ceil(5 * len / 64)` words, so reading index `len` fell outside them whenever
    /// `5 * len % 64 >= 60` — that is `len % 64` in {12, 25, 38, 51} — or whenever `5 * len` was
    /// an exact multiple of 64 (`len % 64 == 0`). The mmap backend never faulted there at all; it
    /// read whatever followed inside the page-rounded mapping, so the two backends could return
    /// opposite verdicts for the same query. Every length below is one of those cases.
    #[test]
    fn test_c_term_check_is_bounded_at_end_of_text() {
        for len in [12usize, 25, 38, 51, 64] {
            let text = format!("{}$", "A".repeat(len - 1));
            assert_eq!(text.len(), len);

            let searcher = searcher_over_text(&text, 1);
            let t = searcher.proteins.text();
            assert_eq!(t.len(), len, "fixture length must survive the builder");

            // One past the last residue is the end of the final protein: a valid C-terminal cut,
            // symmetric with `check_start_of_protein(0)`. Neither argument may panic.
            assert!(searcher.check_tryptic_c_term(t, t.len(), true), "len {len}, last_is_kr");
            assert!(searcher.check_tryptic_c_term(t, t.len(), false), "len {len}, not kr");
        }
    }

    /// The other half of the same bound, on the other entry point.
    ///
    /// `check_tryptic_c_term` guards `match_end >= text.len()`; `check_tryptic_boundaries` carried
    /// an inline copy of that read that did not, so the normal search path panicked at the end of
    /// the text on exactly the lengths the sibling test pins — the extended path had been fixed
    /// and this one had not. Asserting the two agree, rather than just that neither panics, is
    /// what stops them drifting apart again.
    #[test]
    fn test_boundary_check_is_bounded_at_end_of_text() {
        for len in [12usize, 25, 38, 51, 64] {
            let text = format!("{}$", "A".repeat(len - 1));
            assert_eq!(text.len(), len);

            let searcher = searcher_over_text(&text, 1);
            let t = searcher.proteins.text();
            assert_eq!(t.len(), len, "fixture length must survive the builder");

            for last_is_kr in [true, false] {
                for first_not_proline in [true, false] {
                    let query = TrypticQuery::On { first_not_proline, last_is_kr };
                    // `match_start == 0` is a protein start, so the N-terminal conjunct holds and
                    // the verdict is the C-terminal one alone — which is what must not panic.
                    assert_eq!(
                        searcher.check_tryptic_boundaries(t, 0, t.len(), query),
                        searcher.check_tryptic_c_term(t, t.len(), last_is_kr),
                        "len {len}, last_is_kr {last_is_kr}, first_not_proline {first_not_proline}"
                    );
                }
            }
        }
    }

    // Direct test of the protein-boundary checks used by tryptic filtering.
    // Text "AI-CLACVAA-AC-KCRLY$": separators at 2/10/13, termination at 19.
    #[test]
    fn test_check_protein_boundaries() {
        let searcher = example_searcher();

        // start of a protein: index 0, or immediately after a separator
        for i in [0usize, 3, 11, 14] {
            assert!(searcher.check_start_of_protein(i), "expected protein start at {i}");
        }
        for i in [1usize, 5, 12] {
            assert!(!searcher.check_start_of_protein(i), "unexpected protein start at {i}");
        }

        // end of a protein: the position itself holds a separator or termination char
        for i in [2usize, 10, 13, 19] {
            assert!(searcher.check_end_of_protein(i), "expected protein end at {i}");
        }
        for i in [0usize, 1, 5] {
            assert!(!searcher.check_end_of_protein(i), "unexpected protein end at {i}");
        }
    }

    // A tryptic cut is valid iff preceded by K or R and NOT followed by proline (P).
    #[test]
    fn test_check_tryptic_cut() {
        // K(0) A(1) R(2) C(3) K(4) P(5) D(6) $(7)
        let searcher = searcher_over_text("KARCKPD$", 1);

        assert!(searcher.check_tryptic_cut(1)); // after K, C follows
        assert!(searcher.check_tryptic_cut(3)); // after R, C follows
        assert!(!searcher.check_tryptic_cut(5)); // after K but P follows (proline blocks)
        assert!(!searcher.check_tryptic_cut(2)); // preceded by A, not K/R
    }

    // The original four-text-read formulation, kept verbatim as the oracle for the equivalence
    // test below. `check_start_of_protein` MUST stay on the left of the `||`: it is the only
    // thing that stops `check_tryptic_cut(0)` from reading `text[-1]`.
    fn old_tryptic_predicate(searcher: &PreloadedSearcher, match_start: usize, match_end: usize) -> bool {
        (searcher.check_start_of_protein(match_start) || searcher.check_tryptic_cut(match_start))
            && (searcher.check_end_of_protein(match_end) || searcher.check_tryptic_cut(match_end))
    }

    // A text exercising every boundary case the tryptic filter can hit:
    //
    //   index 0..: K  P  A  R  C  K  A  -  R  P  K  R  A  -  P  K  -  A  I  L  K  $
    //
    // - match_start == 0 (the underflow guard)
    // - protein starts right after the separators at 7, 13 and 16
    // - protein ends at the separators 7/13/16 and at the terminator 21
    // - K/R-preceded cuts (positions 1, 4, 6, 11, 12, 16, 21)
    // - P-blocked cuts ("KP" at 0-1, "RP" at 8-9, "KP"-like "K-" vs "AP" at 13-14)
    // - a single-residue protein ("K" between the separators at 13 and 16 is "PK")
    // - I and L adjacent at 18/19, so peptides can differ from the text at L/I positions
    const BOUNDARY_TEXT: &str = "KPARCKA-RPKRA-PK-AILK$";

    // Every candidate the filter can ever see: for each span of the text, the peptide IS the
    // text over that span, which is exactly the invariant the two dropped reads rely on.
    // The new two-read predicate must agree with the old four-read one on every single one.
    #[test]
    fn test_tryptic_boundaries_equivalent_to_original() {
        let searcher = searcher_over_text(BOUNDARY_TEXT, 1);
        let text = searcher.proteins.text();
        let text_len = BOUNDARY_TEXT.len();

        let mut checked = 0usize;
        let mut accepted = 0usize;
        for match_start in 0..text_len {
            // `match_end == text_len` would read one past the end in BOTH formulations, and
            // cannot occur for a real match: the text always ends in the terminator, which no
            // peptide contains.
            for match_end in (match_start + 1)..text_len {
                let peptide: Vec<u8> = (match_start..match_end).map(|i| text.get(i)).collect();
                let query = TrypticQuery::new(true, &peptide);

                let expected = old_tryptic_predicate(&searcher, match_start, match_end);
                let actual = searcher.check_tryptic_boundaries(text, match_start, match_end, query);
                assert_eq!(
                    actual,
                    expected,
                    "mismatch at [{match_start}, {match_end}) for peptide {:?}",
                    String::from_utf8_lossy(&peptide)
                );

                // L/I equating: the index stores L as I, so a peptide may carry an L where the
                // text holds an I. K, R and P survive that normalization untouched, so swapping
                // I <-> L in the peptide must not move the predicate either way.
                let swapped: Vec<u8> = peptide
                    .iter()
                    .map(|&c| match c {
                        b'I' => b'L',
                        b'L' => b'I',
                        other => other
                    })
                    .collect();
                let swapped_actual =
                    searcher.check_tryptic_boundaries(text, match_start, match_end, TrypticQuery::new(true, &swapped));
                assert_eq!(swapped_actual, expected, "I/L swap changed the verdict at [{match_start}, {match_end})");

                checked += 1;
                accepted += expected as usize;
            }
        }

        // Guard against the test silently degenerating: the text must produce a meaningful mix
        // of accepted and rejected candidates, otherwise the equivalence proves nothing.
        assert_eq!(checked, text_len * (text_len - 1) / 2);
        assert!(accepted > 0 && accepted < checked, "degenerate fixture: {accepted}/{checked} accepted");
    }

    // The zero-length fallback keeps reading the text, so it must also match the oracle.
    //
    // `cut == 0` is skipped because it panics identically in both formulations: with
    // `match_start == match_end == 0`, `check_end_of_protein(0)` is false for a text that does
    // not start with a separator, and the `check_tryptic_cut(0)` that follows it underflows on
    // `text[-1]`. That is pre-existing behaviour of the degenerate zero-length path (only the
    // N-terminal `||` has a `match_start == 0` guard, the C-terminal one does not), and it is
    // unreachable in production, where empty peptides are dropped before the search.
    #[test]
    fn test_tryptic_boundaries_zero_length_matches_original() {
        let searcher = searcher_over_text(BOUNDARY_TEXT, 1);
        let text = searcher.proteins.text();

        for cut in 1..BOUNDARY_TEXT.len() {
            assert_eq!(
                searcher.check_tryptic_boundaries(text, cut, cut, TrypticQuery::new(true, b"")),
                old_tryptic_predicate(&searcher, cut, cut),
                "zero-length mismatch at {cut}"
            );
        }
    }

    // `tryptic = false` short-circuits the whole check, without touching the text.
    #[test]
    fn test_tryptic_boundaries_off_accepts_everything() {
        let searcher = searcher_over_text(BOUNDARY_TEXT, 1);
        let text = searcher.proteins.text();

        for match_start in 0..BOUNDARY_TEXT.len() {
            for match_end in match_start..BOUNDARY_TEXT.len() {
                let query = TrypticQuery::new(false, b"PAAP");
                assert!(searcher.check_tryptic_boundaries(text, match_start, match_end, query));
            }
        }
    }

    // The two query-invariant terms are read off the peptide itself, so they must be right for
    // every combination of first/last residue that changes a verdict.
    #[test]
    fn test_tryptic_query_hoisting() {
        assert!(matches!(TrypticQuery::new(false, b"PAAK"), TrypticQuery::Off));
        assert!(matches!(TrypticQuery::new(true, b""), TrypticQuery::ZeroLength));

        let flags = |peptide: &[u8]| match TrypticQuery::new(true, peptide) {
            TrypticQuery::On { first_not_proline, last_is_kr } => (first_not_proline, last_is_kr),
            _ => panic!("expected TrypticQuery::On for {:?}", String::from_utf8_lossy(peptide))
        };

        // (first_not_proline, last_is_kr)
        assert_eq!(flags(b"PAAC"), (false, false)); // starts with P
        assert_eq!(flags(b"PAAK"), (false, true)); // starts with P, ends with K
        assert_eq!(flags(b"AAAP"), (true, false)); // ends with P
        assert_eq!(flags(b"AAAK"), (true, true)); // ends with K
        assert_eq!(flags(b"AAAR"), (true, true)); // ends with R
        assert_eq!(flags(b"AAAC"), (true, false)); // neither
        assert_eq!(flags(b"P"), (false, false)); // single residue: first == last
        assert_eq!(flags(b"K"), (true, true));
        assert_eq!(flags(b"R"), (true, true));
        // K/R are never L/I-normalized, so a lowercase-free peptide of I/L cannot be a cut site.
        assert_eq!(flags(b"ILI"), (true, false));

        // And the flags actually drive the verdict. Text "KPAAK-": the K at 0 precedes the
        // match, so the N-terminal cut hinges entirely on the peptide's own first residue.
        let searcher = searcher_over_text("KPAAK-AAAA$", 1);
        let text = searcher.proteins.text();
        // [1, 5) is "PAAK": preceded by K, but the peptide starts with P -> N-term blocked.
        assert!(!searcher.check_tryptic_boundaries(text, 1, 5, TrypticQuery::new(true, b"PAAK")));
        // [2, 5) is "AAK": preceded by P, which is neither K/R nor a separator -> N-term fails.
        assert!(!searcher.check_tryptic_boundaries(text, 2, 5, TrypticQuery::new(true, b"AAK")));
        // [0, 5) is "KPAAK": match_start == 0 is a protein start, and the match ends on the
        // separator at 5 -> both ends valid.
        assert!(searcher.check_tryptic_boundaries(text, 0, 5, TrypticQuery::new(true, b"KPAAK")));
        // [0, 1) is "K": protein start, but text[1] is P, so the K cut is proline-blocked.
        assert!(!searcher.check_tryptic_boundaries(text, 0, 1, TrypticQuery::new(true, b"K")));
        // [6, 10) is "AAAA": protein start (separator at 5) but the peptide does not end in
        // K/R and text[10] is the terminator -> C-term valid via end-of-protein.
        assert!(searcher.check_tryptic_boundaries(text, 6, 10, TrypticQuery::new(true, b"AAAA")));
        // [6, 9) is "AAA": neither end-of-protein nor a K/R terminus -> rejected.
        assert!(!searcher.check_tryptic_boundaries(text, 6, 9, TrypticQuery::new(true, b"AAA")));
    }
}
