//! Which storage backend each index structure uses in this build.
//!
//! This module is the **only** place in the workspace where a storage feature is read. The
//! libraries compile both backends of every structure unconditionally and are generic over which
//! one they are handed; the features exist here, at the top of the dependency graph, purely to
//! name one concrete type per structure.
//!
//! Four independent choices. `mmap` maps everything; each `preloaded-*` then pulls **one**
//! structure back into owned memory, leaving the rest mapped:
//!
//! | alias | owned | mapped | mapped when |
//! |---|---|---|---|
//! | [`ActiveSa`] | `InMemorySA` | `MmapBackedSA` | `mmap` |
//! | [`ActiveText`] | `InMemoryProteinText` | `MmapBackedProteinText` | `mmap` and not `preloaded-text` |
//! | [`ActiveProteins`] | `InMemoryProteins<T>` | `MmapBackedProteins<T>` | `mmap` and not `preloaded-proteins` |
//! | [`ActiveMapping`] | `InMemorySuffixToProteinMapping` | `MmapBackedSuffixToProteinMapping` | `mmap` and not `preloaded-mapping` |
//!
//! Nine configurations in all: everything preloaded, everything mapped, and the seven mixtures.
//! The point is that the best place for one structure is not the best place for another — the text
//! is the hottest and the metadata table the biggest — so, for instance
//! `--features mmap,preloaded-text` keeps the multi-gigabyte index mapped while the ~190 MB text
//! that search reads once per character compared sits in owned RAM.
//!
//! Two things follow that are easy to trip over:
//!
//! * **No crate declares `default = [...]`**, so a plain `cargo build` gives the fully *preloaded*
//!   configuration. The production server is built `--features mmap`.
//! * The suffix array follows `mmap` and has no override; there is no `preloaded-sa`. A
//!   `preloaded-*` feature without `mmap` is a no-op — everything is preloaded already. Cargo
//!   features are additive and cannot be negated by a dependent crate, so they only ever *remove*
//!   mapping, never add it.
//!
//! The text and the metadata share one file (`proteins.bin`) but are separate axes, which is why
//! [`ActiveProteins`] instantiates the protein struct with [`ActiveText`]: the text axis is already
//! resolved by the time it is substituted, so the two compose without a case per combination.
//! Which reader that pairing needs is not decided here — it is the `LoadIndex` implementation on
//! the pairing itself; see `protein_metadata::mmap`.

// Every type below is named by its full path rather than imported: only one arm of each pair is
// compiled, so an import would be unused in the configuration that does not take it.

/// The suffix-array backend this build uses.
#[cfg(feature = "mmap")]
pub type ActiveSa = sa_index::array::MmapBackedSA;
/// The suffix-array backend this build uses.
#[cfg(not(feature = "mmap"))]
pub type ActiveSa = sa_index::array::InMemorySA;

/// The protein-text backend this build uses.
#[cfg(all(feature = "mmap", not(feature = "preloaded-text")))]
pub type ActiveText = text_compression::MmapBackedProteinText;
/// The protein-text backend this build uses.
#[cfg(any(not(feature = "mmap"), feature = "preloaded-text"))]
pub type ActiveText = text_compression::InMemoryProteinText;

/// The protein-metadata backend this build uses, holding [`ActiveText`].
#[cfg(all(feature = "mmap", not(feature = "preloaded-proteins")))]
pub type ActiveProteins = protein_metadata::MmapBackedProteins<ActiveText>;
/// The protein-metadata backend this build uses, holding [`ActiveText`].
#[cfg(any(not(feature = "mmap"), feature = "preloaded-proteins"))]
pub type ActiveProteins = protein_metadata::InMemoryProteins<ActiveText>;

/// The suffix-to-protein mapping backend this build uses.
#[cfg(all(feature = "mmap", not(feature = "preloaded-mapping")))]
pub type ActiveMapping = sa_index::suffix_to_protein_index::MmapBackedSuffixToProteinMapping;
/// The suffix-to-protein mapping backend this build uses.
#[cfg(any(not(feature = "mmap"), feature = "preloaded-mapping"))]
pub type ActiveMapping = sa_index::suffix_to_protein_index::InMemorySuffixToProteinMapping;

/// The searcher this build serves from — the three aliases above, assembled.
///
/// `Searcher` itself is generic over all three, so this alias is what saves every signature in
/// `sa-server` and `sa-benchmarks` from spelling them out.
pub type ActiveSearcher = sa_index::sa_searcher::Searcher<ActiveSa, ActiveProteins, ActiveMapping>;

/// How the suffix array is stored in this build.
///
/// Baked in at compile time, so it cannot be inspected any other way at runtime — which is why the
/// server logs all four at startup. The configurations have very different memory profiles, and
/// telling them apart otherwise means checking how the binary was built.
pub const SA_BACKEND: &str = if cfg!(feature = "mmap") { "mmap" } else { "preloaded" };
/// How the concatenated protein text is stored in this build.
pub const TEXT_BACKEND: &str =
    if cfg!(all(feature = "mmap", not(feature = "preloaded-text"))) { "mmap" } else { "preloaded" };
/// How the protein metadata table is stored in this build.
pub const PROTEINS_BACKEND: &str =
    if cfg!(all(feature = "mmap", not(feature = "preloaded-proteins"))) { "mmap" } else { "preloaded" };
/// How the suffix-to-protein mapping is stored in this build.
pub const MAPPING_BACKEND: &str =
    if cfg!(all(feature = "mmap", not(feature = "preloaded-mapping"))) { "mmap" } else { "preloaded" };

/// One line naming the storage of every structure, for the startup log.
pub fn backend_summary() -> String {
    format!("sa={SA_BACKEND} text={TEXT_BACKEND} proteins={PROTEINS_BACKEND} mapping={MAPPING_BACKEND}")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Whether `type_name` names a mapped backend, ignoring any type parameters.
    ///
    /// Splitting at `<` is what makes this answer about the *outer* type: `ActiveProteins` may be
    /// `MmapBackedProteins<InMemoryProteinText>` or `InMemoryProteins<MmapBackedProteinText>`, and
    /// the metadata axis is the one outside the angle brackets.
    fn names_a_mapped_type<T>() -> bool {
        std::any::type_name::<T>().split('<').next().unwrap().contains("Mmap")
    }

    /// The four reported strings and the four selected types are separate `#[cfg]` blocks over the
    /// same features, so editing one and not the other would make the server's startup line — the
    /// only runtime evidence of how it was built — quietly lie.
    ///
    /// This checks whichever configuration it is compiled in; CI runs it across all nine.
    #[test]
    fn reported_backends_match_the_types_actually_selected() {
        assert_eq!(SA_BACKEND == "mmap", names_a_mapped_type::<ActiveSa>(), "SA_BACKEND is {SA_BACKEND}");
        assert_eq!(TEXT_BACKEND == "mmap", names_a_mapped_type::<ActiveText>(), "TEXT_BACKEND is {TEXT_BACKEND}");
        assert_eq!(
            PROTEINS_BACKEND == "mmap",
            names_a_mapped_type::<ActiveProteins>(),
            "PROTEINS_BACKEND is {PROTEINS_BACKEND}"
        );
        assert_eq!(
            MAPPING_BACKEND == "mmap",
            names_a_mapped_type::<ActiveMapping>(),
            "MAPPING_BACKEND is {MAPPING_BACKEND}"
        );
    }

    /// The protein struct holds the text, so the two axes are only independent if `ActiveProteins`
    /// is instantiated at `ActiveText` — otherwise a `preloaded-text` build could report a
    /// preloaded text while serving from a mapped one.
    #[test]
    fn the_protein_backend_holds_the_selected_text_backend() {
        assert!(
            std::any::type_name::<ActiveProteins>().contains(std::any::type_name::<ActiveText>()),
            "ActiveProteins is {}, which does not hold ActiveText ({})",
            std::any::type_name::<ActiveProteins>(),
            std::any::type_name::<ActiveText>()
        );
    }

    /// Every alias has to be loadable, which is not automatic: `ActiveProteins` is one of four
    /// pairings and each needs its own `LoadIndex` impl. A missing one is a compile error here
    /// rather than at the call site in `lib.rs`.
    #[test]
    fn every_active_backend_can_be_loaded() {
        fn assert_loadable<T: binary_traits::LoadIndex>() {}
        assert_loadable::<ActiveSa>();
        assert_loadable::<ActiveProteins>();
        assert_loadable::<ActiveMapping>();
    }
}
