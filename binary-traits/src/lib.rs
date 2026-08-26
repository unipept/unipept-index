//! The I/O traits every on-disk index structure is written and read through.
//!
//! An index is built once by `sa-builder` and then read back through one of two backends, so most
//! structures have one writer and two readers. No type implements all four traits: a structure
//! takes either the owned route ([`WriteBinary`] + [`ReadBinary`] + [`LoadIndex`]) or the mapped
//! one ([`ReadBinaryMmap`] + [`LoadIndex`]).
//!
//! * [`WriteBinary`] — serialise. Implemented once per structure, by the preloaded type, and used
//!   by `sa-builder` regardless of which backend will later read the file.
//! * [`ReadBinary`] — deserialise into owned memory. The preloaded backend.
//! * [`ReadBinaryMmap`] — map the file and decode fields in place. The mmap backend.
//! * [`LoadIndex`] — load by whichever of the two routes that concrete type uses, so that code
//!   generic over the backends does not have to know which.
//!
//! Because the writer lives with the preloaded type but its output is also consumed by an mmap
//! reader in a different module, the two are easy to drift apart. Every format is therefore
//! documented at its *writer*, and each reader points back at it. `sa-builder` names only the
//! preloaded types for exactly this reason: one writer per structure, whichever backend reads it
//! later.
//!
//! These traits live in their own crate so that `sa-index`, `protein-metadata` and `text-compression`
//! can implement them for each other's types without a dependency cycle. Each of those crates
//! depends on this one directly and names the traits as `binary_traits::…`, so that the crate a
//! trait is imported from is the crate that defines it.
#![warn(missing_docs)]

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path
};

/// Serialises a structure to a byte stream.
///
/// Takes `self` by value: writing is the last thing done with a built structure, and consuming it
/// lets implementations move large buffers out instead of cloning them.
///
/// Implementations must use [`Write::write_all`] rather than [`Write::write`] — a short write
/// would silently truncate the file and produce an index that loads without complaint.
pub trait WriteBinary {
    /// Writes `self` to `writer`, consuming it.
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>>;
}

/// Reads a structure written by [`WriteBinary`] into owned memory.
///
/// Reads exactly the bytes the corresponding `write_binary` produced and leaves the reader
/// positioned immediately after them, so that composite structures can chain several
/// `read_binary` calls over one stream.
pub trait ReadBinary: Sized {
    /// Reads one structure from `reader`, consuming exactly the bytes `write_binary` wrote.
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>>;
}

/// Maps a file written by [`WriteBinary`] and decodes it in place.
///
/// The counterpart to [`ReadBinary`] for the mmap backend: instead of copying the payload into
/// owned buffers, implementations keep the mapping alive and compute field offsets into it.
///
/// # Contract
///
/// The header is untrusted input. Implementations **must** validate every length and offset
/// against the mapping's actual size before slicing, and return an error rather than panicking on
/// a truncated or corrupt file. A reader that trusts the header turns a damaged index into a
/// panic inside a request handler.
///
/// Every implementation also has an `unsafe { Mmap::map(..) }` call to justify: the mapped file
/// must not be modified or truncated for the lifetime of the mapping, or reads become undefined
/// behaviour. That argument is written out once, at the mapping in `text_compression::mmap`, and
/// covers every `Mmap::map` in the workspace — it turns on how an index file is *replaced*
/// (rename safe, `O_TRUNC` not), which is a deployment property rather than a per-structure one.
/// New implementations should point at that note rather than restate it, and must not weaken it.
pub trait ReadBinaryMmap: Sized {
    /// Maps the file at `path` and decodes its header, borrowing the payload from the mapping.
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>>;
}

/// Loads a structure from an index file, by whichever of the two routes above it uses.
///
/// [`ReadBinary`] and [`ReadBinaryMmap`] describe *how* a structure is decoded; this describes how
/// to get one from a path without knowing which applies. Each concrete type has exactly one route
/// and therefore exactly one implementation: an owned type opens the file and streams it, a mapped
/// one maps it.
///
/// # Why this exists
///
/// A bound cannot say "implements `ReadBinary` *or* `ReadBinaryMmap`", so without this trait no
/// function can be generic over the storage backends — which is what the searcher's test fixtures
/// need in order to exercise every combination in one build.
///
/// It also puts a fact in the type system that used to be a hand-maintained `#[cfg]` predicate:
/// `proteins.bin` holds both the protein text and the metadata table, so it must be *mapped*
/// whenever either section is mapped, not merely when the metadata is. That is now expressed by
/// which of the four `InMemory`/`MmapBacked` protein pairings implements this trait through which
/// route; see `protein_metadata::mmap`.
pub trait LoadIndex: Sized {
    /// Loads the structure stored at `path`.
    fn load(path: &Path) -> Result<Self, Box<dyn Error>>;
}

/// Opens `path` and reads one structure out of it, for [`LoadIndex`] impls that take the owned
/// route.
///
/// Every such impl is a one-line call to this, so the buffering is decided once here rather than
/// per structure.
///
/// That makes the `BufReader` capacity — [`BufReader::new`]'s default 8 KiB — a shared constant
/// three readers are calibrated against: the `READ_CHUNK_BYTES` in `bitarray::binary`,
/// `sa_index::array::preloaded::original` and `sa_index::suffix_to_protein_index::preloaded::bitvec`
/// are all sized *above* it so that `BufReader::read` bypasses its internal buffer instead of
/// copying every byte twice. Switching this to `BufReader::with_capacity` is therefore not a local
/// change: raising the capacity past those constants silently restores the double copy. The
/// measurement behind that is at `bitarray::binary`'s constant.
pub fn load_owned<T: ReadBinary>(path: &Path) -> Result<T, Box<dyn Error>> {
    T::read_binary(&mut BufReader::new(File::open(path)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Toy structure with the same shape as the real ones: a fixed-size header followed by a
    /// variable-length payload.
    #[derive(Debug, PartialEq)]
    struct Toy {
        tag: u8,
        payload: Vec<u8>
    }

    impl WriteBinary for Toy {
        fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
            writer.write_all(&[self.tag])?;
            writer.write_all(&(self.payload.len() as u64).to_le_bytes())?;
            writer.write_all(&self.payload)?;
            Ok(())
        }
    }

    impl ReadBinary for Toy {
        fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
            let mut tag = [0u8; 1];
            reader.read_exact(&mut tag)?;
            let mut len = [0u8; 8];
            reader.read_exact(&mut len)?;
            let mut payload = vec![0u8; u64::from_le_bytes(len) as usize];
            reader.read_exact(&mut payload)?;
            Ok(Self { tag: tag[0], payload })
        }
    }

    #[test]
    fn round_trips_through_a_byte_stream() {
        let original = Toy { tag: 7, payload: b"suffix array".to_vec() };

        let mut buf = Vec::new();
        Toy { tag: 7, payload: b"suffix array".to_vec() }.write_binary(&mut buf).unwrap();

        assert_eq!(Toy::read_binary(&mut buf.as_slice()).unwrap(), original);
    }

    /// The contract that makes composite structures work: `read_binary` consumes exactly what
    /// `write_binary` produced, so several of them chain over a single stream. `InMemoryProteins`
    /// relies on this to read its text and then its protein table from one file.
    #[test]
    fn consumes_exactly_its_own_bytes_so_reads_chain() {
        let mut buf = Vec::new();
        Toy { tag: 1, payload: b"first".to_vec() }.write_binary(&mut buf).unwrap();
        Toy { tag: 2, payload: Vec::new() }.write_binary(&mut buf).unwrap();
        Toy { tag: 3, payload: b"third".to_vec() }.write_binary(&mut buf).unwrap();

        let mut reader = buf.as_slice();
        assert_eq!(Toy::read_binary(&mut reader).unwrap().payload, b"first");
        assert_eq!(Toy::read_binary(&mut reader).unwrap().payload, b"");
        assert_eq!(Toy::read_binary(&mut reader).unwrap().tag, 3);
        assert!(reader.is_empty(), "readers must not over-consume");
    }

    #[test]
    fn truncated_input_errors_rather_than_panicking() {
        let mut buf = Vec::new();
        Toy { tag: 9, payload: b"payload".to_vec() }.write_binary(&mut buf).unwrap();

        for cut in 0..buf.len() {
            assert!(
                Toy::read_binary(&mut &buf[..cut]).is_err(),
                "truncating to {cut} of {} bytes should error",
                buf.len()
            );
        }
    }
}
