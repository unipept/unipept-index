//! Protein text decoded straight out of a memory mapping.
//!
//! The counterpart to [`crate::preloaded`]: same 5-bit packing, same file, but nothing is copied
//! into owned memory. The kernel decides what stays resident, which is what makes serving a
//! multi-gigabyte index in a bounded RSS possible.
//!
//! Only the reading half lives here; the file itself is written by `preloaded`'s `WriteBinary`.

use std::error::Error;
use std::io::Write;
use std::sync::Arc;
use std::{fs::File, path::Path};
use binary_traits::{WriteBinary, ReadBinaryMmap};

use memmap2::Mmap;

use crate::{bit_array_byte_size, BIT5_TO_CHAR, ProteinTextBackend};

// ── MmapBackedProteinText ─────────────────────────────────────────────────────

/// Protein text borrowed from a memory mapping.
///
/// The mapping is shared (`Arc`) because one index file holds several structures — the text and
/// the protein table live in the same `proteins.bin` — and each borrows the same mapping rather
/// than opening the file twice.
pub struct MmapBackedProteinText {
    pub(crate) mmap: Arc<Mmap>,
    pub(crate) data_offset: usize,
    pub(crate) len: usize,
}

impl MmapBackedProteinText {
    /// Borrows `len` residues packed at 5 bits each, starting `data_offset` bytes into `mmap`.
    ///
    /// The caller is responsible for having validated that the mapping is long enough; see
    /// `read_binary_mmap`, or `sa-mappings`, which builds one of these over a shared mapping
    /// after doing its own bounds checks.
    pub fn from_mmap(mmap: Arc<Mmap>, data_offset: usize, len: usize) -> Self {
        Self { mmap, data_offset, len }
    }
}

impl ProteinTextBackend for MmapBackedProteinText {
    /// Decodes the residue at `index`.
    ///
    /// # Why this re-implements `bitarray`'s unpacking instead of reusing it
    ///
    /// This is deliberate, not copy-paste. `bitarray` indexes a `&[u64]`, which requires the
    /// packed data to be an aligned slice of `u64`. Here it is a byte range inside a mapping at
    /// an arbitrary `data_offset`, so reinterpreting it as `&[u64]` would be unaligned — and
    /// making it aligned would mean either copying the data (defeating the point of mmapping it)
    /// or constraining the file layout to keep every structure 8-byte aligned.
    ///
    /// Reading through `u64::from_le_bytes` instead costs nothing: it compiles to the same load
    /// on both supported architectures, and it fuses the `BIT5_TO_CHAR` lookup into the same
    /// function so the hot path stays one call deep.
    ///
    /// The packing must nevertheless stay bit-for-bit identical to `bitarray`'s, since the file
    /// is written through `BitArray<5>`. The parity test in `bitarray` pins the two
    /// implementations there together; `tests::matches_the_preloaded_backend` below pins this one
    /// to them.
    #[inline]
    fn get(&self, index: usize) -> u8 {
        const BITS: usize = 5;
        const MASK: u64 = (1u64 << BITS) - 1;
        let bit_offset = index * BITS;
        let start_block = bit_offset / 64;
        let start_bit = bit_offset % 64;
        let byte_off = self.data_offset + start_block * 8;
        let lo = u64::from_le_bytes(self.mmap[byte_off..byte_off + 8].try_into().unwrap());
        let raw = if start_bit + BITS <= 64 {
            (lo >> (64 - start_bit - BITS)) & MASK
        } else {
            let end_bit = (index + 1) * BITS % 64;
            let hi = u64::from_le_bytes(self.mmap[byte_off + 8..byte_off + 16].try_into().unwrap());
            ((lo << end_bit) | (hi >> (64 - end_bit))) & MASK
        };
        BIT5_TO_CHAR[raw as usize]
    }

    fn len(&self) -> usize { self.len }

    #[inline]
    fn prefetch_at(&self, index: usize) {
        let bit_off = self.data_offset + (index * 5) / 8;
        if bit_off < self.mmap.len() {
            prefetch::prefetch_read(&self.mmap[bit_off] as *const u8);
        }
    }
}

impl WriteBinary for MmapBackedProteinText {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let text_length = self.len as u64;
        writer.write_all(&text_length.to_le_bytes())?;
        let n_bytes = bit_array_byte_size(self.len);
        writer.write_all(&self.mmap[self.data_offset..self.data_offset + n_bytes])?;
        Ok(())
    }
}

impl ReadBinaryMmap for MmapBackedProteinText {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let f = File::open(path)?;
        // SAFETY: `Mmap::map` is unsafe because the mapping aliases a file that another process
        // could modify or truncate underneath us; a shrinking file turns subsequent reads into
        // SIGBUS, and concurrent writes into a data race. Neither applies to an index file: it is
        // written once by `sa-builder` and is thereafter read-only for the lifetime of the
        // server. Deployments must not rebuild an index in place under a running process — swap
        // in a new file and restart instead. This same argument covers every `Mmap::map` in the
        // workspace.
        let mmap = Arc::new(unsafe { Mmap::map(&f)? });

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        if mmap.len() < 8 {
            return Err("File is too small to contain ProteinText header (8 bytes required)".into());
        }

        let text_length = u64::from_le_bytes(mmap[0..8].try_into()
            .map_err(|_| "Failed to parse ProteinText header")?) as usize;

        if mmap.len() < 8 + bit_array_byte_size(text_length) {
            return Err("File is too small to contain ProteinText BitArray data for declared length".into());
        }

        Ok(Self::from_mmap(mmap, 8, text_length))
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use memmap2::Mmap;
    use super::*;
    use crate::preloaded::InMemoryProteinText;

    fn write_protein_text_file(input: &str) -> tempfile::NamedTempFile {
        use bitarray::{Binary, BitArray};
        use std::collections::HashMap;
        let char_to_5bit: HashMap<u8, u8> = BIT5_TO_CHAR.iter().enumerate()
            .map(|(i, &c)| (c, i as u8)).collect();
        let mut ba = BitArray::<5>::with_capacity(input.len());
        for (i, c) in input.bytes().enumerate() {
            ba.set(i, *char_to_5bit.get(&c).unwrap() as u64);
        }
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(&(input.len() as u64).to_le_bytes());
        ba.write_binary(&mut buf).unwrap();
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, &buf).unwrap();
        std::io::Write::flush(&mut tmp).unwrap();
        tmp
    }

    #[test]
    fn test_mmap_roundtrip() {
        let input = "ACACA-CAC$MLPGLALLLL$";
        let tmp = write_protein_text_file(input);
        let f = std::fs::File::open(tmp.path()).unwrap();
        let mmap = Arc::new(unsafe { Mmap::map(&f).unwrap() });
        let text_len = u64::from_le_bytes(mmap[0..8].try_into().unwrap()) as usize;
        let mmap_text = MmapBackedProteinText::from_mmap(Arc::clone(&mmap), 8, text_len);
        assert_eq!(mmap_text.len(), input.len());
        for (i, c) in input.chars().enumerate() {
            assert_eq!(mmap_text.get(i), c as u8, "mismatch at index {}", i);
        }
    }

    #[test]
    fn test_mmap_block_boundary() {
        // 13 characters: 13*5=65 bits, crosses a u64 boundary
        let input = "ABCDEFGHIKLMN";
        let tmp = write_protein_text_file(input);
        let mmap = Arc::new(unsafe {
            let f = std::fs::File::open(tmp.path()).unwrap();
            Mmap::map(&f).unwrap()
        });
        let text_len = u64::from_le_bytes(mmap[0..8].try_into().unwrap()) as usize;
        let mmap_text = MmapBackedProteinText::from_mmap(Arc::clone(&mmap), 8, text_len);
        assert_eq!(mmap_text.len(), input.len());
        for (i, c) in input.chars().enumerate() {
            assert_eq!(mmap_text.get(i), c as u8, "boundary mismatch at index {}", i);
        }
    }

    #[test]
    fn test_mmap_write_and_read_binary_mmap() {
        let input = "ACACA-CAC$";
        let text = InMemoryProteinText::from_string(input);
        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, &buf).unwrap();
        std::io::Write::flush(&mut tmp).unwrap();
        let mmap_text = MmapBackedProteinText::read_binary_mmap(tmp.path()).unwrap();
        assert_eq!(mmap_text.len(), input.len());
        for (i, c) in input.chars().enumerate() {
            assert_eq!(mmap_text.get(i), c as u8);
        }
    }

    /// The two backends decode the same file, so any divergence in the 5-bit unpacking — which
    /// this module deliberately reimplements rather than reusing `bitarray` — silently changes
    /// what the server returns. Nothing else pins them together.
    #[test]
    fn matches_the_preloaded_backend() {
        // Long enough to cross many word boundaries, and covering the whole alphabet including
        // both delimiters, so every 5-bit code and every intra-word alignment is exercised.
        let input: String = std::iter::repeat_n("ABCDEFGHIKLMNOPQRSTUVWXYZ-$", 9).collect();

        let preloaded = InMemoryProteinText::from_string(&input);
        let mut buf: Vec<u8> = Vec::new();
        InMemoryProteinText::from_string(&input).write_binary(&mut buf).unwrap();
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, &buf).unwrap();
        std::io::Write::flush(&mut tmp).unwrap();
        let mapped = MmapBackedProteinText::read_binary_mmap(tmp.path()).unwrap();

        assert_eq!(mapped.len(), preloaded.len());
        for i in 0..input.len() {
            assert_eq!(
                mapped.get(i), preloaded.get(i),
                "backends disagree at index {i} (alignment {} in word {})", (i * 5) % 64, (i * 5) / 64
            );
        }
    }

    /// `read_binary_mmap` takes an untrusted header. Both of its length checks must produce an
    /// error rather than a panic or an out-of-bounds read; neither was covered before.
    #[test]
    fn truncated_files_error_rather_than_panicking() {
        let text = InMemoryProteinText::from_string("ACACA-CAC$MLPGLALLLL$");
        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();

        for cut in 0..buf.len() {
            let mut tmp = tempfile::NamedTempFile::new().unwrap();
            std::io::Write::write_all(&mut tmp, &buf[..cut]).unwrap();
            std::io::Write::flush(&mut tmp).unwrap();

            let err = MmapBackedProteinText::read_binary_mmap(tmp.path())
                .err()
                .unwrap_or_else(|| panic!("{cut} of {} bytes should not load", buf.len()));
            let msg = err.to_string();
            assert!(msg.contains("too small"), "unexpected error at {cut} bytes: {msg}");
        }
    }

    /// A header claiming more text than the file holds must be rejected, not trusted.
    #[test]
    fn overlong_declared_length_is_rejected() {
        let text = InMemoryProteinText::from_string("ACACA-CAC$");
        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();
        buf[0..8].copy_from_slice(&1_000_000_u64.to_le_bytes());

        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, &buf).unwrap();
        std::io::Write::flush(&mut tmp).unwrap();

        assert!(MmapBackedProteinText::read_binary_mmap(tmp.path()).is_err());
    }
}
