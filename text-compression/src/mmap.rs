//! Protein text decoded straight out of a memory mapping.
//!
//! The counterpart to [`crate::preloaded`]: same 5-bit packing, same file, but nothing is copied
//! into owned memory. The kernel decides what stays resident, which is what makes serving a
//! multi-gigabyte index in a bounded RSS possible.
//!
//! Only the reading half lives here; the file itself is written by `preloaded`'s `WriteBinary`.

use std::{error::Error, fs::File, path::Path, sync::Arc};

use binary_traits::ReadBinaryMmap;
use memmap2::Mmap;

use crate::{BIT5_TO_CHAR, ProteinTextBackend, bit_array_byte_size};

/// Page size assumed when warming a mapping. Touching one byte per this many bytes is enough to
/// fault in every page; a larger real page size just means some touches are redundant.
const ASSUMED_PAGE_SIZE: usize = 4096;

/// Reads every page of `mmap[range]` into the page cache.
///
/// Warmup only — never called per query. Serving from a cold mapping means the first requests pay
/// the page faults; `sa-benchmarks` sweeps every structure before it times anything, which is what
/// makes a throughput figure a steady-state one. `sa-server` does **not** call this today, so a
/// freshly started server warms itself on live traffic.
///
/// All three steps matter:
///
/// 1. [`memmap2::Advice::Sequential`] tells the kernel to read far ahead, so the sweep faults in long runs
///    instead of one page at a time.
/// 2. Touching one byte per page forces the fault. The read **must** be laundered through
///    [`std::hint::black_box`]: without it the optimizer deletes a loop whose result is unused,
///    and the warmup silently does nothing.
/// 3. [`memmap2::Advice::Random`] restores the steady-state pattern. The index is probed in an order the
///    kernel cannot predict, so leaving readahead enabled would make every later miss drag in
///    neighbouring pages that will not be used.
///
/// `range` is a byte range into `mmap`; callers pass only their own section, so a structure
/// sharing a file with others does not warm its neighbours.
///
/// Returns the number of bytes swept. Every caller passes it back up so the benchmark harness can
/// divide it by the elapsed time: a sweep running at disk bandwidth and one running at memcpy
/// bandwidth do the same work and take an order of magnitude apart, and without the byte count the
/// two are indistinguishable in a report.
///
/// This lives here because `sa-index` and `sa-mappings` both need it and both already depend on
/// this crate; it previously existed as five near-identical copies.
pub fn touch_all_pages(mmap: &Mmap, range: std::ops::Range<usize>) -> u64 {
    #[cfg(unix)]
    let _ = mmap.advise(memmap2::Advice::Sequential);

    let swept = mmap[range.clone()].len() as u64;
    for chunk in mmap[range].chunks(ASSUMED_PAGE_SIZE) {
        std::hint::black_box(chunk[0]);
    }

    #[cfg(unix)]
    let _ = mmap.advise(memmap2::Advice::Random);

    swept
}

// ── MmapBackedProteinText ─────────────────────────────────────────────────────

/// Protein text borrowed from a memory mapping.
///
/// The mapping is shared (`Arc`) because one index file holds several structures — the text and
/// the protein table live in the same `proteins.bin` — and each borrows the same mapping rather
/// than opening the file twice.
pub struct MmapBackedProteinText {
    pub(crate) mmap: Arc<Mmap>,
    pub(crate) data_offset: usize,
    pub(crate) len: usize
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

    /// `#[inline]` for the same reason as [`Self::get`], and the preloaded backend carries it too:
    /// `Searcher::compare` reads this once per binary-search probe, so without cross-crate LTO it
    /// would be a real call per probe to return a field.
    #[inline]
    fn len(&self) -> usize {
        self.len
    }

    fn touch_all_pages(&self) -> u64 {
        // Only this text's own section of the mapping: `proteins.bin` holds the metadata after it,
        // and whoever owns that warms it — or deliberately does not.
        // Infallible here: `self.len` was validated against the mapping when the file was mapped,
        // so it cannot be the overflowing header this function guards against.
        let end = self.data_offset
            + bit_array_byte_size(self.len).expect("text length was validated when the file was mapped");
        touch_all_pages(&self.mmap, self.data_offset..end)
    }

    #[inline]
    fn prefetch_at(&self, index: usize) {
        let bit_off = self.data_offset + (index * 5) / 8;
        if bit_off < self.mmap.len() {
            prefetch::prefetch_read(&self.mmap[bit_off] as *const u8);
        }
    }
}

// No `WriteBinary` here, deliberately. Serialisation is the preloaded half's job for every
// structure in the index — see `binary_traits` — and `sa-builder` names only preloaded types
// because of it. This module used to carry one anyway, the only mmap type in the workspace that
// did; nothing in the workspace called it, tests included.

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
            return Err("The protein text file is too small to contain the text header".into());
        }

        let text_length =
            u64::from_le_bytes(mmap[0..8].try_into().map_err(|_| "Failed to parse ProteinText header")?) as usize;

        let data_bytes =
            bit_array_byte_size(text_length).ok_or("The protein text header declares an implausible text length")?;
        if mmap.len() < 8 + data_bytes {
            return Err("The protein text file is too small to contain the text data its header declares".into());
        }

        Ok(Self::from_mmap(mmap, 8, text_length))
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use binary_traits::WriteBinary;
    use memmap2::Mmap;

    use super::*;
    use crate::preloaded::InMemoryProteinText;

    fn write_protein_text_file(input: &str) -> tempfile::NamedTempFile {
        use std::collections::HashMap;

        use bitarray::{Binary, BitArray};
        let char_to_5bit: HashMap<u8, u8> = BIT5_TO_CHAR.iter().enumerate().map(|(i, &c)| (c, i as u8)).collect();
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
                mapped.get(i),
                preloaded.get(i),
                "backends disagree at index {i} (alignment {} in word {})",
                (i * 5) % 64,
                (i * 5) / 64
            );
        }
    }

    /// The mapped half of the hint contract; the preloaded half is the twin of this test.
    #[test]
    fn prefetch_hints_are_harmless() {
        let input = "ACACA-CAC$MLPGLALLLL$";
        let mut buf: Vec<u8> = Vec::new();
        InMemoryProteinText::from_string(input).write_binary(&mut buf).unwrap();
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, &buf).unwrap();
        std::io::Write::flush(&mut tmp).unwrap();

        let mapped = MmapBackedProteinText::read_binary_mmap(tmp.path()).unwrap();
        crate::test_utils::assert_prefetch_is_harmless(&mapped, input);
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
