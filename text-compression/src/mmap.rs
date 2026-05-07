// This entire file is mmap-only.
use std::error::Error;
use std::io::Write;
use std::sync::Arc;
use std::{fs::File, path::Path};

use memmap2::Mmap;

use crate::traits::{WriteBinary, ReadBinaryMmap};
use crate::bit_array_byte_size;

const BIT5_TO_CHAR: &[u8; 27] = b"ABCDEFGHIKLMNOPQRSTUVWXYZ-$";

// ── MmapProteinText ───────────────────────────────────────────────────────────

pub struct MmapProteinText {
    pub(crate) mmap: Arc<Mmap>,
    pub(crate) data_offset: usize,
    pub(crate) len: usize,
}

impl MmapProteinText {
    pub fn from_mmap(mmap: Arc<Mmap>, data_offset: usize, len: usize) -> Self {
        Self { mmap, data_offset, len }
    }

    fn get_mmap(mmap: &Mmap, data_offset: usize, index: usize) -> u8 {
        const BITS: usize = 5;
        const MASK: u64 = (1u64 << BITS) - 1;
        let bit_offset = index * BITS;
        let start_block = bit_offset / 64;
        let start_bit = bit_offset % 64;
        let byte_off = data_offset + start_block * 8;
        let lo = u64::from_le_bytes(mmap[byte_off..byte_off + 8].try_into().unwrap());
        let raw = if start_bit + BITS <= 64 {
            (lo >> (64 - start_bit - BITS)) & MASK
        } else {
            let end_bit = (index + 1) * BITS % 64;
            let hi = u64::from_le_bytes(mmap[byte_off + 8..byte_off + 16].try_into().unwrap());
            ((lo << end_bit) | (hi >> (64 - end_bit))) & MASK
        };
        BIT5_TO_CHAR[raw as usize]
    }

    #[inline]
    pub fn get(&self, index: usize) -> u8 {
        Self::get_mmap(&self.mmap, self.data_offset, index)
    }

    pub fn len(&self) -> usize { self.len }
    pub fn is_empty(&self) -> bool { self.len == 0 }

    pub fn iter(&self) -> crate::ProteinTextIterator<'_> {
        crate::ProteinTextIterator { protein_text: self, index: 0 }
    }

    pub fn slice(&self, start: usize, end: usize) -> crate::ProteinTextSlice<'_> {
        crate::ProteinTextSlice::new(self, start, end)
    }

    #[inline]
    pub fn prefetch_at(&self, index: usize) {
        let bit_off = self.data_offset + (index * 5) / 8;
        if bit_off < self.mmap.len() {
            prefetch::prefetch_read(&self.mmap[bit_off] as *const u8);
        }
    }
}

impl WriteBinary for MmapProteinText {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let text_length = self.len as u64;
        writer.write_all(&text_length.to_le_bytes())?;
        let n_bytes = bit_array_byte_size(self.len);
        writer.write_all(&self.mmap[self.data_offset..self.data_offset + n_bytes])?;
        Ok(())
    }
}

impl ReadBinaryMmap for MmapProteinText {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let f = File::open(path)?;
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
    use crate::traits::WriteBinary as _;

    fn write_protein_text_file(input: &str) -> tempfile::NamedTempFile {
        use bitarray::{Binary, BitArray};
        use std::collections::HashMap;
        let char_to_5bit: HashMap<u8, u8> = "ABCDEFGHIKLMNOPQRSTUVWXYZ-$"
            .chars().enumerate().map(|(i, c)| (c as u8, i as u8)).collect();
        let mut ba = BitArray::with_capacity(input.len(), 5);
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
        let mmap_text = MmapProteinText::from_mmap(Arc::clone(&mmap), 8, text_len);
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
        let mmap_text = MmapProteinText::from_mmap(Arc::clone(&mmap), 8, text_len);
        assert_eq!(mmap_text.len(), input.len());
        for (i, c) in input.chars().enumerate() {
            assert_eq!(mmap_text.get(i), c as u8, "boundary mismatch at index {}", i);
        }
    }
}
