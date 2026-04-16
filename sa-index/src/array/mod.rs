use std::{
    fs::File,
    io::{BufRead, Write},
    path::Path
};

use bitarray::{Binary, BitArray};
use memmap2::Mmap;
use text_compression::{WriteBinary, ReadBinary, ReadBinaryMmap};

pub mod original;
pub mod compressed;
pub mod mmap;

pub use original::dump_suffix_array;
pub use compressed::{dump_compressed_suffix_array, load_compressed_suffix_array};

/// Represents a suffix array.
pub enum SuffixArray {
    /// The original suffix array.
    Original(Vec<i64>, u8),
    /// The compressed suffix array.
    Compressed(BitArray, u8),
    /// A suffix array backed by a memory-mapped file. Works for both compressed and uncompressed
    /// formats: bits_per_value == 64 means uncompressed (i64 values), otherwise compressed
    /// (BitArray-style packed bit values with the given bits per element).
    MmapBacked {
        mmap: Mmap,
        data_offset: usize,
        len: usize,
        bits_per_value: usize,
        sample_rate: u8
    }
}

impl SuffixArray {
    /// Returns the length of the suffix array.
    pub fn len(&self) -> usize {
        match self {
            SuffixArray::Original(sa, _) => sa.len(),
            SuffixArray::Compressed(sa, _) => sa.len(),
            SuffixArray::MmapBacked { len, .. } => *len
        }
    }

    /// Returns the number of bits per value in the suffix array.
    pub fn bits_per_value(&self) -> usize {
        match self {
            SuffixArray::Original(_, _) => 64,
            SuffixArray::Compressed(sa, _) => sa.bits_per_value(),
            SuffixArray::MmapBacked { bits_per_value, .. } => *bits_per_value
        }
    }

    /// Returns the sample rate used for the suffix array.
    pub fn sample_rate(&self) -> u8 {
        match self {
            SuffixArray::Original(_, sample_rate) => *sample_rate,
            SuffixArray::Compressed(_, sample_rate) => *sample_rate,
            SuffixArray::MmapBacked { sample_rate, .. } => *sample_rate
        }
    }

    /// Returns the suffix array value at the given index.
    pub fn get(&self, index: usize) -> i64 {
        match self {
            SuffixArray::Original(sa, _) => sa[index],
            SuffixArray::Compressed(sa, _) => sa.get(index) as i64,
            SuffixArray::MmapBacked { mmap, data_offset, bits_per_value, .. } => {
                mmap::get_mmap(mmap, *data_offset, *bits_per_value, index)
            }
        }
    }

    /// Returns whether the suffix array is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Issues an OS prefetch hint (`MADV_WILLNEED`) for the mmap pages covering SA indices
    /// `lo..hi_exclusive`.  No-op for in-memory variants and on non-Unix platforms.
    #[inline]
    pub fn prefetch_sa_range(&self, lo: usize, hi_exclusive: usize) {
        #[cfg(unix)]
        if let SuffixArray::MmapBacked { mmap, data_offset, bits_per_value, .. } = self {
            let byte_lo = data_offset + (lo * bits_per_value) / 8;
            let byte_hi = data_offset + (hi_exclusive * bits_per_value).div_ceil(8);
            let len = byte_hi.saturating_sub(byte_lo);
            if len > 0 && byte_hi <= mmap.len() {
                // SAFETY: MADV_WILLNEED is a non-destructive, read-only prefetch hint.
                let _ = mmap.advise_range(memmap2::Advice::WillNeed, byte_lo, len);
            }
        }
    }
}

impl WriteBinary for SuffixArray {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            SuffixArray::Original(sa, sparseness_factor) => {
                original::dump_suffix_array(sa, sparseness_factor, writer)
            }
            SuffixArray::Compressed(bit_array, sample_rate) => {
                writer.write_all(&[bit_array.bits_per_value() as u8])
                    .map_err(|_| "Could not write the required bits")?;
                writer.write_all(&[sample_rate])
                    .map_err(|_| "Could not write the sparseness factor")?;
                writer.write_all(&(bit_array.len() as u64).to_le_bytes())
                    .map_err(|_| "Could not write the size")?;
                bit_array.write_binary(writer)
                    .map_err(|_| "Could not write the compressed suffix array")?;
                Ok(())
            }
            SuffixArray::MmapBacked { .. } => {
                Err("WriteBinary is not supported for SuffixArray::MmapBacked".into())
            }
        }
    }
}

impl ReadBinary for SuffixArray {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn std::error::Error>> {
        let mut buf1 = [0u8; 1];
        reader.read_exact(&mut buf1).map_err(|_| "Could not read the required bits from the binary file")?;
        let bits_per_value = buf1[0] as usize;

        reader.read_exact(&mut buf1).map_err(|_| "Could not read the sample rate from the binary file")?;
        let sample_rate = buf1[0];

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8).map_err(|_| "Could not read the size of the suffix array from the binary file")?;
        let size = u64::from_le_bytes(buf8) as usize;

        if bits_per_value == 64 {
            let sa = original::load_original(reader, sample_rate, size)?;
            Ok(SuffixArray::Original(sa, sample_rate))
        } else {
            let sa = compressed::load_compressed(reader, bits_per_value, size)?;
            Ok(SuffixArray::Compressed(sa, sample_rate))
        }
    }
}

impl ReadBinaryMmap for SuffixArray {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        if mmap.len() < 10 {
            return Err("The binary file is too small to contain the SA header".into());
        }

        let bits_per_value = mmap[0] as usize;
        let sample_rate = mmap[1];
        let amount_of_items = u64::from_le_bytes(mmap[2..10].try_into()?) as usize;

        let header_bytes = 10usize;
        let total_bits = amount_of_items
            .checked_mul(bits_per_value)
            .ok_or("The SA header declares too many items or bits per value")?;
        let data_bytes = total_bits.div_ceil(64);

        if mmap.len() < header_bytes + data_bytes {
            return Err("The binary file is too small to contain the SA data".into());
        }

        Ok(SuffixArray::MmapBacked {
            mmap,
            data_offset: 10,
            len: amount_of_items,
            bits_per_value,
            sample_rate
        })
    }
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead, Read};

    use bitarray::BitArray;

    use super::*;
    use crate::ReadBinaryMmap;

    struct FailingReader {
        valid_read_count: usize
    }

    impl Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.valid_read_count == 0 {
                return Err(std::io::Error::other("Read failed"));
            }
            self.valid_read_count -= 1;
            Ok(buf.len())
        }
    }

    impl BufRead for FailingReader {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
            Ok(&[])
        }

        fn consume(&mut self, _: usize) {}
    }

    #[test]
    fn test_suffix_array_original() {
        let sa = SuffixArray::Original(vec![1, 2, 3, 4, 5], 1);
        assert_eq!(sa.len(), 5);
        assert_eq!(sa.get(0), 1);
        assert_eq!(sa.get(1), 2);
        assert_eq!(sa.get(2), 3);
        assert_eq!(sa.get(3), 4);
        assert_eq!(sa.get(4), 5);
    }

    #[test]
    fn test_suffix_array_compressed() {
        let mut bitarray = BitArray::with_capacity(5, 40);
        bitarray.set(0, 1_u64);
        bitarray.set(1, 2_u64);
        bitarray.set(2, 3_u64);
        bitarray.set(3, 4_u64);
        bitarray.set(4, 5_u64);

        let sa = SuffixArray::Compressed(bitarray, 1);
        assert_eq!(sa.len(), 5);
        assert_eq!(sa.get(0), 1);
        assert_eq!(sa.get(1), 2);
        assert_eq!(sa.get(2), 3);
        assert_eq!(sa.get(3), 4);
        assert_eq!(sa.get(4), 5);
    }

    #[test]
    fn test_suffix_array_len() {
        let sa = SuffixArray::Original(vec![1, 2, 3, 4, 5], 1);
        assert_eq!(sa.len(), 5);

        let bitarray = BitArray::with_capacity(5, 40);
        let sa = SuffixArray::Compressed(bitarray, 1);
        assert_eq!(sa.len(), 5);
    }

    #[test]
    fn test_suffix_array_bits_per_value() {
        let sa = SuffixArray::Original(vec![1, 2, 3, 4, 5], 1);
        assert_eq!(sa.bits_per_value(), 64);

        let bitarray = BitArray::with_capacity(5, 40);
        let sa = SuffixArray::Compressed(bitarray, 1);
        assert_eq!(sa.bits_per_value(), 40);
    }

    #[test]
    fn test_suffix_array_sample_rate() {
        let sa = SuffixArray::Original(vec![1, 2, 3, 4, 5], 1);
        assert_eq!(sa.sample_rate(), 1);

        let bitarray = BitArray::with_capacity(5, 40);
        let sa = SuffixArray::Compressed(bitarray, 1);
        assert_eq!(sa.sample_rate(), 1);
    }

    #[test]
    fn test_suffix_array_is_empty() {
        let sa = SuffixArray::Original(vec![], 1);
        assert!(sa.is_empty());

        let bitarray = BitArray::with_capacity(0, 0);
        let sa = SuffixArray::Compressed(bitarray, 1);
        assert!(sa.is_empty());
    }

    #[test]
    fn test_load_suffix_array() {
        let buffer = vec![
            64, 1, 5, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
        ];

        let mut reader = buffer.as_slice();
        let sa = SuffixArray::read_binary(&mut reader).unwrap();

        assert_eq!(sa.sample_rate(), 1);
        for i in 0..5 {
            assert_eq!(sa.get(i), i as i64 + 1);
        }
    }

    #[test]
    #[should_panic(expected = "Could not read the required bits from the binary file")]
    fn test_load_suffix_array_fail_sample_rate() {
        let mut reader = FailingReader { valid_read_count: 0 };
        SuffixArray::read_binary(&mut reader).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the sample rate from the binary file")]
    fn test_load_suffix_array_fail_size() {
        let mut reader = FailingReader { valid_read_count: 1 };
        SuffixArray::read_binary(&mut reader).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the size of the suffix array from the binary file")]
    fn test_load_suffix_array_fail_suffix_array() {
        let mut reader = FailingReader { valid_read_count: 2 };
        SuffixArray::read_binary(&mut reader).unwrap();
    }

    #[test]
    fn test_compressed_write_binary_roundtrip() {
        let mut bitarray = BitArray::with_capacity(5, 40);
        bitarray.set(0, 10_u64);
        bitarray.set(1, 20_u64);
        bitarray.set(2, 30_u64);
        bitarray.set(3, 40_u64);
        bitarray.set(4, 50_u64);

        let sa = SuffixArray::Compressed(bitarray, 3);
        let mut buf = Vec::new();
        sa.write_binary(&mut buf).unwrap();

        let mut reader = std::io::BufReader::new(buf.as_slice());
        let restored = SuffixArray::read_binary(&mut reader).unwrap();

        assert_eq!(restored.bits_per_value(), 40);
        assert_eq!(restored.sample_rate(), 3);
        assert_eq!(restored.len(), 5);
        assert_eq!(restored.get(0), 10);
        assert_eq!(restored.get(1), 20);
        assert_eq!(restored.get(2), 30);
        assert_eq!(restored.get(3), 40);
        assert_eq!(restored.get(4), 50);
    }

    #[test]
    fn test_load_suffix_array_mmap_uncompressed() {
        use tempdir::TempDir;

        let tmp = TempDir::new("mmap_test").unwrap();
        let path = tmp.path().join("sa.bin");

        let sa = vec![1_i64, 2, 3, 4, 5];
        let mut file = std::fs::File::create(&path).unwrap();
        dump_suffix_array(sa, 3, &mut file).unwrap();
        drop(file);

        let loaded = SuffixArray::read_binary_mmap(&path).unwrap();

        assert_eq!(loaded.bits_per_value(), 64);
        assert_eq!(loaded.sample_rate(), 3);
        assert_eq!(loaded.len(), 5);
        for i in 0..5 {
            assert_eq!(loaded.get(i), i as i64 + 1);
        }
    }

    #[test]
    fn test_load_suffix_array_mmap_compressed() {
        use tempdir::TempDir;

        let tmp = TempDir::new("mmap_compressed_test").unwrap();
        let path = tmp.path().join("sa_compressed.bin");

        let sa = vec![1_i64, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut file = std::fs::File::create(&path).unwrap();
        dump_compressed_suffix_array(sa, 2, 40, &mut file).unwrap();
        drop(file);

        let loaded = SuffixArray::read_binary_mmap(&path).unwrap();

        assert_eq!(loaded.bits_per_value(), 40);
        assert_eq!(loaded.sample_rate(), 2);
        assert_eq!(loaded.len(), 10);
        for i in 0..10 {
            assert_eq!(loaded.get(i), i as i64 + 1);
        }
    }
}
