use std::{error::Error, io::{BufRead, Write}};

use bitarray::{data_to_writer, Binary, BitArray};

pub fn dump_compressed_suffix_array<const B: usize>(
    sa: Vec<i64>, 
    sparseness_factor: u8, 
    writer: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    // Write the flags to the writer
    // 00000001 indicates that the suffix array is compressed
    writer.write(&[1_u8]).map_err(|_| "Could not write the flags to the writer")?;

    // Write the sparseness factor to the writer
    writer.write(&[sparseness_factor]).map_err(|_| "Could not write the sparseness factor to the writer")?;

    // Write the size of the suffix array to the writer
    writer.write(&(sa.len() as u64).to_le_bytes()).map_err(|_| "Could not write the size of the suffix array to the writer")?;

    // Compress the suffix array and write it to the writer
    data_to_writer::<B>(sa, writer, 8 * 1024).map_err(|_| "Could not write the compressed suffix array to the writer")?;

    Ok(())
}

pub fn load_compressed_suffix_array<const B: usize>(
    reader: &mut impl BufRead,
) -> Result<(u8, BitArray<B>), Box<dyn Error>> {
    // Read the sample rate from the binary file (1 byte)
    let mut sample_rate_buffer = [0_u8; 1];
    reader.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file")?;
    let sample_rate = sample_rate_buffer[0];

    // Read the size of the suffix array from the binary file (8 bytes)
    let mut size_buffer = [0_u8; 8];
    reader.read_exact(&mut size_buffer).map_err(|_| "Could not read the size of the suffix array from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;

    // Read the compressed suffix array from the binary file
    let mut compressed_suffix_array = BitArray::<B>::with_capacity(size);
    compressed_suffix_array.read_binary(reader).map_err(|_| "Could not read the compressed suffix array from the binary file")?;

    Ok((sample_rate, compressed_suffix_array))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dump_compressed_suffix_array() {
        let sa = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

        let mut writer = vec![];
        dump_compressed_suffix_array::<8>(sa, 1, &mut writer).unwrap();

        assert_eq!(writer, vec![
            // flags
            1,
            // sparseness factor
            1,
            // size of the suffix array
            10, 0, 0, 0, 0, 0, 0, 0,
            // compressed suffix array
            8, 7, 6, 5, 4, 3, 2, 1,
            0, 0, 0, 0, 0, 0, 10, 9
        ]);
    }

    #[test]
    fn test_load_compressed_suffix_array() {
        let data = vec![
            // flags
            1,
            // sparseness factor
            1,
            // size of the suffix array
            10, 0, 0, 0, 0, 0, 0, 0,
            // compressed suffix array
            8, 7, 6, 5, 4, 3, 2, 1,
            0, 0, 0, 0, 0, 0, 10, 9
        ];

        let mut reader = std::io::BufReader::new(&data[..]);
        let (sample_rate, compressed_suffix_array) = load_compressed_suffix_array::<8>(&mut reader).unwrap();

        assert_eq!(sample_rate, 1);
        for i in 0..10 {
            assert_eq!(compressed_suffix_array.get(i), i as u64 + 1);
        }
    }
}
