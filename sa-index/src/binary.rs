use std::{error::Error, io::{BufRead, Read, Write}};

/// The `Binary` trait provides methods for reading and writing a struct as binary.
pub trait Binary {
    /// Writes the struct as binary to the given writer.
    ///
    /// # Arguments
    ///
    /// * `writer` - The writer to write the binary data to.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the write operation is successful, or an `Err` if an error occurs.
    fn write_binary<W: Write>(&self, writer: &mut W) -> std::io::Result<()>;

    /// Reads binary data into a struct from the given reader.
    ///
    /// # Arguments
    ///
    /// * `reader` - The reader to read the binary data from.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read operation is successful, or an `Err` if an error occurs.
    fn read_binary<R: BufRead>(&mut self, reader: R) -> std::io::Result<()>;
}

impl Binary for Vec<i64> {
    fn write_binary<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        for value in self {
            writer.write_all(&value.to_le_bytes())?;
        }

        Ok(())
    }

    fn read_binary<R: BufRead>(&mut self, mut reader: R) -> std::io::Result<()> {
        self.clear();

        let mut buffer = vec![0; 8 * 1024];
 
        loop {
            let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer);
            for buffer_slice in buffer[..bytes_read].chunks_exact(8) {
                self.push(i64::from_le_bytes(buffer_slice.try_into().unwrap()));
            }

            if finished {
                break;
            }
        }

        Ok(())
    }
}

/// Loads the suffix array from the file with the given `filename`
///
/// # Arguments
/// * `filename` - The filename of the file where the suffix array is stored
///
/// # Returns
///
/// Returns the sample rate of the suffix array, together with the suffix array
///
/// # Errors
///
/// Returns any error from opening the file or reading the file
pub fn load_suffix_array(reader: &mut impl BufRead) -> Result<(u8, Vec<i64>), Box<dyn Error>> {
    // Read the sample rate from the binary file (1 byte)
    let mut sample_rate_buffer = [0_u8; 1];
    reader.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file")?;
    let sample_rate = sample_rate_buffer[0];

    // Read the size of the suffix array from the binary file (8 bytes)
    let mut size_buffer = [0_u8; 8];
    reader.read_exact(&mut size_buffer).map_err(|_| "Could not read the size of the suffix array from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;

    let mut sa = Vec::with_capacity(size);
    sa.read_binary(reader).map_err(|_| "Could not read the suffix array from the binary file")?;

    Ok((sample_rate, sa))
}

pub fn dump_suffix_array(
    sa: &Vec<i64>,
    sparseness_factor: u8,
    writer: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    // Write the flags to the writer
    // 00000000 indicates that the suffix array is not compressed
    writer.write(&[64_u8]).map_err(|_| "Could not write the flags to the writer")?;

    // Write the sparseness factor to the writer
    writer.write(&[sparseness_factor]).map_err(|_| "Could not write the sparseness factor to the writer")?;

    // Write the size of the suffix array to the writer
    let sa_len = sa.len();
    writer.write(&(sa_len).to_le_bytes()).map_err(|_| "Could not write the size of the suffix array to the writer")?;

    // Write the suffix array to the writer
    let sa = Vec::with_capacity(sa_len);
    sa.write_binary(writer).map_err(|_| "Could not write the suffix array to the writer")?;

    Ok(())
}

/// Fills the buffer with data read from the input.
///
/// # Arguments
///
/// * `input` - The input source to read data from.
/// * `buffer` - The buffer to fill with data.
///
/// # Returns
///
/// Returns a tuple `(finished, bytes_read)` where `finished` indicates whether the end of the input is reached,
/// and `bytes_read` is the number of bytes read into the buffer.
fn fill_buffer<T: Read>(input: &mut T, buffer: &mut Vec<u8>) -> (bool, usize) {
    // Store the buffer size in advance, because rust will complain
    // about the buffer being borrowed mutably while it's borrowed
    let buffer_size = buffer.len();

    let mut writable_buffer_space = buffer.as_mut();

    loop {
        match input.read(writable_buffer_space) {
            // No bytes written, which means we've completely filled the buffer
            // or we've reached the end of the file
            Ok(0) => {
                return (
                    !writable_buffer_space.is_empty(),
                    buffer_size - writable_buffer_space.len()
                );
            }

            // We've read {bytes_read} bytes
            Ok(bytes_read) => {
                // Shrink the writable buffer slice
                writable_buffer_space = writable_buffer_space[bytes_read..].as_mut();
            }

            Err(err) => {
                panic!("Error while reading input: {}", err);
            }
        }
    }
}
