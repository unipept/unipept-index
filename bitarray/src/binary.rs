use std::io::{BufRead, Read, Result, Write};

/// Trait for reading and writing a struct as packed binary data.
pub trait Binary {
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<()>;
    fn read_binary<R: BufRead>(&mut self, reader: R) -> Result<()>;
}

// ── Shared helpers ────────────────────────────────────────────────────────────

/// Writes each `u64` word as little-endian bytes.
pub(crate) fn write_words<W: Write>(data: &[u64], writer: &mut W) -> Result<()> {
    for value in data {
        writer.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

/// Clears `data` and refills it by reading little-endian `u64` words from `reader`.
pub(crate) fn read_words_into<R: BufRead>(data: &mut Vec<u64>, mut reader: R) -> Result<()> {
    data.clear();
    let mut buffer = vec![0u8; 8 * 1024];
    loop {
        let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer)?;
        for chunk in buffer[..bytes_read].chunks_exact(8) {
            data.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        if finished {
            break;
        }
    }
    Ok(())
}

// ── fill_buffer ───────────────────────────────────────────────────────────────

/// Fills `buffer` as fully as possible from `input`.
///
/// Returns `(finished, bytes_read)` where `finished` is true when the reader
/// returned EOF before the buffer was full.
pub(crate) fn fill_buffer<T: Read>(input: &mut T, buffer: &mut Vec<u8>) -> Result<(bool, usize)> {
    let buffer_size = buffer.len();
    let mut writable = buffer.as_mut();
    loop {
        match input.read(writable) {
            Ok(0) => return Ok((!writable.is_empty(), buffer_size - writable.len())),
            Ok(n) => writable = writable[n..].as_mut(),
            Err(e) => return Err(e)
        }
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    struct ErrorInput;

    impl Read for ErrorInput {
        fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("read error"))
        }
    }

    #[test]
    fn test_fill_buffer() {
        let input_str = "a".repeat(8_000);
        let mut input = input_str.as_bytes();
        let mut buffer = vec![0; 800];
        loop {
            let (finished, bytes_read) = fill_buffer(&mut input, &mut buffer).unwrap();
            if finished {
                assert!(bytes_read < 800);
                break;
            } else {
                assert_eq!(bytes_read, 800);
            }
        }
    }

    #[test]
    fn test_fill_buffer_read_error() {
        let mut input = ErrorInput;
        let mut buffer = vec![0; 800];
        assert!(fill_buffer(&mut input, &mut buffer).is_err());
    }
}
