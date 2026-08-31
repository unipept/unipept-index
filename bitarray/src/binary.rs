//! Reading and writing a bit array's backing words as raw little-endian `u64`s.
//!
//! Deliberately headerless: the width and the value count live in the *containing* structure's
//! header, which the caller reads and writes itself. See [`Binary`] for what that leaves the
//! caller responsible for.

use std::io::{BufRead, Read, Result, Write};

/// Reads and writes a bit array's backing words as packed binary data.
///
/// The payload is words only — no length, no width, no marker. Everything needed to interpret it
/// lives in a header the caller owns, which is what makes a bit array embeddable in a larger file
/// alongside other sections.
///
/// # This is not `binary_traits::ReadBinary`
///
/// The workspace has a second pair of traits with the same method names, in `binary-traits`, for
/// whole index structures. They are not interchangeable, and the contracts are opposite where it
/// matters most: `binary_traits::ReadBinary` consumes *exactly* the bytes its writer produced, so
/// several structures chain over one stream, whereas [`read_binary`](Self::read_binary) below
/// consumes the reader to EOF. A structure implements that trait; its bit array implements this
/// one, usually from inside that implementation.
pub trait Binary {
    /// Writes every backing word as little-endian bytes, and nothing else.
    ///
    /// The output is `word_len() * 8` bytes. A caller embedding this in a larger file must have
    /// written whatever header is needed to interpret it — at minimum the width and the value
    /// count, neither of which is recoverable from the payload.
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<()>;

    /// Refills the backing store from `reader`, **consuming it to EOF**.
    ///
    /// Two consequences, both of which every caller in this workspace has had to handle:
    ///
    /// * *It reads to the end of the stream.* Since the payload carries no length, there is
    ///   nothing to stop at. A bit array is therefore either the last section of its file — as in
    ///   the compressed suffix array — or the caller must hand over a bounded reader, as
    ///   `protein-text` does with `Read::take` when the protein metadata follows the text in
    ///   the same file. Passing an unbounded reader for a non-final section silently swallows the
    ///   sections after it.
    /// * *It validates nothing.* The store ends up holding however many words the reader happened
    ///   to yield, which says nothing about the count the caller's header declared. A truncated
    ///   file therefore loads cleanly and panics later, on the first lookup past the real data,
    ///   unless the caller compares `word_len()` against `required_words()` afterwards. Both
    ///   implementations expose that pair for exactly this check.
    ///
    /// Errors only on an I/O error from `reader`; a short stream is not one.
    fn read_binary<R: BufRead>(&mut self, reader: R) -> Result<()>;
}

/// Writes each `u64` word as little-endian bytes.
pub(crate) fn write_words<W: Write>(data: &[u64], writer: &mut W) -> Result<()> {
    for value in data {
        writer.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

/// Bytes pulled out of the reader per iteration of the refill loop below.
///
/// Chosen to be **larger than the caller's `BufReader` capacity**, not for its own sake. Every
/// preloaded structure is read through `binary_traits::load_owned`, whose `BufReader` uses the
/// default 8 KiB; `BufReader::read` only bypasses that internal buffer when the destination slice
/// is at least as large as its capacity. At 8 KiB it never was, so each request was served partly
/// from the internal buffer and partly from a fresh refill — one syscall *and two memcpys* per
/// 8 KiB, over the whole suffix array and the whole protein text. Above the capacity the bypass
/// applies: one syscall and one copy straight into this buffer.
///
/// So this is not a knob to tidy downwards. Anything below `load_owned`'s `BufReader` capacity
/// silently restores the double copy.
const READ_CHUNK_BYTES: usize = 4 << 20;

/// Clears `data` and refills it by reading little-endian `u64` words from `reader`.
pub(crate) fn read_words_into<R: BufRead>(data: &mut Vec<u64>, mut reader: R) -> Result<()> {
    data.clear();
    let mut buffer = vec![0u8; READ_CHUNK_BYTES];
    loop {
        let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer)?;
        for chunk in buffer[..bytes_read].as_chunks::<8>().0 {
            data.push(u64::from_le_bytes(*chunk));
        }
        if finished {
            break;
        }
    }
    Ok(())
}

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
