//! Fixtures for the preloaded suffix-array tests, on top of the backend-agnostic ones in
//! [`array::test_utils`](crate::array::test_utils).
//!
//! Both are failing halves of an I/O pair. The preloaded backend is the only one that streams
//! through `Read`/`Write` — the mmap backend is handed a whole file — so it is the only one whose
//! error paths need a reader or writer that gives out partway.

use std::io::{BufRead, Read, Write};

/// A writer that reports success for `valid_write_count` calls and then fails. Each call claims to
/// have written a single byte, so a `write_all` over an n-byte field consumes n of the count.
pub struct FailingWriter {
    pub valid_write_count: usize
}

impl Write for FailingWriter {
    fn write(&mut self, _: &[u8]) -> Result<usize, std::io::Error> {
        if self.valid_write_count == 0 {
            return Err(std::io::Error::other("Write failed"));
        }
        self.valid_write_count -= 1;
        Ok(1)
    }

    fn flush(&mut self) -> Result<(), std::io::Error> {
        Ok(())
    }
}

/// A reader that fills `valid_read_count` buffers and then fails.
pub struct FailingReader {
    pub valid_read_count: usize
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
