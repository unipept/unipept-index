use std::{error::Error, io::{BufRead, Write}, path::Path};

pub trait WriteBinary {
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>>;
}

pub trait ReadBinary: Sized {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>>;
}

pub trait ReadBinaryMmap: Sized {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>>;
}
