//! Independent gzip members encoded by libdeflate on the existing gzp workers.
//! gzp owns bounded queues, ordered writes and worker shutdown; each member
//! carries its own header, CRC32 and uncompressed length.
use bytes::Bytes;
use gzp::{check::PassThroughCheck, Compression, FormatSpec, GzpError};
use libdeflater::{CompressionLvl, Compressor};
use std::io;

#[derive(Clone, Copy, Debug)]
pub struct Gzip;

pub(crate) fn threads_per_output(outputs: usize) -> usize {
    // Share the compression budget across files, leaving room for decoding and
    // processing. A single output can use the whole budget.
    num_cpus::get()
        .div_ceil(2)
        .clamp(1, 8)
        .div_ceil(outputs.max(1))
}

impl FormatSpec for Gzip {
    type C = PassThroughCheck;
    type Compressor = Compressor;
    const DEFAULT_BUFSIZE: usize = 256 * 1024;

    fn new() -> Self {
        Self
    }

    fn create_compressor(&self, level: Compression) -> Result<Compressor, GzpError> {
        let level = CompressionLvl::new(level.level() as i32)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{e:?}")))?;
        Ok(Compressor::new(level))
    }

    fn needs_dict(&self) -> bool {
        false
    }

    fn encode(
        &self,
        input: &[u8],
        encoder: &mut Compressor,
        _level: Compression,
        _dict: Option<&Bytes>,
        _is_last: bool,
    ) -> Result<Vec<u8>, GzpError> {
        let mut output = vec![0; encoder.gzip_compress_bound(input.len())];
        let size = encoder
            .gzip_compress(input, &mut output)
            .map_err(|e| io::Error::other(format!("gzip compression failed: {e:?}")))?;
        output.truncate(size);
        Ok(output)
    }

    fn header(&self, _level: Compression) -> Vec<u8> {
        Vec::new()
    }

    fn footer(&self, _check: &PassThroughCheck) -> Vec<u8> {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gzp::{par::compress::ParCompressBuilder, ZWriter};
    use std::io::{Read, Write};

    #[test]
    fn parallel_members_roundtrip_levels_boundaries_and_empty_output() {
        let data: Vec<_> = (0..Gzip::DEFAULT_BUFSIZE * 3 + 31)
            .map(|i| ((i * 71 + i / 239) % 256) as u8)
            .collect();
        let temp = tempfile::TempDir::new().unwrap();
        for level in [1, 4, 9] {
            for len in [
                0,
                1,
                Gzip::DEFAULT_BUFSIZE - 1,
                Gzip::DEFAULT_BUFSIZE,
                data.len(),
            ] {
                let path = temp.path().join("out.gz");
                let mut writer = ParCompressBuilder::<Gzip>::new()
                    .compression_level(Compression::new(level))
                    .num_threads(2)
                    .unwrap()
                    .from_writer(std::fs::File::create(&path).unwrap());
                for block in data[..len].chunks(27013) {
                    writer.write_all(block).unwrap();
                }
                writer.flush().unwrap();
                writer.finish().unwrap();
                let bytes = std::fs::read(path).unwrap();
                assert!(
                    bytes.starts_with(&[31, 139]),
                    "empty output must still be gzip"
                );
                let mut actual = Vec::new();
                flate2::read::MultiGzDecoder::new(&bytes[..])
                    .read_to_end(&mut actual)
                    .unwrap();
                assert_eq!(actual, data[..len]);
            }
        }
    }

    #[test]
    fn compressed_output_errors_are_returned_by_finish() {
        struct Broken;
        impl Write for Broken {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed output"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let mut writer = ParCompressBuilder::<Gzip>::new()
            .num_threads(2)
            .unwrap()
            .from_writer(Broken);
        writer.write_all(b"some input").unwrap();
        assert!(writer.finish().is_err());
    }
}
