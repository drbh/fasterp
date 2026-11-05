//! Internals exposed for the fuzz targets in `fuzz/` and for integration
//! tests. Not a stable API.
use std::io::{self, BufReader, Read};

/// `input` as one gzip member written by the FASTQ-aware encoder, whether or
/// not the output path would choose it for this block.
pub fn fq_gzip_member(input: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    crate::fqdeflate::gzip_member(input, &mut out);
    out
}

/// Whether gzip output uses the FASTQ-aware encoder for this block.
pub fn fq_suits(input: &[u8]) -> bool {
    crate::fqdeflate::suits(input)
}

/// Decompress `gz` with the parallel single-member decoder, in regions of
/// `region` compressed bytes on `threads` workers.
pub fn pgzip(gz: &[u8], region: usize, threads: usize) -> io::Result<Vec<u8>> {
    let (sender, receiver) = crossbeam_channel::unbounded();
    std::thread::scope(|scope| {
        scope.spawn(|| {
            crate::pgzip::decompress_source(
                &gz,
                region.max(1),
                threads.max(1),
                sender,
                crossbeam_channel::never(),
            )
        });
        let mut out = Vec::new();
        for part in receiver {
            let part = part?;
            if part.is_empty() {
                return Ok(out);
            }
            out.extend_from_slice(&part);
        }
        Err(io::Error::other("decoder stopped without an end marker"))
    })
}

/// Decompress `gz` as local gzip files are read when they cannot be split:
/// member by member, with each member's checksum verified by the reader.
pub fn sequential_gzip(gz: &[u8]) -> io::Result<Vec<u8>> {
    let input: Box<dyn Read + Send> = Box::new(io::Cursor::new(gz.to_vec()));
    let mut out = Vec::new();
    crate::io::read_gzip_ahead(BufReader::new(input)).read_to_end(&mut out)?;
    Ok(out)
}

/// The decompressed contents of `path`, read as fasterp reads its inputs.
pub fn read_input(path: &str) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    crate::io::open_input(path)
        .map_err(io::Error::other)?
        .read_to_end(&mut out)?;
    Ok(out)
}

/// Decode every single-member gzip input in parallel, in regions of `region`
/// compressed bytes; 0 restores the normal choice.
pub fn force_split_decoding(region: usize) {
    crate::io::force_split_decoding(region);
}
