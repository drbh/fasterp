//! Turns fuzzer bytes into inputs shaped like the data fasterp handles, so
//! the fuzzer reaches the FASTQ-specific paths rather than only random bytes.
#![allow(dead_code)]
use std::io::{Read, Write};

/// Reads choices from the fuzzer's bytes; zeros once they run out.
pub struct Choices<'a> {
    data: &'a [u8],
    at: usize,
}

impl<'a> Choices<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, at: 0 }
    }
    pub fn byte(&mut self) -> u8 {
        let b = self.data.get(self.at).copied().unwrap_or(0);
        self.at += 1;
        b
    }
    pub fn below(&mut self, n: usize) -> usize {
        let v = u16::from_le_bytes([self.byte(), self.byte()]) as usize;
        v % n.max(1)
    }
    pub fn rest(&self) -> &'a [u8] {
        &self.data[self.at.min(self.data.len())..]
    }
}

/// Bytes up to about `limit` long: FASTQ records (binned, full-range or
/// arbitrary qualities, duplicate reads, CRLF, truncated lines) mixed with
/// raw fuzzer bytes and long runs.
pub fn payload(c: &mut Choices, limit: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let quals: Vec<u8> = match c.byte() % 4 {
        0 => b"F:,#".to_vec(),
        1 => b"FFFFFFF:,".to_vec(),
        2 => (b'!'..=b'J').collect(),
        _ => (0..=255).collect(),
    };
    let target = c.below(limit + 1);
    let mut last: Vec<u8> = Vec::new();
    let mut i = 0usize;
    while out.len() < target {
        i += 1;
        let before = out.len();
        match c.byte() % 32 {
            0 => {
                let n = c.below(512);
                let rest = c.rest();
                out.extend_from_slice(&rest[..n.min(rest.len())]);
                for _ in 0..n.min(rest.len()) {
                    c.byte();
                }
            }
            1 => out.extend(std::iter::repeat_n(b"AF@+\n"[c.below(5)], c.below(40_000))),
            2 => out.extend_from_slice(b"@\n+\n@@\n\n+\r\n"),
            _ => {
                let nl: &[u8] = if c.byte() == 0 { b"\r\n" } else { b"\n" };
                write!(out, "@R.{i} {}:{}/1", c.byte() % 3, i * 7).unwrap();
                out.extend_from_slice(nl);
                let seq: Vec<u8> = if c.byte() % 5 == 0 && !last.is_empty() {
                    last.clone()
                } else {
                    let len = 20 + c.below(140);
                    let mut s = Vec::with_capacity(len);
                    while s.len() < len {
                        let b = c.byte();
                        for k in 0..4 {
                            s.push(b"ACGT"[usize::from(b >> (2 * k)) & 3]);
                        }
                        if b == 0xff {
                            s.push(b'N');
                        }
                    }
                    s.truncate(len);
                    s
                };
                out.extend_from_slice(&seq);
                out.extend_from_slice(nl);
                out.push(b'+');
                out.extend_from_slice(nl);
                let n = if c.byte() % 40 == 0 { c.below(seq.len() + 1) } else { seq.len() };
                for _ in 0..n {
                    out.push(quals[c.below(quals.len())]);
                }
                out.extend_from_slice(nl);
                last = seq;
            }
        }
        // Every step adds something, also once the fuzzer's bytes run out.
        if out.len() == before {
            out.push(b'\n');
        }
    }
    out.truncate(limit);
    out
}

/// One gzip member of `data` from one of several compressors, so decoders
/// meet stored, fixed and dynamic blocks of many sizes.
pub fn gzip(c: &mut Choices, data: &[u8]) -> Vec<u8> {
    match c.byte() % 4 {
        0 => {
            let level = u32::from(c.byte() % 10);
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
            e.write_all(data).unwrap();
            e.finish().unwrap()
        }
        1 => {
            let level = i32::from(c.byte() % 12) + 1;
            let mut d = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
            let mut out = vec![0; d.gzip_compress_bound(data.len())];
            let n = d.gzip_compress(data, &mut out).unwrap();
            out.truncate(n);
            out
        }
        _ => fasterp::testing::fq_gzip_member(data),
    }
}

/// What flate2 makes of a (possibly multi-member) gzip stream.
pub fn flate2_decode(gz: &[u8]) -> std::io::Result<Vec<u8>> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(gz).read_to_end(&mut out)?;
    Ok(out)
}
