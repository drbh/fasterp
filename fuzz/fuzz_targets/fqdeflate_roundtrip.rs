//! The FASTQ-aware encoder must write gzip that every inflater decodes back
//! to the input, for any input.
#![no_main]
#[path = "common.rs"]
mod common;
use libfuzzer_sys::fuzz_target;
use std::io::Read;

fuzz_target!(|data: &[u8]| {
    let mut c = common::Choices::new(data);
    let input = if c.byte() % 4 == 0 {
        c.rest().to_vec()
    } else {
        common::payload(&mut c, 300_000)
    };
    let gz = fasterp::testing::fq_gzip_member(&input);
    let _ = fasterp::testing::fq_suits(&input);

    let mut back = Vec::with_capacity(input.len());
    flate2::read::GzDecoder::new(&gz[..])
        .read_to_end(&mut back)
        .expect("flate2 rejects the encoder's output");
    assert!(back == input, "flate2 decodes different bytes");

    let mut d = libdeflater::Decompressor::new();
    let mut buf = vec![0; input.len() + 1];
    let n = d
        .gzip_decompress(&gz, &mut buf)
        .expect("libdeflate rejects the encoder's output");
    assert!(buf[..n] == input[..], "libdeflate decodes different bytes");
});
