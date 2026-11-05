//! The parallel decoder must give exactly what flate2 gives, for streams
//! from every compressor, with regions small enough that block searches,
//! false starts and region boundaries land everywhere.
#![no_main]
#[path = "common.rs"]
mod common;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut c = common::Choices::new(data);
    let region = 64 + c.below(16 << 10);
    let threads = 1 + c.below(4);
    let input = common::payload(&mut c, 400_000);
    let mut gz = common::gzip(&mut c, &input);
    let mut expected = input;
    if c.byte() % 4 == 0 {
        // Later members are decoded sequentially after the first.
        let more = common::payload(&mut c, 20_000);
        gz.extend(common::gzip(&mut c, &more));
        expected.extend(more);
    }
    let got = fasterp::testing::pgzip(&gz, region, threads)
        .unwrap_or_else(|e| panic!("valid stream rejected: {e} (region {region}, threads {threads})"));
    assert!(got == expected, "different bytes (region {region}, threads {threads})");
});
