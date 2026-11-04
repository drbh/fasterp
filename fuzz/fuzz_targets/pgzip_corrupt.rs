//! Damaged gzip must make the parallel decoder fail cleanly: no panic, no
//! hang, and never output that flate2 would not also produce.
#![no_main]
#[path = "common.rs"]
mod common;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut c = common::Choices::new(data);
    let region = 64 + c.below(8 << 10);
    let threads = 1 + c.below(4);
    let input = common::payload(&mut c, 100_000);
    let mut gz = common::gzip(&mut c, &input);
    for _ in 0..1 + c.below(4) {
        match c.byte() % 3 {
            0 => {
                let at = c.below(gz.len());
                gz[at] ^= 1 << (c.byte() % 8);
            }
            1 => {
                let at = c.below(gz.len());
                gz[at] = c.byte();
            }
            _ => {
                let keep = c.below(gz.len() + 1);
                gz.truncate(keep.max(1));
            }
        }
    }
    if let Ok(path) = std::env::var("FUZZ_DUMP") {
        std::fs::write(&path, &gz).unwrap();
        std::fs::write(format!("{path}.args"), format!("{region} {threads}")).unwrap();
    }
    if let Ok(got) = fasterp::testing::pgzip(&gz, region, threads) {
        let expected = common::flate2_decode(&gz)
            .expect("parallel decoder accepted a stream flate2 rejects");
        assert!(got == expected, "accepted with different bytes than flate2");
    }
});
