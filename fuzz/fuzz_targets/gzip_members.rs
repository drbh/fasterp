//! The sequential reader checks each member's trailer on the reading side.
//! It must accept exactly the multi-member streams flate2 accepts, with the
//! same bytes, whatever the member sizes and damage.
#![no_main]
#[path = "common.rs"]
mod common;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut c = common::Choices::new(data);
    let mut gz = Vec::new();
    for _ in 0..1 + c.below(6) {
        let part = common::payload(&mut c, 60_000);
        gz.extend(common::gzip(&mut c, &part));
    }
    if c.byte() % 2 == 0 {
        for _ in 0..1 + c.below(3) {
            if gz.is_empty() {
                break;
            }
            let at = c.below(gz.len());
            gz[at] ^= 1 << (c.byte() % 8);
        }
    }
    if let Ok(path) = std::env::var("FUZZ_DUMP") {
        std::fs::write(path, &gz).unwrap();
    }
    let ours = fasterp::testing::sequential_gzip(&gz);
    let theirs = common::flate2_decode(&gz);
    match (ours, theirs) {
        (Ok(a), Ok(b)) => assert!(a == b, "different bytes than flate2"),
        (Err(_), Err(_)) => {}
        (Ok(_), Err(e)) => panic!("accepted a stream flate2 rejects: {e}"),
        (Err(e), Ok(_)) => panic!("rejected a stream flate2 accepts: {e}"),
    }
});
