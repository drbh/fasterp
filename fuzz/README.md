# Fuzzing

Coverage-guided fuzz targets for the code that writes and reads gzip. Each
compares against an independent implementation (flate2, libdeflate), so a
target fails on wrong bytes, not only on crashes.

| Target | Checks |
|---|---|
| `fqdeflate_roundtrip` | The FASTQ-aware encoder's output decodes back to its input with flate2 and libdeflate |
| `pgzip_differential` | The parallel single-member decoder gives flate2's bytes, for streams from zlib, libdeflate and the FASTQ-aware encoder, with regions small enough that boundaries land everywhere |
| `pgzip_corrupt` | On damaged streams the parallel decoder errors, or gives flate2's bytes, and never panics, hangs or runs out of memory |
| `gzip_members` | The sequential multi-member reader accepts exactly the streams flate2 accepts, with the same bytes |

Inputs are built from the fuzzer's bytes by `fuzz_targets/common.rs`, which
shapes them like FASTQ (binned and full-range qualities, duplicate reads,
CRLF, truncated records) so the targets reach the FASTQ-specific paths.

Run a target with nightly Rust and [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz):

```bash
cargo install cargo-fuzz
cd fuzz
cargo +nightly fuzz run -O pgzip_differential -- -max_total_time=600 -max_len=65536 -timeout=20 \
  -rss_limit_mb=8192 -malloc_limit_mb=2048
```

`-malloc_limit_mb` catches a single runaway allocation. The resident-set limit is
loose because the decoder targets start threads on every input, and on macOS
AddressSanitizer's per-thread bookkeeping makes the fuzzer's memory grow by
about a gigabyte a minute, although the decoder itself does not leak.

CI runs every target for a minute on each push and pull request; a manual
run of the CI workflow takes a longer `fuzz_seconds`. A failing input is
saved under `artifacts/<target>/` and uploaded by CI. To see the gzip stream a
`pgzip_corrupt` or `gzip_members` input produced, replay it with
`FUZZ_DUMP=stream.gz`:

```bash
FUZZ_DUMP=/tmp/stream.gz cargo +nightly fuzz run -O pgzip_corrupt artifacts/pgzip_corrupt/<file>
```
