# Contributing

The native CLI targets the read-processing behavior of **fastp 1.4.0**. Preserve ordered FASTQ output and meaningful report values when changing an algorithm. [COMPATIBILITY.md](COMPATIBILITY.md) lists intentional differences and test coverage.

## Build and test

```sh
cargo build --locked --release
cargo test --locked --release --lib --bin fasterp --test compatibility_v140
python3 -m unittest discover -s scripts -p 'test_*.py'
```

The focused tests generate their own data. The complete `cargo test --locked --release` suite also requires the fixtures from `test_data/getdata.sh`, a reference `fastp` on `PATH`, and sequencing datasets used by the older integration tests.

Build the compatibility reference from its pinned source commit using `scripts/build_fastp_reference.sh`. Install a C++ compiler and the ISA-L, libdeflate and Highway development libraries first. Use a new destination directory.

```sh
bash scripts/build_fastp_reference.sh target/fastp-reference
cargo build --locked --release
python3 scripts/compare_fastp.py --fastp target/fastp-reference/fastp \
  --fasterp target/release/fasterp --workers 1 4 8 --workdir target/fastp-parity
```

On Apple silicon with Homebrew libraries, pass `FASTP_INCLUDE_DIRS=/opt/homebrew/include FASTP_LIBRARY_DIRS=/opt/homebrew/lib CXX=clang++` before the helper command. Other installations can supply their own include/library paths. The reference commit is `ee125a52fddae9d8189e13cd8a33ad04988a1351`.

The differential runner checks decompressed FASTQ bytes, counters and report curves, including adapter discovery, merging and gzip output. Keep scalar reference tests for SIMD and search optimizations.

## Performance changes

Use [BENCHMARKS.md](BENCHMARKS.md) for the datasets and reproduction commands. Compare fasterp with the reference in one session, with equivalent options at 1, 4 and 8 processing workers. Include elapsed time, process CPU time, peak RSS and compressed output sizes. Check output and report agreement before using timing results.

Freeze the normal release executable before running tests, which can replace `target/release/fasterp`. Run benchmarks without concurrent builds, tests or downloads. The current benchmark runner uses macOS `/usr/bin/time -lp`; Linux needs its own measurement adapter and evidence.

The local Cargo configuration enables `target-cpu=native`. Use `RUSTFLAGS=''` when compiling portable distribution binaries or WebAssembly. To check the shared WASM code:

```sh
cargo test --locked --no-default-features --features wasm --lib
RUSTFLAGS='' cargo build --locked --release --target wasm32-unknown-unknown \
  --no-default-features --features wasm --lib
```

## Credits

fastp 1.4.0 compatibility and performance contributions by [@cheezeeerd](https://github.com/cheezeeerd) and [@deminden](https://github.com/deminden).

Preserve the project's MIT license, contributor credit, and [OpenGene's MIT notice](LICENSE.fastp) for upstream-derived algorithms and adapter data. Describe the bug or workload, the change, and how you tested it.
