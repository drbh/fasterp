# fasterp


**fasterp** — a Rust reimplementation of [**fastp**](https://github.com/OpenGene/fastp).
Same interface and behavior, often **>5× faster** but your mileage may vary depending on dataset and CPU.

| Feature              | fastp | fasterp |
| :------------------- | :---: | :-----: |
| Filtering            |   ✓   |    ✓    |
| Quality trimming     |   ✓   |    ✓    |
| PolyG/PolyX trimming |   ✓   |    ✓    |
| JSON report          |   ✓   |    ✓    |
| Multi-threading      |   ✓   |    ✓    |
| gzip I/O             |   ✓   |    ✓    |
| Kmer counting        |   ✓   |    ✓    |
| Adapter trimming     |   ✓   |    -    |
| Rust safety          |   -   |    ✓    |
| SIMD acceleration    |   -   |    ✓    |

## Installation

```bash
cargo install --git https://github.com/drbh/fasterp.git
fasterp --version
```

## Usage

```bash
fasterp -i input.fq -o output.fq
```

Uses the same CLI and options as `fastp`.

## Performance

Run benchmarks with:

```bash
cargo run --example bench_chart
```

Example:

```txt
Small dataset (1k reads):
  fastp    ████████████████████████████ 147ms
  fasterp  ██ 12ms  ⚡ 12.2× faster

Medium dataset (10k reads):
  fastp    ████████████████████████████ 151ms
  fasterp  ████ 23ms  ⚡ 6.5× faster

With quality trimming (10k reads):
  fastp    ████████████████████████████ 152ms
  fasterp  ████ 24ms  ⚡ 6.4× faster
```

*results from [recent run](https://github.com/drbh/fasterp/actions/runs/19215579022/job/54924377219#step:7:123) on github actions

## Correctness

`fasterp` produces identical output to `fastp` and can replace it directly in existing workflows.
Integration tests confirm equivalence across datasets (`cargo test`), and all checks run automatically in CI.

See [integration tests](tests/integration_tests.rs) for ~40 examples comparing the inputs and outputs of both tools, these can be run locally via `cargo test` for verification.

## Sanity Check

Quick verification that `fasterp` and `fastp` produce identical results:

```bash
# Download test data
wget -q https://raw.githubusercontent.com/OpenGene/fastp/master/testdata/R1.fq

# Run fastp
fastp -i R1.fq -o /tmp/fastp_out.fq -j /tmp/fastp.json 2>/dev/null

# Run fasterp
fasterp -i R1.fq -o /tmp/fasterp_out.fq -j /tmp/fasterp.json

# Compare outputs - hashes should match
sha256sum /tmp/fastp_out.fq /tmp/fasterp_out.fq
# b3bed5be02bb6ffad48800c4befaaae02c0baea5349cd675a3206efddf9b8912  /tmp/fastp_out.fq
# b3bed5be02bb6ffad48800c4befaaae02c0baea5349cd675a3206efddf9b8912  /tmp/fasterp_out.fq

# Compare kmer counts from JSON reports
jq '.read1_before_filtering.kmer_count | {AAAAA, TTTTT, CCCCC, GGGGG}' /tmp/fastp.json
jq '.read1_before_filtering.kmer_count | {AAAAA, TTTTT, CCCCC, GGGGG}' /tmp/fasterp.json
# {
#   "AAAAA": 43,
#   "TTTTT": 3,
#   "CCCCC": 6,
#   "GGGGG": 0
# }
# {
#   "AAAAA": 43,
#   "TTTTT": 3,
#   "CCCCC": 6,
#   "GGGGG": 0
# }
``` 