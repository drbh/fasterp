# fasterp

**fasterp** — a Rust reimplementation of [**fastp**](https://github.com/OpenGene/fastp).
Same interface and behavior, typically **3–10× faster**.

| Feature              | fastp | fasterp |
| :------------------- | :---: | :-----: |
| Filtering            |   ✓   |    ✓    |
| Quality trimming     |   ✓   |    ✓    |
| PolyG/PolyX trimming |   ✓   |    ✓    |
| JSON report          |   ✓   |    ✓    |
| Multi-threading      |   ✓   |    ✓    |
| gzip I/O             |   ✓   |    ✓    |
| Kmer counting        |   ✓   |    ✓    |
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
  fastp    ████████████████████████████ 141ms
  fasterp  ██ 12ms  ⚡ 11.5× faster

Medium dataset (10k reads):
  fastp    ████████████████████████████ 145ms
  fasterp  ████ 23ms  ⚡ 6.2× faster

With quality trimming (10k reads):
  fastp    ████████████████████████████ 146ms
  fasterp  ████ 23ms  ⚡ 6.1× faster
```

*results from [recent run](https://github.com/drbh/fasterp/actions/runs/19214620217/job/54922026577#step:7:114) on github actions

## Correctness

`fasterp` produces identical output to `fastp` and can replace it directly in existing workflows.
Integration tests confirm equivalence across datasets (`cargo test`), and all checks run automatically in CI.

See [integration tests](tests/integration_tests.rs) for ~40 examples comparing the inputs and outputs of both tools, these can be run locally via `cargo test` for verification. 