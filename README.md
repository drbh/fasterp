# fasterp

**fasterp** — a fast, modern, and simplified reimplementation of [fastp](https://github.com/OpenGene/fastp) in **Rust**.

> 2–4× faster than fastp, with full feature parity and exactly the same command-line interface.

| Feature              | fastp | fasterp |
| :------------------- | :---: | :-----: |
| Filtering            |   ✅   |    ✅    |
| Quality trimming     |   ✅   |    ✅    |
| PolyG/PolyX trimming |   ✅   |    ✅    |
| JSON report          |   ✅   |    ✅    |
| Multi-threading      |   ✅   |    ✅    |
| gzip I/O             |   ✅   |    ✅    |
| Kmer counting        |   ✅   |    ✅    |
| Rust safety          |   ❌   |    ✅    |
| Speed                |   🐢   |    ⚡    |

## Usage

```bash
fasterp -i input.fq -o output.fq
```

Same arguments as `fastp`, but faster.

## Performance

`fasterp` is often much faster than `fastp`, depending on the dataset and hardware used. Check on your machine with the following command:

```bash
cargo run --example bench_chart
```

results

```txt
Performance Comparison: fasterp vs fastp

Small dataset (1k reads):
  fastp    ████████████████████████████ 67ms
  fasterp  ████████ 19ms  ⚡ 3.4× faster

Medium dataset (10k reads):
  fastp    ████████████████████████████ 100ms
  fasterp  ██████████ 37ms  ⚡ 2.7× faster

With quality trimming (10k reads):
  fastp    ████████████████████████████ 100ms
  fasterp  ██████████ 39ms  ⚡ 2.6× faster
```

## Correctness

`fasterp` produces output identical to `fastp` and aims to maintain full feature parity and interface compatibility.

You should be able to replace `fastp` with `fasterp` in your pipelines without any changes, except for the speedup.

We include a large number of integration tests to ensure correctness in `tests/integration_tests.rs`, each of these test compares the output of `fasterp` against `fastp` on various datasets and options. 

These test can be run with: `cargo test` and are run automatically on every commit via GitHub Actions.