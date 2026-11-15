
```
    ░████                          ░██                                   
   ░██                             ░██                                   
░████████  ░██████    ░███████  ░████████  ░███████  ░██░████ ░████████  
   ░██          ░██  ░██           ░██    ░██    ░██ ░███     ░██    ░██ 
   ░██     ░███████   ░███████     ░██    ░█████████ ░██      ░██    ░██ 
   ░██    ░██   ░██         ░██    ░██    ░██        ░██      ░███   ░██ 
   ░██     ░█████░██  ░███████      ░████  ░███████  ░██      ░██░█████  
                                                              ░██        
                                                              ░██        
```
fasterp - A faster [fastp](https://github.com/OpenGene/fastp) — same interface, same output. The [benchmarks](https://drbh.github.io/fasterp/benchmarks.html) show how much faster on each dataset, workload and worker count.

```
10M read pairs:  fastp 10s → fasterp 2.9s
```

Built with SIMD acceleration (AVX2/NEON), efficient multi-threading, and minimal memory allocation.

## Install

```bash
cargo install --git https://github.com/drbh/fasterp.git
```

## Usage

```bash
# Single-end
fasterp -i input.fq -o output.fq

# Paired-end
fasterp -i R1.fq -I R2.fq -o out1.fq -O out2.fq
```

Drop-in replacement for fastp — all the same flags work.

## Using fasterp as a library

`fasterp::Run` runs the command line's engine on readers or on readkit
batch sources and returns columnar batches instead of writing files: the
passing reads are what `--out1`/`--out2` would hold, in the same order, and
`finish()` returns the `-j` report. Downstream code reads the sequence and quality columns in
place. The examples:

| example | what it shows | run it |
| --- | --- | --- |
| `zero_copy` | a stand-in decoder fills batch columns, fasterp trims them, a k-mer counter reads the output columns, buffers are recycled | `cargo run --release --example zero_copy` |
| `arrow_export` | passing reads as an arrow-rs `RecordBatch`, without copying | `cargo run --release --features arrow --example arrow_export` |
| `sra_pipeline` | an `.sra` run decoded in process by readkit-sra (NCBI VDB) straight into fasterp, k-mers and stats from the output columns | `cargo run --release --manifest-path examples/sra-pipeline/Cargo.toml -- RUN.sra [THREADS]` |

`sra_pipeline` is a crate of its own (`examples/sra-pipeline`), so its
dependencies stay out of fasterp's builds. It needs the SRA Toolkit's
`libncbi-ngs` (set `NGS_PREFIX` if it is not under `/opt/homebrew`). With
`--fastq OUT` it also writes the passing reads, which match
`fasterp -o` run on `fasterq-dump --stdout --split-spot RUN.sra` apart from
the `+` line.

## Documentation

**[Full documentation and tutorial](https://drbh.github.io/fasterp)**
