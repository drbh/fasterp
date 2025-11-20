# fasterp

A faster [fastp](https://github.com/OpenGene/fastp) — same interface, same output, **2-10x faster**.

```
10M read pairs:  fastp 11s → fasterp 6s
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

## Documentation

**[Full documentation and tutorial](https://drbh.github.io/fasterp)**
