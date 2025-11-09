# fasterp

A fast and simple FASTQ preprocessor written in Rust, implementing core functionality compatible with [fastp](https://github.com/OpenGene/fastp).

## Features

- **Length filtering**: Filter reads below minimum length threshold
- **Quality filtering**: Filter reads by mean quality score (Phred+33 encoding)
- **N-base filtering**: Filter reads with too many ambiguous bases
- **JSON reporting**: Generates fastp-compatible JSON statistics report
- **Kmer counting**: Counts all 1024 possible 5-mers (AAAAA through TTTTT)
- **FASTQ validation**: Detects and skips invalid records (mismatched sequence/quality lengths)
- **Fast performance**: 5-16x faster than fastp on typical datasets
- **Clean codebase**: Simple, well-documented Rust implementation

## Installation

```bash
cargo build --release
```

The binary will be available at `target/release/fasterp`

## Usage

```bash
# Basic usage (default parameters)
fasterp -i input.fq -o output.fq

# Custom filtering
fasterp -i input.fq -o output.fq -l 20 -q 30 -n 3

# View all options
fasterp --help
```

### Options

- `-i, --input <INPUT>`: Input FASTQ file (required)
- `-o, --output <OUTPUT>`: Output FASTQ file (required)
- `-l, --length-required <LENGTH>`: Minimum read length (default: 15)
- `-q, --qualified-quality-phred <QUALITY>`: Mean quality score threshold (default: 0, disabled)
- `-n, --n-base-limit <N>`: Maximum number of N bases allowed (default: 5)
- `-j, --json <JSON>`: JSON report file (default: fastp.json)

