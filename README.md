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
fasterp - Fast FASTQ preprocessing in Rust, with familiar [fastp](https://github.com/OpenGene/fastp) options.

Trim adapters, filter reads, merge pairs, and generate HTML and JSON quality reports.

Faster than fastp 1.4.0 in **32 of 33 benchmark cases** on an Apple M4 Pro. See the [results and reproduction steps](BENCHMARKS.md).

## Install

```bash
cargo install --git https://github.com/drbh/fasterp.git --locked
```

## Usage

```bash
# Single-end
fasterp -i input.fq.gz -o output.fq.gz

# Paired-end
fasterp -i R1.fq.gz -I R2.fq.gz -o out1.fq.gz -O out2.fq.gz
```

Open `fasterp.html` for the quality report. Use `--help` for options and see [fastp compatibility](COMPATIBILITY.md) for supported behavior and differences.

## Documentation

[Full documentation and tutorial](https://drbh.github.io/fasterp) · [Build and test](CONTRIBUTING.md)

MIT. See [LICENSE](LICENSE) and [LICENSE.fastp](LICENSE.fastp).
