# Benchmarks

Wall time in seconds for fastp and fasterp on the same inputs and options, measured with [hyperfine](https://github.com/sharkdp/hyperfine). Lower is better. The label next to each pair is fastp time divided by fasterp time, and the whiskers span the fastest and slowest run. The CSV also records each program's CPU time, peak memory, total gzip output size and a SHA-256 of its decompressed output, so a faster setting cannot hide larger files and any difference from fastp's records shows up.

<div class="bench-legend"><span class="fastp">fastp</span><span class="fasterp">fasterp</span></div>
<div id="bench-chart"></div>
<p id="bench-meta" class="bench-meta"></p>

fasterp compresses gzip output with a FASTQ-aware encoder at levels 2 to 4. The stream is standard gzip, and the records are the same as fastp's, as the digests in the CSV show. With binned qualities, as NovaSeq writes them, the files are as small as fastp's. With full-range qualities they are larger: about 4% on the 148 bp reads of SRR1766443 and 5 to 8% on the 50 bp reads of SRR5808766 and SRR891269. Reads sorted by alignment position, as in ERR3239279, overlap their neighbours, and level 4 compresses them to half the size the FASTQ-aware encoder does, so their blocks stay with level 4.

The ERR3239279-400M rows were measured with the `0.3.0` release and have not been rerun since.

## Workloads

All workloads read gzip input and write gzip output, with the default filters and an HTML and JSON report.

| Workload | Options |
|---|---|
| single-end | `-i R1 -o out.fq.gz` |
| paired-end | `-i R1 -I R2 -o out1.fq.gz -O out2.fq.gz` |
| paired-correct-merge | paired-end plus `-c -m --merged_out merged.fq.gz` |

## Datasets

The files are mirrored at `genedata.dholtz.com` so they download quickly. `.github/scripts/download_test_data.sh` fetches them by set name and checks each file against a pinned SHA-256.

| Set | Data | Prepared as |
|---|---|---|
| `ci` (SRR5808766) | Mouse ATAC-seq run [SRR5808766](https://www.ebi.ac.uk/ena/browser/view/SRR5808766) used in the [tutorial](./tutorial.md), 9.8 million pairs of 51 bp | Unchanged ENA files |
| `srr1766443-2m` | Human WGS run [SRR1766443](https://www.ebi.ac.uk/ena/browser/view/SRR1766443), 2x148 bp | First 2,000,000 records per mate from ENA, recompressed with `gzip -4 -n` |
| `srr1766443-2m-bgzf` | The same records as `srr1766443-2m` | Compressed with `bgzip -l 4`, to exercise BGZF input |
| `err3534523-2m` | 1000 Genomes 30x run [ERR3534523](https://www.ebi.ac.uk/ena/browser/view/ERR3534523), NovaSeq 6000, 2x150 bp with binned qualities | First 2,000,000 records per mate, recompressed with `gzip -4 -n` |
| `err3534523` | The full ERR3534523 run, 33.1 million pairs | Recompressed with `gzip -4 -n` |
| `srr891269-2m` | GM12878 ATAC-seq run [SRR891269](https://www.ebi.ac.uk/ena/browser/view/SRR891269) from Buenrostro et al. 2013, 2x50 bp with short inserts and Nextera adapters | First 2,000,000 records per mate from ENA, recompressed with `gzip -4 -n` |
| `err3239279-400m` | Human WGS 30x run [ERR3239279](https://www.ebi.ac.uk/ena/browser/view/ERR3239279), NovaSeq 6000, 2x150 bp with binned qualities, 25 GB gzipped. Run only at 8 workers, once, because fastp takes over four minutes per run | First 400,000,000 pairs, converted from NCBI's copy with `fastq-dump --split-spot` and compressed with `pigz -4 -n`. Not in the `all` set because of its size |

ERR3534523 was converted from NCBI's copy with `fasterq-dump`, and ERR3239279 with `fastq-dump`. Neither keeps the original Illumina read names, so the headers look like `@ERR3534523.1 1/1`. For ERR3534523, the sequences and qualities of the first 2,000,000 records of mate 2 match ENA's file. Read names starting with `@ERR` also mean fastp does not detect a two-color instrument, so poly-G trimming has to be enabled with `-g` (`ARGS=-g` below).

## Results

<div id="bench-table" class="bench-table"></div>

## Reproduce

Build fasterp and install fastp 1.4.0 and hyperfine, then run on an otherwise idle machine:

```bash
cargo build --release
CSV=book/src/benchmarks.csv bash scripts/bench.sh fastp target/release/fasterp
```

The results are written to `book/src/benchmarks.csv`, which this page reads. Rows for the same dataset are replaced and other datasets are kept, so set `R1`, `R2` and `DATASET` to add or update another dataset.

Programs that never change are measured once and their rows are kept as snapshots. That covers fastp, which takes most of the time, and any fasterp release whose binary is named in `RELEASES`. When the CSV already has a row for the same program version, dataset, workload and worker count, measured on the same machine, the script keeps that row and does not run the program. It still compares fasterp's output digest with the fastp snapshot's, and times fastp again if they differ, which is what happens when the inputs or `ARGS` have changed. `REFERENCE=rerun` times everything regardless.

The build under test is always timed. Its row is labelled with the version and the commit, as in `0.3.0+13331c6`, unless the commit is that version's tag. A release gets a row only where its output digest equals the build under test's, so every row of a group did the same work. Release rows stay in the CSV when their binary is not given, so the history can be plotted later. This page charts fastp against the newest fasterp row.

Each row also records the CPU time (user plus system) and the peak resident memory, taken from `/usr/bin/time`. `scripts/bench_plot.py` draws wall time against peak memory for every program in the CSV:

```bash
RELEASES="path/to/fasterp-0.2.1 path/to/fasterp-0.3.0" CSV=book/src/benchmarks.csv bash scripts/bench.sh
python3 scripts/bench_plot.py book/src/benchmarks.csv memory.png --dataset SRR5808766
```

`WORKERS`, `RUNS` and `WARMUP` change the worker counts, the number of timed runs and the number of untimed runs before them, and `ARGS` adds options to every command. For example, for the NovaSeq prefix:

```bash
bash .github/scripts/download_test_data.sh test_data err3534523-2m
R1=test_data/ERR3534523_2M_1.fastq.gz R2=test_data/ERR3534523_2M_2.fastq.gz \
  DATASET=ERR3534523-2M ARGS=-g CSV=book/src/benchmarks.csv bash scripts/bench.sh
```

<script src="benchmarks.js"></script>
