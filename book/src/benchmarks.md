# Benchmarks

Wall time in seconds for fastp and fasterp on the same inputs and options, measured with [hyperfine](https://github.com/sharkdp/hyperfine). Lower is better. The label next to each pair is fastp time divided by fasterp time, and the whiskers span the fastest and slowest run.

<div class="bench-legend"><span class="fastp">fastp</span><span class="fasterp">fasterp</span></div>
<div id="bench-chart"></div>
<p id="bench-meta" class="bench-meta"></p>

fasterp `0.2.1-dev` is the development build ahead of the next release. Its gzip output is new: blocks with binned qualities, as NovaSeq writes them, use a FASTQ-aware encoder that is as small as level 4 at about 2.4 times its speed, and other blocks are compressed as before.

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

ERR3534523 was converted from NCBI's copy with `fasterq-dump`, which does not keep the original Illumina read names, so its headers look like `@ERR3534523.1 1/1`. The sequences and qualities of the first 2,000,000 records of mate 2 match ENA's file. Read names starting with `@ERR` also mean fastp does not detect a two-color instrument, so poly-G trimming has to be enabled with `-g` (`ARGS=-g` below).

## Results

<div id="bench-table" class="bench-table"></div>

## Reproduce

Build fasterp and install fastp 1.4.0 and hyperfine, then run on an otherwise idle machine:

```bash
cargo build --release
CSV=book/src/benchmarks.csv bash scripts/bench.sh fastp target/release/fasterp
```

The results are written to `book/src/benchmarks.csv`, which this page reads. Rows for the same dataset are replaced and other datasets are kept, so set `R1`, `R2` and `DATASET` to add or update another dataset. `WORKERS` and `RUNS` change the worker counts and the number of timed runs, and `ARGS` adds options to every command. For example, for the NovaSeq prefix:

```bash
bash .github/scripts/download_test_data.sh test_data err3534523-2m
R1=test_data/ERR3534523_2M_1.fastq.gz R2=test_data/ERR3534523_2M_2.fastq.gz \
  DATASET=ERR3534523-2M ARGS=-g CSV=book/src/benchmarks.csv bash scripts/bench.sh
```

<script src="benchmarks.js"></script>
