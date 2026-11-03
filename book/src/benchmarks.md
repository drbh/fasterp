# Benchmarks

Wall time in seconds for fastp and fasterp on the same inputs and options, measured with [hyperfine](https://github.com/sharkdp/hyperfine). Lower is better. The label next to each pair is fastp time divided by fasterp time, and the whiskers span the fastest and slowest run.

<div class="bench-legend"><span class="fastp">fastp</span><span class="fasterp">fasterp</span></div>
<div id="bench-chart"></div>
<p id="bench-meta" class="bench-meta"></p>

## Workloads

All workloads read gzip input and write gzip output, with the default filters and an HTML and JSON report.

| Workload | Options |
|---|---|
| single-end | `-i R1 -o out.fq.gz` |
| paired-end | `-i R1 -I R2 -o out1.fq.gz -O out2.fq.gz` |
| paired-correct-merge | paired-end plus `-c -m --merged_out merged.fq.gz` |

## Datasets

- **SRR1766443-2M**: the first 2,000,000 records per mate of human WGS run [SRR1766443](https://www.ebi.ac.uk/ena/browser/view/SRR1766443), 2x148 bp, recompressed with gzip level 4
- **SRR5808766**: the full mouse ATAC-seq run used in the [tutorial](./tutorial.md), 9.8 million pairs of 51 bp

## Results

<div id="bench-table" class="bench-table"></div>

## Reproduce

Build fasterp and install fastp 1.4.0 and hyperfine, then run on an otherwise idle machine:

```bash
cargo build --release
CSV=book/src/benchmarks.csv bash scripts/bench.sh fastp target/release/fasterp
```

The results are written to `book/src/benchmarks.csv`, which this page reads. Rows for the same dataset are replaced and other datasets are kept, so set `R1`, `R2` and `DATASET` to add or update another dataset. `WORKERS` and `RUNS` change the worker counts and the number of timed runs.

<script src="benchmarks.js"></script>
