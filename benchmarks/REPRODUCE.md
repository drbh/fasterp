# Reproduce the benchmarks

Build fastp 1.4.0 using the [reference-build instructions](../CONTRIBUTING.md). The timing runner uses macOS `/usr/bin/time -lp`.

```sh
cargo build --locked --release
mkdir -p target/benchmarks
cp target/release/fasterp target/benchmarks/fasterp
python3 scripts/prepare_ena_subset.py --accession SRR1766443 --records 500000 --out target/benchmarks/wgs
python3 scripts/prepare_ena_subset.py --accession SRR13493714 --records 500000 --out target/benchmarks/rna
python3 scripts/prepare_ena_subset.py --accession SRR5808766 --records 2000000 --out target/benchmarks/short
python3 scripts/prepare_common_benchmark.py --wgs target/benchmarks/wgs/plan.json \
  --rna target/benchmarks/rna/plan.json --short target/benchmarks/short/plan.json \
  --out target/benchmarks/common
python3 scripts/prepare_ena_subset.py --accession SRR1766443 --records 2000000 --out target/benchmarks/wgs2m
python3 scripts/prepare_stream_benchmark.py --common target/benchmarks/common/plan.json \
  --large-wgs target/benchmarks/wgs2m/plan.json --out target/benchmarks/streams
python3 scripts/benchmark_fastp.py --plan target/benchmarks/streams/plan.json \
  --program fastp=target/fastp-reference/fastp --program fasterp=target/benchmarks/fasterp \
  --repeats 3 --out target/benchmarks/results
```

The common-workload helper also prepares a short-read control; the 33-case stream matrix uses only WGS and RNA-seq. Helpers check input hashes and require empty output directories. Run without other builds or benchmarks competing for resources.

The published JSON contains every fastp and fasterp run from this matrix, with local paths replaced by portable labels. Timing values, input/output hashes and measured binary hashes are unchanged.
