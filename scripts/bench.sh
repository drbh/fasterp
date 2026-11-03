#!/usr/bin/env bash
# Time fastp and fasterp on the same inputs with hyperfine.
#
#   bash scripts/bench.sh [fastp] [fasterp]
#
# Defaults to `fastp` on PATH and target/release/fasterp, and the SRR5808766
# pair from .github/scripts/download_test_data.sh. Writes one row per program,
# workload and worker count to $CSV (default target/bench/results.csv). Rows
# already in $CSV for the same DATASET are replaced, so
#
#   CSV=book/src/benchmarks.csv bash scripts/bench.sh
#
# updates that dataset in the book's benchmarks chapter. Override with R1, R2,
# DATASET, WORKERS (comma separated), RUNS and OUT. ARGS is added to every
# command for both programs, for example ARGS=-g to trim poly-G tails.
set -euo pipefail

fastp=${1:-fastp}
fasterp=${2:-target/release/fasterp}
r1=${R1:-test_data/SRR5808766_1.fastq.gz}
r2=${R2:-test_data/SRR5808766_2.fastq.gz}
dataset=${DATASET:-$(basename "$r1" | sed 's/_1\..*//')}
workers=${WORKERS:-1,4,8}
runs=${RUNS:-3}
out=${OUT:-target/bench}
extra=${ARGS:-}
csv=${CSV:-$out/results.csv}

command -v hyperfine >/dev/null || { echo "hyperfine is required (brew install hyperfine)" >&2; exit 1; }
tmp=$out/tmp
mkdir -p "$tmp"
cpus=$(getconf _NPROCESSORS_ONLN)
load() { uptime | sed 's/.*load averages*: //; s/,* .*//'; }
# Timings are only comparable on an idle machine.
if [ "$(load | cut -d. -f1)" -ge $((cpus / 2)) ]; then
  echo "warning: load average $(load) on $cpus cpus, results will be noisy" >&2
fi
cpu=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || sed -n 's/^model name\s*: //p' /proc/cpuinfo | head -1)
machine="$cpu $cpus cpus $(uname -s)"
fastp_version=$("$fastp" --version 2>&1 | awk '{print $2}')
fasterp_version=$("$fasterp" --version | awk '{print $2}')
date=$(date +%F)

rows=$tmp/rows.csv
: >"$rows"

bench() {
  local name=$1 args=$2 w
  for w in ${workers//,/ }; do
    hyperfine --warmup 1 --runs "$runs" --export-csv "$tmp/run.csv" \
      -n "fastp $w" "$fastp $args $extra -j $tmp/fastp.json -h $tmp/fastp.html -w $w" \
      -n "fasterp $w" "$fasterp $args $extra -j $tmp/fasterp.json -h $tmp/fasterp.html -w $w"
    # hyperfine columns: command,mean,stddev,median,user,system,min,max
    tail -n +2 "$tmp/run.csv" | while IFS=, read -r command mean stddev _ _ _ min max; do
      program=${command% *}
      version=$fastp_version
      [ "$program" = fasterp ] && version=$fasterp_version
      printf '%s,%s,%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%s,%s,%s\n' "$dataset" "$name" "$w" "$program" \
        "$version" "$mean" "$stddev" "$min" "$max" "$runs" "$machine" "$date" >>"$rows"
    done
  done
}

bench single-end "-i $r1 -o $tmp/se.fq.gz"
bench paired-end "-i $r1 -I $r2 -o $tmp/r1.fq.gz -O $tmp/r2.fq.gz"
bench paired-correct-merge "-i $r1 -I $r2 -c -m --merged_out $tmp/merged.fq.gz -o $tmp/r1.fq.gz -O $tmp/r2.fq.gz"

# keep other datasets' rows, replace this one's
{
  echo dataset,workload,workers,program,version,mean,stddev,min,max,runs,machine,date
  if [ -f "$csv" ]; then tail -n +2 "$csv" | awk -F, -v d="$dataset" '$1 != d'; fi
  cat "$rows"
} >"$tmp/merged.csv"
mkdir -p "$(dirname "$csv")"
mv "$tmp/merged.csv" "$csv"
echo "wrote $(wc -l <"$rows" | tr -d ' ') rows for $dataset to $csv"
