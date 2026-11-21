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
# DATASET, WORKERS (comma separated), RUNS, WARMUP and OUT. ARGS is added to every
# command for both programs, for example ARGS=-g to trim poly-G tails.
#
# Rows of programs that never change are snapshots. That is fastp, the
# reference, which takes most of the time, and any fasterp release named in
# RELEASES (binaries, space separated). A row measured with the same program
# version on the same machine, for the same dataset, workload and worker count,
# is kept and the program is not run again. Release snapshots stay in $CSV when
# their binary is not given. A release gets a row only where its output equals
# the build under test's. fasterp's output must still match the fastp
# snapshot's digest; if it does not, the snapshot describes other inputs or
# options and fastp is timed again. REFERENCE=rerun times everything regardless.
#
# The build under test is always timed. Its row is labelled with the version
# and, unless HEAD is that version's tag, the commit, as in 0.3.0+13331c6. A
# version with a + is never kept as a snapshot. LABEL overrides the label.
#
# Each row also has the CPU time (user plus system, mean per run) and the peak
# resident memory (largest of the runs), the latter from /usr/bin/time.
set -euo pipefail

fastp=${1:-fastp}
fasterp=${2:-target/release/fasterp}
r1=${R1:-test_data/SRR5808766_1.fastq.gz}
r2=${R2:-test_data/SRR5808766_2.fastq.gz}
dataset=${DATASET:-$(basename "$r1" | sed 's/_1\..*//')}
workers=${WORKERS:-1,4,8}
runs=${RUNS:-3}
warmup=${WARMUP:-1}
out=${OUT:-target/bench}
extra=${ARGS:-}
csv=${CSV:-$out/results.csv}
reference=${REFERENCE:-snapshot}
releases=${RELEASES:-}
header=dataset,workload,workers,program,version,mean,stddev,min,max,runs,machine,date,output_bytes,output_sha256,cpu_seconds,peak_mb

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
label=${LABEL:-}
if [ -z "$label" ]; then
  label=$fasterp_version
  commit=$(git rev-parse --short HEAD 2>/dev/null || true)
  if [ -n "$commit" ] && [ "$(git describe --tags --exact-match 2>/dev/null || true)" != "v$fasterp_version" ]; then
    label="$fasterp_version+$commit"
  fi
fi

# BSD time reports the peak resident set in bytes with -l, GNU time in
# kilobytes with -v. Without either the memory column stays empty.
timer= per_mb=1
if /usr/bin/time -l true >/dev/null 2>&1; then
  timer="/usr/bin/time -l" per_mb=1048576
elif /usr/bin/time -v true >/dev/null 2>&1; then
  timer="/usr/bin/time -v" per_mb=1024
fi

rows=$tmp/rows.csv
: >"$rows"

# Total size of the gzip files a program wrote, so output size is recorded
# alongside time and a faster encoder cannot hide larger files.
output_bytes() {
  find "$1" -name '*.gz' -exec wc -c {} + | awk '$2 != "total" {s += $1} END {print s + 0}'
}

sha256() {
  if command -v sha256sum >/dev/null; then sha256sum | cut -c1-64; else shasum -a 256 | cut -c1-64; fi
}

# Digest of the decompressed outputs: each file's name and SHA-256, in name
# order, hashed together. Equal digests mean both programs wrote the same
# records byte for byte, whatever their compression.
output_sha256() {
  local f
  for f in $(cd "$1" && ls *.gz 2>/dev/null | sort); do
    echo "$f $(gzip -dc "$1/$f" | sha256)"
  done | sha256 | cut -c1-16
}

# The row already in $csv for a program version, workload and worker count,
# measured on this machine with memory recorded, or nothing.
snapshot() {
  [ "$reference" = snapshot ] && [ -f "$csv" ] || return 0
  awk -F, -v d="$dataset" -v p="$1" -v v="$2" -v n="$3" -v w="$4" -v m="$machine" \
    '$1 == d && $2 == n && $3 == w && $4 == p && $5 == v && $11 == m && $16 != ""' "$csv" | tail -1
}

# Largest resident set of the runs logged by $timer, in megabytes.
peak_mb() {
  [ -n "$timer" ] || return 0
  awk -v per="$per_mb" '/[Mm]aximum resident set size/ {
    v = ($1 ~ /^[0-9]+$/) ? $1 : $NF; if (v > m) m = v
  } END { printf "%.0f", m / per }' "$tmp/time.log"
}

# Time one program and print its row. Arguments name output files as @/file,
# and every program writes to the same directory, emptied first. hyperfine
# columns: command,mean,stddev,median,user,system,min,max
time_program() {
  local program=$1 version=$2 name=$3 w=$4 command=${5//@/$tmp/out}
  rm -rf "$tmp/out"
  mkdir -p "$tmp/out"
  : >"$tmp/time.log"
  hyperfine --warmup "$warmup" --runs "$runs" --export-csv "$tmp/run.csv" -n "$program $version -w $w" \
    "${timer:+$timer -a -o $tmp/time.log }$command $extra --json $tmp/report.json --html $tmp/report.html -w $w" >&2 || return 1
  [ "$(output_bytes "$tmp/out")" -gt 0 ] || return 1
  tail -n +2 "$tmp/run.csv" | while IFS=, read -r _ mean stddev _ user system min max; do
    printf '%s,%s,%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%s,%s,%s,%s,%s,%.2f,%s\n' "$dataset" "$name" "$w" "$program" \
      "$version" "$mean" "$stddev" "$min" "$max" "$runs" "$machine" "$date" \
      "$(output_bytes "$tmp/out")" "$(output_sha256 "$tmp/out")" \
      "$(echo "$user + $system" | bc -l)" "$(peak_mb)"
  done
}

bench() {
  local name=$1 args=$2 w ours kept release version given
  for w in ${workers//,/ }; do
    ours=$(time_program fasterp "$label" "$name" "$w" "$fasterp $args")
    kept=$(snapshot fastp "$fastp_version" "$name" "$w")
    if [ -n "$kept" ] && [ "$(echo "$kept" | cut -d, -f14)" != "$(echo "$ours" | cut -d, -f14)" ]; then
      echo "$dataset $name -w $w: output differs from the fastp snapshot, timing fastp again" >&2
      kept=
    fi
    if [ -n "$kept" ]; then
      reused=$((reused + 1))
    else
      kept=$(time_program fastp "$fastp_version" "$name" "$w" "$fastp $args")
    fi
    echo "$kept" >>"$rows"
    given=" $label "
    for release in $releases; do
      version=$("$release" --version | awk '{print $2}')
      case $given in *" $version "*) continue ;; esac
      given="$given$version "
      kept=$(snapshot fasterp "$version" "$name" "$w")
      if [ -n "$kept" ]; then
        reused=$((reused + 1))
      else
        kept=$(time_program fasterp "$version" "$name" "$w" "$release $args") || {
          echo "$dataset $name -w $w: fasterp $version does not run this workload, no row written" >&2
          continue
        }
      fi
      # a release that writes other records did other work, and is not comparable
      if [ "$(echo "$kept" | cut -d, -f14)" != "$(echo "$ours" | cut -d, -f14)" ]; then
        echo "$dataset $name -w $w: fasterp $version writes different output, no row written" >&2
        continue
      fi
      echo "$kept" >>"$rows"
    done
    # snapshots of releases whose binary was not given
    if [ -f "$csv" ]; then
      awk -F, -v d="$dataset" -v n="$name" -v w="$w" -v given="$given" \
        '$1 == d && $2 == n && $3 == w && $4 == "fasterp" && $5 !~ /\+/ && $16 != "" && !index(given, " " $5 " ")' \
        "$csv" >>"$rows"
    fi
    echo "$ours" >>"$rows"
  done
}

reused=0
bench single-end "-i $r1 -o @/se.fq.gz"
bench paired-end "-i $r1 -I $r2 -o @/r1.fq.gz -O @/r2.fq.gz"
bench paired-correct-merge "-i $r1 -I $r2 -c -m --merged_out @/merged.fq.gz -o @/r1.fq.gz -O @/r2.fq.gz"

# keep other datasets' rows, replace this one's
{
  echo "$header"
  if [ -f "$csv" ]; then
    tail -n +2 "$csv" | awk -F, -v OFS=, -v d="$dataset" '$1 != d { if (NF == 14) $16 = ""; print }'
  fi
  cat "$rows"
} >"$tmp/merged.csv"
mkdir -p "$(dirname "$csv")"
mv "$tmp/merged.csv" "$csv"
echo "wrote $(wc -l <"$rows" | tr -d ' ') rows for $dataset to $csv ($reused snapshots reused)"
