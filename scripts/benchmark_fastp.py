#!/usr/bin/env python3
"""Reproducible native benchmark with child-only RSS and FASTQ fingerprints.

Example:
  python3 scripts/benchmark_fastp.py --plan target/benchmarks/plan.json \
    --program fastp=../fastp-upstream/fastp \
    --program fasterp=target/benchmarks/baseline/fasterp --repeats 3

Runs are serial, program order rotates, and inputs are read into the OS cache
before each case. Validation and file hashing run outside the timed child.
"""
import argparse
from contextlib import ExitStack
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import time
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path):
    ordered = hashlib.sha256()
    total = 0
    xor = 0
    count = 0
    bases = 0
    modulus = 1 << 256
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as stream:
        while True:
            header = stream.readline()
            if not header:
                break
            seq = stream.readline()
            plus = stream.readline()
            qual = stream.readline()
            if not header.startswith(b"@") or not plus.startswith(b"+") or len(seq.rstrip(b"\r\n")) != len(qual.rstrip(b"\r\n")):
                raise ValueError(f"Malformed output in {path}, record {count+1}")
            record = header + seq + plus + qual
            ordered.update(record)
            item = int.from_bytes(hashlib.blake2b(record, digest_size=32).digest(), "big")
            total = (total + item) % modulus
            xor ^= item
            count += 1
            bases += len(seq.rstrip(b"\r\n"))
    return {"records": count, "bases": bases, "ordered_sha256": ordered.hexdigest(),
            "multiset_sum_blake2b256": f"{total:064x}", "multiset_xor_blake2b256": f"{xor:064x}"}


def compare_reports(reference, actual):
    differences = []
    for stage in ("before_filtering", "after_filtering"):
        for key, expected in reference["summary"][stage].items():
            observed = actual["summary"][stage].get(key)
            if isinstance(expected, (int, float)):
                matches = observed is not None and (expected == observed if isinstance(expected, int) else math.isclose(expected, observed, rel_tol=1e-5, abs_tol=1e-6))
                if not matches:
                    differences.append(f"summary.{stage}.{key}: {expected} != {observed}")
    for section in ("filtering_result", "adapter_cutting", "duplication"):
        for key, expected in reference.get(section, {}).items():
            observed = actual.get(section, {}).get(key)
            if isinstance(expected, (int, float)):
                matches = observed is not None and (expected == observed if isinstance(expected, int) else math.isclose(expected, observed, rel_tol=1e-5, abs_tol=1e-6))
                if not matches:
                    differences.append(f"{section}.{key}: {expected} != {observed}")
    for section in ("read1_before_filtering", "read2_before_filtering", "read1_after_filtering", "read2_after_filtering", "merged_and_filtered"):
        if section not in reference:
            continue
        for key, expected in reference[section].items():
            observed = actual.get(section, {}).get(key)
            if key in ("quality_curves", "content_curves"):
                for base, curve in expected.items():
                    other = (observed or {}).get(base, [])
                    if len(curve) != len(other) or any(not math.isclose(x, y, rel_tol=1e-5, abs_tol=1e-6) for x, y in zip(curve, other)):
                        differences.append(f"{section}.{key}.{base} differs")
            elif expected != observed:
                differences.append(f"{section}.{key} differs")
    return differences


def parse_time(path):
    text = path.read_text()
    result = {}
    for name in ("real", "user", "sys"):
        match = re.search(rf"^{name}\s+([0-9.]+)$", text, re.M)
        if not match:
            raise ValueError(f"Missing {name} in {path}: {text}")
        result[name + "_seconds"] = float(match.group(1))
    match = re.search(r"^\s*(\d+)\s+maximum resident set size", text, re.M)
    if not match:
        raise ValueError(f"Missing RSS in {path}")
    result["peak_rss_bytes"] = int(match.group(1))
    match = re.search(r"^\s*(\d+)\s+peak memory footprint", text, re.M)
    if match:
        result["peak_footprint_bytes"] = int(match.group(1))
    return result


def warm_inputs(case):
    for path in case["inputs"]:
        with open(path, "rb") as stream:
            while stream.read(4 * 1024 * 1024):
                pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--program", action="append", required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cases", help="Comma-separated case names")
    parser.add_argument("--out", type=Path, default=Path("target/benchmarks/results"))
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if platform.system() != "Darwin":
        parser.error("This runner uses macOS /usr/bin/time -lp (RSS in bytes).")
    programs = [(name, str(Path(path).resolve())) for name, path in (item.split("=", 1) for item in args.program)]
    if len({name for name, _ in programs}) != len(programs):
        parser.error("Program names must be unique")
    plan = json.loads(args.plan.read_text())
    if args.cases and set(args.cases.split(",")) - {c["name"] for c in plan["cases"]}:
        parser.error("--cases contains an unknown case name")
    cases = [c for c in plan["cases"] if not args.cases or c["name"] in args.cases.split(",")]
    if not cases:
        parser.error("No benchmark cases selected")
    root = args.out.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Use an empty output directory to preserve earlier measurements")
    root.mkdir(parents=True, exist_ok=True)
    metadata = {"platform": platform.platform(), "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "repeats": args.repeats, "input_cache": "warm", "program_order": "rotating positions; reversed on alternate cases", "measurement": "macOS /usr/bin/time -lp, child maximum RSS in bytes",
                "programs": [{"name": name, "binary": path, "sha256": sha256(path), "version": subprocess.check_output([path, "--version"], text=True).strip()} for name, path in programs],
                "plan": plan, "machine": "Apple silicon"}
    (root / "metadata.json").write_text(json.dumps(metadata, indent=2))
    runs = []
    accuracy = {}
    for case_index, case in enumerate(cases):
        print(f"CASE {case['name']}", flush=True)
        warm_inputs(case)
        for repeat in range(args.repeats):
            offset = (repeat + case_index) % len(programs)
            order = programs[offset:] + programs[:offset]
            if case_index % 2:
                order = list(reversed(order))
            for name, binary in order:
                folder = root / case["name"] / f"{name}-{repeat}"
                folder.mkdir(parents=True, exist_ok=True)
                command = [binary, "-w", str(case["threads"]), "-j", "report.json", "-h", "report.html"]
                command += ["--stdin"] if case.get("stdin") else ["-i", case["inputs"][0]]
                if len(case["inputs"]) == 2:
                    command += ["-I", case["inputs"][1]]
                for flag, filename in case["outputs"]:
                    command += [flag, filename]
                command += case.get("options", [])
                command += case.get("program_options", {}).get(name, [])
                output_files = [filename for _, filename in case["outputs"]]
                if case.get("stdout"):
                    command += ["--stdout"]
                    output_files.append(case["stdout"])
                (folder / "command.json").write_text(json.dumps(command))
                # Redirection is opened outside the timed child. This measures
                # stdin/stdout handling, not the cost of an external producer.
                with ExitStack() as stack:
                    stderr = stack.enter_context(open(folder / "stderr.log", "wb"))
                    stdin = stack.enter_context(open(case["inputs"][0], "rb")) if case.get("stdin") else subprocess.DEVNULL
                    stdout = stack.enter_context(open(folder / case["stdout"], "wb")) if case.get("stdout") else subprocess.DEVNULL
                    result = subprocess.run(["/usr/bin/time", "-lp", "-o", str(folder / "time.txt"), *command],
                                            cwd=folder, stdin=stdin, stdout=stdout, stderr=stderr, timeout=600)
                if result.returncode:
                    raise RuntimeError(f"{name} failed in {case['name']}: {(folder/'stderr.log').read_text()[-2000:]}")
                measurement = parse_time(folder / "time.txt")
                measurement.update({"case": case["name"], "program": name, "repeat": repeat,
                                    "threads": case["threads"], "command": command,
                                    "output_bytes": {filename: (folder / filename).stat().st_size for filename in output_files}})
                runs.append(measurement)
                print(f"  {name} #{repeat+1}: {measurement['real_seconds']:.2f}s, {measurement['peak_rss_bytes']/2**20:.1f} MiB", flush=True)
                if repeat == 0:
                    accuracy.setdefault(case["name"], {})[name] = {
                        "outputs": {filename: fingerprint(folder / filename) for filename in output_files},
                        "report": json.loads((folder / "report.json").read_text())}
                for filename in output_files:
                    (folder / filename).unlink()
                (root / "runs.json").write_text(json.dumps(runs, indent=2))
        (root / "accuracy.json").write_text(json.dumps(accuracy, indent=2))
    summary = []
    for case in cases:
        reference = accuracy[case["name"]][programs[0][0]]
        for name, _ in programs:
            values = [r for r in runs if r["case"] == case["name"] and r["program"] == name]
            actual = accuracy[case["name"]][name]
            same_records = all(all(ref[k] == actual["outputs"][filename][k] for k in ("records", "bases", "multiset_sum_blake2b256", "multiset_xor_blake2b256")) for filename, ref in reference["outputs"].items())
            same_order = all(ref["ordered_sha256"] == actual["outputs"][filename]["ordered_sha256"] for filename, ref in reference["outputs"].items())
            differences = compare_reports(reference["report"], actual["report"])
            summary.append({"case": case["name"], "program": name, "threads": case["threads"],
                            "median_seconds": statistics.median(r["real_seconds"] for r in values),
                            "min_seconds": min(r["real_seconds"] for r in values), "max_seconds": max(r["real_seconds"] for r in values),
                            "median_peak_rss_mib": statistics.median(r["peak_rss_bytes"] / 2**20 for r in values),
                            "max_peak_rss_mib": max(r["peak_rss_bytes"] / 2**20 for r in values),
                            "median_cpu_seconds": statistics.median(r["user_seconds"] + r["sys_seconds"] for r in values),
                            "median_output_bytes": statistics.median(sum(r["output_bytes"].values()) for r in values),
                            "records_match": same_records, "order_matches": same_order,
                            "reports_match": not differences, "report_differences": differences})
    (root / "summary.json").write_text(json.dumps(summary, indent=2))
    with open(root / "summary.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(json.dumps(summary, indent=2), flush=True)
    if any(not (row["records_match"] and row["order_matches"] and row["reports_match"]) for row in summary):
        raise SystemExit("Benchmark validation failed; inspect summary.json before using timings.")


if __name__ == "__main__":
    main()
