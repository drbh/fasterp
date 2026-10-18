#!/usr/bin/env python3
"""Add reproducible stream cases using a prepare_common_benchmark.py manifest.

Existing unchanged WGS mates and gzip copies are reused after checksum checks.
Stdin/stdout are redirected files; producer/consumer timing is not included.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--common", type=Path, required=True)
    parser.add_argument("--large-wgs", type=Path,
                        help="Optional two-million-pair SRR1766443 manifest; adds paired file and RNA-seq controls")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.common.read_text())
    large = json.loads(args.large_wgs.read_text()) if args.large_wgs else None
    if large and (large.get("accession") != "SRR1766443" or large.get("records_per_mate") != 2_000_000):
        parser.error("Expected two million SRR1766443 pairs for --large-wgs")
    root = args.out.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Use an empty output directory")
    cases = {c["name"]: c for c in source["cases"]}
    plain = cases["srr1766443_paired_t1"]["inputs"]
    compressed = cases["human_gzip_io_t1"]["inputs"]
    records = cases["srr1766443_paired_t1"]["records_per_mate"]
    files = {}
    for path in plain + compressed:
        info = source["files"][path]
        if digest(path) != info["sha256"]:
            parser.error(f"Input checksum differs from manifest: {path}")
        files[path] = info
    if large:
        for path in large["cases"][0]["inputs"]:
            info = large["files"][path]
            if digest(path) != info["sha256"]:
                parser.error(f"Input checksum differs from manifest: {path}")
            files[path] = info
    root.mkdir(parents=True, exist_ok=True)
    interleaved = root / "interleaved.fq"
    with open(plain[0], "rb") as first, open(plain[1], "rb") as second, interleaved.open("wb") as output:
        for index in range(records):
            for mate in (first, second):
                lines = [mate.readline() for _ in range(4)]
                if not lines[0].startswith(b"@") or not lines[2].startswith(b"+") or not lines[3]:
                    raise ValueError(f"Invalid input record {index + 1}")
                output.write(b"".join(lines))
        if first.read(1) or second.read(1):
            raise ValueError("Manifest record count does not cover the input")
    files[str(interleaved)] = dict(sources=plain, records=records * 2,
                                  bytes=interleaved.stat().st_size, sha256=digest(interleaved),
                                  transformation="Unchanged alternating R1/R2 records")
    matrix = []
    for workers in (1, 4, 8):
        for name, inputs, stdin, stdout, options in [
            ("stdin_file", plain[:1], True, False, ["-A", "-G"]),
            ("stdin_stdout", plain[:1], True, True, ["-A", "-G"]),
            ("paired_stdout", plain, False, True, []),
            ("gzip_paired_stdout", compressed, False, True, []),
            ("interleaved_stdin_stdout", [str(interleaved)], True, True, ["--interleaved_in", "-A", "-G"]),
            ("merge_stdout", plain, False, True, ["-m", "--include_unmerged"]),
        ]:
            case = dict(name=f"{name}_t{workers}", inputs=inputs, threads=workers,
                        records_per_mate=records, stdin=stdin, options=options,
                        outputs=[] if stdout else [["-o", "out.fq"]])
            if stdout:
                case["stdout"] = "stdout.fq"
            matrix.append(case)
    datasets = [d for d in source["datasets"] if d.get("accession") == "SRR1766443"]
    if large:
        datasets += [d for d in source["datasets"] if d.get("accession") == "SRR13493714"]
        datasets.append({k: v for k, v in large.items() if k not in ("files", "cases")})
        large_plain = large["cases"][0]["inputs"]
        large_gzip = []
        for mate, path in enumerate(large_plain, 1):
            target = root / f"large_R{mate}.fq.gz"
            with open(path, "rb") as raw, target.open("wb") as packed:
                with gzip.GzipFile(filename="", mode="wb", fileobj=packed, compresslevel=4, mtime=0) as output:
                    for block in iter(lambda: raw.read(4 * 1024 * 1024), b""):
                        output.write(block)
            files[str(target)] = dict(source_uncompressed=path, records=large["records_per_mate"],
                                      bytes=target.stat().st_size, sha256=digest(target), compression_level=4)
            large_gzip.append(str(target))
        for workers in (1, 4, 8):
            for name, inputs, gz in [("wgs2m_plain", large_plain, False),
                                     ("wgs2m_gzip_input", large_gzip, False),
                                     ("wgs2m_gzip_io", large_gzip, True)]:
                ext = ".fq.gz" if gz else ".fq"
                matrix.append(dict(name=f"{name}_t{workers}", inputs=inputs, threads=workers,
                                   records_per_mate=large["records_per_mate"], options=[],
                                   outputs=[["-o", "out1" + ext], ["-O", "out2" + ext]]))
            matrix.append(dict(name=f"wgs2m_single_gzip_io_t{workers}", inputs=large_gzip[:1],
                               threads=workers, records_per_mate=large["records_per_mate"],
                               options=[], outputs=[["-o", "out1.fq.gz"]]))
            control = cases[f"rna_trim_t{workers}"]
            matrix.append(control)
            for path in control["inputs"]:
                if path not in files and digest(path) != source["files"][path]["sha256"]:
                    parser.error(f"Input checksum differs from manifest: {path}")
                files[path] = source["files"][path]
    result = dict(purpose=__doc__, datasets=datasets, files=files, cases=matrix)
    (root / "plan.json").write_text(json.dumps(result, indent=2) + "\n")
    print(root / "plan.json")


if __name__ == "__main__":
    main()
