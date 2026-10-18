#!/usr/bin/env python3
"""Build the common-workload matrix from three prepare_ena_subset.py manifests.

The input subsets remain unchanged. Only the gzip copies use new compression.
The output directory must be empty; downloads are handled by the subset helper.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wgs", type=Path, required=True, help="SRR1766443 plan.json")
    parser.add_argument("--rna", type=Path, required=True, help="SRR13493714 plan.json")
    parser.add_argument("--short", type=Path, required=True, help="SRR5808766 plan.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.out.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Use an empty output directory")
    plans = [json.loads(path.read_text()) for path in (args.wgs, args.rna, args.short)]
    accessions = ("SRR1766443", "SRR13493714", "SRR5808766")
    inputs = []
    files = {}
    for plan, accession in zip(plans, accessions):
        if accession not in (plan.get("accession", "") + plan.get("source", "")):
            parser.error(f"Expected {accession} input provenance")
        mates = [str(Path(p).resolve()) for p in plan["cases"][0]["inputs"]]
        if len(mates) != 2 or any(not Path(p).is_file() for p in mates):
            parser.error("Each manifest must describe two existing FASTQ mates")
        for path in mates:
            info = plan["files"][path]
            digest = hashlib.sha256()
            with open(path, "rb") as source:
                for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest() != info["sha256"]:
                parser.error(f"Input checksum differs from manifest: {path}")
            files[path] = info
        inputs.append(mates)
    root.mkdir(parents=True, exist_ok=True)
    compressed = []
    for label, mates, plan in zip(("human", "rna"), inputs, plans):
        pair = []
        for mate, source_path in enumerate(mates, 1):
            path = root / f"{label}_R{mate}.fq.gz"
            with open(source_path, "rb") as source, open(path, "wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=4, mtime=0) as target:
                    for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                        target.write(block)
            digest = hashlib.sha256()
            with open(path, "rb") as source:
                for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                    digest.update(block)
            files[str(path)] = dict(source_uncompressed=source_path,
                                    records=plan["records_per_mate"], bytes=path.stat().st_size,
                                    sha256=digest.hexdigest(), compression_level=4)
            pair.append(str(path))
        compressed.append(pair)
    cases = []

    def add(name, paths, workers, records, gz=False, options=()):
        ext = ".fq.gz" if gz else ".fq"
        outputs = [["-o", "out1" + ext]]
        if len(paths) == 2:
            outputs.append(["-O", "out2" + ext])
        cases.append(dict(name=f"{name}_t{workers}", inputs=paths, outputs=outputs,
                          threads=workers, records_per_mate=records, options=list(options)))

    for workers in (1, 4, 8):
        wgs_records, rna_records, short_records = [p["records_per_mate"] for p in plans]
        add("human_gzip_io", compressed[0], workers, wgs_records, gz=True)
        add("human_gzip_input", compressed[0], workers, wgs_records)
        add("human_gzip_output", inputs[0], workers, wgs_records, gz=True)
        add("rna_plain", inputs[1], workers, rna_records)
        add("rna_gzip", compressed[1], workers, rna_records, gz=True)
        add("rna_trim", compressed[1], workers, rna_records, gz=True,
            options=("--detect_adapter_for_pe", "-5", "-3", "-x"))
        add("srr1766443_paired", inputs[0], workers, wgs_records)
        add("srr1766443_single", inputs[0][:1], workers, wgs_records)
        add("real_pe_51bp", inputs[2], workers, short_records)
    result = dict(purpose="Common Illumina paired WGS and RNA-seq, gzip and trimming, with plain and single-end controls",
                  datasets=[{k: v for k, v in p.items() if k not in ("files", "cases")} for p in plans],
                  files=files, cases=cases)
    (root / "plan.json").write_text(json.dumps(result, indent=2) + "\n")
    print(root / "plan.json")


if __name__ == "__main__":
    main()
