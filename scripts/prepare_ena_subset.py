#!/usr/bin/env python3
"""Download a bounded, unchanged paired FASTQ subset and record its provenance."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--records", type=int, default=500_000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.records <= 0:
        parser.error("--records must be positive")
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        parser.error("Use an empty output directory to preserve previous inputs and manifests")
    query = urllib.parse.urlencode(dict(
        accession=args.accession, result="read_run", format="json",
        fields="run_accession,scientific_name,instrument_model,library_strategy,library_layout,read_count,base_count,fastq_ftp,fastq_bytes"))
    metadata_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query
    with urllib.request.urlopen(metadata_url, timeout=60) as response:
        records = json.load(response)
    if len(records) != 1 or records[0]["library_layout"] != "PAIRED":
        parser.error("Expected one paired sequencing run")
    metadata = records[0]
    urls = metadata["fastq_ftp"].split(";")
    mates = []
    for mate in (1, 2):
        matches = [u for u in urls if u.endswith(f"_{mate}.fastq.gz")]
        if len(matches) != 1:
            parser.error(f"Could not identify mate {mate} unambiguously")
        mates.append((mate, "https://" + matches[0]))

    def download(item):
        mate, url = item
        target = root / f"real_R{mate}.fq"
        partial = target.with_suffix(".part")
        digest = hashlib.sha256()
        bases, minimum, maximum = 0, None, 0
        request = urllib.request.Request(url, headers={"Accept-Encoding": "identity"})
        with urllib.request.urlopen(request, timeout=60) as response:
            with gzip.GzipFile(fileobj=response) as source, partial.open("wb") as output:
                for index in range(args.records):
                    lines = [source.readline() for _ in range(4)]
                    length = len(lines[1].rstrip(b"\r\n"))
                    if not lines[0].startswith(b"@") or not lines[2].startswith(b"+") or not lines[3] or length != len(lines[3].rstrip(b"\r\n")):
                        raise ValueError(f"Invalid or incomplete mate {mate}, record {index + 1}")
                    data = b"".join(lines)
                    output.write(data)
                    digest.update(data)
                    bases += length
                    minimum = length if minimum is None else min(minimum, length)
                    maximum = max(maximum, length)
        partial.replace(target)
        result = dict(source=url, records=args.records, bytes=target.stat().st_size,
                      sha256=digest.hexdigest(), bases=bases,
                      min_length=minimum, max_length=maximum)
        print(f"Prepared mate {mate}: {args.records} records, lengths {minimum}-{maximum}", flush=True)
        return str(target), result

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        files = dict(pool.map(download, mates))
    cases = []
    for workers in (1, 4, 8):
        cases.append(dict(name=f"{args.accession.lower()}_paired_t{workers}",
                          inputs=list(files), threads=workers, options=[],
                          outputs=[["-o", "out1.fq"], ["-O", "out2.fq"]]))
    plan = dict(accession=args.accession, records_per_mate=args.records,
                selection=f"First {args.records} unchanged records per mate; download streams closed after subset",
                metadata_url=metadata_url, metadata=metadata, files=files, cases=cases)
    (root / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    print(root / "plan.json")


if __name__ == "__main__":
    main()
