#!/usr/bin/env python3
"""Compare deterministic synthetic FASTQ cases with the pinned fastp release.

No sequencing downloads or Python packages are needed. Example:
  python3 scripts/compare_fastp.py --fastp ../fastp-upstream/fastp \
      --fasterp target/release/fasterp
"""
import argparse
import gzip
import json
import math
import random
import subprocess
import tempfile
from pathlib import Path


def fixtures(root):
    rng = random.Random(140)
    rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    records = [[], []]
    for i in range(400):
        insert = rng.randrange(15, 250)
        dna = "".join(rng.choices("ACGT", k=insert))
        for mate, seq in enumerate((dna, rc(dna))):
            adapter = ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")[mate]
            seq = (seq + adapter * 5)[:150]
            quality = "I" * len(seq)
            if i % 7 == 0:
                quality = "!" * 8 + quality[8:-12] + "!" * 12
            if i % 11 == 0 and mate == 1:
                quality = "!" * len(seq)
            if i % 13 == 0:
                seq = seq[:20] + "N" * 8 + seq[28:]
            records[mate].append(f"@SEQ:{i} {mate+1}:N:0:ACGT+TGCA\n{seq}\n+\n{quality}\n")
    for i, reads in enumerate(records):
        (root / f"r{i+1}.fq").write_text("".join(reads))
    (root / "duplicate1.fq").write_text("".join(records[0])*2)
    (root / "duplicate2.fq").write_text("".join(records[1])*2)
    (root / "interleaved.fq").write_text("".join(a+b for a,b in zip(*records)))
    raw = "".join(records[0]).encode()
    (root / "multi.fq.gz").write_bytes(gzip.compress(raw[:len(raw)//2]) + gzip.compress(b"") + gzip.compress(raw[len(raw)//2:]))
    # Also test a gzip stream without an extension and with Phred+64 encoding.
    (root / "gzip.data").write_bytes(gzip.compress(raw))
    lines = raw.decode().splitlines()
    for i in range(3, len(lines), 4):
        lines[i] = "".join(chr(ord(q)+31) for q in lines[i])
    (root / "phred64.fq").write_text("\n".join(lines)+"\n")
    (root / "adapters.fa").write_text(">r1\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n>r2\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n")
    (root / "indexes.txt").write_text("ACGT\n")
    (root / "empty.fq").write_bytes(b"")
    for mate in (1, 2):
        (root / f"no_overlap{mate}.fq").write_text(f"@pair/{mate}\n{'A'*80}\n+\n{'I'*80}\n")
    for i in range(2):
        (root / f"large{i+1}.fq").write_text("".join(records[i])*25)
    custom=[]
    for i in range(10000):
        dna="".join(rng.choices("ACGT",k=rng.randrange(25,120)))
        seq=(dna+"TCGATGCTACGATTCGAGTCATGCATGTAG"*6)[:150]
        custom.append(f"@CUSTOM:{i}\n{seq}\n+\n{'I'*len(seq)}\n")
    (root / "custom.fq").write_text("".join(custom))


def cases():
    return [
        ("gzip_out_1", True, ["-z", "1"]),
        ("gzip_out_4", True, ["-z", "4"]),
        ("gzip_out_9", True, ["-z", "9"]),
        ("gzip_out_empty", False, []),
        ("dedup_se",False,["-D","--dup_calc_accuracy","1"]),
        ("dedup_pe",True,["-D","--dup_calc_accuracy","1"]),
        ("large_se_auto",False,[]), ("large_pe_auto",True,["-2"]),
        ("large_discovery",False,[]),
        ("gzip_members",False,[]),
        ("phred64",False,["-6"]), ("empty",False,[]),
        ("interleaved",True,["--interleaved_in"]),
        ("stdout",True,["--stdout"]),
        ("stdout_se", False, ["--stdout"]),
        ("stdout_merge", True, ["--stdout", "-m"]),
        ("stdout_unmerged", True, ["--stdout", "-m", "--include_unmerged"]),
        ("stdout_no_overlap", True, ["--stdout", "-m", "-A", "-G"]),
        ("stdin_se", False, ["--stdin", "--stdout", "-A", "-G"]),
        ("stdin_interleaved", True, ["--stdin", "--stdout", "--interleaved_in", "-A", "-G"]),
        ("stdin_interleaved_trim", True, ["--stdin", "--stdout", "--interleaved_in", "-5", "-3", "-x", "-U", "--umi_loc", "per_read", "--umi_len", "4"]),
        ("stdin_interleaved_merge", True, ["--stdin", "--stdout", "--interleaved_in", "-m", "--include_unmerged"]),
        ("se_default", False, []), ("pe_default", True, []),
        ("se_raw", False, ["-A", "-G"]), ("pe_raw", True, ["-A", "-G"]),
        ("quality_off", True, ["-Q", "-A", "-G"]),
        ("length_limit", True, ["--length_limit", "100"]),
        ("length_off", True, ["-L", "--length_limit", "100"]),
        ("complexity", True, ["-y", "-Y", "60"]),
        ("average_quality", True, ["-e", "35"]),
        ("qualified_quality", True, ["-q", "35", "-u", "10", "-n", "2"]),
        ("fixed_trim", True, ["-f", "8", "-t", "4"]),
        ("explicit_r2_zero", True, ["-f", "8", "-F", "0", "-t", "4", "-T", "0", "-b", "90", "-B", "0"]),
        ("cut_front", True, ["-5"]), ("cut_tail", True, ["-3"]),
        ("cut_right", True, ["-r"]),
        ("cut_both", True, ["-5", "-3", "-W", "5", "-M", "25"]),
        ("cut_overrides", True, ["-5", "-3", "--cut_front_window_size", "2", "--cut_front_mean_quality", "30", "--cut_tail_window_size", "6", "--cut_tail_mean_quality", "25"]),
        ("max_length", True, ["-b", "80"]),
        ("poly_g", True, ["-g", "--poly_g_min_len", "4"]),
        ("poly_x", True, ["-x", "--poly_x_min_len", "4"]),
        ("poly_g_and_x", True, ["-g", "-x", "--poly_g_min_len", "4", "--poly_x_min_len", "6"]),
        ("manual_adapter", True, ["-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "--adapter_sequence_r2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]),
        ("adapter_fasta", True, ["--adapter_fasta", "{adapters.fa}"]),
        ("adapter_dimer", True, ["--adapter_fasta", "{adapters.fa}", "--dimer_max_len", "50", "-L", "-Q"]),
        ("pe_detection", True, ["-2"]),
        ("gap_overlap", True, ["--allow_gap_overlap_trimming"]),
        ("correction", True, ["-c"]),
        ("overlap_thresholds", True, ["--overlap_len_require", "15", "--overlap_diff_limit", "2", "--overlap_diff_percent_limit", "10"]),
        ("merge", True, ["-m", "--merged_out", "merged.fq"]),
        ("merge_after_trim", True, ["-m", "--merged_out", "merged.fq", "-f", "8", "-b", "100", "-3"]),
        ("merge_corrected", True, ["-m", "--merged_out", "merged.fq", "-c"]),
        ("merge_poly_x", True, ["-m", "--merged_out", "merged.fq", "-x", "--poly_x_min_len", "4"]),
        ("merge_manual_adapter", True, ["-m", "--merged_out", "merged.fq", "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "--adapter_sequence_r2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]),
        ("merge_gap_max_length", True, ["-m", "--merged_out", "merged.fq", "--allow_gap_overlap_trimming", "-b", "65", "-B", "80"]),
        ("include_unmerged", True, ["-m", "--merged_out", "merged.fq", "--include_unmerged"]),
        ("failed_after_cut",True,["--unpaired1","unpaired.fq","--failed_out","failed.fq","-5"]),
        ("failed_and_unpaired", True, ["--unpaired1", "unpaired.fq", "--failed_out", "failed.fq"]),
        ("failed_and_two_unpaired", True, ["--unpaired1", "unpaired1.fq", "--unpaired2", "unpaired2.fq", "--failed_out", "failed.fq"]),
        ("overlapped", True, ["--overlapped_out", "overlapped.fq"]),
        ("read_limit", True, ["--reads_to_process", "27"]),
        ("umi_read", True, ["-U", "--umi_loc", "read1", "--umi_len", "6", "--umi_skip", "2"]),
        ("umi_per_read", True, ["-U", "--umi_loc", "per_read", "--umi_len", "6", "--umi_prefix", "UMI", "--umi_delim", "_"]),
        ("umi_index", True, ["-U", "--umi_loc", "index1"]),
        ("umi_per_index", True, ["-U", "--umi_loc", "per_index"]),
        ("index_filter", True, ["--filter_by_index1", "{indexes.txt}"]),
        ("overrepresentation", True, ["-p", "-P", "1"]),
    ]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fastp", required=True, type=Path)
    p.add_argument("--fasterp", required=True, type=Path)
    p.add_argument("--workdir", type=Path)
    p.add_argument("--workers", type=int, nargs="+", default=[1, 3],
                   help="fasterp worker counts to compare with the single-worker reference")
    args = p.parse_args()
    bins = [args.fastp.resolve(), args.fasterp.resolve()]
    version = subprocess.check_output([str(bins[0]), "--version"], text=True).strip()
    if version != "fastp 1.4.0":
        p.error(f"Expected fastp 1.4.0, got {version}")
    temporary = tempfile.TemporaryDirectory(prefix="fasterp-parity-")
    root = args.workdir.resolve() if args.workdir else Path(temporary.name)
    root.mkdir(parents=True, exist_ok=True)
    fixtures(root)
    failed = []
    for name, paired, options in cases():
        for workers in args.workers:
            folders = []
            reports = []
            for program, binary in enumerate(bins):
                folder = Path(tempfile.mkdtemp(prefix=f"{name}-{workers}-{program}-",dir=root))
                folders.append(folder)
                input1={"dedup_se":"duplicate1.fq","dedup_pe":"duplicate1.fq","large_se_auto":"large1.fq","large_pe_auto":"large1.fq","large_discovery":"custom.fq","gzip_members":"multi.fq.gz","gzip_magic":"gzip.data","phred64":"phred64.fq","empty":"empty.fq","interleaved":"interleaved.fq"}.get(name,"r1.fq")
                if name == "gzip_out_empty":
                    input1 = "empty.fq"
                elif name.startswith("gzip_out_"):
                    input1 = "large1.fq"
                elif name.startswith("stdin_interleaved"):
                    input1 = "interleaved.fq"
                elif name == "stdout_no_overlap":
                    input1 = "no_overlap1.fq"
                suffix = ".fq.gz" if name.startswith("gzip_out_") else ".fq"
                command = [str(binary), "-i", str(root / input1), "-o", "out1" + suffix, "-j", "report.json", "-h", "report.html", "--dont_eval_duplication", "-w", str(1 if program == 0 else workers)]
                if paired:
                    command += ["-O","out2" + suffix]
                    if "--interleaved_in" not in options:
                        command += ["-I",str(root / ("large2.fq" if name=="large_pe_auto" or name.startswith("gzip_out_") else "duplicate2.fq" if name=="dedup_pe" else "no_overlap2.fq" if name=="stdout_no_overlap" else "r2.fq"))]
                command += [str(root / opt[1:-1]) if opt.startswith("{") else opt for opt in options]
                with open(root / input1, "rb") as source:
                    result = subprocess.run(command, cwd=folder, capture_output=True, timeout=60,
                                            stdin=source if "--stdin" in options else subprocess.DEVNULL)
                (folder / "stderr.log").write_bytes(result.stderr)
                if result.returncode:
                    failed.append((name, workers, f"{binary.name} failed: {result.stderr.decode()[-500:]}"))
                    break
                reports.append(json.loads((folder / "report.json").read_text()))
                (folder / "stdout.fq").write_bytes(result.stdout)
            else:
                differences = []
                for path in sorted(set(f.name for d in folders for f in d.glob("*.fq*"))):
                    opener = gzip.open if path.endswith(".gz") else open
                    content = []
                    for folder in folders:
                        if (folder / path).exists():
                            with opener(folder / path, "rb") as stream:
                                content.append(stream.read())
                        else:
                            content.append(b"")
                    if content[0] != content[1]:
                        differences.append(f"{path} differs ({len(content[0])} vs {len(content[1])} bytes)")
                for key in ("filtering_result", "adapter_cutting"):
                    for field, expected in reports[0].get(key, {}).items():
                        if isinstance(expected, (int, float)) and expected != reports[1].get(key, {}).get(field):
                            differences.append(f"{key}.{field}: {expected} vs {reports[1].get(key, {}).get(field)}")
                for stage in ("before_filtering", "after_filtering"):
                    for field in ("total_reads", "total_bases", "q20_bases", "q30_bases", "read1_mean_length", "read2_mean_length"):
                        expected = reports[0]["summary"][stage].get(field)
                        actual = reports[1]["summary"][stage].get(field)
                        if expected != actual:
                            differences.append(f"summary.{stage}.{field}: {expected} vs {actual}")
                for section in ("read1_before_filtering","read1_after_filtering","read2_before_filtering","read2_after_filtering","merged_and_filtered"):
                    if section not in reports[0]:
                        continue
                    for field in ("total_cycles","q40_bases","kmer_count","overrepresented_sequences"):
                        if reports[0][section].get(field)!=reports[1].get(section,{}).get(field):
                            differences.append(f"{section}.{field} differs")
                    for curves in ("quality_curves","content_curves"):
                        for base,expected in reports[0][section].get(curves,{}).items():
                            actual=reports[1].get(section,{}).get(curves,{}).get(base,[])
                            if len(expected)!=len(actual) or any(not math.isclose(x,y,rel_tol=1e-5,abs_tol=1e-6) for x,y in zip(expected,actual)):
                                differences.append(f"{section}.{curves}.{base} differs")
                if differences:
                    failed.append((name, workers, "; ".join(differences)))
    for name, workers, reason in failed:
        print(f"FAIL {name} (threads={workers}): {reason}")
    total = len(cases()) * len(args.workers)
    print(f"{total-len(failed)}/{total} comparisons passed against {version}")
    if args.workdir:
        print(f"Results: {root}")
    raise SystemExit(bool(failed))


if __name__ == "__main__":
    main()
