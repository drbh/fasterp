//! Differential tests against fastp 1.4.0.
//!
//! Each case runs `fastp` (on PATH, must report `fastp 1.4.0`) once with one
//! worker and fasterp at 1, 4 and 8 workers on the same synthetic input, then
//! compares every FASTQ output (decompressed), the filtering and adapter
//! counters, the summary and the per-read quality and content curves.
//!
//! Run one case with `cargo test --release --test fastp_parity merge`.
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde_json::Value;
use std::{
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::OnceLock,
};

const WORKERS: [usize; 3] = [1, 4, 8];
const ADAPTER_R1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const ADAPTER_R2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

struct Case {
    name: &'static str,
    /// R1 input, or the interleaved input.
    input1: &'static str,
    /// R2 input. `None` with `paired` means interleaved input.
    input2: Option<&'static str>,
    paired: bool,
    gzip_out: bool,
    args: &'static [&'static str],
}

impl Case {
    const fn se(name: &'static str, args: &'static [&'static str]) -> Self {
        Case {
            name,
            input1: "r1.fq",
            input2: None,
            paired: false,
            gzip_out: false,
            args,
        }
    }
    const fn pe(name: &'static str, args: &'static [&'static str]) -> Self {
        Case {
            name,
            input1: "r1.fq",
            input2: Some("r2.fq"),
            paired: true,
            gzip_out: false,
            args,
        }
    }
    const fn inputs(mut self, input1: &'static str, input2: Option<&'static str>) -> Self {
        self.input1 = input1;
        self.input2 = input2;
        self
    }
    const fn interleaved(mut self, input: &'static str) -> Self {
        self.input1 = input;
        self.input2 = None;
        self.paired = true;
        self
    }
    const fn gzip_out(mut self) -> Self {
        self.gzip_out = true;
        self
    }
}

/// Synthetic inputs shared by every case, regenerated once per test run.
fn fixtures() -> &'static Path {
    static ROOT: OnceLock<PathBuf> = OnceLock::new();
    ROOT.get_or_init(|| {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("fastp-parity");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        write_fixtures(&dir);
        dir
    })
}

fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn random_dna(rng: &mut StdRng, len: usize) -> String {
    (0..len)
        .map(|_| b"ACGT"[rng.gen_range(0..4)] as char)
        .collect()
}

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn write_fixtures(root: &Path) {
    let mut rng = StdRng::seed_from_u64(140);
    let mut mates = [String::new(), String::new()];
    // 400 pairs with 15-250bp inserts in 150bp reads, so most carry adapter.
    // Every 7th has low-quality ends, every 11th R2 is all low quality and
    // every 13th has an N run.
    for i in 0..400 {
        let len = rng.gen_range(15..250);
        let insert = random_dna(&mut rng, len);
        for (mate, (seq, adapter)) in [
            (insert.clone(), ADAPTER_R1),
            (reverse_complement(&insert), ADAPTER_R2),
        ]
        .into_iter()
        .enumerate()
        {
            let mut seq: String = (seq + &adapter.repeat(5)).chars().take(150).collect();
            let mut qual = "I".repeat(seq.len());
            if i % 7 == 0 {
                qual = format!(
                    "{}{}{}",
                    "!".repeat(8),
                    &qual[8..qual.len() - 12],
                    "!".repeat(12)
                );
            }
            if i % 11 == 0 && mate == 1 {
                qual = "!".repeat(seq.len());
            }
            if i % 13 == 0 {
                seq.replace_range(20..28, "NNNNNNNN");
            }
            mates[mate] += &format!("@SEQ:{i} {}:N:0:ACGT+TGCA\n{seq}\n+\n{qual}\n", mate + 1);
        }
    }
    let [r1, r2] = &mates;
    let write = |name: &str, data: &[u8]| fs::write(root.join(name), data).unwrap();
    write("r1.fq", r1.as_bytes());
    write("r2.fq", r2.as_bytes());
    write("duplicate1.fq", r1.repeat(2).as_bytes());
    write("duplicate2.fq", r2.repeat(2).as_bytes());
    write("large1.fq", r1.repeat(25).as_bytes());
    write("large2.fq", r2.repeat(25).as_bytes());
    let interleaved: String = r1
        .split_inclusive('\n')
        .collect::<Vec<_>>()
        .chunks(4)
        .zip(r2.split_inclusive('\n').collect::<Vec<_>>().chunks(4))
        .map(|(a, b)| a.concat() + &b.concat())
        .collect();
    write("interleaved.fq", interleaved.as_bytes());
    // Three gzip members (one empty), split mid-record.
    let half = r1.len() / 2;
    let mut multi = gzip(&r1.as_bytes()[..half]);
    multi.extend(gzip(b""));
    multi.extend(gzip(&r1.as_bytes()[half..]));
    write("multi.fq.gz", &multi);
    let phred64: String = r1
        .lines()
        .enumerate()
        .map(|(i, line)| {
            let line = if i % 4 == 3 {
                line.bytes().map(|q| (q + 31) as char).collect()
            } else {
                line.to_string()
            };
            line + "\n"
        })
        .collect();
    write("phred64.fq", phred64.as_bytes());
    write(
        "adapters.fa",
        format!(">r1\n{ADAPTER_R1}\n>r2\n{ADAPTER_R2}\n").as_bytes(),
    );
    write("indexes.txt", b"ACGT\n");
    write("empty.fq", b"");
    for mate in 1..=2 {
        write(
            &format!("no_overlap{mate}.fq"),
            format!("@pair/{mate}\n{}\n+\n{}\n", "A".repeat(80), "I".repeat(80)).as_bytes(),
        );
    }
    // An adapter fastp does not know, for auto detection.
    let custom: String = (0..10000)
        .map(|i| {
            let len = rng.gen_range(25..120);
            let seq: String = (random_dna(&mut rng, len)
                + &"TCGATGCTACGATTCGAGTCATGCATGTAG".repeat(6))
                .chars()
                .take(150)
                .collect();
            format!("@CUSTOM:{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
        })
        .collect();
    write("custom.fq", custom.as_bytes());
}

/// Run one program in a fresh directory and return that directory.
fn run(program: &Path, case: &Case, workers: usize) -> tempfile::TempDir {
    let root = fixtures();
    let dir = tempfile::tempdir_in(root).unwrap();
    let suffix = if case.gzip_out { ".fq.gz" } else { ".fq" };
    let mut cmd = Command::new(program);
    cmd.current_dir(dir.path())
        .arg("-i")
        .arg(root.join(case.input1))
        .args([
            "-o",
            &format!("out1{suffix}"),
            "-j",
            "report.json",
            "-h",
            "report.html",
        ])
        .args(["--dont_eval_duplication", "-w", &workers.to_string()]);
    if case.paired {
        cmd.args(["-O", &format!("out2{suffix}")]);
    }
    if let Some(input2) = case.input2 {
        cmd.arg("-I").arg(root.join(input2));
    }
    for arg in case.args {
        // `@file` names a fixture.
        match arg.strip_prefix('@') {
            Some(file) => cmd.arg(root.join(file)),
            None => cmd.arg(arg),
        };
    }
    let stdin = if case.args.contains(&"--stdin") {
        Stdio::from(fs::File::open(root.join(case.input1)).unwrap())
    } else {
        Stdio::null()
    };
    let output = cmd
        .stdin(stdin)
        .output()
        .unwrap_or_else(|e| panic!("cannot run {}: {e}", program.display()));
    assert!(
        output.status.success(),
        "{} failed on {} (-w {workers}): {}",
        program.display(),
        case.name,
        String::from_utf8_lossy(&output.stderr)
    );
    fs::write(dir.path().join("stdout.fq"), &output.stdout).unwrap();
    dir
}

/// Decompressed contents. A missing or 0-byte file reads as empty.
fn read_fastq(path: &Path) -> Vec<u8> {
    let raw = fs::read(path).unwrap_or_default();
    if raw.is_empty() || path.extension().is_none_or(|e| e != "gz") {
        return raw;
    }
    let mut data = Vec::new();
    MultiGzDecoder::new(&raw[..])
        .read_to_end(&mut data)
        .unwrap();
    data
}

fn fastq_outputs(dirs: [&Path; 2]) -> Vec<String> {
    let mut names: Vec<String> = dirs
        .iter()
        .flat_map(|d| fs::read_dir(d).unwrap())
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .filter(|n| n.contains(".fq"))
        .collect();
    names.sort();
    names.dedup();
    names
}

fn close(a: &Value, b: &Value) -> bool {
    match (a.as_f64(), b.as_f64()) {
        (Some(x), Some(y)) => (x - y).abs() <= 1e-6 + 1e-5 * y.abs(),
        _ => a == b,
    }
}

fn differences(expected: &Path, actual: &Path) -> Vec<String> {
    let mut diffs = Vec::new();
    for name in fastq_outputs([expected, actual]) {
        let (a, b) = (
            read_fastq(&expected.join(&name)),
            read_fastq(&actual.join(&name)),
        );
        if a != b {
            diffs.push(format!("{name} differs ({} vs {} bytes)", a.len(), b.len()));
        }
    }
    let report = |dir: &Path| -> Value {
        serde_json::from_slice(&fs::read(dir.join("report.json")).unwrap()).unwrap()
    };
    let (a, b) = (report(expected), report(actual));
    for key in ["filtering_result", "adapter_cutting"] {
        for (field, value) in a[key].as_object().into_iter().flatten() {
            if value.is_number() && !close(value, &b[key][field]) {
                diffs.push(format!("{key}.{field}: {value} vs {}", b[key][field]));
            }
        }
    }
    for stage in ["before_filtering", "after_filtering"] {
        for field in [
            "total_reads",
            "total_bases",
            "q20_bases",
            "q30_bases",
            "read1_mean_length",
            "read2_mean_length",
        ] {
            let (x, y) = (&a["summary"][stage][field], &b["summary"][stage][field]);
            if x != y {
                diffs.push(format!("summary.{stage}.{field}: {x} vs {y}"));
            }
        }
    }
    for section in [
        "read1_before_filtering",
        "read1_after_filtering",
        "read2_before_filtering",
        "read2_after_filtering",
        "merged_and_filtered",
    ] {
        if a.get(section).is_none() {
            continue;
        }
        for field in [
            "total_cycles",
            "q40_bases",
            "kmer_count",
            "overrepresented_sequences",
        ] {
            if a[section].get(field) != b[section].get(field) {
                diffs.push(format!("{section}.{field} differs"));
            }
        }
        for curves in ["quality_curves", "content_curves"] {
            for (base, x) in a[section][curves].as_object().into_iter().flatten() {
                let (x, y) = (x.as_array().unwrap(), b[section][curves][base].as_array());
                if y.is_none_or(|y| {
                    y.len() != x.len() || x.iter().zip(y).any(|(x, y)| !close(x, y))
                }) {
                    diffs.push(format!("{section}.{curves}.{base} differs"));
                }
            }
        }
    }
    diffs
}

fn fastp() -> &'static Path {
    static FASTP: OnceLock<PathBuf> = OnceLock::new();
    FASTP.get_or_init(|| {
        let path = PathBuf::from("fastp");
        let version = Command::new(&path)
            .arg("--version")
            .output()
            .expect("fastp 1.4.0 must be on PATH");
        let version = String::from_utf8_lossy(&[version.stdout, version.stderr].concat())
            .trim()
            .to_string();
        assert_eq!(version, "fastp 1.4.0", "expected fastp 1.4.0 on PATH");
        path
    })
}

fn check(case: Case) {
    let reference = run(fastp(), &case, 1);
    let fasterp = PathBuf::from(assert_cmd::cargo::cargo_bin!("fasterp"));
    let failures: Vec<String> = WORKERS
        .iter()
        .flat_map(|&w| {
            let out = run(&fasterp, &case, w);
            differences(reference.path(), out.path())
                .into_iter()
                .map(move |d| format!("-w {w}: {d}"))
        })
        .collect();
    assert!(
        failures.is_empty(),
        "{} differs from fastp 1.4.0:\n{}",
        case.name,
        failures.join("\n")
    );
}

macro_rules! cases {
    ($($name:ident: $case:expr;)*) => {
        $(#[test] fn $name() { check($case) })*
    };
}

const LARGE: (&str, Option<&str>) = ("large1.fq", Some("large2.fq"));
const MERGE: &str = "merged.fq";

cases! {
    gzip_out_1: Case::pe("gzip_out_1", &["-z", "1"]).inputs(LARGE.0, LARGE.1).gzip_out();
    gzip_out_4: Case::pe("gzip_out_4", &["-z", "4"]).inputs(LARGE.0, LARGE.1).gzip_out();
    gzip_out_9: Case::pe("gzip_out_9", &["-z", "9"]).inputs(LARGE.0, LARGE.1).gzip_out();
    gzip_out_empty: Case::se("gzip_out_empty", &[]).inputs("empty.fq", None).gzip_out();
    dedup_se: Case::se("dedup_se", &["-D", "--dup_calc_accuracy", "1"]).inputs("duplicate1.fq", None);
    dedup_pe: Case::pe("dedup_pe", &["-D", "--dup_calc_accuracy", "1"]).inputs("duplicate1.fq", Some("duplicate2.fq"));
    large_se_auto: Case::se("large_se_auto", &[]).inputs("large1.fq", None);
    large_pe_auto: Case::pe("large_pe_auto", &["-2"]).inputs(LARGE.0, LARGE.1);
    large_discovery: Case::se("large_discovery", &[]).inputs("custom.fq", None);
    gzip_members: Case::se("gzip_members", &[]).inputs("multi.fq.gz", None);
    phred64: Case::se("phred64", &["-6"]).inputs("phred64.fq", None);
    empty: Case::se("empty", &[]).inputs("empty.fq", None);
    interleaved: Case::pe("interleaved", &["--interleaved_in"]).interleaved("interleaved.fq");
    stdout_pe: Case::pe("stdout_pe", &["--stdout"]);
    stdout_se: Case::se("stdout_se", &["--stdout"]);
    stdout_merge: Case::pe("stdout_merge", &["--stdout", "-m"]);
    stdout_unmerged: Case::pe("stdout_unmerged", &["--stdout", "-m", "--include_unmerged"]);
    stdout_no_overlap: Case::pe("stdout_no_overlap", &["--stdout", "-m", "-A", "-G"]).inputs("no_overlap1.fq", Some("no_overlap2.fq"));
    stdin_se: Case::se("stdin_se", &["--stdin", "--stdout", "-A", "-G"]);
    stdin_interleaved: Case::pe("stdin_interleaved", &["--stdin", "--stdout", "--interleaved_in", "-A", "-G"]).interleaved("interleaved.fq");
    stdin_interleaved_trim: Case::pe("stdin_interleaved_trim", &["--stdin", "--stdout", "--interleaved_in", "-5", "-3", "-x", "-U", "--umi_loc", "per_read", "--umi_len", "4"]).interleaved("interleaved.fq");
    stdin_interleaved_merge: Case::pe("stdin_interleaved_merge", &["--stdin", "--stdout", "--interleaved_in", "-m", "--include_unmerged"]).interleaved("interleaved.fq");
    se_default: Case::se("se_default", &[]);
    pe_default: Case::pe("pe_default", &[]);
    se_raw: Case::se("se_raw", &["-A", "-G"]);
    pe_raw: Case::pe("pe_raw", &["-A", "-G"]);
    quality_off: Case::pe("quality_off", &["-Q", "-A", "-G"]);
    length_limit: Case::pe("length_limit", &["--length_limit", "100"]);
    length_off: Case::pe("length_off", &["-L", "--length_limit", "100"]);
    complexity: Case::pe("complexity", &["-y", "-Y", "60"]);
    average_quality: Case::pe("average_quality", &["-e", "35"]);
    qualified_quality: Case::pe("qualified_quality", &["-q", "35", "-u", "10", "-n", "2"]);
    fixed_trim: Case::pe("fixed_trim", &["-f", "8", "-t", "4"]);
    explicit_r2_zero: Case::pe("explicit_r2_zero", &["-f", "8", "-F", "0", "-t", "4", "-T", "0", "-b", "90", "-B", "0"]);
    cut_front: Case::pe("cut_front", &["-5"]);
    cut_tail: Case::pe("cut_tail", &["-3"]);
    cut_right: Case::pe("cut_right", &["-r"]);
    cut_both: Case::pe("cut_both", &["-5", "-3", "-W", "5", "-M", "25"]);
    cut_overrides: Case::pe("cut_overrides", &["-5", "-3", "--cut_front_window_size", "2", "--cut_front_mean_quality", "30", "--cut_tail_window_size", "6", "--cut_tail_mean_quality", "25"]);
    max_length: Case::pe("max_length", &["-b", "80"]);
    poly_g: Case::pe("poly_g", &["-g", "--poly_g_min_len", "4"]);
    poly_x: Case::pe("poly_x", &["-x", "--poly_x_min_len", "4"]);
    poly_g_and_x: Case::pe("poly_g_and_x", &["-g", "-x", "--poly_g_min_len", "4", "--poly_x_min_len", "6"]);
    manual_adapter: Case::pe("manual_adapter", &["-a", ADAPTER_R1, "--adapter_sequence_r2", ADAPTER_R2]);
    adapter_fasta: Case::pe("adapter_fasta", &["--adapter_fasta", "@adapters.fa"]);
    adapter_dimer: Case::pe("adapter_dimer", &["--adapter_fasta", "@adapters.fa", "--dimer_max_len", "50", "-L", "-Q"]);
    pe_detection: Case::pe("pe_detection", &["-2"]);
    gap_overlap: Case::pe("gap_overlap", &["--allow_gap_overlap_trimming"]);
    correction: Case::pe("correction", &["-c"]);
    overlap_thresholds: Case::pe("overlap_thresholds", &["--overlap_len_require", "15", "--overlap_diff_limit", "2", "--overlap_diff_percent_limit", "10"]);
    merge: Case::pe("merge", &["-m", "--merged_out", MERGE]);
    merge_after_trim: Case::pe("merge_after_trim", &["-m", "--merged_out", MERGE, "-f", "8", "-b", "100", "-3"]);
    merge_corrected: Case::pe("merge_corrected", &["-m", "--merged_out", MERGE, "-c"]);
    merge_poly_x: Case::pe("merge_poly_x", &["-m", "--merged_out", MERGE, "-x", "--poly_x_min_len", "4"]);
    merge_manual_adapter: Case::pe("merge_manual_adapter", &["-m", "--merged_out", MERGE, "-a", ADAPTER_R1, "--adapter_sequence_r2", ADAPTER_R2]);
    merge_gap_max_length: Case::pe("merge_gap_max_length", &["-m", "--merged_out", MERGE, "--allow_gap_overlap_trimming", "-b", "65", "-B", "80"]);
    include_unmerged: Case::pe("include_unmerged", &["-m", "--merged_out", MERGE, "--include_unmerged"]);
    failed_after_cut: Case::pe("failed_after_cut", &["--unpaired1", "unpaired.fq", "--failed_out", "failed.fq", "-5"]);
    failed_and_unpaired: Case::pe("failed_and_unpaired", &["--unpaired1", "unpaired.fq", "--failed_out", "failed.fq"]);
    failed_and_two_unpaired: Case::pe("failed_and_two_unpaired", &["--unpaired1", "unpaired1.fq", "--unpaired2", "unpaired2.fq", "--failed_out", "failed.fq"]);
    overlapped: Case::pe("overlapped", &["--overlapped_out", "overlapped.fq"]);
    read_limit: Case::pe("read_limit", &["--reads_to_process", "27"]);
    umi_read: Case::pe("umi_read", &["-U", "--umi_loc", "read1", "--umi_len", "6", "--umi_skip", "2"]);
    umi_per_read: Case::pe("umi_per_read", &["-U", "--umi_loc", "per_read", "--umi_len", "6", "--umi_prefix", "UMI", "--umi_delim", "_"]);
    umi_index: Case::pe("umi_index", &["-U", "--umi_loc", "index1"]);
    umi_per_index: Case::pe("umi_per_index", &["-U", "--umi_loc", "per_index"]);
    index_filter: Case::pe("index_filter", &["--filter_by_index1", "@indexes.txt"]);
    overrepresentation: Case::pe("overrepresentation", &["-p", "-P", "1"]);
}
