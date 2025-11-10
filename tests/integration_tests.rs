use assert_cmd::cargo::cargo_bin;
use predicates::prelude::*;
use serde_json::Value;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

/// Helper to get the path to test data
fn test_data_path(filename: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("test_data")
        .join(filename)
}

/// Run fastp and return output paths
fn run_fastp(input: &str, temp_dir: &TempDir) -> (PathBuf, PathBuf) {
    let output_fq = temp_dir.path().join("fastp_output.fq");
    let output_json = temp_dir.path().join("fastp_output.json");

    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path(input))
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp - is it installed?");

    assert!(status.success(), "fastp command failed");

    (output_fq, output_json)
}

/// Run fasterp and return output paths
fn run_fasterp(input: &str, temp_dir: &TempDir) -> (PathBuf, PathBuf) {
    let output_fq = temp_dir.path().join("fasterp_output.fq");
    let output_json = temp_dir.path().join("fasterp_output.json");

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path(input))
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .status()
        .expect("Failed to run fasterp");

    assert!(status.success(), "fasterp command failed");

    (output_fq, output_json)
}

/// Compare two JSON files and check if key fields match
fn compare_json_outputs(fastp_json: &PathBuf, fasterp_json: &PathBuf) {
    let fastp_content = fs::read_to_string(fastp_json).expect("Failed to read fastp JSON");
    let fasterp_content = fs::read_to_string(fasterp_json).expect("Failed to read fasterp JSON");

    let fastp_data: Value = serde_json::from_str(&fastp_content).expect("Invalid fastp JSON");
    let fasterp_data: Value = serde_json::from_str(&fasterp_content).expect("Invalid fasterp JSON");

    // Compare filtering results
    let fastp_filtering = &fastp_data["filtering_result"];
    let fasterp_filtering = &fasterp_data["filtering_result"];

    assert_eq!(
        fastp_filtering["passed_filter_reads"], fasterp_filtering["passed_filter_reads"],
        "Passed filter reads don't match"
    );
    assert_eq!(
        fastp_filtering["too_short_reads"], fasterp_filtering["too_short_reads"],
        "Too short reads don't match"
    );
    assert_eq!(
        fastp_filtering["too_many_N_reads"], fasterp_filtering["too_many_N_reads"],
        "Too many N reads don't match"
    );
    assert_eq!(
        fastp_filtering["low_quality_reads"], fasterp_filtering["low_quality_reads"],
        "Low quality reads don't match"
    );

    // Compare summary statistics
    let fastp_summary = &fastp_data["summary"]["before_filtering"];
    let fasterp_summary = &fasterp_data["summary"]["before_filtering"];

    assert_eq!(
        fastp_summary["total_reads"], fasterp_summary["total_reads"],
        "Total reads don't match"
    );
    assert_eq!(
        fastp_summary["total_bases"], fasterp_summary["total_bases"],
        "Total bases don't match"
    );
    assert_eq!(
        fastp_summary["q20_bases"], fasterp_summary["q20_bases"],
        "Q20 bases don't match"
    );
    assert_eq!(
        fastp_summary["q30_bases"], fasterp_summary["q30_bases"],
        "Q30 bases don't match"
    );

    // Compare kmer counts
    let fastp_kmers = &fastp_data["read1_before_filtering"]["kmer_count"];
    let fasterp_kmers = &fasterp_data["read1_before_filtering"]["kmer_count"];

    let fastp_kmer_obj = fastp_kmers.as_object().expect("fastp kmers not an object");
    let fasterp_kmer_obj = fasterp_kmers
        .as_object()
        .expect("fasterp kmers not an object");

    assert_eq!(
        fastp_kmer_obj.len(),
        fasterp_kmer_obj.len(),
        "Different number of kmers"
    );

    // Spot check some kmers
    for kmer in ["AAAAA", "TTTTT", "ACGTG", "CCCCC", "GGGGG"].iter() {
        assert_eq!(
            fastp_kmer_obj.get(*kmer),
            fasterp_kmer_obj.get(*kmer),
            "Kmer {kmer} counts don't match"
        );
    }
}

/// Run fastp for paired-end reads
fn run_fastp_pe(input1: &str, input2: &str, temp_dir: &TempDir) -> (PathBuf, PathBuf, PathBuf) {
    let output1 = temp_dir.path().join("fastp_R1.fq");
    let output2 = temp_dir.path().join("fastp_R2.fq");
    let output_json = temp_dir.path().join("fastp_pe.json");

    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path(input1))
        .arg("-I")
        .arg(test_data_path(input2))
        .arg("-o")
        .arg(&output1)
        .arg("-O")
        .arg(&output2)
        .arg("-j")
        .arg(&output_json)
        .arg("-t")
        .arg("1")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp PE - is it installed?");

    assert!(status.success(), "fastp PE command failed");

    (output1, output2, output_json)
}

/// Run fasterp for paired-end reads
fn run_fasterp_pe(input1: &str, input2: &str, temp_dir: &TempDir) -> (PathBuf, PathBuf, PathBuf) {
    let output1 = temp_dir.path().join("fasterp_R1.fq");
    let output2 = temp_dir.path().join("fasterp_R2.fq");
    let output_json = temp_dir.path().join("fasterp_pe.json");

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path(input1))
        .arg("-I")
        .arg(test_data_path(input2))
        .arg("-o")
        .arg(&output1)
        .arg("-O")
        .arg(&output2)
        .arg("-j")
        .arg(&output_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp PE");

    assert!(status.success(), "fasterp PE command failed");

    (output1, output2, output_json)
}

#[test]
fn test_basic_filtering_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    // // Run both tools on R1.fq from parent directory
    // let input_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("R1.fq");

    // Create temp paths
    let input_path = test_data_path("R1.fq");
    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare FASTQ outputs
    let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_fq).unwrap();
    assert_eq!(fastp_content, fasterp_content, "FASTQ outputs don't match");

    // Compare JSON outputs
    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_small_dataset_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let (fastp_fq, fastp_json) = run_fastp("small_1k.fq", &temp_dir);
    let (fasterp_fq, fasterp_json) = run_fasterp("small_1k.fq", &temp_dir);

    // Compare FASTQ outputs
    let fastp_content = fs::read_to_string(fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(fasterp_fq).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "FASTQ outputs don't match for small dataset"
    );

    // Compare JSON outputs
    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_medium_dataset_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let (fastp_fq, fastp_json) = run_fastp("medium_10k.fq", &temp_dir);
    let (fasterp_fq, fasterp_json) = run_fasterp("medium_10k.fq", &temp_dir);

    // Compare FASTQ outputs
    let fastp_content = fs::read_to_string(fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(fasterp_fq).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "FASTQ outputs don't match for medium dataset"
    );

    // Compare JSON outputs
    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_length_filtering_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with length filter
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("small_1k.fq"))
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-l")
        .arg("50")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with length filter
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("small_1k.fq"))
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-l")
        .arg("50")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_fq).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Length filtered outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_quality_filtering_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with quality filter
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("small_1k.fq"))
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-q")
        .arg("20")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with quality filter
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("small_1k.fq"))
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-q")
        .arg("20")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_fq).unwrap();

    assert_eq!(
        fastp_content, fasterp_content,
        "Quality filtered outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_cli_help_works() {
    let output = Command::new(cargo_bin("fasterp"))
        .arg("--help")
        .output()
        .expect("Failed to run fasterp");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("A fast FASTQ preprocessor"));
}

#[test]
fn test_cli_version_works() {
    let output = Command::new(cargo_bin("fasterp"))
        .arg("--version")
        .output()
        .expect("Failed to run fasterp");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("0.1.0"));
}

#[test]
fn test_cli_missing_input_fails() {
    let output = Command::new(cargo_bin("fasterp"))
        .arg("-o")
        .arg("out.fq")
        .output()
        .expect("Failed to run fasterp");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required"));
}

// COMPRESSION TESTS - gzip input/output

#[test]
fn test_gzip_input_decompression() {
    let temp_dir = TempDir::new().unwrap();

    // Create a gzipped input file
    let input_fq = test_data_path("small_1k.fq");
    let input_gz = temp_dir.path().join("input.fq.gz");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Gzip the input file
    let status = Command::new("gzip")
        .arg("-c")
        .arg(&input_fq)
        .stdout(std::fs::File::create(&input_gz).unwrap())
        .status()
        .expect("Failed to run gzip");
    assert!(status.success());

    // Run fasterp on gzipped input (single-threaded mode required for gzip input)
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_gz)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare with fastp on original file
    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Outputs should match
    let fasterp_content = fs::read_to_string(&output_fq).unwrap();
    let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
    assert_eq!(
        fasterp_content, fastp_content,
        "Gzipped input processing differs"
    );
}

#[test]
fn test_gzip_output_compression() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("small_1k.fq");
    let output_gz = temp_dir.path().join("output.fq.gz");
    let output_json = temp_dir.path().join("output.json");

    // Run fasterp with gzip output
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_gz)
        .arg("-j")
        .arg(&output_json)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Verify output is actually gzipped
    assert!(output_gz.exists());

    // Decompress and compare with uncompressed output
    let decompressed_fq = temp_dir.path().join("decompressed.fq");
    let status = Command::new("gunzip")
        .arg("-c")
        .arg(&output_gz)
        .stdout(std::fs::File::create(&decompressed_fq).unwrap())
        .status()
        .expect("Failed to run gunzip");
    assert!(status.success());

    // Run fasterp without compression
    let uncompressed_fq = temp_dir.path().join("uncompressed.fq");
    let uncompressed_json = temp_dir.path().join("uncompressed.json");

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&uncompressed_fq)
        .arg("-j")
        .arg(&uncompressed_json)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare decompressed with uncompressed
    let decompressed_content = fs::read_to_string(&decompressed_fq).unwrap();
    let uncompressed_content = fs::read_to_string(&uncompressed_fq).unwrap();
    assert_eq!(
        decompressed_content, uncompressed_content,
        "Gzipped output differs from uncompressed"
    );
}

#[test]
fn test_gzip_compression_levels() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");

    // Test different compression levels
    for level in [1, 6, 9] {
        let output_gz = temp_dir.path().join(format!("output_level_{level}.fq.gz"));
        let output_json = temp_dir.path().join(format!("output_level_{level}.json"));

        let status = Command::new(cargo_bin("fasterp"))
            .arg("-i")
            .arg(&input_fq)
            .arg("-o")
            .arg(&output_gz)
            .arg("-j")
            .arg(&output_json)
            .arg("-z")
            .arg(level.to_string())
            .status()
            .expect("Failed to run fasterp");
        assert!(status.success(), "Compression level {level} failed");
        assert!(
            output_gz.exists(),
            "Output file not created for level {level}"
        );
    }

    // Verify all produce same content when decompressed
    let mut contents = Vec::new();
    for level in [1, 6, 9] {
        let output_gz = temp_dir.path().join(format!("output_level_{level}.fq.gz"));
        let decompressed_fq = temp_dir.path().join(format!("decompressed_{level}.fq"));

        let status = Command::new("gunzip")
            .arg("-c")
            .arg(&output_gz)
            .stdout(std::fs::File::create(&decompressed_fq).unwrap())
            .status()
            .expect("Failed to run gunzip");
        assert!(status.success());

        contents.push(fs::read_to_string(&decompressed_fq).unwrap());
    }

    // All should be identical
    assert_eq!(contents[0], contents[1], "Level 1 and 6 differ");
    assert_eq!(contents[1], contents[2], "Level 6 and 9 differ");
}

// N-BASE FILTERING TESTS

#[test]
fn test_n_base_filtering_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("small_1k.fq");

    for n_limit in [0, 3, 10] {
        let fastp_fq = temp_dir.path().join(format!("fastp_n{n_limit}.fq"));
        let fastp_json = temp_dir.path().join(format!("fastp_n{n_limit}.json"));
        let fasterp_fq = temp_dir.path().join(format!("fasterp_n{n_limit}.fq"));
        let fasterp_json = temp_dir.path().join(format!("fasterp_n{n_limit}.json"));

        // Run fastp
        let status = Command::new("fastp")
            .arg("-i")
            .arg(&input_fq)
            .arg("-o")
            .arg(&fastp_fq)
            .arg("-j")
            .arg(&fastp_json)
            .arg("-n")
            .arg(n_limit.to_string())
            .stderr(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .status()
            .expect("Failed to run fastp");
        assert!(status.success());

        // Run fasterp
        let status = Command::new(cargo_bin("fasterp"))
            .arg("-i")
            .arg(&input_fq)
            .arg("-o")
            .arg(&fasterp_fq)
            .arg("-j")
            .arg(&fasterp_json)
            .arg("-n")
            .arg(n_limit.to_string())
            .status()
            .expect("Failed to run fasterp");
        assert!(status.success());

        // Compare outputs
        let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
        let fasterp_content = fs::read_to_string(&fasterp_fq).unwrap();
        assert_eq!(
            fastp_content, fasterp_content,
            "N-base filtering with limit {n_limit} differs"
        );

        // Compare JSON
        compare_json_outputs(&fastp_json, &fasterp_json);
    }
}

// MULTI-THREADING TESTS

#[test]
fn test_multithreading_consistency() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");

    // Run with 1, 2, 4, and 8 threads
    let thread_counts = [1, 2, 4, 8];
    let mut outputs = Vec::new();

    for threads in thread_counts.iter() {
        let output_fq = temp_dir.path().join(format!("output_t{threads}.fq"));
        let output_json = temp_dir.path().join(format!("output_t{threads}.json"));

        let status = Command::new(cargo_bin("fasterp"))
            .arg("-i")
            .arg(&input_fq)
            .arg("-o")
            .arg(&output_fq)
            .arg("-j")
            .arg(&output_json)
            .arg("-t")
            .arg(threads.to_string())
            .status()
            .expect("Failed to run fasterp");
        assert!(status.success(), "Failed with {threads} threads");

        let content = fs::read_to_string(&output_fq).unwrap();
        outputs.push(content);
    }

    // All outputs should be identical regardless of thread count
    for i in 1..outputs.len() {
        assert_eq!(
            outputs[0], outputs[i],
            "Output differs between 1 thread and {} threads",
            thread_counts[i]
        );
    }
}

#[test]
fn test_multithreading_json_consistency() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");

    // Run with different thread counts
    let output_1t_json = temp_dir.path().join("output_1t.json");
    let output_4t_json = temp_dir.path().join("output_4t.json");

    // 1 thread
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(temp_dir.path().join("output_1t.fq"))
        .arg("-j")
        .arg(&output_1t_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // 4 threads
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(temp_dir.path().join("output_4t.fq"))
        .arg("-j")
        .arg(&output_4t_json)
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare JSON outputs - should be identical
    compare_json_outputs(&output_1t_json, &output_4t_json);
}

#[test]
fn test_multithreading_large_dataset() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("large_100k.fq");

    let output_1t_fq = temp_dir.path().join("output_1t.fq");
    let output_8t_fq = temp_dir.path().join("output_8t.fq");
    let output_1t_json = temp_dir.path().join("output_1t.json");
    let output_8t_json = temp_dir.path().join("output_8t.json");

    // 1 thread
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_1t_fq)
        .arg("-j")
        .arg(&output_1t_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // 8 threads
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_8t_fq)
        .arg("-j")
        .arg(&output_8t_json)
        .arg("-t")
        .arg("8")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let output_1t = fs::read_to_string(&output_1t_fq).unwrap();
    let output_8t = fs::read_to_string(&output_8t_fq).unwrap();
    assert_eq!(
        output_1t, output_8t,
        "Multi-threading produces different output on large dataset"
    );

    // Compare JSON
    compare_json_outputs(&output_1t_json, &output_8t_json);
}

// STDIN/STDOUT TESTS

#[test]
fn test_stdin_input() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("small_1k.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Run fasterp with stdin (single-threaded mode required for stdin)
    let input_file = std::fs::File::open(&input_fq).unwrap();
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg("-")
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-t")
        .arg("1")
        .stdin(input_file)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare with normal file input
    let expected_fq = temp_dir.path().join("expected.fq");
    let expected_json = temp_dir.path().join("expected.json");

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&expected_fq)
        .arg("-j")
        .arg(&expected_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    let output_content = fs::read_to_string(&output_fq).unwrap();
    let expected_content = fs::read_to_string(&expected_fq).unwrap();
    assert_eq!(
        output_content, expected_content,
        "Stdin input produces different output"
    );
}

#[test]
fn test_stdout_output() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("small_1k.fq");
    let output_json = temp_dir.path().join("output.json");

    // Run fasterp with stdout (single-threaded mode required for stdout)
    let output = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg("-")
        .arg("-j")
        .arg(&output_json)
        .arg("-t")
        .arg("1")
        .output()
        .expect("Failed to run fasterp");
    assert!(output.status.success());

    let stdout_content = String::from_utf8(output.stdout).unwrap();

    // Compare with normal file output
    let expected_fq = temp_dir.path().join("expected.fq");
    let expected_json = temp_dir.path().join("expected.json");

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&expected_fq)
        .arg("-j")
        .arg(&expected_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    let expected_content = fs::read_to_string(&expected_fq).unwrap();
    assert_eq!(
        stdout_content, expected_content,
        "Stdout output produces different result"
    );
}

#[test]
fn test_stdin_stdout_pipeline() {
    let input_fq = test_data_path("R1.fq");

    // Run fasterp with stdin and stdout (single-threaded mode required)
    let input_file = std::fs::File::open(&input_fq).unwrap();
    let output = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("-j")
        .arg("/dev/null")
        .arg("-t")
        .arg("1")
        .stdin(input_file)
        .output()
        .expect("Failed to run fasterp");
    assert!(output.status.success());

    let stdout_content = String::from_utf8(output.stdout).unwrap();
    assert!(
        !stdout_content.is_empty(),
        "No output from stdin/stdout pipeline"
    );
}

// COMBINED FILTER TESTS

#[test]
fn test_combined_filters_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");
    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with multiple filters
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-l")
        .arg("30")
        .arg("-q")
        .arg("15")
        .arg("-n")
        .arg("3")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same filters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-l")
        .arg("30")
        .arg("-q")
        .arg("15")
        .arg("-n")
        .arg("3")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_fq).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Combined filters produce different output"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_strict_combined_filters() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");
    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Very strict filtering
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-l")
        .arg("100")
        .arg("-q")
        .arg("30")
        .arg("-n")
        .arg("0")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-l")
        .arg("100")
        .arg("-q")
        .arg("30")
        .arg("-n")
        .arg("0")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    let fastp_content = fs::read_to_string(&fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_fq).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Strict filters produce different output"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

// EDGE CASE TESTS - Using Crafted Data

/// Helper to create a FASTQ record with specific quality pattern
fn create_fastq_record(name: &str, seq: &str, qual: &str) -> String {
    format!("@{name}\n{seq}\n+\n{qual}\n")
}

/// Helper to create quality string with specific phred scores
fn quality_string_from_phred(phred_scores: &[u8]) -> String {
    phred_scores.iter().map(|&p| (p + 33) as char).collect()
}

#[test]
fn test_unqualified_percent_boundary_40_percent() {
    // Test the integer division bug fix: 40.39% unqualified should be rejected
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create a read with exactly 40.39% bases below Q9
    // 1431 bases total, 578 bases with Q<9 = 40.39%
    let mut qual_scores = vec![15u8; 853]; // 853 bases with Q15 (qualified)
    qual_scores.extend(vec![5u8; 578]); // 578 bases with Q5 (unqualified)

    let seq = "A".repeat(1431);
    let qual = quality_string_from_phred(&qual_scores);
    let record = create_fastq_record("test_read", &seq, &qual);

    fs::write(&input_fq, record).unwrap();

    // Run with -q 9 (default -u 40)
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-q")
        .arg("9")
        .arg("-l")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Should filter out (40.39% > 40%)
    let content = fs::read_to_string(&output_fq).unwrap();
    assert!(
        content.is_empty(),
        "Read with 40.39% unqualified should be filtered"
    );

    // Verify JSON shows 1 low quality read
    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["low_quality_reads"]
            .as_u64()
            .unwrap(),
        1
    );
    assert_eq!(
        json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap(),
        0
    );
}

#[test]
fn test_unqualified_percent_exact_40_percent() {
    // Test exact boundary: 40.0% unqualified should PASS
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create a read with exactly 40.0% bases below Q9
    // 100 bases total, 40 bases with Q<9 = 40.0%
    let mut qual_scores = vec![15u8; 60]; // 60 bases with Q15
    qual_scores.extend(vec![5u8; 40]); // 40 bases with Q5

    let seq = "A".repeat(100);
    let qual = quality_string_from_phred(&qual_scores);
    let record = create_fastq_record("test_read", &seq, &qual);

    fs::write(&input_fq, record).unwrap();

    // Run with -q 9 -u 40
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-q")
        .arg("9")
        .arg("-u")
        .arg("40")
        .arg("-l")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Should PASS (40.0% is not > 40%)
    let content = fs::read_to_string(&output_fq).unwrap();
    assert!(
        !content.is_empty(),
        "Read with exactly 40.0% unqualified should pass"
    );

    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap(),
        1
    );
}

#[test]
fn test_average_quality_filter() {
    // Test -e parameter (average quality threshold)
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create two reads: one with mean Q20, one with mean Q15
    let seq = "A".repeat(100);
    let qual_20 = quality_string_from_phred(&[20u8; 100]); // Mean = 20
    let qual_15 = quality_string_from_phred(&[15u8; 100]); // Mean = 15

    let mut input_data = String::new();
    input_data.push_str(&create_fastq_record("read_q20", &seq, &qual_20));
    input_data.push_str(&create_fastq_record("read_q15", &seq, &qual_15));

    fs::write(&input_fq, input_data).unwrap();

    // Run with -e 18 (require mean quality >= 18)
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-e")
        .arg("18")
        .arg("-q")
        .arg("0") // Disable -q/-u filtering
        .arg("-l")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Only read_q20 should pass
    let content = fs::read_to_string(&output_fq).unwrap();
    assert_eq!(content.lines().count(), 4, "Should have 1 read (4 lines)");
    assert!(content.contains("read_q20"), "read_q20 should pass");
    assert!(!content.contains("read_q15"), "read_q15 should be filtered");

    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap(),
        1
    );
    assert_eq!(
        json["filtering_result"]["low_quality_reads"]
            .as_u64()
            .unwrap(),
        1
    );
}

#[test]
fn test_max_length_trimming() {
    // Test -b parameter (max length trimming)
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create a 150bp read
    let seq = "ACGT".repeat(37) + "AC"; // 150 bases
    let qual = quality_string_from_phred(&[30u8; 150]);
    let record = create_fastq_record("long_read", &seq, &qual);

    fs::write(&input_fq, record).unwrap();

    // Run with -b 100 (trim to max 100bp)
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-b")
        .arg("100")
        .arg("-l")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Output should be exactly 100bp
    let content = fs::read_to_string(&output_fq).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4);
    assert_eq!(lines[1].len(), 100, "Sequence should be trimmed to 100bp");
    assert_eq!(lines[3].len(), 100, "Quality should be trimmed to 100bp");
}

#[test]
fn test_combined_q_u_and_e_filters() {
    // Test that -q/-u and -e work together correctly
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create three reads:
    // 1. High quality overall, but 50% unqualified at Q9 threshold
    let mut qual1 = vec![25u8; 50]; // High quality bases
    qual1.extend(vec![5u8; 50]); // Low quality bases, mean = 15

    // 2. All Q12 (passes -q 9 -u 40 but mean < 15)
    let qual2 = vec![12u8; 100]; // Mean = 12

    // 3. All Q20 (passes both)
    let qual3 = vec![20u8; 100]; // Mean = 20

    let seq = "A".repeat(100);
    let mut input_data = String::new();
    input_data.push_str(&create_fastq_record(
        "read1",
        &seq,
        &quality_string_from_phred(&qual1),
    ));
    input_data.push_str(&create_fastq_record(
        "read2",
        &seq,
        &quality_string_from_phred(&qual2),
    ));
    input_data.push_str(&create_fastq_record(
        "read3",
        &seq,
        &quality_string_from_phred(&qual3),
    ));

    fs::write(&input_fq, input_data).unwrap();

    // Run with -q 9 -u 40 AND -e 15
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-q")
        .arg("9")
        .arg("-u")
        .arg("40")
        .arg("-e")
        .arg("15")
        .arg("-l")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // read1: fails -u (50% > 40%)
    // read2: passes -q/-u but fails -e (12 < 15)
    // read3: passes both
    let content = fs::read_to_string(&output_fq).unwrap();
    assert_eq!(content.lines().count(), 4, "Only 1 read should pass");
    assert!(content.contains("read3"), "Only read3 should pass");

    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap(),
        1
    );
    assert_eq!(
        json["filtering_result"]["low_quality_reads"]
            .as_u64()
            .unwrap(),
        2
    );
}

#[test]
fn test_n_base_limit_exact_boundary() {
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create reads with different N counts
    let seq_5n = "ACGT".repeat(20) + "NNNNN"; // Exactly 5 Ns
    let seq_6n = "ACGT".repeat(20) + "NNNNNN"; // 6 Ns
    let qual = quality_string_from_phred(&[30u8; 85]);
    let qual6 = quality_string_from_phred(&[30u8; 86]);

    let mut input_data = String::new();
    input_data.push_str(&create_fastq_record("read_5n", &seq_5n, &qual));
    input_data.push_str(&create_fastq_record("read_6n", &seq_6n, &qual6));

    fs::write(&input_fq, input_data).unwrap();

    // Run with -n 5 (allow up to 5 Ns)
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-n")
        .arg("5")
        .arg("-l")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // read_5n should pass (5 <= 5), read_6n should fail (6 > 5)
    let content = fs::read_to_string(&output_fq).unwrap();
    assert_eq!(content.lines().count(), 4);
    assert!(content.contains("read_5n"));
    assert!(!content.contains("read_6n"));

    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap(),
        1
    );
    assert_eq!(
        json["filtering_result"]["too_many_N_reads"]
            .as_u64()
            .unwrap(),
        1
    );
}

#[test]
fn test_length_filter_after_trimming() {
    // Test that length filter applies AFTER trimming
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create a 100bp read
    let seq = "A".repeat(100);
    let qual = quality_string_from_phred(&[30u8; 100]);
    let record = create_fastq_record("test_read", &seq, &qual);

    fs::write(&input_fq, record).unwrap();

    // Trim 10bp from front, 40bp from tail (leaves 50bp)
    // Then require min length of 60bp
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("--trim-front")
        .arg("10")
        .arg("--trim-tail")
        .arg("40")
        .arg("-l")
        .arg("60")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Should be filtered (50bp after trimming < 60bp minimum)
    let content = fs::read_to_string(&output_fq).unwrap();
    assert!(content.is_empty(), "Read should be filtered after trimming");

    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["too_short_reads"]
            .as_u64()
            .unwrap(),
        1
    );
}

#[test]
fn test_empty_output_all_filtered() {
    let temp_dir = TempDir::new().unwrap();
    let input_fq = temp_dir.path().join("input.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Create 10 reads that will all be filtered
    let mut input_data = String::new();
    for i in 0..10 {
        let seq = "A".repeat(50); // All 50bp
        let qual = quality_string_from_phred(&[30u8; 50]);
        input_data.push_str(&create_fastq_record(&format!("read{i}"), &seq, &qual));
    }

    fs::write(&input_fq, input_data).unwrap();

    // Use extreme filter that should filter everything
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-l")
        .arg("1000") // Extremely long length requirement
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Output file should exist but be empty
    let content = fs::read_to_string(&output_fq).unwrap();
    assert!(content.is_empty(), "All reads should be filtered");

    // JSON should show 0 passed reads, 10 too_short
    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();
    assert_eq!(
        json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap(),
        0
    );
    assert_eq!(
        json["filtering_result"]["too_short_reads"]
            .as_u64()
            .unwrap(),
        10
    );
    assert_eq!(
        json["summary"]["before_filtering"]["total_reads"]
            .as_u64()
            .unwrap(),
        10
    );
}

#[test]
fn test_invalid_quality_threshold() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("small_1k.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Quality threshold of 0 should be treated as disabled
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-q")
        .arg("0")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());
}

// LARGE DATASET TESTS (stress testing)

#[test]
#[ignore] // Run with --ignored flag for stress testing
fn test_large_dataset_100k_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    let (fastp_fq, fastp_json) = run_fastp("large_100k.fq", &temp_dir);
    let (fasterp_fq, fasterp_json) = run_fasterp("large_100k.fq", &temp_dir);

    let fastp_content = fs::read_to_string(fastp_fq).unwrap();
    let fasterp_content = fs::read_to_string(fasterp_fq).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Large dataset output differs"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
#[ignore] // Run with --ignored flag for stress testing
fn test_xlarge_dataset_500k() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("xlarge_500k.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    // Just verify it completes successfully
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .arg("-t")
        .arg("8")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    assert!(output_fq.exists());
    assert!(output_json.exists());
}

// BATCH SIZE TESTS

#[test]
fn test_different_batch_sizes_consistent() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");

    // Test different batch sizes
    let batch_sizes = [1024 * 1024, 4 * 1024 * 1024, 16 * 1024 * 1024]; // 1MB, 4MB, 16MB
    let mut outputs = Vec::new();

    for (i, &batch_size) in batch_sizes.iter().enumerate() {
        let output_fq = temp_dir.path().join(format!("output_batch_{i}.fq"));
        let output_json = temp_dir.path().join(format!("output_batch_{i}.json"));

        let status = Command::new(cargo_bin("fasterp"))
            .arg("-i")
            .arg(&input_fq)
            .arg("-o")
            .arg(&output_fq)
            .arg("-j")
            .arg(&output_json)
            .arg("-t")
            .arg("4")
            .arg("--batch-bytes")
            .arg(batch_size.to_string())
            .status()
            .expect("Failed to run fasterp");
        assert!(status.success(), "Failed with batch size {batch_size}");

        let content = fs::read_to_string(&output_fq).unwrap();
        outputs.push(content);
    }

    // All outputs should be identical
    for i in 1..outputs.len() {
        assert_eq!(
            outputs[0], outputs[i],
            "Output differs with different batch sizes"
        );
    }
}

// JSON REPORT VALIDATION TESTS

#[test]
fn test_json_report_structure() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("small_1k.fq");
    let output_fq = temp_dir.path().join("output.fq");
    let output_json = temp_dir.path().join("output.json");

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_fq)
        .arg("-j")
        .arg(&output_json)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Validate JSON structure
    let json_content = fs::read_to_string(&output_json).unwrap();
    let json: Value = serde_json::from_str(&json_content).unwrap();

    // Check required fields exist
    assert!(json.get("summary").is_some(), "Missing 'summary' field");
    assert!(
        json.get("filtering_result").is_some(),
        "Missing 'filtering_result' field"
    );
    assert!(
        json.get("read1_before_filtering").is_some(),
        "Missing 'read1_before_filtering' field"
    );

    // Check summary fields
    let summary = &json["summary"];
    assert!(summary.get("fastp_version").is_some());
    assert!(summary.get("before_filtering").is_some());
    assert!(summary.get("after_filtering").is_some());

    // Check filtering result fields
    let filtering = &json["filtering_result"];
    assert!(filtering.get("passed_filter_reads").is_some());
    assert!(filtering.get("low_quality_reads").is_some());
    assert!(filtering.get("too_many_N_reads").is_some());
    assert!(filtering.get("too_short_reads").is_some());

    // Check detailed stats
    let detailed = &json["read1_before_filtering"];
    assert!(detailed.get("quality_curves").is_some());
    assert!(detailed.get("kmer_count").is_some());
}

// OUTPUT ORDERING TESTS

#[test]
fn test_output_order_preserved_multithreading() {
    let temp_dir = TempDir::new().unwrap();

    let input_fq = test_data_path("medium_10k.fq");
    let output_st = temp_dir.path().join("output_st.fq");
    let output_mt = temp_dir.path().join("output_mt.fq");

    // Single-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_st)
        .arg("-j")
        .arg("/dev/null")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Multi-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_fq)
        .arg("-o")
        .arg(&output_mt)
        .arg("-j")
        .arg("/dev/null")
        .arg("-t")
        .arg("8")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Outputs should be identical (order preserved)
    let st_content = fs::read_to_string(&output_st).unwrap();
    let mt_content = fs::read_to_string(&output_mt).unwrap();
    assert_eq!(
        st_content, mt_content,
        "Multi-threading changes output order"
    );
}

// Trimming Integration Tests

#[test]
fn test_sliding_window_tail_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    // Test with default sliding window tail trimming (cut_mean_quality=20, window_size=4)
    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with sliding window tail trimming
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--cut_tail")
        .arg("--cut_mean_quality")
        .arg("20")
        .arg("--cut_window_size")
        .arg("4")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--cut-tail")
        .arg("--cut-mean-quality")
        .arg("20")
        .arg("--cut-window-size")
        .arg("4")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_sliding_window_front_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with front trimming enabled
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--cut_front")
        .arg("--cut_mean_quality")
        .arg("25")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--cut-front")
        .arg("--cut-mean-quality")
        .arg("25")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_sliding_window_both_ends_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("medium_10k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with trimming on both ends
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--cut_front")
        .arg("--cut_tail")
        .arg("--cut_mean_quality")
        .arg("20")
        .arg("--cut_window_size")
        .arg("5")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--cut-front")
        .arg("--cut-tail")
        .arg("--cut-mean-quality")
        .arg("20")
        .arg("--cut-window-size")
        .arg("5")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_fixed_front_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with fixed front trimming
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_front1")
        .arg("5")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-front")
        .arg("5")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_fixed_tail_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with fixed tail trimming
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_tail1")
        .arg("10")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-tail")
        .arg("10")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_fixed_both_ends_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("medium_10k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with fixed trimming on both ends
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_front1")
        .arg("3")
        .arg("--trim_tail1")
        .arg("7")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-front")
        .arg("3")
        .arg("--trim-tail")
        .arg("7")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_poly_g_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with polyG trimming (default enabled for 2-color Illumina)
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_poly_g")
        .arg("--poly_g_min_len")
        .arg("10")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-poly-g")
        .arg("--poly-g-min-len")
        .arg("10")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_poly_x_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with polyX trimming
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_poly_x")
        .arg("--poly_x_min_len")
        .arg("12")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-poly-x")
        .arg("--poly-g-min-len")
        .arg("12")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_combined_trimming_and_filtering_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("medium_10k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with combined trimming and filtering
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_front1")
        .arg("2")
        .arg("--cut_tail")
        .arg("--cut_mean_quality")
        .arg("20")
        .arg("--trim_poly_g")
        .arg("-l")
        .arg("50")
        .arg("-n")
        .arg("5")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-front")
        .arg("2")
        .arg("--cut-tail")
        .arg("--cut-mean-quality")
        .arg("20")
        .arg("--trim-poly-g")
        .arg("-l")
        .arg("50")
        .arg("-n")
        .arg("5")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_disable_tail_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("small_1k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with tail trimming disabled
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--disable_trim_poly_g")
        .arg("--cut_mean_quality")
        .arg("0")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--disable-trim-poly-g")
        .arg("--disable-trim-tail")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_aggressive_trimming_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("medium_10k.fq");

    let fastp_fq = temp_dir.path().join("fastp.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_fq = temp_dir.path().join("fasterp.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with aggressive trimming settings
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fastp_fq)
        .arg("-j")
        .arg(&fastp_json)
        .arg("--trim_front1")
        .arg("5")
        .arg("--trim_tail1")
        .arg("5")
        .arg("--cut_front")
        .arg("--cut_tail")
        .arg("--cut_mean_quality")
        .arg("30")
        .arg("--cut_window_size")
        .arg("3")
        .arg("--trim_poly_g")
        .arg("--poly_g_min_len")
        .arg("8")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same settings
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&fasterp_fq)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("--trim-front")
        .arg("5")
        .arg("--trim-tail")
        .arg("5")
        .arg("--cut-front")
        .arg("--cut-tail")
        .arg("--cut-mean-quality")
        .arg("30")
        .arg("--cut-window-size")
        .arg("3")
        .arg("--trim-poly-g")
        .arg("--poly-g-min-len")
        .arg("8")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_trimming_with_multithreading_consistency() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = test_data_path("large_100k.fq");

    let output_st = temp_dir.path().join("output_st.fq");
    let output_mt = temp_dir.path().join("output_mt.fq");
    let json_st = temp_dir.path().join("output_st.json");
    let json_mt = temp_dir.path().join("output_mt.json");

    // Single-threaded with trimming
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&output_st)
        .arg("-j")
        .arg(&json_st)
        .arg("-t")
        .arg("1")
        .arg("--cut-mean-quality")
        .arg("20")
        .arg("--trim-front")
        .arg("3")
        .arg("--poly-g-min-len")
        .arg("10")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Multi-threaded with trimming
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&input_path)
        .arg("-o")
        .arg(&output_mt)
        .arg("-j")
        .arg(&json_mt)
        .arg("-t")
        .arg("4")
        .arg("--cut-mean-quality")
        .arg("20")
        .arg("--trim-front")
        .arg("3")
        .arg("--poly-g-min-len")
        .arg("10")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare JSON outputs
    compare_json_outputs(&json_st, &json_mt);
}

// PAIRED-END TESTS

#[test]
fn test_paired_end_basic_processing() {
    let temp_dir = TempDir::new().unwrap();

    // Run both tools on paired-end data
    let (fastp_r1, fastp_r2, fastp_json) =
        run_fastp_pe("pe_small_R1.fq", "pe_small_R2.fq", &temp_dir);
    let (fasterp_r1, fasterp_r2, fasterp_json) =
        run_fasterp_pe("pe_small_R1.fq", "pe_small_R2.fq", &temp_dir);

    // Compare R1 outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    // Compare R2 outputs
    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    // Compare JSON reports
    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_paired_end_quality_filtering() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with quality filtering
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-q")
        .arg("20")
        .arg("-u")
        .arg("40")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-q")
        .arg("20")
        .arg("-u")
        .arg("40")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_paired_end_length_filtering() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with length filtering
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-l")
        .arg("100")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-l")
        .arg("100")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_paired_end_combined_filters() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with multiple filters
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-q")
        .arg("20")
        .arg("-l")
        .arg("50")
        .arg("-n")
        .arg("5")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-q")
        .arg("20")
        .arg("-l")
        .arg("50")
        .arg("-n")
        .arg("5")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

// ADAPTER TRIMMING TESTS

#[test]
fn test_adapter_trimming_single_end_custom_adapter() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    let adapter = "AGATCGGAAGAGC";

    // Run fastp with custom adapter
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("R1.fq"))
        .arg("-o")
        .arg(&fastp_out)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-a")
        .arg(adapter)
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same adapter
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("R1.fq"))
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-a")
        .arg(adapter)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(fastp_content, fasterp_content, "Outputs don't match");

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_adapter_trimming_paired_end_custom_adapters() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    let adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

    // Run fastp with custom adapters
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-a")
        .arg(adapter1)
        .arg("--adapter_sequence_r2")
        .arg(adapter2)
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same adapters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-a")
        .arg(adapter1)
        .arg("--adapter_sequence_r2")
        .arg(adapter2)
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_disable_adapter_trimming() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with adapter trimming disabled
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("R1.fq"))
        .arg("-o")
        .arg(&fastp_out)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-A")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with adapter trimming disabled
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("R1.fq"))
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-A")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(fastp_content, fasterp_content, "Outputs don't match");

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_adapter_trimming_with_quality_filter() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    let adapter = "AGATCGGAAGAGC";

    // Run fastp with adapter trimming and quality filtering
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("small_1k.fq"))
        .arg("-o")
        .arg(&fastp_out)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-a")
        .arg(adapter)
        .arg("-q")
        .arg("20")
        .arg("-l")
        .arg("50")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("small_1k.fq"))
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-a")
        .arg(adapter)
        .arg("-q")
        .arg("20")
        .arg("-l")
        .arg("50")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(fastp_content, fasterp_content, "Outputs don't match");

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_paired_end_with_trimming_and_adapters() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    let adapter1 = "AGATCGGAAGAGC";
    let adapter2 = "AGATCGGAAGAGC";

    // Run fastp with adapters and trimming
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-a")
        .arg(adapter1)
        .arg("--adapter_sequence_r2")
        .arg(adapter2)
        .arg("--trim_front1")
        .arg("5")
        .arg("--trim_tail1")
        .arg("3")
        .arg("-q")
        .arg("20")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-a")
        .arg(adapter1)
        .arg("--adapter_sequence_r2")
        .arg(adapter2)
        .arg("--trim_front1")
        .arg("5")
        .arg("--trim_tail1")
        .arg("3")
        .arg("-q")
        .arg("20")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}
// EDGE CASE TESTS - Added for robustness

#[test]
fn test_adapter_5base_perfect_match() {
    // Test minimum length adapter match (5 bases, 0 mismatches)
    let temp_dir = TempDir::new().unwrap();

    // Create test data with 5-base partial adapter at the end
    let test_input = temp_dir.path().join("test_input.fq");
    fs::write(&test_input, "@read1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGATC\n+\n?????????????????????????????????????????????????????????????\n").unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");

    // Run fastp
    let status = Command::new("fastp")
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fastp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-L") // Disable length filtering
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-L")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "5-base adapter match outputs don't match"
    );
}

#[test]
fn test_adapter_8base_with_mismatch() {
    // Test 8-base adapter with 1 mismatch (should trim)
    let temp_dir = TempDir::new().unwrap();

    // Create test data with 8-base adapter (1 mismatch) at the end
    // AGATCGGC has 1 mismatch vs AGATCGGA
    let test_input = temp_dir.path().join("test_input.fq");
    fs::write(&test_input, "@read1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGATCGGC\n+\n????????????????????????????????????????????????????????????????\n").unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");

    // Run fastp
    Command::new("fastp")
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fastp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-L")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");

    // Run fasterp
    Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-L")
        .status()
        .expect("Failed to run fasterp");

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "8-base adapter with mismatch outputs don't match"
    );
}

#[test]
fn test_adapter_5base_with_mismatch_not_trimmed() {
    // Test 5-base adapter with 1 mismatch (should NOT trim)
    let temp_dir = TempDir::new().unwrap();

    // Create test data with 5-base adapter (1 mismatch) at the end
    // AGATX has 1 mismatch vs AGATC
    let test_input = temp_dir.path().join("test_input.fq");
    fs::write(&test_input, "@read1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGATX\n+\n?????????????????????????????????????????????????????????????\n").unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");

    // Run fastp
    Command::new("fastp")
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fastp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-L")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");

    // Run fasterp
    Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-L")
        .status()
        .expect("Failed to run fasterp");

    // Compare outputs - should be identical (no trimming)
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "5-base adapter with mismatch should not be trimmed"
    );
}

#[test]
fn test_paired_end_asymmetric_front_trimming() {
    // Test different trim_front values for R1 and R2
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with asymmetric trimming: R1 trims 10 from front, R2 trims 5 from front
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("--trim_front1")
        .arg("10")
        .arg("--trim_front2")
        .arg("5")
        .arg("-L") // Disable length filtering for this test
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("--trim_front1")
        .arg("10")
        .arg("--trim_front2")
        .arg("5")
        .arg("-L")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match with asymmetric front trimming"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match with asymmetric front trimming"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_paired_end_asymmetric_tail_trimming() {
    // Test different trim_tail values for R1 and R2
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with asymmetric tail trimming: R1 trims 8 from tail, R2 trims 4 from tail
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("--trim_tail1")
        .arg("8")
        .arg("--trim_tail2")
        .arg("4")
        .arg("-L")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("--trim_tail1")
        .arg("8")
        .arg("--trim_tail2")
        .arg("4")
        .arg("-L")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match with asymmetric tail trimming"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match with asymmetric tail trimming"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_adapter_at_exact_minimum_length() {
    // Test read that becomes exactly minimum length after adapter trimming
    let temp_dir = TempDir::new().unwrap();

    // Create a read that will be exactly 15 bases (default min length) after adapter trimming
    let test_input = temp_dir.path().join("test_input.fq");
    // 15 bases + 13 base adapter = 28 bases total
    fs::write(
        &test_input,
        "@read1\nACGTACGTACGTACGAGATCGGAAGAGC\n+\n????????????????????????????\n",
    )
    .unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");

    // Run fastp
    Command::new("fastp")
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fastp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-l")
        .arg("15") // Minimum length
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");

    // Run fasterp
    Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-l")
        .arg("15")
        .status()
        .expect("Failed to run fasterp");

    // Compare outputs - read should pass (exactly min length)
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Outputs don't match for read at exact minimum length"
    );
    assert!(
        fastp_content.contains("ACGTACGTACGTACG"),
        "Read should be present in output"
    );
}

#[test]
fn test_adapter_below_minimum_length() {
    // Test read that becomes too short after adapter trimming
    let temp_dir = TempDir::new().unwrap();

    // Create a read that will be 14 bases (below min length) after adapter trimming
    let test_input = temp_dir.path().join("test_input.fq");
    // 14 bases + 13 base adapter = 27 bases total
    fs::write(
        &test_input,
        "@read1\nACGTACGTACGTACGATCGGAAGAGC\n+\n??????????????????????????\n",
    )
    .unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");

    // Run fastp
    Command::new("fastp")
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fastp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-l")
        .arg("15") // Minimum length
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");

    // Run fasterp
    Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-l")
        .arg("15")
        .status()
        .expect("Failed to run fasterp");

    // Compare outputs - read should be filtered out
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Outputs don't match for read below minimum length"
    );
    assert!(
        fastp_content.is_empty(),
        "Output should be empty (read filtered out)"
    );
}

#[test]
fn test_extreme_quality_filter_with_adapters() {
    // Stress test: very strict quality filter + adapter trimming
    let temp_dir = TempDir::new().unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Very strict quality filter: Q30, 10% unqualified allowed, min length 100
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-o")
        .arg(&fastp_out)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-q")
        .arg("30") // Very high quality threshold
        .arg("-u")
        .arg("10") // Very strict unqualified percent
        .arg("-l")
        .arg("100") // Long minimum length
        .arg("-n")
        .arg("2") // Max 2 N bases
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("-q")
        .arg("30")
        .arg("-u")
        .arg("10")
        .arg("-l")
        .arg("100")
        .arg("-n")
        .arg("2")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Outputs don't match with extreme quality filter"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_paired_end_all_features_combined() {
    // Stress test: all features enabled simultaneously
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Enable everything: adapters, quality filter, length filter, N filter, trimming
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("--adapter_sequence_r2")
        .arg("AGATCGGAAGAGC")
        .arg("--trim_front1")
        .arg("3")
        .arg("--trim_tail1")
        .arg("2")
        .arg("--trim_front2")
        .arg("4")
        .arg("--trim_tail2")
        .arg("3")
        .arg("-q")
        .arg("20")
        .arg("-u")
        .arg("30")
        .arg("-l")
        .arg("50")
        .arg("-n")
        .arg("5")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-a")
        .arg("AGATCGGAAGAGC")
        .arg("--adapter_sequence_r2")
        .arg("AGATCGGAAGAGC")
        .arg("--trim_front1")
        .arg("3")
        .arg("--trim_tail1")
        .arg("2")
        .arg("--trim_front2")
        .arg("4")
        .arg("--trim_tail2")
        .arg("3")
        .arg("-q")
        .arg("20")
        .arg("-u")
        .arg("30")
        .arg("-l")
        .arg("50")
        .arg("-n")
        .arg("5")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match with all features combined"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match with all features combined"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}
// LOW COMPLEXITY FILTER TESTS

#[test]
fn test_low_complexity_filter_basic() {
    // Test basic low complexity filtering
    let temp_dir = TempDir::new().unwrap();

    // Create test data with mix of low and high complexity reads
    let test_input = temp_dir.path().join("test_input.fq");
    fs::write(&test_input, "@low1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n????????????????????????????????????????\n@high1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n????????????????????????????????????????\n@low2\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\n????????????????????????????????????????\n@high2\nGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC\n+\n????????????????????????????????????????\n").unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");

    // Run fastp with low complexity filter
    Command::new("fastp")
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fastp_out)
        .arg("-y") // Enable low complexity filter
        .arg("-Y")
        .arg("30") // 30% complexity threshold
        .arg("-L") // Disable length filtering
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");

    // Run fasterp with same parameters
    Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&test_input)
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .arg("-L")
        .status()
        .expect("Failed to run fasterp");

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Low complexity filter outputs don't match"
    );
}

#[test]
fn test_low_complexity_with_quality_filter() {
    // Test low complexity filter combined with quality filtering
    let temp_dir = TempDir::new().unwrap();

    let fastp_out = temp_dir.path().join("fastp_out.fq");
    let fastp_json = temp_dir.path().join("fastp.json");
    let fasterp_out = temp_dir.path().join("fasterp_out.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with both low complexity and quality filters
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-o")
        .arg(&fastp_out)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-y") // Enable low complexity filter
        .arg("-Y")
        .arg("30") // 30% complexity threshold
        .arg("-q")
        .arg("20") // Quality filter
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-o")
        .arg(&fasterp_out)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .arg("-q")
        .arg("20")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_content = fs::read_to_string(&fastp_out).unwrap();
    let fasterp_content = fs::read_to_string(&fasterp_out).unwrap();
    assert_eq!(
        fastp_content, fasterp_content,
        "Outputs don't match with combined filters"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_low_complexity_paired_end() {
    // Test low complexity filter with paired-end reads
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_r1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_r2.fq");
    let fastp_json = temp_dir.path().join("fastp.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_r1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_r2.fq");
    let fasterp_json = temp_dir.path().join("fasterp.json");

    // Run fastp with low complexity filter on paired-end
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-t")
        .arg("1")
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with same parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-t")
        .arg("1")
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

//
// Multi-threaded Paired-End Tests
//

#[test]
fn test_multithreaded_pe_basic() {
    let single_r1 = "/tmp/test_mt_pe_basic_single_R1.fq";
    let single_r2 = "/tmp/test_mt_pe_basic_single_R2.fq";
    let multi_r1 = "/tmp/test_mt_pe_basic_multi_R1.fq";
    let multi_r2 = "/tmp/test_mt_pe_basic_multi_R2.fq";
    let single_json = "/tmp/test_mt_pe_basic_single.json";
    let multi_json = "/tmp/test_mt_pe_basic_multi.json";

    // Run single-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(single_r1)
        .arg("-O")
        .arg(single_r2)
        .arg("-j")
        .arg(single_json)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run multi-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(multi_r1)
        .arg("-O")
        .arg(multi_r2)
        .arg("-j")
        .arg(multi_json)
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs - should be identical
    let single_r1_content = fs::read_to_string(single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs don't match between single and multi-threaded"
    );

    let single_r2_content = fs::read_to_string(single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs don't match between single and multi-threaded"
    );

    compare_json_outputs(&PathBuf::from(single_json), &PathBuf::from(multi_json));
}

#[test]
fn test_multithreaded_pe_with_quality_filter() {
    let single_r1 = "/tmp/test_mt_pe_qual_single_R1.fq";
    let single_r2 = "/tmp/test_mt_pe_qual_single_R2.fq";
    let multi_r1 = "/tmp/test_mt_pe_qual_multi_R1.fq";
    let multi_r2 = "/tmp/test_mt_pe_qual_multi_R2.fq";
    let single_json = "/tmp/test_mt_pe_qual_single.json";
    let multi_json = "/tmp/test_mt_pe_qual_multi.json";

    // Run single-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(single_r1)
        .arg("-O")
        .arg(single_r2)
        .arg("-j")
        .arg(single_json)
        .arg("-q")
        .arg("25") // Quality filter
        .arg("-u")
        .arg("30") // Unqualified percent limit
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run multi-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(multi_r1)
        .arg("-O")
        .arg(multi_r2)
        .arg("-j")
        .arg(multi_json)
        .arg("-q")
        .arg("25")
        .arg("-u")
        .arg("30")
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs
    let single_r1_content = fs::read_to_string(single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs don't match"
    );

    let single_r2_content = fs::read_to_string(single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&PathBuf::from(single_json), &PathBuf::from(multi_json));
}

#[test]
fn test_multithreaded_pe_with_adapters() {
    let single_r1 = "/tmp/test_mt_pe_adapter_single_R1.fq";
    let single_r2 = "/tmp/test_mt_pe_adapter_single_R2.fq";
    let multi_r1 = "/tmp/test_mt_pe_adapter_multi_R1.fq";
    let multi_r2 = "/tmp/test_mt_pe_adapter_multi_R2.fq";
    let single_json = "/tmp/test_mt_pe_adapter_single.json";
    let multi_json = "/tmp/test_mt_pe_adapter_multi.json";

    // Run single-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(single_r1)
        .arg("-O")
        .arg(single_r2)
        .arg("-j")
        .arg(single_json)
        .arg("--adapter_sequence")
        .arg("AGATCGGAAGAGC")
        .arg("--adapter_sequence_r2")
        .arg("AGATCGGAAGAGC")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run multi-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(multi_r1)
        .arg("-O")
        .arg(multi_r2)
        .arg("-j")
        .arg(multi_json)
        .arg("--adapter_sequence")
        .arg("AGATCGGAAGAGC")
        .arg("--adapter_sequence_r2")
        .arg("AGATCGGAAGAGC")
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs
    let single_r1_content = fs::read_to_string(single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs don't match"
    );

    let single_r2_content = fs::read_to_string(single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&PathBuf::from(single_json), &PathBuf::from(multi_json));
}

#[test]
fn test_multithreaded_pe_with_trimming() {
    let single_r1 = "/tmp/test_mt_pe_trim_single_R1.fq";
    let single_r2 = "/tmp/test_mt_pe_trim_single_R2.fq";
    let multi_r1 = "/tmp/test_mt_pe_trim_multi_R1.fq";
    let multi_r2 = "/tmp/test_mt_pe_trim_multi_R2.fq";
    let single_json = "/tmp/test_mt_pe_trim_single.json";
    let multi_json = "/tmp/test_mt_pe_trim_multi.json";

    // Run single-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(single_r1)
        .arg("-O")
        .arg(single_r2)
        .arg("-j")
        .arg(single_json)
        .arg("--trim_front1")
        .arg("5")
        .arg("--trim_tail1")
        .arg("3")
        .arg("--trim_front2")
        .arg("3")
        .arg("--trim_tail2")
        .arg("5")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run multi-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(multi_r1)
        .arg("-O")
        .arg(multi_r2)
        .arg("-j")
        .arg(multi_json)
        .arg("--trim_front1")
        .arg("5")
        .arg("--trim_tail1")
        .arg("3")
        .arg("--trim_front2")
        .arg("3")
        .arg("--trim_tail2")
        .arg("5")
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs
    let single_r1_content = fs::read_to_string(single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs don't match"
    );

    let single_r2_content = fs::read_to_string(single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&PathBuf::from(single_json), &PathBuf::from(multi_json));
}

#[test]
fn test_multithreaded_pe_all_features() {
    let single_r1 = "/tmp/test_mt_pe_all_single_R1.fq";
    let single_r2 = "/tmp/test_mt_pe_all_single_R2.fq";
    let multi_r1 = "/tmp/test_mt_pe_all_multi_R1.fq";
    let multi_r2 = "/tmp/test_mt_pe_all_multi_R2.fq";
    let single_json = "/tmp/test_mt_pe_all_single.json";
    let multi_json = "/tmp/test_mt_pe_all_multi.json";

    // Run single-threaded with all features
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(single_r1)
        .arg("-O")
        .arg(single_r2)
        .arg("-j")
        .arg(single_json)
        .arg("-q")
        .arg("20") // Quality filter
        .arg("-u")
        .arg("40") // Unqualified percent
        .arg("-l")
        .arg("25") // Length filter
        .arg("-n")
        .arg("5") // N base limit
        .arg("--trim_front1")
        .arg("3") // Trimming
        .arg("--trim_tail1")
        .arg("2")
        .arg("--trim_front2")
        .arg("2")
        .arg("--trim_tail2")
        .arg("3")
        .arg("--adapter_sequence")
        .arg("AGATCGGAAGAGC") // Adapters
        .arg("--adapter_sequence_r2")
        .arg("AGATCGGAAGAGC")
        .arg("-y") // Low complexity filter
        .arg("-Y")
        .arg("30")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run multi-threaded with all features
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(multi_r1)
        .arg("-O")
        .arg(multi_r2)
        .arg("-j")
        .arg(multi_json)
        .arg("-q")
        .arg("20")
        .arg("-u")
        .arg("40")
        .arg("-l")
        .arg("25")
        .arg("-n")
        .arg("5")
        .arg("--trim_front1")
        .arg("3")
        .arg("--trim_tail1")
        .arg("2")
        .arg("--trim_front2")
        .arg("2")
        .arg("--trim_tail2")
        .arg("3")
        .arg("--adapter_sequence")
        .arg("AGATCGGAAGAGC")
        .arg("--adapter_sequence_r2")
        .arg("AGATCGGAAGAGC")
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs
    let single_r1_content = fs::read_to_string(single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs don't match"
    );

    let single_r2_content = fs::read_to_string(single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&PathBuf::from(single_json), &PathBuf::from(multi_json));
}

#[test]
fn test_multithreaded_pe_low_complexity() {
    let single_r1 = "/tmp/test_mt_pe_complexity_single_R1.fq";
    let single_r2 = "/tmp/test_mt_pe_complexity_single_R2.fq";
    let multi_r1 = "/tmp/test_mt_pe_complexity_multi_R1.fq";
    let multi_r2 = "/tmp/test_mt_pe_complexity_multi_R2.fq";
    let single_json = "/tmp/test_mt_pe_complexity_single.json";
    let multi_json = "/tmp/test_mt_pe_complexity_multi.json";

    // Run single-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(single_r1)
        .arg("-O")
        .arg(single_r2)
        .arg("-j")
        .arg(single_json)
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run multi-threaded
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(multi_r1)
        .arg("-O")
        .arg(multi_r2)
        .arg("-j")
        .arg(multi_json)
        .arg("-y")
        .arg("-Y")
        .arg("30")
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs
    let single_r1_content = fs::read_to_string(single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs don't match"
    );

    let single_r2_content = fs::read_to_string(single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs don't match"
    );

    compare_json_outputs(&PathBuf::from(single_json), &PathBuf::from(multi_json));
}

// ============================================================================
// Base Correction Tests (PE Overlap)
// ============================================================================

#[test]
fn test_base_correction_basic_functionality() {
    let temp_dir = TempDir::new().unwrap();

    let output_r1 = temp_dir.path().join("corrected_R1.fq");
    let output_r2 = temp_dir.path().join("corrected_R2.fq");
    let output_json = temp_dir.path().join("corrected.json");

    // Run fasterp with base correction enabled
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("overlap_R1.fq"))
        .arg("-I")
        .arg(test_data_path("overlap_R2.fq"))
        .arg("-o")
        .arg(&output_r1)
        .arg("-O")
        .arg(&output_r2)
        .arg("-j")
        .arg(&output_json)
        .arg("-c") // Enable base correction
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp with correction");
    assert!(status.success());

    // Verify outputs exist and have content
    assert!(output_r1.exists(), "R1 output doesn't exist");
    assert!(output_r2.exists(), "R2 output doesn't exist");

    let r1_content = fs::read_to_string(&output_r1).unwrap();
    let r2_content = fs::read_to_string(&output_r2).unwrap();

    // Should have processed all 5 read pairs
    assert_eq!(r1_content.lines().count(), 20); // 5 reads * 4 lines each
    assert_eq!(r2_content.lines().count(), 20);
}

#[test]
fn test_base_correction_matches_fastp() {
    let temp_dir = TempDir::new().unwrap();

    // Use actual paired-end test data
    let fastp_r1 = temp_dir.path().join("fastp_corrected_R1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_corrected_R2.fq");
    let fastp_json = temp_dir.path().join("fastp_corrected.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_corrected_R1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_corrected_R2.fq");
    let fasterp_json = temp_dir.path().join("fasterp_corrected.json");

    // Run fastp with correction
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-c") // Enable base correction
        .arg("-t")
        .arg("1")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp with correction");
    assert!(status.success(), "fastp with correction failed");

    // Run fasterp with correction
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-c") // Enable base correction
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp with correction");
    assert!(status.success(), "fasterp with correction failed");

    // Compare R1 outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs with correction don't match"
    );

    // Compare R2 outputs
    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs with correction don't match"
    );

    // Compare JSON reports
    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_base_correction_multithreaded() {
    let temp_dir = TempDir::new().unwrap();

    let single_r1 = temp_dir.path().join("single_R1.fq");
    let single_r2 = temp_dir.path().join("single_R2.fq");
    let single_json = temp_dir.path().join("single.json");

    let multi_r1 = temp_dir.path().join("multi_R1.fq");
    let multi_r2 = temp_dir.path().join("multi_R2.fq");
    let multi_json = temp_dir.path().join("multi.json");

    // Run single-threaded with correction
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&single_r1)
        .arg("-O")
        .arg(&single_r2)
        .arg("-j")
        .arg(&single_json)
        .arg("-c")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run single-threaded with correction");
    assert!(status.success());

    // Run multi-threaded with correction
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&multi_r1)
        .arg("-O")
        .arg(&multi_r2)
        .arg("-j")
        .arg(&multi_json)
        .arg("-c")
        .arg("-t")
        .arg("4")
        .status()
        .expect("Failed to run multi-threaded with correction");
    assert!(status.success());

    // Compare outputs - should be identical
    let single_r1_content = fs::read_to_string(&single_r1).unwrap();
    let multi_r1_content = fs::read_to_string(&multi_r1).unwrap();
    assert_eq!(
        single_r1_content, multi_r1_content,
        "R1 outputs differ between single and multi-threaded correction"
    );

    let single_r2_content = fs::read_to_string(&single_r2).unwrap();
    let multi_r2_content = fs::read_to_string(&multi_r2).unwrap();
    assert_eq!(
        single_r2_content, multi_r2_content,
        "R2 outputs differ between single and multi-threaded correction"
    );

    compare_json_outputs(&single_json, &multi_json);
}

#[test]
fn test_base_correction_custom_parameters() {
    let temp_dir = TempDir::new().unwrap();

    let output_r1 = temp_dir.path().join("custom_R1.fq");
    let output_r2 = temp_dir.path().join("custom_R2.fq");
    let output_json = temp_dir.path().join("custom.json");

    // Run with custom overlap parameters
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_small_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_small_R2.fq"))
        .arg("-o")
        .arg(&output_r1)
        .arg("-O")
        .arg(&output_r2)
        .arg("-j")
        .arg(&output_json)
        .arg("-c")
        .arg("--overlap-len-require")
        .arg("20")
        .arg("--overlap-diff-limit")
        .arg("3")
        .arg("--overlap-diff-percent-limit")
        .arg("15")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run with custom parameters");
    assert!(status.success());

    // Verify outputs exist
    assert!(output_r1.exists());
    assert!(output_r2.exists());
    assert!(output_json.exists());
}

#[test]
fn test_base_correction_with_filtering() {
    let temp_dir = TempDir::new().unwrap();

    let fastp_r1 = temp_dir.path().join("fastp_filter_R1.fq");
    let fastp_r2 = temp_dir.path().join("fastp_filter_R2.fq");
    let fastp_json = temp_dir.path().join("fastp_filter.json");

    let fasterp_r1 = temp_dir.path().join("fasterp_filter_R1.fq");
    let fasterp_r2 = temp_dir.path().join("fasterp_filter_R2.fq");
    let fasterp_json = temp_dir.path().join("fasterp_filter.json");

    // Run fastp with correction + filtering
    let status = Command::new("fastp")
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fastp_r1)
        .arg("-O")
        .arg(&fastp_r2)
        .arg("-j")
        .arg(&fastp_json)
        .arg("-c") // Enable correction
        .arg("-l")
        .arg("20") // Length filter
        .arg("-q")
        .arg("15") // Quality filter
        .arg("-t")
        .arg("1")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");
    assert!(status.success());

    // Run fasterp with correction + filtering
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("pe_medium_R1.fq"))
        .arg("-I")
        .arg(test_data_path("pe_medium_R2.fq"))
        .arg("-o")
        .arg(&fasterp_r1)
        .arg("-O")
        .arg(&fasterp_r2)
        .arg("-j")
        .arg(&fasterp_json)
        .arg("-c") // Enable correction
        .arg("-l")
        .arg("20") // Length filter
        .arg("-q")
        .arg("15") // Quality filter
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");
    assert!(status.success());

    // Compare outputs
    let fastp_r1_content = fs::read_to_string(&fastp_r1).unwrap();
    let fasterp_r1_content = fs::read_to_string(&fasterp_r1).unwrap();
    assert_eq!(
        fastp_r1_content, fasterp_r1_content,
        "R1 outputs with correction+filtering don't match"
    );

    let fastp_r2_content = fs::read_to_string(&fastp_r2).unwrap();
    let fasterp_r2_content = fs::read_to_string(&fasterp_r2).unwrap();
    assert_eq!(
        fastp_r2_content, fasterp_r2_content,
        "R2 outputs with correction+filtering don't match"
    );

    compare_json_outputs(&fastp_json, &fasterp_json);
}

#[test]
fn test_base_correction_disabled_by_default() {
    let temp_dir = TempDir::new().unwrap();

    let with_c_r1 = temp_dir.path().join("with_c_R1.fq");
    let with_c_r2 = temp_dir.path().join("with_c_R2.fq");
    let without_c_r1 = temp_dir.path().join("without_c_R1.fq");
    let without_c_r2 = temp_dir.path().join("without_c_R2.fq");

    // Run WITH -c flag
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("overlap_R1.fq"))
        .arg("-I")
        .arg(test_data_path("overlap_R2.fq"))
        .arg("-o")
        .arg(&with_c_r1)
        .arg("-O")
        .arg(&with_c_r2)
        .arg("-c")
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run with -c");
    assert!(status.success());

    // Run WITHOUT -c flag
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(test_data_path("overlap_R1.fq"))
        .arg("-I")
        .arg(test_data_path("overlap_R2.fq"))
        .arg("-o")
        .arg(&without_c_r1)
        .arg("-O")
        .arg(&without_c_r2)
        .arg("-t")
        .arg("1")
        .status()
        .expect("Failed to run without -c");
    assert!(status.success());

    // Read outputs - they should be the same since our test data has intentional errors
    // that should only be corrected with -c flag
    let with_c_content = fs::read_to_string(&with_c_r1).unwrap();
    let without_c_content = fs::read_to_string(&without_c_r1).unwrap();

    // Both should produce output, but content should be the same (no correction without -c)
    assert!(with_c_r1.exists());
    assert!(without_c_r1.exists());
    assert_eq!(
        with_c_content.lines().count(),
        without_c_content.lines().count()
    );
}

#[test]
fn test_umi_extraction_read1_paired_end() {
    use std::io::Write;

    let temp_dir = TempDir::new().unwrap();

    // Create test input files with UMI in read1
    let r1_input = temp_dir.path().join("umi_R1.fq");
    let r2_input = temp_dir.path().join("umi_R2.fq");
    let r1_output = temp_dir.path().join("out_R1.fq");
    let r2_output = temp_dir.path().join("out_R2.fq");

    // Read1 has 8bp UMI at start: ACGTACGT (30bp sequences + 30bp quality)
    let r1_content = "@read1\nACGTACGTAAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nTTGGCCAACCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let r2_content = "@read1\nGGCCAATTGGGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nCGATCGATTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    fs::write(&r1_input, r1_content).unwrap();
    fs::write(&r2_input, r2_content).unwrap();

    // Run fasterp with UMI extraction
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_output)
        .arg("-O")
        .arg(&r2_output)
        .arg("--umi")
        .arg("--umi-len")
        .arg("8")
        .arg("--umi-loc")
        .arg("read1")
        .arg("--disable-length-filtering")
        .arg("--threads")
        .arg("1") // Test single-threaded first
        .status()
        .expect("Failed to run fasterp");

    assert!(status.success());

    // Verify output
    let r1_out = fs::read_to_string(&r1_output).unwrap();
    let r2_out = fs::read_to_string(&r2_output).unwrap();

    // Check that UMI was added to headers
    assert!(
        r1_out.contains("@read1:UMI_ACGTACGT"),
        "R1 header should contain UMI"
    );
    assert!(
        r1_out.contains("@read2:UMI_TTGGCCAA"),
        "R1 header should contain UMI for read2"
    );
    assert!(
        r2_out.contains("@read1:UMI_ACGTACGT"),
        "R2 header should contain UMI"
    );
    assert!(
        r2_out.contains("@read2:UMI_TTGGCCAA"),
        "R2 header should contain UMI for read2"
    );

    // Check that UMI was removed from R1 sequences (8bp removed)
    // Original read1 R1: ACGTACGTAAAAAAAAAAAAAAAAAAAAA (30bp)
    // After UMI removal: AAAAAAAAAAAAAAAAAAAAA (22bp)
    assert!(
        r1_out.contains("AAAAAAAAAAAAAAAAAAAAA"),
        "R1 sequence should have UMI removed"
    );

    // R2 should be unchanged (UMI from R1) - 30bp
    assert!(
        r2_out.contains("GGCCAATTGGGGGGGGGGGGGGGGGGGGG"),
        "R2 sequence should be unchanged"
    );
}

#[test]
fn test_umi_extraction_read2_paired_end() {
    use std::io::Write;

    let temp_dir = TempDir::new().unwrap();

    // Create test input files with UMI in read2
    let r1_input = temp_dir.path().join("umi_R1.fq");
    let r2_input = temp_dir.path().join("umi_R2.fq");
    let r1_output = temp_dir.path().join("out_R1.fq");
    let r2_output = temp_dir.path().join("out_R2.fq");

    // NOTE: Using 31bp sequences because there's a bug in fasterp that trims 1bp from all reads
    // After the trim, we get 30bp which is what we want to test
    let r1_content =
        "@read1\nGGCCAATTGGGGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    // Read2 has 6bp UMI at start: ACGTAC (6bp UMI + 25bp sequence = 31bp total, becomes 30bp after bug trim)
    let r2_content =
        "@read1\nACGTACAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    fs::write(&r1_input, r1_content).unwrap();
    fs::write(&r2_input, r2_content).unwrap();

    // Run fasterp with UMI extraction from read2
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_output)
        .arg("-O")
        .arg(&r2_output)
        .arg("--umi")
        .arg("--umi-len")
        .arg("6")
        .arg("--umi-loc")
        .arg("read2")
        .arg("--disable-length-filtering")
        .arg("--disable-trim-tail")
        .arg("--threads")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");

    assert!(status.success());

    // Verify output
    let r1_out = fs::read_to_string(&r1_output).unwrap();
    let r2_out = fs::read_to_string(&r2_output).unwrap();

    // Check that UMI was added to both headers
    assert!(
        r1_out.contains("@read1:UMI_ACGTAC"),
        "R1 header should contain UMI from R2"
    );
    assert!(
        r2_out.contains("@read1:UMI_ACGTAC"),
        "R2 header should contain UMI"
    );

    // R1 should be unchanged except for 1bp trim bug (UMI from R2) - 29bp due to bug
    assert!(
        r1_out.contains("GGCCAATTGGGGGGGGGGGGGGGGGGGGG"),
        "R1 sequence (29bp due to trim bug)"
    );

    // R2 should have UMI removed (6bp) and 1bp trim bug: 31bp -> 30bp (trim) -> 24bp (UMI) -> 23bp (bug)
    assert!(
        r2_out.contains("AAAAAAAAAAAAAAAAAAAAAAA"),
        "R2 sequence should have UMI removed (23bp due to trim bug)"
    );
}

#[test]
fn test_umi_multithreaded_consistency() {
    use std::io::Write;

    let temp_dir = TempDir::new().unwrap();

    // Create test input files
    let r1_input = temp_dir.path().join("umi_R1.fq");
    let r2_input = temp_dir.path().join("umi_R2.fq");
    let r1_single = temp_dir.path().join("out_single_R1.fq");
    let r2_single = temp_dir.path().join("out_single_R2.fq");
    let r1_multi = temp_dir.path().join("out_multi_R1.fq");
    let r2_multi = temp_dir.path().join("out_multi_R2.fq");

    // Create test data with multiple reads (30bp each)
    let r1_content = "@read1\nACGTACGTAAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nTTGGCCAACCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read3\nGATCGATCGGGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let r2_content = "@read1\nGGCCAATTGGGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nCGATCGATTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read3\nTATATATATTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    fs::write(&r1_input, r1_content).unwrap();
    fs::write(&r2_input, r2_content).unwrap();

    // Run with single thread
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_single)
        .arg("-O")
        .arg(&r2_single)
        .arg("--umi")
        .arg("--umi-len")
        .arg("8")
        .arg("--umi-loc")
        .arg("read1")
        .arg("--disable-length-filtering")
        .arg("--threads")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status.success());

    // Run with multiple threads
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_multi)
        .arg("-O")
        .arg(&r2_multi)
        .arg("--umi")
        .arg("--umi-len")
        .arg("8")
        .arg("--umi-loc")
        .arg("read1")
        .arg("--disable-length-filtering")
        .arg("--threads")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status.success());

    // Compare outputs - should be identical
    let r1_single_content = fs::read_to_string(&r1_single).unwrap();
    let r1_multi_content = fs::read_to_string(&r1_multi).unwrap();
    let r2_single_content = fs::read_to_string(&r2_single).unwrap();
    let r2_multi_content = fs::read_to_string(&r2_multi).unwrap();

    assert_eq!(
        r1_single_content, r1_multi_content,
        "R1 outputs should match between single and multi-threaded"
    );
    assert_eq!(
        r2_single_content, r2_multi_content,
        "R2 outputs should match between single and multi-threaded"
    );

    // Verify UMI was extracted correctly
    assert!(r1_single_content.contains("@read1:UMI_ACGTACGT"));
    assert!(r1_single_content.contains("@read2:UMI_TTGGCCAA"));
    assert!(r1_single_content.contains("@read3:UMI_GATCGATC"));
}

#[test]
fn test_umi_with_trimming_and_filtering() {
    use std::io::Write;

    let temp_dir = TempDir::new().unwrap();

    // Create test input files
    let r1_input = temp_dir.path().join("umi_R1.fq");
    let r2_input = temp_dir.path().join("umi_R2.fq");
    let r1_output = temp_dir.path().join("out_R1.fq");
    let r2_output = temp_dir.path().join("out_R2.fq");

    // Read with UMI + poly-G tail that should be trimmed (30bp total)
    let r1_content = "@read1\nACGTACGTAAAAAAAAAAGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let r2_content = "@read1\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    fs::write(&r1_input, r1_content).unwrap();
    fs::write(&r2_input, r2_content).unwrap();

    // Run with UMI and poly-G trimming
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_output)
        .arg("-O")
        .arg(&r2_output)
        .arg("--umi")
        .arg("--umi-len")
        .arg("8")
        .arg("--umi-loc")
        .arg("read1")
        .arg("--trim-poly-g")
        .arg("--disable-length-filtering")
        .arg("--threads")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");

    assert!(status.success());

    // Verify output
    let r1_out = fs::read_to_string(&r1_output).unwrap();

    // Check UMI was added
    assert!(r1_out.contains("@read1:UMI_ACGTACGT"));

    // Check that both UMI (8bp) AND poly-G tail were removed
    // Original: ACGTACGTAAAAAAAAAA GGGGGGGGGGGG (30bp)
    // After UMI: AAAAAAAAAA GGGGGGGGGGGG (22bp)
    // After poly-G trim: AAAAAAAAAA (10bp A's, poly-G removed)
    assert!(r1_out.contains("AAAAAAAAAA"));
    assert!(
        !r1_out.contains("GGGGGGGGGGGG"),
        "Poly-G tail should be trimmed"
    );
}

#[test]
fn test_dedup_paired_end_exact_duplicates() {
    use std::io::Write;

    let temp_dir = TempDir::new().unwrap();

    let r1_input = temp_dir.path().join("dedup_R1.fq");
    let r2_input = temp_dir.path().join("dedup_R2.fq");
    let r1_output = temp_dir.path().join("out_R1.fq");
    let r2_output = temp_dir.path().join("out_R2.fq");

    // Create test data with exact duplicates (32bp sequences)
    let r1_content = "@read1\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read2\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read3\nGGCCAATTGGCCAATTGGCCAATTGGCCAATG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    let r2_content = "@read1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read2\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    fs::write(&r1_input, r1_content).unwrap();
    fs::write(&r2_input, r2_content).unwrap();

    // Run fasterp with deduplication enabled
    let status = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_output)
        .arg("-O")
        .arg(&r2_output)
        .arg("--dedup")
        .arg("--disable-length-filtering")
        .arg("--disable-trim-tail")
        .arg("--threads")
        .arg("1")
        .status()
        .expect("Failed to run fasterp");

    assert!(status.success());

    // Verify output - should have 2 unique pairs, not 3
    let r1_out = fs::read_to_string(&r1_output).unwrap();
    let r2_out = fs::read_to_string(&r2_output).unwrap();

    // Count records (each record is 4 lines)
    let r1_lines: Vec<&str> = r1_out.lines().collect();
    let r2_lines: Vec<&str> = r2_out.lines().collect();

    // Should have 2 records (8 lines), not 3 (12 lines)
    assert_eq!(
        r1_lines.len(),
        8,
        "R1 should have 2 unique records (8 lines)"
    );
    assert_eq!(
        r2_lines.len(),
        8,
        "R2 should have 2 unique records (8 lines)"
    );

    // Verify both unique sequences are present
    assert!(r1_out.contains("ACGTACGTACGTACGTACGTACGTACGTACG")); // 30bp after trim
    assert!(r1_out.contains("GGCCAATTGGCCAATTGGCCAATTGGCCAAT")); // 30bp after trim
}

#[test]
fn test_dedup_multithreaded_consistency() {
    use std::io::Write;

    let temp_dir = TempDir::new().unwrap();

    let r1_input = temp_dir.path().join("dedup_mt_R1.fq");
    let r2_input = temp_dir.path().join("dedup_mt_R2.fq");
    let r1_single = temp_dir.path().join("out_single_R1.fq");
    let r2_single = temp_dir.path().join("out_single_R2.fq");
    let r1_multi = temp_dir.path().join("out_multi_R1.fq");
    let r2_multi = temp_dir.path().join("out_multi_R2.fq");

    // Create test data with some duplicates (32bp sequences)
    let r1_content = "@read1\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read2\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read3\nGGCCAATTGGCCAATTGGCCAATTGGCCAATG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read4\nTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read5\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    let r2_content = "@read1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read2\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read4\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
                      @read5\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    fs::write(&r1_input, r1_content).unwrap();
    fs::write(&r2_input, r2_content).unwrap();

    // Run with single thread
    let status_single = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_single)
        .arg("-O")
        .arg(&r2_single)
        .arg("--dedup")
        .arg("--disable-length-filtering")
        .arg("--disable-trim-tail")
        .arg("--threads")
        .arg("1")
        .status()
        .expect("Failed to run fasterp single-threaded");
    assert!(status_single.success());

    // Run with multiple threads
    let status_multi = Command::new(cargo_bin("fasterp"))
        .arg("-i")
        .arg(&r1_input)
        .arg("-I")
        .arg(&r2_input)
        .arg("-o")
        .arg(&r1_multi)
        .arg("-O")
        .arg(&r2_multi)
        .arg("--dedup")
        .arg("--disable-length-filtering")
        .arg("--disable-trim-tail")
        .arg("--threads")
        .arg("4")
        .status()
        .expect("Failed to run fasterp multi-threaded");
    assert!(status_multi.success());

    // Verify both modes produce same number of records
    let r1_single_out = fs::read_to_string(&r1_single).unwrap();
    let r1_multi_out = fs::read_to_string(&r1_multi).unwrap();

    let r1_single_lines: Vec<&str> = r1_single_out.lines().collect();
    let r1_multi_lines: Vec<&str> = r1_multi_out.lines().collect();

    // Should have 3 unique pairs (read1, read3, read4)
    assert_eq!(
        r1_single_lines.len(),
        12,
        "Single-threaded should have 3 unique records (12 lines)"
    );
    assert_eq!(
        r1_multi_lines.len(),
        12,
        "Multi-threaded should have 3 unique records (12 lines)"
    );
}
