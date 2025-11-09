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
            "Kmer {} counts don't match",
            kmer
        );
    }
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
