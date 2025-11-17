use criterion::{Criterion, criterion_group, criterion_main};
use std::process::Command;
use tempfile::NamedTempFile;

/// Helper to run fasterp with given arguments
fn run_fasterp(input: &str, args: &[&str]) -> std::time::Duration {
    let output = NamedTempFile::new().unwrap();
    let json = NamedTempFile::new().unwrap();

    let start = std::time::Instant::now();

    let status = Command::new("./target/release/fasterp")
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output.path())
        .arg("-j")
        .arg(json.path())
        .args(args)
        .status()
        .expect("Failed to run fasterp");

    let duration = start.elapsed();
    assert!(status.success(), "fasterp failed");
    duration
}

/// Helper to run fastp with given arguments
fn run_fastp(input: &str, args: &[&str]) -> std::time::Duration {
    let output = NamedTempFile::new().unwrap();
    let json = NamedTempFile::new().unwrap();

    let start = std::time::Instant::now();

    let status = Command::new("fastp")
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output.path())
        .arg("-j")
        .arg(json.path())
        .args(args)
        .stderr(std::process::Stdio::null())
        .status()
        .expect("Failed to run fastp");

    let duration = start.elapsed();
    assert!(status.success(), "fastp failed");
    duration
}

fn bench_basic_filtering(c: &mut Criterion) {
    let mut group = c.benchmark_group("basic_filtering");
    group.sample_size(10); // Reduce sample size for large dataset

    group.bench_function("fasterp", |b| {
        b.iter(|| run_fasterp("test_data/huge_5m.fq", &[]));
    });

    group.bench_function("fastp", |b| {
        b.iter(|| run_fastp("test_data/huge_5m.fq", &[]));
    });

    group.finish();
}

fn bench_quality_trimming(c: &mut Criterion) {
    let mut group = c.benchmark_group("quality_trimming");
    group.sample_size(10); // Reduce sample size for large dataset

    group.bench_function("fasterp", |b| {
        b.iter(|| {
            run_fasterp(
                "test_data/huge_5m.fq",
                &["--cut-tail", "--cut-mean-quality", "20"],
            )
        });
    });

    group.bench_function("fastp", |b| {
        b.iter(|| {
            run_fastp(
                "test_data/huge_5m.fq",
                &["--cut_tail", "--cut_mean_quality", "20"],
            )
        });
    });

    group.finish();
}

fn bench_aggressive_trimming(c: &mut Criterion) {
    let mut group = c.benchmark_group("aggressive_trimming");
    group.sample_size(10); // Reduce sample size for large dataset

    group.bench_function("fasterp", |b| {
        b.iter(|| {
            run_fasterp(
                "test_data/huge_5m.fq",
                &[
                    "--trim-front",
                    "5",
                    "--trim-tail",
                    "5",
                    "--cut-front",
                    "--cut-tail",
                    "--cut-mean-quality",
                    "20",
                    "--trim-poly-g",
                ],
            )
        });
    });

    group.bench_function("fastp", |b| {
        b.iter(|| {
            run_fastp(
                "test_data/huge_5m.fq",
                &[
                    "--trim_front1",
                    "5",
                    "--trim_tail1",
                    "5",
                    "--cut_front",
                    "--cut_tail",
                    "--cut_mean_quality",
                    "20",
                    "--trim_poly_g",
                ],
            )
        });
    });

    group.finish();
}

fn bench_multithreading(c: &mut Criterion) {
    let mut group = c.benchmark_group("multithreading");
    group.sample_size(10); // Reduce sample size for large dataset

    // fasterp single-threaded
    group.bench_function("fasterp_1thread", |b| {
        b.iter(|| run_fasterp("test_data/huge_5m.fq", &["-t", "1"]));
    });

    // fasterp multi-threaded
    group.bench_function("fasterp_4threads", |b| {
        b.iter(|| run_fasterp("test_data/huge_5m.fq", &["-t", "4"]));
    });

    // fastp single-threaded
    group.bench_function("fastp_1thread", |b| {
        b.iter(|| run_fastp("test_data/huge_5m.fq", &["-w", "1"]));
    });

    // fastp multi-threaded
    group.bench_function("fastp_4threads", |b| {
        b.iter(|| run_fastp("test_data/huge_5m.fq", &["-w", "4"]));
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_basic_filtering,
    bench_quality_trimming,
    bench_aggressive_trimming,
    bench_multithreading
);
criterion_main!(benches);
