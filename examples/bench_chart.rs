/// Example: Generate ASCII bar charts comparing fasterp vs fastp performance
///
/// Usage: cargo run --release --example bench_chart

use std::process::Command;
use std::time::Instant;

struct BenchResult {
    name: String,
    fastp_ms: f64,
    fasterp_ms: f64,
}

fn main() {
    println!("Running benchmarks...\n");

    // Build release binary if needed
    println!("Building fasterp...");
    Command::new("cargo")
        .args(["build", "--release"])
        .output()
        .expect("Failed to build fasterp");

    let mut results = Vec::new();

    // Benchmark 1: Small dataset (1k reads)
    if let Some(result) = run_benchmark("test_data/small_1k.fq", "Small dataset (1k reads)", &[]) {
        results.push(result);
    }

    // Benchmark 2: Medium dataset (10k reads)
    if let Some(result) = run_benchmark("test_data/medium_10k.fq", "Medium dataset (10k reads)", &[]) {
        results.push(result);
    }

    // Benchmark 3: Quality trimming
    if let Some(result) = run_benchmark(
        "test_data/medium_10k.fq",
        "With quality trimming (10k reads)",
        &["--cut-tail", "--cut-mean-quality", "20"],
    ) {
        results.push(result);
    }

    println!("\nPerformance Comparison: fasterp vs fastp\n");

    for result in results {
        print_benchmark_chart(&result);
        println!();
    }
}

fn run_benchmark(input: &str, name: &str, extra_args: &[&str]) -> Option<BenchResult> {
    // Check if input file exists
    if !std::path::Path::new(input).exists() {
        eprintln!("Skipping {}: {} not found", name, input);
        return None;
    }

    println!("Running {}...", name);

    // Create temp files
    let output_dir = std::env::temp_dir();
    let fastp_out = output_dir.join("fastp_bench.fq");
    let fastp_json = output_dir.join("fastp_bench.json");
    let fasterp_out = output_dir.join("fasterp_bench.fq");
    let fasterp_json = output_dir.join("fasterp_bench.json");

    // Run fastp 3 times and take median
    let mut fastp_times = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();

        // Convert args for fastp (all dashes to underscores after --)
        let fastp_args: Vec<String> = extra_args.iter()
            .map(|s| {
                if s.starts_with("--") {
                    let rest = &s[2..];
                    format!("--{}", rest.replace('-', "_"))
                } else {
                    s.to_string()
                }
            })
            .collect();

        let status = Command::new("fastp")
            .arg("-i").arg(input)
            .arg("-o").arg(&fastp_out)
            .arg("-j").arg(&fastp_json)
            .args(&fastp_args)
            .stderr(std::process::Stdio::null())
            .status();

        let duration = start.elapsed();
        if let Ok(status) = status {
            if status.success() {
                fastp_times.push(duration.as_secs_f64() * 1000.0);
            }
        }
    }

    if fastp_times.is_empty() {
        eprintln!("  Warning: Failed to run fastp for {}", name);
        return None;
    }
    fastp_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let fastp_ms = if fastp_times.len() >= 3 {
        fastp_times[1] // median
    } else {
        fastp_times[0] // just use first if we don't have 3
    };

    // Run fasterp 3 times and take median
    let mut fasterp_times = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        let status = Command::new("./target/release/fasterp")
            .arg("-i").arg(input)
            .arg("-o").arg(&fasterp_out)
            .arg("-j").arg(&fasterp_json)
            .args(extra_args)
            .status();

        let duration = start.elapsed();
        if let Ok(status) = status {
            if status.success() {
                fasterp_times.push(duration.as_secs_f64() * 1000.0);
            }
        }
    }

    if fasterp_times.is_empty() {
        eprintln!("  Warning: Failed to run fasterp for {}", name);
        return None;
    }
    fasterp_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let fasterp_ms = if fasterp_times.len() >= 3 {
        fasterp_times[1] // median
    } else {
        fasterp_times[0] // just use first if we don't have 3
    };

    // Cleanup
    let _ = std::fs::remove_file(fastp_out);
    let _ = std::fs::remove_file(fastp_json);
    let _ = std::fs::remove_file(fasterp_out);
    let _ = std::fs::remove_file(fasterp_json);

    Some(BenchResult {
        name: name.to_string(),
        fastp_ms,
        fasterp_ms,
    })
}

fn print_benchmark_chart(result: &BenchResult) {
    let speedup = result.fastp_ms / result.fasterp_ms;

    // Calculate bar lengths (fastp gets full width, fasterp is proportional)
    let max_bar_width = 28;
    let fastp_bar_len = max_bar_width;
    let fasterp_bar_len = (result.fasterp_ms / result.fastp_ms * max_bar_width as f64) as usize;

    println!("{}:", result.name);
    println!(
        "  fastp    {} {}ms",
        "█".repeat(fastp_bar_len),
        result.fastp_ms as u32
    );
    println!(
        "  fasterp  {} {}ms  ⚡ {:.1}× faster",
        "█".repeat(fasterp_bar_len),
        result.fasterp_ms as u32,
        speedup
    );
}
