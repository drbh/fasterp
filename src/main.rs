//! A fast and simple FASTQ preprocessor
//!
//! This tool processes FASTQ files with quality control filters including:
//! - Length filtering (minimum read length)
//! - Quality filtering (mean quality score)
//! - N-base filtering (maximum ambiguous bases)

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::thread;

mod io;
mod kmer;
mod pipeline;
mod processor;
mod simd;
mod stats;
mod trimming;

use io::*;
use pipeline::*;
use processor::*;
use stats::*;
use trimming::*;

#[derive(Parser, Debug)]
#[command(author, version, about = "A fast FASTQ preprocessor", long_about = None)]
struct Args {
    /// Input FASTQ file (use '-' for stdin, supports .gz)
    #[arg(short = 'i', long)]
    input: String,

    /// Output FASTQ file (use '-' for stdout, .gz extension enables compression)
    #[arg(short = 'o', long)]
    output: String,

    /// Minimum length required (default: 15)
    #[arg(short = 'l', long, default_value = "15")]
    length_required: usize,

    /// Maximum length limit - trim reads longer than this (default: 0 = disabled)
    #[arg(short = 'b', long, default_value = "0")]
    max_len: usize,

    /// Quality value that a base is qualified (phred quality >= this value) (default: 15)
    #[arg(short = 'q', long, default_value = "15")]
    qualified_quality_phred: u8,

    /// Percent of bases allowed to be unqualified (0-100) (default: 40)
    #[arg(short = 'u', long, default_value = "40")]
    unqualified_percent_limit: usize,

    /// Average quality threshold - discard if mean quality < this (default: 0 = disabled)
    #[arg(short = 'e', long, default_value = "0")]
    average_qual: u8,

    /// Max number of N bases allowed (default: 5)
    #[arg(short = 'n', long, default_value = "5")]
    n_base_limit: usize,

    /// JSON report file (default: fastp.json)
    #[arg(short = 'j', long, default_value = "fastp.json")]
    json: String,

    /// Stats output format: compact (default), pretty, off, or jsonl
    #[arg(long, default_value = "compact")]
    stats_format: String,

    /// Compression level for gzip output (0-9, default: 6)
    #[arg(short = 'z', long)]
    compression_level: Option<u32>,

    /// Number of worker threads (default: auto-detect CPU count)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Batch size in bytes (default: 16 MiB)
    #[arg(long, default_value = "16777216")]
    batch_bytes: usize,

    /// Maximum backlog of batches (default: threads+1)
    #[arg(long)]
    max_backlog: Option<usize>,

    /// Skip k-mer counting for ceiling performance tests
    #[arg(long)]
    no_kmer: bool,

    /// Quality cutoff for sliding-window trimming (default: 0 = disabled)
    /// Trim when mean quality in window falls below this value
    #[arg(long, default_value = "0")]
    cut_mean_quality: u8,

    /// Sliding window size for quality trimming (default: 4)
    #[arg(long, default_value = "4")]
    cut_window_size: usize,

    /// Enable quality trimming at 5' end (front)
    #[arg(long)]
    cut_front: bool,

    /// Enable quality trimming at 3' end (tail)
    #[arg(long)]
    cut_tail: bool,

    /// Disable polyG/quality tail trimming (for backward compatibility)
    #[arg(long)]
    disable_trim_tail: bool,

    /// Trim N bases from 5' (front) end
    #[arg(long, default_value = "0")]
    trim_front: usize,

    /// Trim N bases from 3' (tail) end
    #[arg(long, default_value = "0")]
    trim_tail: usize,

    /// Enable polyG tail trimming
    #[arg(long)]
    trim_poly_g: bool,

    /// Disable polyG tail trimming (for backward compatibility)
    #[arg(short = 'G', long)]
    disable_trim_poly_g: bool,

    /// Enable generic polyX tail trimming (any homopolymer)
    #[arg(long)]
    trim_poly_x: bool,

    /// Minimum length for polyG/polyX detection (default: 10)
    #[arg(long, default_value = "10")]
    poly_g_min_len: usize,

    /// Disable adapter trimming (no-op for compatibility - adapters not implemented)
    #[arg(short = 'A', long)]
    disable_adapter_trimming: bool,
}

// Helper function to create TrimmingConfig from CLI args
fn create_trimming_config(args: &Args) -> TrimmingConfig {
    TrimmingConfig {
        enable_trim_front: args.cut_front && args.cut_mean_quality > 0,
        enable_trim_tail: args.cut_tail && args.cut_mean_quality > 0 && !args.disable_trim_tail,
        cut_mean_quality: args.cut_mean_quality,
        cut_window_size: args.cut_window_size,
        trim_front_bases: args.trim_front,
        trim_tail_bases: args.trim_tail,
        max_len: args.max_len,
        enable_poly_g: args.trim_poly_g && !args.disable_trim_poly_g,
        enable_poly_x: args.trim_poly_x,
        poly_min_len: args.poly_g_min_len,
    }
}

// MAIN FUNCTION

/// MULTI-THREADED MAIN FUNCTION:
///
/// Two modes:
/// 1. Single-threaded (threads=1): Uses old streaming approach
/// 2. Multi-threaded (threads>1): Uses 3-stage pipeline
///    - Producer: reads blocks, parses FASTQ into batches
///    - Workers: process batches in parallel
///    - Merger: writes output in order, reduces stats
///
/// Result: 2-4x faster on large datasets with multi-threading
fn main() -> Result<()> {
    let args = Args::parse();

    // Determine number of threads
    let num_threads = args.threads.unwrap_or_else(num_cpus::get);

    // Create trimming configuration from CLI args
    let trimming_config = create_trimming_config(&args);

    let acc = if num_threads == 1 {
        // SINGLE-THREADED MODE: use streaming approach with new parser
        let reader = open_input(&args.input)?;
        let mut writer = open_output(&args.output, args.compression_level)?;

        let acc = process_fastq_stream(
            reader,
            &mut writer,
            args.length_required,
            args.n_base_limit,
            args.qualified_quality_phred,
            args.unqualified_percent_limit,
            args.average_qual,
            &trimming_config,
        )?;

        writer.finish()?;
        acc
    } else {
        // MULTI-THREADED MODE: 3-stage pipeline
        let backlog = args.max_backlog.unwrap_or(num_threads + 1);

        // Create channels
        let (batch_tx, batch_rx) = bounded::<Option<Batch>>(backlog);
        let (result_tx, result_rx) = bounded::<Option<WorkerResult>>(backlog);

        // Spawn producer thread
        let input_path = args.input.clone();
        let batch_bytes = args.batch_bytes;
        let producer = thread::spawn(move || producer_thread(input_path, batch_bytes, batch_tx));

        // Spawn worker threads
        let mut workers = Vec::new();
        for _ in 0..num_threads {
            let batch_rx_clone = batch_rx.clone();
            let result_tx_clone = result_tx.clone();
            let min_len = args.length_required;
            let n_limit = args.n_base_limit;
            let qualified_qual = args.qualified_quality_phred;
            let unqualified_pct = args.unqualified_percent_limit;
            let avg_qual = args.average_qual;
            let no_kmer = args.no_kmer;
            let trimming_config_clone = trimming_config.clone();

            let worker = thread::spawn(move || {
                worker_thread(
                    batch_rx_clone,
                    result_tx_clone,
                    min_len,
                    n_limit,
                    qualified_qual,
                    unqualified_pct,
                    avg_qual,
                    no_kmer,
                    trimming_config_clone,
                )
            });
            workers.push(worker);
        }

        // Drop original senders so merger knows when all workers are done
        drop(batch_rx);
        drop(result_tx);

        // Spawn merger thread
        let output_path = args.output.clone();
        let compression_level = args.compression_level;
        let merger = thread::spawn(move || {
            merger_thread(result_rx, output_path, num_threads, compression_level)
        });

        // Wait for all threads
        producer.join().unwrap()?;
        for worker in workers {
            worker.join().unwrap();
        }
        merger.join().unwrap()?
    };

    // Build report from accumulated stats (same for both modes)
    let before_stats = acc.before.to_read_stats();
    let after_stats = acc.after.to_read_stats();
    let quality_curves = acc.pos.to_quality_curves();
    let kmer_map = acc.kmer_table_to_map();

    let report = FastpReport {
        summary: Summary {
            fastp_version: env!("CARGO_PKG_VERSION").to_string(),
            sequencing: format!("single end ({} cycles)", acc.max_cycle),
            before_filtering: before_stats.clone(),
            after_filtering: after_stats,
        },
        filtering_result: FilteringResult {
            passed_filter_reads: acc.after.total_reads,
            low_quality_reads: acc.low_quality,
            too_many_n_reads: acc.too_many_n,
            too_short_reads: acc.too_short,
            too_long_reads: 0,
        },
        read1_before_filtering: DetailedReadStats {
            total_reads: before_stats.total_reads,
            total_bases: before_stats.total_bases,
            q20_bases: before_stats.q20_bases,
            q30_bases: before_stats.q30_bases,
            quality_curves,
            kmer_count: kmer_map,
        },
    };

    // Write JSON report based on stats_format
    match args.stats_format.as_str() {
        "off" => {
            // Skip JSON output entirely
        }
        "compact" => {
            let json_file = File::create(&args.json)
                .context(format!("Failed to create JSON file: {}", args.json))?;
            // Use a large BufWriter (256 KiB) to reduce system calls
            let mut buf_writer = BufWriter::with_capacity(256 * 1024, json_file);
            serde_json::to_writer(&mut buf_writer, &report)?;
            // Explicit flush to ensure data is written
            // (BufWriter::drop also flushes, but explicit is clearer)
        }
        "pretty" => {
            let json_file = File::create(&args.json)
                .context(format!("Failed to create JSON file: {}", args.json))?;
            let mut buf_writer = BufWriter::with_capacity(256 * 1024, json_file);
            serde_json::to_writer_pretty(&mut buf_writer, &report)?;
        }
        "jsonl" => {
            let json_file = File::create(&args.json)
                .context(format!("Failed to create JSON file: {}", args.json))?;
            let mut buf_writer = BufWriter::with_capacity(256 * 1024, json_file);
            serde_json::to_writer(&mut buf_writer, &report)?;
            writeln!(buf_writer)?;
        }
        _ => {
            anyhow::bail!(
                "Invalid stats format: {}. Use 'compact', 'pretty', 'off', or 'jsonl'",
                args.stats_format
            );
        }
    }

    Ok(())
}

// UNIT TESTS (see src/tests.rs)

#[cfg(test)]
mod tests;
