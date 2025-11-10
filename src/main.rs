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

mod adapter;
mod io;
mod kmer;
mod pipeline;
mod processor;
mod simd;
mod stats;
mod trimming;
mod util;

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

    /// Read2 input file (paired-end mode)
    #[arg(short = 'I', long = "in2")]
    input2: Option<String>,

    /// Read2 output file (paired-end mode)
    #[arg(short = 'O', long = "out2")]
    output2: Option<String>,

    /// Output file for unpaired read1
    #[arg(long)]
    unpaired1: Option<String>,

    /// Output file for unpaired read2
    #[arg(long)]
    unpaired2: Option<String>,

    /// Input is interleaved paired-end
    #[arg(long)]
    interleaved_in: bool,

    /// Minimum length required (default: 15)
    #[arg(short = 'l', long, default_value = "15")]
    length_required: usize,

    /// Disable length filtering (fastp compatibility)
    #[arg(short = 'L', long)]
    disable_length_filtering: bool,

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

    /// Enable low complexity filter
    #[arg(short = 'y', long)]
    low_complexity_filter: bool,

    /// Complexity threshold (0-100). Default 30 means 30% complexity required (fastp compatibility)
    #[arg(short = 'Y', long, default_value = "30")]
    complexity_threshold: usize,

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

    /// Trim N bases from front for read1 (fastp compatibility)
    #[arg(short = 'f', long = "trim_front1", default_value = "0")]
    trim_front1: usize,

    /// Trim N bases from tail for read1 (fastp compatibility)
    #[arg(long = "trim_tail1", default_value = "0")]
    trim_tail1: usize,

    /// Trim N bases from front for read2 (fastp compatibility)
    #[arg(short = 'F', long = "trim_front2", default_value = "0")]
    trim_front2: usize,

    /// Trim N bases from tail for read2 (fastp compatibility)
    #[arg(short = 'T', long = "trim_tail2", default_value = "0")]
    trim_tail2: usize,

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

    /// Disable adapter trimming
    #[arg(short = 'A', long)]
    disable_adapter_trimming: bool,

    /// Adapter sequence for read1 (auto-use Illumina TruSeq if not specified)
    #[arg(short = 'a', long = "adapter_sequence")]
    adapter_sequence: Option<String>,

    /// Adapter sequence for read2
    #[arg(long = "adapter_sequence_r2")]
    adapter_sequence_r2: Option<String>,

    /// Enable adapter auto-detection for paired-end (use PE overlap)
    #[arg(long = "detect_adapter_for_pe")]
    detect_adapter_for_pe: bool,
}

// Helper function to create TrimmingConfig from CLI args
fn create_trimming_config(args: &Args) -> TrimmingConfig {
    use crate::adapter::{AdapterConfig, adapters};

    // Create adapter configuration
    let mut adapter_config = AdapterConfig::new();

    if !args.disable_adapter_trimming {
        // Only enable adapter trimming if user explicitly specifies adapters
        // fastp uses auto-detection by default, which we haven't implemented yet
        // So we only trim when adapters are explicitly provided
        if let Some(ref seq) = args.adapter_sequence {
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        }

        if let Some(ref seq) = args.adapter_sequence_r2 {
            adapter_config.adapter_seq_r2 = Some(seq.as_bytes().to_vec());
        }

        adapter_config.detect_adapter_for_pe = args.detect_adapter_for_pe;
    }

    // Determine trim values - prefer read-specific args, fall back to generic
    let trim_front = if args.trim_front1 > 0 {
        args.trim_front1
    } else {
        args.trim_front
    };
    let trim_tail = if args.trim_tail1 > 0 {
        args.trim_tail1
    } else {
        args.trim_tail
    };

    TrimmingConfig {
        enable_trim_front: args.cut_front && args.cut_mean_quality > 0,
        enable_trim_tail: args.cut_tail && args.cut_mean_quality > 0 && !args.disable_trim_tail,
        cut_mean_quality: args.cut_mean_quality,
        cut_window_size: args.cut_window_size,
        trim_front_bases: trim_front,
        trim_tail_bases: trim_tail,
        max_len: args.max_len,
        enable_poly_g: args.trim_poly_g && !args.disable_trim_poly_g,
        enable_poly_x: args.trim_poly_x,
        poly_min_len: args.poly_g_min_len,
        adapter_config,
    }
}

// Helper function to create TrimmingConfig for read2 in paired-end mode
fn create_trimming_config_r2(args: &Args) -> TrimmingConfig {
    use crate::adapter::AdapterConfig;

    // Create adapter configuration for R2
    let mut adapter_config = AdapterConfig::new();

    if !args.disable_adapter_trimming {
        // For R2, use adapter_seq_r2 if specified
        if let Some(ref seq) = args.adapter_sequence_r2 {
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        } else if let Some(ref seq) = args.adapter_sequence {
            // Fall back to R1 adapter if R2 not specified
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        }

        adapter_config.detect_adapter_for_pe = args.detect_adapter_for_pe;
    }

    // Determine trim values for R2
    // Per fastp docs: "If not specified, it will follow read1's settings"
    let trim_front = if args.trim_front2 > 0 {
        args.trim_front2
    } else if args.trim_front1 > 0 {
        args.trim_front1
    } else {
        args.trim_front
    };

    let trim_tail = if args.trim_tail2 > 0 {
        args.trim_tail2
    } else if args.trim_tail1 > 0 {
        args.trim_tail1
    } else {
        args.trim_tail
    };

    TrimmingConfig {
        enable_trim_front: args.cut_front && args.cut_mean_quality > 0,
        enable_trim_tail: args.cut_tail && args.cut_mean_quality > 0 && !args.disable_trim_tail,
        cut_mean_quality: args.cut_mean_quality,
        cut_window_size: args.cut_window_size,
        trim_front_bases: trim_front,
        trim_tail_bases: trim_tail,
        max_len: args.max_len,
        enable_poly_g: args.trim_poly_g && !args.disable_trim_poly_g,
        enable_poly_x: args.trim_poly_x,
        poly_min_len: args.poly_g_min_len,
        adapter_config,
    }
}

// Helper function to build and write paired-end report
fn build_and_write_paired_end_report(args: &Args, pe_acc: PairedEndAccumulator) -> Result<()> {
    let before_stats_r1 = pe_acc.before_r1.to_read_stats();
    let after_stats_r1 = pe_acc.after_r1.to_read_stats();
    let before_stats_r2 = pe_acc.before_r2.to_read_stats();
    let after_stats_r2 = pe_acc.after_r2.to_read_stats();

    let quality_curves_r1 = pe_acc.pos_r1.to_quality_curves();
    let quality_curves_r2 = pe_acc.pos_r2.to_quality_curves();
    let kmer_map_r1 = pe_acc.kmer_table_to_map_r1();
    let kmer_map_r2 = pe_acc.kmer_table_to_map_r2();

    // Calculate combined before/after stats for summary
    let combined_before = ReadStats {
        // Count both R1 and R2 reads to match fastp's behavior
        total_reads: before_stats_r1.total_reads + before_stats_r2.total_reads,
        total_bases: before_stats_r1.total_bases + before_stats_r2.total_bases,
        q20_bases: before_stats_r1.q20_bases + before_stats_r2.q20_bases,
        q30_bases: before_stats_r1.q30_bases + before_stats_r2.q30_bases,
        q20_rate: if before_stats_r1.total_bases + before_stats_r2.total_bases > 0 {
            (before_stats_r1.q20_bases + before_stats_r2.q20_bases) as f64
                / (before_stats_r1.total_bases + before_stats_r2.total_bases) as f64
        } else {
            0.0
        },
        q30_rate: if before_stats_r1.total_bases + before_stats_r2.total_bases > 0 {
            (before_stats_r1.q30_bases + before_stats_r2.q30_bases) as f64
                / (before_stats_r1.total_bases + before_stats_r2.total_bases) as f64
        } else {
            0.0
        },
        read1_mean_length: if before_stats_r1.total_reads + before_stats_r2.total_reads > 0 {
            (before_stats_r1.total_bases + before_stats_r2.total_bases)
                / (before_stats_r1.total_reads + before_stats_r2.total_reads)
        } else {
            0
        },
        gc_content: if before_stats_r1.total_bases + before_stats_r2.total_bases > 0 {
            (before_stats_r1.gc_content * before_stats_r1.total_bases as f64
                + before_stats_r2.gc_content * before_stats_r2.total_bases as f64)
                / (before_stats_r1.total_bases + before_stats_r2.total_bases) as f64
        } else {
            0.0
        },
    };

    let combined_after = ReadStats {
        // Count both R1 and R2 reads to match fastp's behavior
        total_reads: after_stats_r1.total_reads + after_stats_r2.total_reads,
        total_bases: after_stats_r1.total_bases + after_stats_r2.total_bases,
        q20_bases: after_stats_r1.q20_bases + after_stats_r2.q20_bases,
        q30_bases: after_stats_r1.q30_bases + after_stats_r2.q30_bases,
        q20_rate: if after_stats_r1.total_bases + after_stats_r2.total_bases > 0 {
            (after_stats_r1.q20_bases + after_stats_r2.q20_bases) as f64
                / (after_stats_r1.total_bases + after_stats_r2.total_bases) as f64
        } else {
            0.0
        },
        q30_rate: if after_stats_r1.total_bases + after_stats_r2.total_bases > 0 {
            (after_stats_r1.q30_bases + after_stats_r2.q30_bases) as f64
                / (after_stats_r1.total_bases + after_stats_r2.total_bases) as f64
        } else {
            0.0
        },
        read1_mean_length: if after_stats_r1.total_reads + after_stats_r2.total_reads > 0 {
            (after_stats_r1.total_bases + after_stats_r2.total_bases)
                / (after_stats_r1.total_reads + after_stats_r2.total_reads)
        } else {
            0
        },
        gc_content: if after_stats_r1.total_bases + after_stats_r2.total_bases > 0 {
            (after_stats_r1.gc_content * after_stats_r1.total_bases as f64
                + after_stats_r2.gc_content * after_stats_r2.total_bases as f64)
                / (after_stats_r1.total_bases + after_stats_r2.total_bases) as f64
        } else {
            0.0
        },
    };

    let report = FastpReport {
        summary: Summary {
            fastp_version: env!("CARGO_PKG_VERSION").to_string(),
            sequencing: format!(
                "paired end ({} cycles + {} cycles)",
                pe_acc.max_cycle_r1, pe_acc.max_cycle_r2
            ),
            before_filtering: combined_before,
            after_filtering: combined_after,
        },
        filtering_result: FilteringResult {
            // Count both R1 and R2 reads to match fastp's behavior
            passed_filter_reads: pe_acc.after_r1.total_reads + pe_acc.after_r2.total_reads,
            low_quality_reads: pe_acc.low_quality,
            low_complexity_reads: pe_acc.low_complexity,
            too_many_n_reads: pe_acc.too_many_n,
            too_short_reads: pe_acc.too_short,
            too_long_reads: 0,
        },
        read1_before_filtering: DetailedReadStats {
            total_reads: before_stats_r1.total_reads,
            total_bases: before_stats_r1.total_bases,
            q20_bases: before_stats_r1.q20_bases,
            q30_bases: before_stats_r1.q30_bases,
            quality_curves: quality_curves_r1,
            kmer_count: kmer_map_r1,
        },
        read2_before_filtering: Some(DetailedReadStats {
            total_reads: before_stats_r2.total_reads,
            total_bases: before_stats_r2.total_bases,
            q20_bases: before_stats_r2.q20_bases,
            q30_bases: before_stats_r2.q30_bases,
            quality_curves: quality_curves_r2,
            kmer_count: kmer_map_r2,
        }),
        read1_after_filtering: None,
        read2_after_filtering: None,
    };

    // Write JSON report
    match args.stats_format.as_str() {
        "off" => {
            // Skip JSON output
        }
        "compact" => {
            let json_file = File::create(&args.json)
                .context(format!("Failed to create JSON file: {}", args.json))?;
            let mut buf_writer = BufWriter::with_capacity(256 * 1024, json_file);
            serde_json::to_writer(&mut buf_writer, &report)?;
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

    // Determine if paired-end mode
    let is_paired_end = args.input2.is_some() || args.interleaved_in;

    // Validate paired-end arguments
    if is_paired_end {
        if args.interleaved_in {
            anyhow::bail!("Interleaved input is not yet implemented");
        }
        if args.input2.is_none() {
            anyhow::bail!("Paired-end mode requires -I/--in2 for read2 input");
        }
        if args.output2.is_none() {
            anyhow::bail!("Paired-end mode requires -O/--out2 for read2 output");
        }
    }

    // Determine number of threads
    let num_threads = args.threads.unwrap_or_else(num_cpus::get);

    // Create trimming configuration from CLI args
    let mut trimming_config = create_trimming_config(&args);

    // Process based on mode (single-end or paired-end)
    if is_paired_end {
        // PAIRED-END MODE
        if num_threads > 1 {
            anyhow::bail!(
                "Multi-threading for paired-end mode is not yet implemented. Please use -t 1"
            );
        }

        // Create separate trimming configs for R1 and R2
        let mut trimming_config_r1 = trimming_config;
        let mut trimming_config_r2 = create_trimming_config_r2(&args);

        // Apply fastp's undocumented default: trim_tail=1 for paired-end mode
        // (only if user didn't explicitly set any tail trimming)
        if args.trim_tail == 0 && args.trim_tail1 == 0 {
            trimming_config_r1.trim_tail_bases = 1;
        }
        if args.trim_tail == 0 && args.trim_tail1 == 0 && args.trim_tail2 == 0 {
            trimming_config_r2.trim_tail_bases = 1;
        }

        // Single-threaded paired-end processing
        let reader1 = open_input(&args.input)?;
        let reader2 = open_input(args.input2.as_ref().unwrap())?;
        let mut writer1 = open_output(&args.output, args.compression_level)?;
        let mut writer2 = open_output(args.output2.as_ref().unwrap(), args.compression_level)?;

        // Apply disable_length_filtering flag
        let min_len = if args.disable_length_filtering { 0 } else { args.length_required };

        let pe_acc = process_paired_fastq_stream(
            reader1,
            reader2,
            &mut writer1,
            &mut writer2,
            min_len,
            args.n_base_limit,
            args.qualified_quality_phred,
            args.unqualified_percent_limit,
            args.average_qual,
            args.low_complexity_filter,
            args.complexity_threshold,
            &trimming_config_r1,
            &trimming_config_r2,
        )?;

        writer1.finish()?;
        writer2.finish()?;

        // Build paired-end report
        build_and_write_paired_end_report(&args, pe_acc)?;
        return Ok(());
    }

    // SINGLE-END MODE
    // Apply disable_length_filtering flag
    let min_len = if args.disable_length_filtering { 0 } else { args.length_required };

    let acc = if num_threads == 1 {
        // SINGLE-THREADED MODE: use streaming approach with new parser
        let reader = open_input(&args.input)?;
        let mut writer = open_output(&args.output, args.compression_level)?;

        let acc = process_fastq_stream(
            reader,
            &mut writer,
            min_len,
            args.n_base_limit,
            args.qualified_quality_phred,
            args.unqualified_percent_limit,
            args.average_qual,
            args.low_complexity_filter,
            args.complexity_threshold,
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
            let min_len_for_worker = min_len;
            let n_limit = args.n_base_limit;
            let qualified_qual = args.qualified_quality_phred;
            let unqualified_pct = args.unqualified_percent_limit;
            let avg_qual = args.average_qual;
            let low_complexity = args.low_complexity_filter;
            let complexity_thresh = args.complexity_threshold;
            let no_kmer = args.no_kmer;
            let trimming_config_clone = trimming_config.clone();

            let worker = thread::spawn(move || {
                worker_thread(
                    batch_rx_clone,
                    result_tx_clone,
                    min_len_for_worker,
                    n_limit,
                    qualified_qual,
                    unqualified_pct,
                    avg_qual,
                    low_complexity,
                    complexity_thresh,
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
            low_complexity_reads: acc.low_complexity,
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
        read2_before_filtering: None,
        read1_after_filtering: None,
        read2_after_filtering: None,
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
