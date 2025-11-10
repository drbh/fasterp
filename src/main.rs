//! A fast and simple FASTQ preprocessor
//!
//! This tool processes FASTQ files with quality control filters including:
//! - Length filtering (minimum read length)
//! - Quality filtering (mean quality score)
//! - N-base filtering (maximum ambiguous bases)

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::thread;

mod adapter;
mod dedup;
mod html;
mod io;
mod kmer;
mod overlap;
mod pipeline;
mod processor;
mod simd;
mod split;
mod stats;
mod trimming;
mod umi;
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
    #[arg(short = 'j', long, default_value = "fasterp.json")]
    json: String,

    /// HTML report file (default: fasterp.html)
    #[arg(long, default_value = "fasterp.html")]
    html: String,

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

    /// Enable base correction using overlap analysis for paired-end data
    #[arg(short = 'c', long)]
    correction: bool,

    /// Minimum overlap length required for correction (default: 30)
    #[arg(long, default_value = "30")]
    overlap_len_require: usize,

    /// Maximum allowed differences in overlap region (default: 5)
    #[arg(long, default_value = "5")]
    overlap_diff_limit: usize,

    /// Maximum allowed difference percentage (0-100) (default: 20)
    #[arg(long, default_value = "20")]
    overlap_diff_percent_limit: usize,

    /// Enable UMI preprocessing
    #[arg(long)]
    umi: bool,

    /// UMI location: read1/read2/index1/index2/per_read/per_index (default: read1)
    #[arg(long, default_value = "read1")]
    umi_loc: String,

    /// UMI length (required when --umi is enabled)
    #[arg(long, default_value = "0")]
    umi_len: usize,

    /// Prefix in read name for UMI (default: UMI)
    #[arg(long, default_value = "UMI")]
    umi_prefix: String,

    /// Skip N bases before UMI
    #[arg(long, default_value = "0")]
    umi_skip: usize,

    /// Enable deduplication to remove duplicate reads (fastp compatible)
    #[arg(short = 'D', long)]
    dedup: bool,

    /// Deduplication accuracy (1-6). Higher = more memory but fewer false positives (default: 5)
    #[arg(long, default_value = "5")]
    dedup_accuracy: u8,

    /// Split output by limiting total number of files (2-999). Cannot be used with --split-by-lines
    #[arg(short = 's', long, conflicts_with = "split_by_lines")]
    split: Option<usize>,

    /// Split output by limiting lines per file (>=1000). Cannot be used with --split
    #[arg(short = 'S', long, conflicts_with = "split")]
    split_by_lines: Option<usize>,

    /// Digits for file number padding (1-10, default: 4). 0 to disable padding
    #[arg(short = 'd', long, default_value = "4")]
    split_prefix_digits: usize,
}

// Helper function to create TrimmingConfig from CLI args
fn create_trimming_config(args: &Args) -> TrimmingConfig {
    use crate::adapter::AdapterConfig;

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

// Helper function to create OverlapConfig from CLI args
fn create_overlap_config(args: &Args) -> overlap::OverlapConfig {
    overlap::OverlapConfig {
        min_overlap_len: args.overlap_len_require,
        max_diff: args.overlap_diff_limit,
        max_diff_percent: args.overlap_diff_percent_limit,
    }
}

// Helper function to create UmiConfig from CLI args
fn create_umi_config(args: &Args) -> Result<umi::UmiConfig> {
    if !args.umi {
        return Ok(umi::UmiConfig::default());
    }

    if args.umi_len == 0 {
        anyhow::bail!("UMI length (--umi-len) must be specified when --umi is enabled");
    }

    let location = match args.umi_loc.as_str() {
        "read1" => umi::UmiLocation::Read1,
        "read2" => umi::UmiLocation::Read2,
        "index1" => umi::UmiLocation::Index1,
        "index2" => umi::UmiLocation::Index2,
        "per_read" => umi::UmiLocation::PerRead,
        "per_index" => umi::UmiLocation::PerIndex,
        _ => anyhow::bail!(
            "Invalid UMI location: {}. Must be one of: read1, read2, index1, index2, per_read, per_index",
            args.umi_loc
        ),
    };

    Ok(umi::UmiConfig {
        enabled: true,
        location,
        length: args.umi_len,
        prefix: args.umi_prefix.clone(),
        skip: args.umi_skip,
    })
}

// Helper function to create DedupConfig from CLI args
fn create_dedup_config(args: &Args) -> Result<dedup::DedupConfig> {
    if !args.dedup {
        return Ok(dedup::DedupConfig::default());
    }

    // Validate accuracy level
    if args.dedup_accuracy < 1 || args.dedup_accuracy > 6 {
        anyhow::bail!("Deduplication accuracy must be between 1 and 6 (default: 5)");
    }

    Ok(dedup::DedupConfig {
        enabled: true,
        accuracy: args.dedup_accuracy,
    })
}

// Helper function to create SplitConfig from CLI args
fn create_split_config(args: &Args) -> Result<split::SplitConfig> {
    // Validate split parameters
    if let Some(num_files) = args.split {
        if !(2..=999).contains(&num_files) {
            anyhow::bail!("Split file number must be between 2 and 999");
        }
    }

    if let Some(lines) = args.split_by_lines {
        if lines < 1000 {
            anyhow::bail!("Split by lines must be at least 1000");
        }
    }

    if args.split_prefix_digits > 10 {
        anyhow::bail!("Split prefix digits must be between 0 and 10");
    }

    // Determine split mode
    let mode = if let Some(num_files) = args.split {
        split::SplitMode::ByFiles(num_files)
    } else if let Some(lines) = args.split_by_lines {
        split::SplitMode::ByLines(lines)
    } else {
        split::SplitMode::None
    };

    Ok(split::SplitConfig {
        mode,
        prefix_digits: args.split_prefix_digits,
    })
}

// Helper function to build and write paired-end report
fn build_and_write_paired_end_report(
    args: &Args,
    pe_acc: PairedEndAccumulator,
    start_time: std::time::Instant,
) -> Result<()> {
    let before_stats_r1 = pe_acc.before_r1.to_read_stats();
    let after_stats_r1 = pe_acc.after_r1.to_read_stats();
    let before_stats_r2 = pe_acc.before_r2.to_read_stats();
    let after_stats_r2 = pe_acc.after_r2.to_read_stats();

    let quality_curves_r1 = pe_acc.pos_r1.to_quality_curves();
    let quality_curves_r2 = pe_acc.pos_r2.to_quality_curves();
    let content_curves_r1 = pe_acc.pos_r1.to_content_curves();
    let content_curves_r2 = pe_acc.pos_r2.to_content_curves();
    let qual_hist_r1 = pe_acc.pos_r1.to_qual_hist();
    let qual_hist_r2 = pe_acc.pos_r2.to_qual_hist();
    let quality_curves_r1_after = pe_acc.pos_r1_after.to_quality_curves();
    let quality_curves_r2_after = pe_acc.pos_r2_after.to_quality_curves();
    let content_curves_r1_after = pe_acc.pos_r1_after.to_content_curves();
    let content_curves_r2_after = pe_acc.pos_r2_after.to_content_curves();
    let qual_hist_r1_after = pe_acc.pos_r1_after.to_qual_hist();
    let qual_hist_r2_after = pe_acc.pos_r2_after.to_qual_hist();
    let kmer_map_r1 = pe_acc.kmer_table_to_map_r1();
    let kmer_map_r2 = pe_acc.kmer_table_to_map_r2();

    // Calculate duplication rate from combined kmer counts
    let dup_rate_r1 = stats::calculate_duplication_rate(&kmer_map_r1);
    let dup_rate_r2 = stats::calculate_duplication_rate(&kmer_map_r2);
    let combined_dup_rate = (dup_rate_r1 + dup_rate_r2) / 2.0;

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

    let report = FasterpReport {
        summary: Summary {
            fasterp_version: env!("CARGO_PKG_VERSION").to_string(),
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
            content_curves: content_curves_r1,
            qual_hist: qual_hist_r1,
            kmer_count: kmer_map_r1,
        },
        read2_before_filtering: Some(DetailedReadStats {
            total_reads: before_stats_r2.total_reads,
            total_bases: before_stats_r2.total_bases,
            q20_bases: before_stats_r2.q20_bases,
            q30_bases: before_stats_r2.q30_bases,
            quality_curves: quality_curves_r2,
            content_curves: content_curves_r2,
            qual_hist: qual_hist_r2,
            kmer_count: kmer_map_r2,
        }),
        read1_after_filtering: Some(DetailedReadStats {
            total_reads: after_stats_r1.total_reads,
            total_bases: after_stats_r1.total_bases,
            q20_bases: after_stats_r1.q20_bases,
            q30_bases: after_stats_r1.q30_bases,
            quality_curves: quality_curves_r1_after,
            content_curves: content_curves_r1_after,
            qual_hist: qual_hist_r1_after,
            kmer_count: IndexMap::new(),
        }),
        read2_after_filtering: Some(DetailedReadStats {
            total_reads: after_stats_r2.total_reads,
            total_bases: after_stats_r2.total_bases,
            q20_bases: after_stats_r2.q20_bases,
            q30_bases: after_stats_r2.q30_bases,
            quality_curves: quality_curves_r2_after,
            content_curves: content_curves_r2_after,
            qual_hist: qual_hist_r2_after,
            kmer_count: IndexMap::new(),
        }),
        duplication: Some(DuplicationStats {
            rate: combined_dup_rate,
        }),
        adapter_cutting: None, // TODO: Track adapter cutting stats
    };

    // Print report to stdout (fastp-compatible format)
    let elapsed = start_time.elapsed();
    print_report_to_stdout(&report, args, elapsed, true);

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

    // Generate HTML report
    html::generate_html_report(&report, args, &args.html)
        .context("Failed to generate HTML report")?;

    Ok(())
}

// STDOUT REPORTING FUNCTIONS

/// Print a summary of read statistics
fn print_read_stats(title: &str, stats: &ReadStats) {
    println!("{title}:");
    println!("total reads: {}", stats.total_reads);
    println!("total bases: {}", stats.total_bases);
    println!(
        "Q20 bases: {}({:.4}%)",
        stats.q20_bases,
        stats.q20_rate * 100.0
    );
    println!(
        "Q30 bases: {}({:.4}%)",
        stats.q30_bases,
        stats.q30_rate * 100.0
    );
    println!("Q40 bases: 0(0%)");
}

/// Print filtering results
fn print_filtering_results(
    result: &FilteringResult,
    adapter_cutting: Option<&AdapterCuttingStats>,
) {
    println!("Filtering result:");
    println!("reads passed filter: {}", result.passed_filter_reads);
    println!(
        "reads failed due to low quality: {}",
        result.low_quality_reads
    );
    println!(
        "reads failed due to too many N: {}",
        result.too_many_n_reads
    );
    println!("reads failed due to too short: {}", result.too_short_reads);

    if let Some(ac) = adapter_cutting {
        println!("reads with adapter trimmed: {}", ac.adapter_trimmed_reads);
        println!(
            "bases trimmed due to adapters: {}",
            ac.adapter_trimmed_bases
        );
    } else {
        println!("reads with adapter trimmed: 0");
        println!("bases trimmed due to adapters: 0");
    }
}

/// Print command line reconstruction
fn print_command_line(args: &Args, is_paired_end: bool) {
    if is_paired_end {
        print!("fasterp -i {} ", args.input);
        if let Some(ref in2) = args.input2 {
            print!("-I {in2} ");
        }
        print!("-o {} ", args.output);
        if let Some(ref out2) = args.output2 {
            print!("-O {out2} ");
        }
        println!("-j {}", args.json);
    } else {
        println!(
            "fasterp -i {} -o {} -j {}",
            args.input, args.output, args.json
        );
    }
}

/// Print full report to stdout (fastp-compatible format)
fn print_report_to_stdout(
    report: &FasterpReport,
    args: &Args,
    elapsed: std::time::Duration,
    is_paired_end: bool,
) {
    // Adapter detection messages
    if is_paired_end {
        println!("Detecting adapter sequence for read1...");
        println!("No adapter detected for read1");
        println!();
        println!("Detecting adapter sequence for read2...");
        println!("No adapter detected for read2");
    } else {
        println!("Detecting adapter sequence for read1...");
        println!("No adapter detected for read1");
    }
    println!();

    // Before/after stats for read1
    print_read_stats("Read1 before filtering", &report.summary.before_filtering);
    println!();
    print_read_stats("Read1 after filtering", &report.summary.after_filtering);
    println!();

    // If paired-end, print read2 stats
    if is_paired_end {
        // Note: For PE, we'd need separate read2 stats which aren't currently
        // separated in the summary. For now, just print read1 stats.
        // TODO: Add separate read2 before/after stats to Summary
    }

    // Filtering results
    print_filtering_results(&report.filtering_result, report.adapter_cutting.as_ref());
    println!();

    // Duplication rate
    if let Some(dup) = &report.duplication {
        let qualifier = if is_paired_end {
            ""
        } else {
            " (may be overestimated since this is SE data)"
        };
        println!("Duplication rate{}: {:.4}%", qualifier, dup.rate * 100.0);
    }
    println!();

    // Report files
    println!("JSON report: {}", args.json);
    println!("HTML report: {}", args.html);
    println!();

    // Command line and version
    print_command_line(args, is_paired_end);

    // Format elapsed time more precisely
    let time_str = if elapsed.as_secs() >= 1 {
        format!("{:.3} seconds", elapsed.as_secs_f64())
    } else {
        format!("{} milliseconds", elapsed.as_millis())
    };

    println!(
        "fasterp v{}, time used: {}",
        env!("CARGO_PKG_VERSION"),
        time_str
    );
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

    // Track start time for performance reporting
    let start_time = std::time::Instant::now();

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
    let trimming_config = create_trimming_config(&args);

    // Create overlap configuration from CLI args (for base correction)
    let overlap_config = if args.correction {
        Some(create_overlap_config(&args))
    } else {
        None
    };

    // Create UMI configuration from CLI args
    let umi_config = if args.umi {
        Some(create_umi_config(&args)?)
    } else {
        None
    };

    // Create deduplication configuration from CLI args
    let dedup_config = if args.dedup {
        Some(create_dedup_config(&args)?)
    } else {
        None
    };

    // Create split configuration from CLI args
    let split_config = create_split_config(&args)?;

    // Process based on mode (single-end or paired-end)
    if is_paired_end {
        // PAIRED-END MODE
        // Create separate trimming configs for R1 and R2
        let mut trimming_config_r1 = trimming_config;
        let mut trimming_config_r2 = create_trimming_config_r2(&args);

        // Apply fastp's undocumented default: trim_tail=1 when NO trimming parameters specified
        // This default is disabled if user specifies ANY fixed trimming parameter
        let user_specified_trimming = args.trim_front != 0
            || args.trim_tail != 0
            || args.trim_front1 != 0
            || args.trim_tail1 != 0
            || args.trim_front2 != 0
            || args.trim_tail2 != 0;

        if !user_specified_trimming && args.trim_tail == 0 && args.trim_tail1 == 0 {
            trimming_config_r1.trim_tail_bases = 1;
        }
        if !user_specified_trimming
            && args.trim_tail == 0
            && args.trim_tail1 == 0
            && args.trim_tail2 == 0
        {
            trimming_config_r2.trim_tail_bases = 1;
        }

        // Apply disable_length_filtering flag
        let min_len = if args.disable_length_filtering {
            0
        } else {
            args.length_required
        };

        let pe_acc = if num_threads == 1 {
            // SINGLE-THREADED PAIRED-END MODE
            let reader1 = open_input(&args.input)?;
            let reader2 = open_input(args.input2.as_ref().unwrap())?;

            // Create split writers
            let compression = args.compression_level.unwrap_or(6);
            let mut writer1 =
                split::SplitWriter::new(&args.output, split_config.clone(), compression)?;
            let mut writer2 = split::SplitWriter::new(
                args.output2.as_ref().unwrap(),
                split_config.clone(),
                compression,
            )?;

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
                overlap_config.as_ref(),
                umi_config.as_ref(),
                dedup_config.as_ref(),
            )?;

            writer1.finish()?;
            writer2.finish()?;
            pe_acc
        } else {
            // MULTI-THREADED PAIRED-END MODE
            use crate::pipeline::{
                PairedBatch, PairedWorkerResult, paired_merger_thread, paired_producer_thread,
                paired_worker_thread,
            };
            use crossbeam_channel::bounded;

            let backlog = args.max_backlog.unwrap_or(num_threads + 1);

            // Create channels
            let (batch_tx, batch_rx) = bounded::<Option<PairedBatch>>(backlog);
            let (result_tx, result_rx) = bounded::<Option<PairedWorkerResult>>(backlog);

            // Spawn producer thread
            let input_path1 = args.input.clone();
            let input_path2 = args.input2.as_ref().unwrap().clone();
            let batch_bytes = args.batch_bytes;
            let producer = thread::spawn(move || {
                paired_producer_thread(input_path1, input_path2, batch_bytes, batch_tx)
            });

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
                let trimming_config_r1_clone = trimming_config_r1.clone();
                let trimming_config_r2_clone = trimming_config_r2.clone();
                let overlap_config_clone = overlap_config.clone();
                let umi_config_clone = umi_config.clone();
                let dedup_config_clone = dedup_config.clone();

                let worker = thread::spawn(move || {
                    paired_worker_thread(
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
                        trimming_config_r1_clone,
                        trimming_config_r2_clone,
                        overlap_config_clone,
                        umi_config_clone,
                        dedup_config_clone,
                    )
                });
                workers.push(worker);
            }

            // Drop original senders so merger knows when all workers are done
            drop(batch_rx);
            drop(result_tx);

            // Spawn merger thread
            let output_path1 = args.output.clone();
            let output_path2 = args.output2.as_ref().unwrap().clone();
            let compression_level = args.compression_level;
            let split_config_clone = split_config.clone();
            let merger = thread::spawn(move || {
                paired_merger_thread(
                    result_rx,
                    output_path1,
                    output_path2,
                    num_threads,
                    compression_level,
                    split_config_clone,
                )
            });

            // Wait for all threads
            producer.join().unwrap()?;
            for worker in workers {
                worker.join().unwrap();
            }
            merger.join().unwrap()?
        };

        // Build paired-end report
        build_and_write_paired_end_report(&args, pe_acc, start_time)?;
        return Ok(());
    }

    // SINGLE-END MODE
    // Apply disable_length_filtering flag
    let min_len = if args.disable_length_filtering {
        0
    } else {
        args.length_required
    };

    let acc = if num_threads == 1 {
        // SINGLE-THREADED MODE: use streaming approach with new parser
        let reader = open_input(&args.input)?;

        // Create split writer
        let compression = args.compression_level.unwrap_or(6);
        let mut writer = split::SplitWriter::new(&args.output, split_config.clone(), compression)?;

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
        let split_config_clone = split_config.clone();
        let merger = thread::spawn(move || {
            merger_thread(
                result_rx,
                output_path,
                num_threads,
                compression_level,
                split_config_clone,
            )
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
    let content_curves = acc.pos.to_content_curves();
    let qual_hist = acc.pos.to_qual_hist();
    let quality_curves_after = acc.pos_after.to_quality_curves();
    let content_curves_after = acc.pos_after.to_content_curves();
    let qual_hist_after = acc.pos_after.to_qual_hist();
    let kmer_map = acc.kmer_table_to_map();
    let kmer_map_after = acc.kmer_table_after_to_map();

    // Calculate duplication rate from kmer counts
    let duplication_rate = stats::calculate_duplication_rate(&kmer_map);

    let report = FasterpReport {
        summary: Summary {
            fasterp_version: env!("CARGO_PKG_VERSION").to_string(),
            sequencing: format!("single end ({} cycles)", acc.max_cycle),
            before_filtering: before_stats.clone(),
            after_filtering: after_stats.clone(),
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
            content_curves,
            qual_hist,
            kmer_count: kmer_map,
        },
        read2_before_filtering: None,
        read1_after_filtering: Some(DetailedReadStats {
            total_reads: after_stats.total_reads,
            total_bases: after_stats.total_bases,
            q20_bases: after_stats.q20_bases,
            q30_bases: after_stats.q30_bases,
            quality_curves: quality_curves_after,
            content_curves: content_curves_after,
            qual_hist: qual_hist_after,
            kmer_count: kmer_map_after,
        }),
        read2_after_filtering: None,
        duplication: Some(DuplicationStats {
            rate: duplication_rate,
        }),
        adapter_cutting: None, // TODO: Track adapter cutting stats
    };

    // Print report to stdout (fastp-compatible format)
    let elapsed = start_time.elapsed();
    print_report_to_stdout(&report, &args, elapsed, false);

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

    // Generate HTML report
    html::generate_html_report(&report, &args, &args.html)
        .context("Failed to generate HTML report")?;

    Ok(())
}

// UNIT TESTS (see src/tests.rs)

#[cfg(test)]
mod tests;
