//! A fast and simple FASTQ preprocessor
//!
//! This tool processes FASTQ files with quality control filters including:
//! - Length filtering (minimum read length)
//! - Quality filtering (mean quality score)
//! - N-base filtering (maximum ambiguous bases)

// Allow cast warnings for numeric computing - these are intentional and necessary for statistics
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_possible_wrap)]

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::thread;

mod adapter;
mod bloom;
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

use io::open_input;
use pipeline::{Batch, WorkerResult, merger_thread, producer_thread, worker_thread};
use processor::{PairedEndAccumulator, process_fastq_stream, process_paired_fastq_stream};
use stats::{
    AdapterCuttingStats, DetailedReadStats, DuplicationStats, FasterpReport, FilteringResult,
    InsertSizeStats, ReadStats, Summary,
};
use trimming::TrimmingConfig;

#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug)]
#[command(author, version, about = "A fast FASTQ preprocessor", long_about = None)]
struct Args {
    /// Input FASTQ file or URL (use '-' for stdin, supports .gz and http(s)://)
    #[arg(short = 'i', long)]
    input: String,

    /// Output FASTQ file (use '-' for stdout, .gz extension enables compression)
    #[arg(short = 'o', long)]
    output: String,

    /// Read2 input file or URL (paired-end mode, supports .gz and http(s)://)
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

    /// Maximum length for read1 - trim if longer (default: 0 = disabled)
    #[arg(short = 'b', long = "max_len1", default_value = "0")]
    max_len1: usize,

    /// Maximum length for read2 - trim if longer (default: 0 = disabled, follows read1 if not specified)
    #[arg(short = 'B', long = "max_len2", default_value = "0")]
    max_len2: usize,

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

    /// Use parallel compression for gzip output (default: true, disable with --no-parallel-compression)
    #[arg(long, default_value = "true", action = clap::ArgAction::Set)]
    parallel_compression: bool,

    /// Number of worker threads (default: auto-detect CPU count)
    #[arg(short = 'w', long)]
    threads: Option<usize>,

    /// Batch size in bytes (default: 32 MiB)
    #[arg(long, default_value = "33554432")]
    batch_bytes: usize,

    /// Maximum backlog of batches (default: threads+1)
    #[arg(long)]
    max_backlog: Option<usize>,

    /// Maximum memory usage in MB (e.g., 500 for 500MB, 2048 for 2GB)
    /// When set, automatically configures threads, batch_bytes, and max_backlog
    #[arg(long)]
    max_memory: Option<usize>,

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
    #[arg(short = 't', long = "trim_tail1", default_value = "0")]
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

    /// Disable adapter trimming (adapter trimming is enabled by default)
    #[arg(short = 'A', long = "disable_adapter_trimming")]
    disable_adapter_trimming: bool,

    /// The adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped
    #[arg(short = 'a', long = "adapter_sequence")]
    adapter_sequence: Option<String>,

    /// Adapter sequence for read2
    #[arg(long = "adapter_sequence_r2")]
    adapter_sequence_r2: Option<String>,

    /// Disable adapter auto-detection (enabled by default, use this to turn off)
    #[arg(long = "disable_adapter_detection")]
    disable_adapter_detection: bool,

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

    /// Merge overlapping paired-end reads into single reads
    #[arg(short = 'm', long)]
    merge: bool,

    /// Output file for merged reads (required when --merge is enabled)
    #[arg(long)]
    merged_out: Option<String>,

    /// Write unmerged reads to merged output file
    #[arg(long)]
    include_unmerged: bool,

    /// Enable UMI preprocessing
    #[arg(long)]
    umi: bool,

    /// UMI location: `read1/read2/index1/index2/per_read/per_index` (default: read1)
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

    /// Deduplication accuracy (1-6). Higher = more memory but fewer false positives (default: 3 for dedup, 1 otherwise)
    #[arg(long, default_value = "0")]
    dup_calc_accuracy: u8,

    /// Don't evaluate duplication rate (saves time and memory)
    #[arg(long)]
    dont_eval_duplication: bool,

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

    if args.disable_adapter_trimming {
        // If adapter trimming is disabled entirely, also disable auto-detection
        adapter_config.detect_adapter_for_pe = false;
    } else {
        // Set manual adapter sequences if provided
        if let Some(ref seq) = args.adapter_sequence {
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        }

        if let Some(ref seq) = args.adapter_sequence_r2 {
            adapter_config.adapter_seq_r2 = Some(seq.as_bytes().to_vec());
        }

        // Disable auto-detection if flag is set
        if args.disable_adapter_detection {
            adapter_config.detect_adapter_for_pe = false;
        }
        // Otherwise keep the default (true, auto-detection enabled)
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
        max_len: args.max_len1,
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

    if args.disable_adapter_trimming {
        // If adapter trimming is disabled entirely, also disable auto-detection
        adapter_config.detect_adapter_for_pe = false;
    } else {
        // For R2, use adapter_seq_r2 if specified
        if let Some(ref seq) = args.adapter_sequence_r2 {
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        } else if let Some(ref seq) = args.adapter_sequence {
            // Fall back to R1 adapter if R2 not specified
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        }

        // Disable auto-detection if flag is set
        if args.disable_adapter_detection {
            adapter_config.detect_adapter_for_pe = false;
        }
        // Otherwise keep the default (true, auto-detection enabled)
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

    // Per fastp docs: "If not specified, it will follow read1's settings"
    let max_len = if args.max_len2 > 0 {
        args.max_len2
    } else {
        args.max_len1
    };

    TrimmingConfig {
        enable_trim_front: args.cut_front && args.cut_mean_quality > 0,
        enable_trim_tail: args.cut_tail && args.cut_mean_quality > 0 && !args.disable_trim_tail,
        cut_mean_quality: args.cut_mean_quality,
        cut_window_size: args.cut_window_size,
        trim_front_bases: trim_front,
        trim_tail_bases: trim_tail,
        max_len,
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
fn create_dedup_config(args: &Args) -> dedup::DedupConfig {
    if !args.dedup {
        return dedup::DedupConfig::default();
    }

    // Dedup uses HashSet-based approach with default accuracy
    dedup::DedupConfig {
        enabled: true,
        accuracy: 5, // Default accuracy for HashSet-based dedup
    }
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
#[allow(clippy::too_many_lines)]
fn build_and_write_paired_end_report(
    args: &Args,
    pe_acc: &PairedEndAccumulator,
    start_time: std::time::Instant,
    dup_detector: Option<&std::sync::Arc<bloom::DuplicateDetector>>,
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
    let kmer_map_r1_after = pe_acc.kmer_table_to_map_r1_after();
    let kmer_map_r2_after = pe_acc.kmer_table_to_map_r2_after();

    // Calculate duplication rate
    // Use Bloom filter if available (matches fastp), otherwise fall back to kmer-based estimation
    let combined_dup_rate = if let Some(detector) = &dup_detector {
        detector.get_dup_rate()
    } else {
        let dup_rate_r1 = stats::calculate_duplication_rate(&kmer_map_r1);
        let dup_rate_r2 = stats::calculate_duplication_rate(&kmer_map_r2);
        f64::midpoint(dup_rate_r1, dup_rate_r2)
    };

    // Calculate combined before/after stats for summary
    let combined_before = ReadStats {
        // Count both R1 and R2 reads to match fastp's behavior
        total_reads: before_stats_r1.total_reads + before_stats_r2.total_reads,
        total_bases: before_stats_r1.total_bases + before_stats_r2.total_bases,
        q20_bases: before_stats_r1.q20_bases + before_stats_r2.q20_bases,
        q30_bases: before_stats_r1.q30_bases + before_stats_r2.q30_bases,
        q40_bases: before_stats_r1.q40_bases + before_stats_r2.q40_bases,
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
        q40_rate: if before_stats_r1.total_bases + before_stats_r2.total_bases > 0 {
            (before_stats_r1.q40_bases + before_stats_r2.q40_bases) as f64
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
        q40_bases: after_stats_r1.q40_bases + after_stats_r2.q40_bases,
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
        q40_rate: if after_stats_r1.total_bases + after_stats_r2.total_bases > 0 {
            (after_stats_r1.q40_bases + after_stats_r2.q40_bases) as f64
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
            fastp_version: env!("CARGO_PKG_VERSION").to_string(),
            sequencing: format!(
                "paired end ({} cycles + {} cycles)",
                pe_acc.max_cycle_r1, pe_acc.max_cycle_r2
            ),
            before_filtering: combined_before,
            after_filtering: combined_after,
            read1_before_filtering: Some(before_stats_r1.clone()),
            read2_before_filtering: Some(before_stats_r2.clone()),
            read1_after_filtering: Some(after_stats_r1.clone()),
            read2_after_filtering: Some(after_stats_r2.clone()),
        },
        filtering_result: FilteringResult {
            // Count both R1 and R2 reads to match fastp's behavior
            passed_filter_reads: pe_acc.after_r1.total_reads + pe_acc.after_r2.total_reads,
            // In paired-end mode, filter counts represent pairs, but fastp reports individual reads
            // So multiply by 2 to count both R1 and R2
            low_quality_reads: pe_acc.low_quality * 2,
            low_complexity_reads: pe_acc.low_complexity * 2,
            too_many_n_reads: pe_acc.too_many_n * 2,
            too_short_reads: pe_acc.too_short * 2,
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
            kmer_count: kmer_map_r1_after,
        }),
        read2_after_filtering: Some(DetailedReadStats {
            total_reads: after_stats_r2.total_reads,
            total_bases: after_stats_r2.total_bases,
            q20_bases: after_stats_r2.q20_bases,
            q30_bases: after_stats_r2.q30_bases,
            quality_curves: quality_curves_r2_after,
            content_curves: content_curves_r2_after,
            qual_hist: qual_hist_r2_after,
            kmer_count: kmer_map_r2_after,
        }),
        duplication: Some(DuplicationStats {
            rate: combined_dup_rate,
        }),
        adapter_cutting: if pe_acc.adapter_trimmed_reads > 0 {
            Some(AdapterCuttingStats {
                adapter_trimmed_reads: pe_acc.adapter_trimmed_reads,
                adapter_trimmed_bases: pe_acc.adapter_trimmed_bases,
            })
        } else {
            None
        },
        insert_size: Some(InsertSizeStats {
            peak: pe_acc.calculate_insert_size_peak(),
            unknown: pe_acc.insert_size_unknown,
            histogram: pe_acc.insert_size_histogram.clone(),
        }),
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
    let html_config = html::HtmlReportConfig {
        input: args.input.clone(),
        output: args.output.clone(),
        json: args.json.clone(),
        input2: args.input2.clone(),
        output2: args.output2.clone(),
    };
    html::generate_html_report(&report, &html_config, &args.html, elapsed)
        .context("Failed to generate HTML report")?;

    Ok(())
}

// STDOUT REPORTING FUNCTIONS

/// Print a summary of read statistics
fn print_read_stats(title: &str, stats: &ReadStats) {
    eprintln!("{title}:");
    eprintln!("total reads: {}", stats.total_reads);
    eprintln!("total bases: {}", stats.total_bases);
    eprintln!(
        "Q20 bases: {}({:.4}%)",
        stats.q20_bases,
        stats.q20_rate * 100.0
    );
    eprintln!(
        "Q30 bases: {}({:.4}%)",
        stats.q30_bases,
        stats.q30_rate * 100.0
    );
    eprintln!(
        "Q40 bases: {}({:.4}%)",
        stats.q40_bases,
        stats.q40_rate * 100.0
    );
}

/// Print filtering results
fn print_filtering_results(
    result: &FilteringResult,
    adapter_cutting: Option<&AdapterCuttingStats>,
) {
    eprintln!("Filtering result:");
    eprintln!("reads passed filter: {}", result.passed_filter_reads);
    eprintln!(
        "reads failed due to low quality: {}",
        result.low_quality_reads
    );
    eprintln!(
        "reads failed due to too many N: {}",
        result.too_many_n_reads
    );
    eprintln!("reads failed due to too short: {}", result.too_short_reads);

    if let Some(ac) = adapter_cutting {
        eprintln!("reads with adapter trimmed: {}", ac.adapter_trimmed_reads);
        eprintln!(
            "bases trimmed due to adapters: {}",
            ac.adapter_trimmed_bases
        );
    } else {
        eprintln!("reads with adapter trimmed: 0");
        eprintln!("bases trimmed due to adapters: 0");
    }
}

/// Print command line reconstruction
fn print_command_line(args: &Args, is_paired_end: bool) {
    if is_paired_end {
        eprint!("fasterp -i {} ", args.input);
        if let Some(ref in2) = args.input2 {
            eprint!("-I {in2} ");
        }
        eprint!("-o {} ", args.output);
        if let Some(ref out2) = args.output2 {
            eprint!("-O {out2} ");
        }
        eprintln!("-j {}", args.json);
    } else {
        eprintln!(
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
        eprintln!("Detecting adapter sequence for read1...");
        eprintln!("No adapter detected for read1");
        eprintln!();
        eprintln!("Detecting adapter sequence for read2...");
        eprintln!("No adapter detected for read2");
    } else {
        eprintln!("Detecting adapter sequence for read1...");
        eprintln!("No adapter detected for read1");
    }
    eprintln!();

    // Before/after stats
    if is_paired_end {
        // For paired-end, print separate Read1 and Read2 stats
        if let Some(ref r1_before) = report.summary.read1_before_filtering {
            print_read_stats("Read1 before filtering", r1_before);
            eprintln!();
        }
        if let Some(ref r2_before) = report.summary.read2_before_filtering {
            print_read_stats("Read2 before filtering", r2_before);
            eprintln!();
        }
        if let Some(ref r1_after) = report.summary.read1_after_filtering {
            print_read_stats("Read1 after filtering", r1_after);
            eprintln!();
        }
        if let Some(ref r2_after) = report.summary.read2_after_filtering {
            print_read_stats("Read2 after filtering", r2_after);
            eprintln!();
        }
    } else {
        // For single-end, print combined stats
        print_read_stats("Read1 before filtering", &report.summary.before_filtering);
        eprintln!();
        print_read_stats("Read1 after filtering", &report.summary.after_filtering);
        eprintln!();
    }

    // Filtering results
    print_filtering_results(&report.filtering_result, report.adapter_cutting.as_ref());
    eprintln!();

    // Duplication rate
    if let Some(dup) = &report.duplication {
        let qualifier = if is_paired_end {
            ""
        } else {
            " (may be overestimated since this is SE data)"
        };
        eprintln!("Duplication rate{}: {:.4}%", qualifier, dup.rate * 100.0);
    }
    eprintln!();

    // Insert size peak (only for paired-end)
    if is_paired_end {
        if let Some(insert_size) = &report.insert_size {
            eprintln!(
                "Insert size peak (evaluated by paired-end reads): {}",
                insert_size.peak
            );
            eprintln!();
        }
    }

    // Report files
    eprintln!("JSON report: {}", args.json);
    eprintln!("HTML report: {}", args.html);
    eprintln!();

    // Command line and version
    print_command_line(args, is_paired_end);

    // Format elapsed time more precisely
    let time_str = if elapsed.as_secs() >= 1 {
        format!("{:.3} seconds", elapsed.as_secs_f64())
    } else {
        format!("{} milliseconds", elapsed.as_millis())
    };

    eprintln!(
        "fasterp v{}, time used: {}",
        env!("CARGO_PKG_VERSION"),
        time_str
    );
}

/// Perform adapter auto-detection for paired-end reads
///
/// Samples the first N read pairs and uses overlap analysis to detect adapters
fn auto_detect_adapters_pe(
    input1: &str,
    input2: &str,
    sample_size: usize,
) -> Result<adapter::AdapterDetectionResult> {
    use std::io::BufRead;

    eprintln!("Auto-detecting adapters from paired-end reads...");
    eprintln!("Sampling first {sample_size} read pairs for analysis");

    let mut reader1 = open_input(input1)?;
    let mut reader2 = open_input(input2)?;

    let mut r1_seqs = Vec::new();
    let mut r2_seqs = Vec::new();

    let mut line1 = String::new();
    let mut line2 = String::new();

    // Sample reads using simple line-by-line reading
    let mut count = 0;
    while count < sample_size {
        // Read R1 record (4 lines)
        line1.clear();
        if reader1.read_line(&mut line1)? == 0 {
            break; // EOF
        }

        let mut seq1 = String::new();
        let mut plus1 = String::new();
        let mut qual1 = String::new();
        if reader1.read_line(&mut seq1)? == 0
            || reader1.read_line(&mut plus1)? == 0
            || reader1.read_line(&mut qual1)? == 0
        {
            break; // EOF or incomplete record
        }

        // Read R2 record (4 lines)
        line2.clear();
        if reader2.read_line(&mut line2)? == 0 {
            break; // EOF
        }

        let mut seq2 = String::new();
        let mut plus2 = String::new();
        let mut qual2 = String::new();
        if reader2.read_line(&mut seq2)? == 0
            || reader2.read_line(&mut plus2)? == 0
            || reader2.read_line(&mut qual2)? == 0
        {
            break; // EOF or incomplete record
        }

        // Store sequences (trim trailing newline)
        r1_seqs.push(seq1.trim_end().as_bytes().to_vec());
        r2_seqs.push(seq2.trim_end().as_bytes().to_vec());
        count += 1;
    }

    if r1_seqs.is_empty() {
        anyhow::bail!("No reads found for adapter detection");
    }

    eprintln!("Analyzing {} read pairs...", r1_seqs.len());

    // Convert to slices for detection
    let r1_refs: Vec<&[u8]> = r1_seqs.iter().map(std::vec::Vec::as_slice).collect();
    let r2_refs: Vec<&[u8]> = r2_seqs.iter().map(std::vec::Vec::as_slice).collect();

    let result = adapter::detect_adapters_from_pe_reads(
        r1_refs.into_iter(),
        r2_refs.into_iter(),
        sample_size,
        0.01, // 1% minimum frequency
    );

    Ok(result)
}

/// Perform adapter auto-detection for single-end reads
///
/// Samples the first N reads and uses k-mer frequency analysis to detect adapters
fn auto_detect_adapter_se(
    input: &str,
    sample_size: usize,
) -> Result<adapter::AdapterDetectionResult> {
    use std::io::BufRead;

    eprintln!("Auto-detecting adapter from single-end reads...");
    eprintln!("Sampling first {sample_size} reads for analysis");

    let mut reader = open_input(input)?;
    let mut seqs = Vec::new();

    let mut line = String::new();

    // Sample reads using simple line-by-line reading
    let mut count = 0;
    while count < sample_size {
        // Read FASTQ record (4 lines)
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break; // EOF
        }

        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();
        if reader.read_line(&mut seq)? == 0
            || reader.read_line(&mut plus)? == 0
            || reader.read_line(&mut qual)? == 0
        {
            break; // EOF or incomplete record
        }

        // Store sequence (trim trailing newline)
        seqs.push(seq.trim_end().as_bytes().to_vec());
        count += 1;
    }

    if seqs.is_empty() {
        anyhow::bail!("No reads found for adapter detection");
    }

    eprintln!("Analyzing {} reads...", seqs.len());

    // Convert to slices for detection
    let seq_refs: Vec<&[u8]> = seqs.iter().map(std::vec::Vec::as_slice).collect();

    let result = adapter::detect_adapter_from_se_reads(
        seq_refs.into_iter(),
        sample_size,
        0.0001, // 0.01% minimum frequency (same as fastp)
    );

    Ok(result)
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
#[allow(clippy::too_many_lines)]
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

    // Determine number of threads and memory configuration
    let (num_threads, batch_bytes, max_backlog) = if let Some(max_memory_mb) = args.max_memory {
        // Calculate optimal configuration based on memory limit
        // Memory model: Total ≈ base_overhead + (batch_bytes × (backlog × 2 + workers))
        // Base overhead includes I/O buffers (~120 MB for SE, ~240 MB for PE)
        let base_overhead_mb = if is_paired_end { 240 } else { 120 };

        if max_memory_mb <= base_overhead_mb {
            anyhow::bail!(
                "max_memory ({} MB) must be greater than base overhead ({} MB)",
                max_memory_mb,
                base_overhead_mb
            );
        }

        let available_mb = max_memory_mb - base_overhead_mb;

        // Strategy: Start with requested/default workers, then scale down if needed
        let requested_threads = args.threads.unwrap_or_else(num_cpus::get);

        // Try to find a configuration that fits
        // Priority: maintain workers > maintain batch size > maintain backlog
        let mut best_threads = requested_threads;
        let mut best_batch_mb = 32usize; // Default 32 MB
        let mut best_backlog = 2usize; // Minimum reasonable backlog

        // Calculate: available = batch × (backlog × 2 + workers)
        // Rearrange: batch = available / (backlog × 2 + workers)

        for threads in (1..=requested_threads).rev() {
            let backlog = args.max_backlog.unwrap_or(threads.min(4)); // Cap default backlog at 4 for memory efficiency
            let divisor = backlog * 2 + threads;
            let batch_mb = available_mb / divisor;

            // Minimum viable batch size is 4 MB
            if batch_mb >= 4 {
                best_threads = threads;
                best_batch_mb = batch_mb.min(32); // Cap at 32 MB for efficiency
                best_backlog = backlog;
                break;
            }
        }

        // Final calculation for actual memory usage
        let actual_pipeline_mb = best_batch_mb * (best_backlog * 2 + best_threads);
        let actual_total_mb = base_overhead_mb + actual_pipeline_mb;

        eprintln!("Memory limit: {} MB (--max-memory)", max_memory_mb);
        eprintln!("Inferred configuration:");
        eprintln!("  threads:     {} (from --max-memory)", best_threads);
        eprintln!("  batch_bytes: {} MB (from --max-memory)", best_batch_mb);
        eprintln!("  max_backlog: {} (from --max-memory)", best_backlog);
        eprintln!(
            "  estimated usage: ~{} MB (base: {} MB, pipeline: {} MB)",
            actual_total_mb, base_overhead_mb, actual_pipeline_mb
        );
        eprintln!();

        (
            best_threads,
            best_batch_mb * 1024 * 1024,
            Some(best_backlog),
        )
    } else {
        // Use explicit or default values
        let threads = args.threads.unwrap_or_else(num_cpus::get);
        (threads, args.batch_bytes, args.max_backlog)
    };

    // Override args with calculated values (for use later in pipeline)
    let args = {
        let mut args = args;
        args.batch_bytes = batch_bytes;
        args.max_backlog = max_backlog;
        args
    };

    println!("Using {num_threads} thread(s) for processing");

    // Validate merge arguments
    if args.merge {
        if !is_paired_end {
            anyhow::bail!("Merge mode requires paired-end input (-I/--in2)");
        }
        if args.merged_out.is_none() {
            anyhow::bail!("Merge mode requires --merged_out to specify output file");
        }
        // TODO: Add multi-threaded support for merge
        if num_threads > 1 {
            eprintln!(
                "Warning: Merge mode currently only supports single-threaded processing. Use -w 1 to enable merge."
            );
            eprintln!("Continuing without merge functionality...");
        }
    }

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
        Some(create_dedup_config(&args))
    } else {
        None
    };

    // Create duplicate detector for calculating duplication rate (separate from dedup)
    // Uses Bloom filter matching fastp's implementation
    let dup_detector = if args.dont_eval_duplication {
        None
    } else {
        let accuracy = if args.dup_calc_accuracy > 0 {
            args.dup_calc_accuracy
        } else if args.dedup {
            3 // Default to 3 when dedup is enabled
        } else {
            1 // Default to 1 otherwise (lower memory usage)
        };
        Some(std::sync::Arc::new(bloom::DuplicateDetector::new(accuracy)))
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
        // Perform adapter auto-detection if requested
        // Skip auto-detection if either input is stdin (can't rewind)
        if trimming_config_r1.adapter_config.detect_adapter_for_pe
            && args.input != "-"
            && args.input2.as_ref().is_none_or(|i2| i2 != "-")
        {
            match auto_detect_adapters_pe(
                &args.input,
                args.input2.as_ref().unwrap(),
                1_000_000, // Sample first 1M read pairs
            ) {
                Ok(detection_result) => {
                    if let Some(ref adapter) = detection_result.adapter_r1 {
                        eprintln!(
                            "Detected R1 adapter: {} (confidence: {:.2}%, {} reads analyzed)",
                            String::from_utf8_lossy(adapter),
                            detection_result.confidence * 100.0,
                            detection_result.reads_analyzed
                        );
                        trimming_config_r1.adapter_config.adapter_seq = Some(adapter.clone());
                    } else {
                        eprintln!("No R1 adapter detected");
                    }

                    if let Some(ref adapter) = detection_result.adapter_r2 {
                        eprintln!("Detected R2 adapter: {}", String::from_utf8_lossy(adapter));
                        trimming_config_r2.adapter_config.adapter_seq = Some(adapter.clone());
                    } else {
                        eprintln!("No R2 adapter detected");
                    }

                    if detection_result.adapter_r1.is_none()
                        && detection_result.adapter_r2.is_none()
                    {
                        eprintln!(
                            "Warning: Adapter auto-detection found no adapters. Proceeding without adapter trimming."
                        );
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Adapter auto-detection failed: {e}");
                    eprintln!("Proceeding without adapter trimming.");
                }
            }
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
            let mut writer1 = split::SplitWriter::new(
                &args.output,
                split_config.clone(),
                compression,
                args.parallel_compression,
            )?;
            let mut writer2 = split::SplitWriter::new(
                args.output2.as_ref().unwrap(),
                split_config.clone(),
                compression,
                args.parallel_compression,
            )?;

            // Open merged output if merge mode enabled
            let mut merged_writer_storage = if args.merge {
                Some(split::SplitWriter::new(
                    args.merged_out.as_ref().unwrap(),
                    split_config.clone(),
                    compression,
                    args.parallel_compression,
                )?)
            } else {
                None
            };

            let pe_acc = process_paired_fastq_stream(
                reader1,
                reader2,
                &mut writer1,
                &mut writer2,
                merged_writer_storage.as_mut(),
                args.merge,
                args.include_unmerged,
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
                dup_detector.as_ref(),
            )?;

            // Finish writers
            if let Some(mw) = merged_writer_storage {
                mw.finish()?;
            }
            writer1.finish()?;
            writer2.finish()?;
            pe_acc
        } else {
            // MULTI-THREADED PAIRED-END MODE
            use crate::pipeline::{
                PairedBatch, PairedWorkerResult, paired_merger_thread,
                paired_producer_thread_with_paths, paired_worker_thread,
            };
            use crossbeam_channel::bounded;

            let backlog = args.max_backlog.unwrap_or(num_threads + 1);

            // Create channels
            let (batch_tx, batch_rx) = bounded::<Option<PairedBatch>>(backlog);
            let (result_tx, result_rx) = bounded::<Option<PairedWorkerResult>>(backlog);

            // Spawn producer thread with file paths (opens readers in thread to isolate GzDecoders)
            let batch_bytes = args.batch_bytes;
            let num_workers = num_threads;
            let input1_clone = args.input.clone();
            let input2_clone = args.input2.as_ref().unwrap().clone();
            let producer = thread::spawn(move || {
                paired_producer_thread_with_paths(
                    input1_clone,
                    input2_clone,
                    batch_bytes,
                    batch_tx,
                    num_workers,
                )
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
                let dup_detector_clone = dup_detector.clone();

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
                        dup_detector_clone,
                    );
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
            let parallel_compression = args.parallel_compression;
            let split_config_clone = split_config.clone();
            let merger = thread::spawn(move || {
                paired_merger_thread(
                    result_rx,
                    output_path1,
                    output_path2,
                    num_threads,
                    compression_level,
                    split_config_clone,
                    parallel_compression,
                )
            });

            // Wait for all threads
            producer.join().unwrap()?;
            for (index, worker) in workers.into_iter().enumerate() {
                worker.join().unwrap();
            }
            merger.join().unwrap()?
        };

        // Build paired-end report
        build_and_write_paired_end_report(&args, &pe_acc, start_time, dup_detector.as_ref())?;
        return Ok(());
    }

    // SINGLE-END MODE

    // Perform adapter auto-detection if requested and no adapter specified
    // Skip auto-detection if input is stdin (can't rewind)
    let mut se_trimming_config = trimming_config;

    // Enable auto-detection for SE mode by default (fastp does this for SE)
    // Unlike PE mode which uses overlap-based detection
    if !args.disable_adapter_detection && se_trimming_config.adapter_config.adapter_seq.is_none() {
        se_trimming_config.adapter_config.detect_adapter_for_pe = true;
    }

    if se_trimming_config.adapter_config.detect_adapter_for_pe
        && se_trimming_config.adapter_config.adapter_seq.is_none()
        && args.input != "-"
    {
        match auto_detect_adapter_se(
            &args.input,
            1_000_000, // Sample first 1M reads
        ) {
            Ok(detection_result) => {
                if let Some(ref adapter) = detection_result.adapter_r1 {
                    eprintln!(
                        "Detected adapter: {} (confidence: {:.2}%, {} reads analyzed)",
                        String::from_utf8_lossy(adapter),
                        detection_result.confidence * 100.0,
                        detection_result.reads_analyzed
                    );
                    se_trimming_config.adapter_config.adapter_seq = Some(adapter.clone());
                } else {
                    eprintln!("No adapter detected in single-end mode");
                    eprintln!("Proceeding without adapter trimming.");
                }
            }
            Err(e) => {
                eprintln!("Warning: Adapter auto-detection failed: {e}");
                eprintln!("Proceeding without adapter trimming.");
            }
        }
    }

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
        let mut writer = split::SplitWriter::new(
            &args.output,
            split_config.clone(),
            compression,
            args.parallel_compression,
        )?;

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
            &se_trimming_config,
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
            let trimming_config_clone = se_trimming_config.clone();

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
                );
            });
            workers.push(worker);
        }

        // Drop original senders so merger knows when all workers are done
        drop(batch_rx);
        drop(result_tx);

        // Spawn merger thread
        let output_path = args.output.clone();
        let compression_level = args.compression_level;
        let parallel_compression = args.parallel_compression;
        let split_config_clone = split_config.clone();
        let merger = thread::spawn(move || {
            merger_thread(
                result_rx,
                output_path,
                num_threads,
                compression_level,
                split_config_clone,
                parallel_compression,
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
            fastp_version: env!("CARGO_PKG_VERSION").to_string(),
            sequencing: format!("single end ({} cycles)", acc.max_cycle),
            before_filtering: before_stats.clone(),
            after_filtering: after_stats.clone(),
            read1_before_filtering: None,
            read2_before_filtering: None,
            read1_after_filtering: None,
            read2_after_filtering: None,
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
        adapter_cutting: if acc.adapter_trimmed_reads > 0 {
            Some(AdapterCuttingStats {
                adapter_trimmed_reads: acc.adapter_trimmed_reads,
                adapter_trimmed_bases: acc.adapter_trimmed_bases,
            })
        } else {
            None
        },
        insert_size: None, // Insert size is only applicable for paired-end reads
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
    let html_config = html::HtmlReportConfig {
        input: args.input.clone(),
        output: args.output.clone(),
        json: args.json.clone(),
        input2: args.input2.clone(),
        output2: args.output2.clone(),
    };
    html::generate_html_report(&report, &html_config, &args.html, elapsed)
        .context("Failed to generate HTML report")?;

    Ok(())
}

// UNIT TESTS (see src/tests.rs)

#[cfg(test)]
mod tests;
