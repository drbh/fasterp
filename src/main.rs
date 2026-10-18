//! Native command-line interface. Compatibility target: OpenGene/fastp v1.4.0.
use anyhow::Result;
use clap::Parser;

mod adapter;
mod auto_adapter;
mod compat;
mod gzip;
mod html;
mod io;
mod kmer;
mod overlap;
mod simd;
mod split;
mod stats;

#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "A fast FASTQ preprocessor (fastp 1.4.0 interface)",
    rename_all = "snake_case",
    disable_help_flag = true,
    disable_version_flag = true,
    args_override_self = true
)]
struct Args {
    /// Input FASTQ file or URL (use '-' for stdin, supports .gz and http(s)://)
    #[arg(short = 'i', long = "in1", alias = "input", default_value = "")]
    input: String,

    /// Output FASTQ file (use '-' for stdout, .gz extension enables compression)
    #[arg(short = 'o', long = "out1", alias = "output", default_value = "")]
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
    #[arg(long, alias = "interleaved-in")]
    interleaved_in: bool,

    /// Minimum length required (default: 15)
    #[arg(short = 'l', long, default_value = "15", alias = "length-required")]
    length_required: usize,

    /// Disable length filtering (fastp compatibility)
    #[arg(short = 'L', long, alias = "disable-length-filtering")]
    disable_length_filtering: bool,

    /// Maximum length for read1 - trim if longer (default: 0 = disabled)
    #[arg(
        short = 'b',
        long = "max_len1",
        default_value = "0",
        alias = "max-len1"
    )]
    max_len1: usize,

    /// Maximum length for read2 - trim if longer (default: 0 = disabled, follows read1 if not specified)
    #[arg(short = 'B', long = "max_len2", alias = "max-len2")]
    max_len2: Option<usize>,

    /// Quality value that a base is qualified (phred quality >= this value) (default: 15)
    #[arg(
        short = 'q',
        long,
        default_value = "15",
        alias = "qualified-quality-phred"
    )]
    qualified_quality_phred: u8,

    /// Percent of bases allowed to be unqualified (0-100) (default: 40)
    #[arg(
        short = 'u',
        long,
        default_value = "40",
        alias = "unqualified-percent-limit"
    )]
    unqualified_percent_limit: usize,

    /// Average quality threshold - discard if mean quality < this (default: 0 = disabled)
    #[arg(short = 'e', long, default_value = "0", alias = "average-qual")]
    average_qual: u8,

    /// Max number of N bases allowed (default: 5)
    #[arg(short = 'n', long, default_value = "5", alias = "n-base-limit")]
    n_base_limit: usize,

    /// Enable low complexity filter
    #[arg(short = 'y', long, alias = "low-complexity-filter")]
    low_complexity_filter: bool,

    /// Complexity threshold (0-100). Default 30 means 30% complexity required (fastp compatibility)
    #[arg(
        short = 'Y',
        long,
        default_value = "30",
        alias = "complexity-threshold"
    )]
    complexity_threshold: usize,

    /// JSON report file (default: fasterp.json)
    #[arg(short = 'j', long, default_value = "fasterp.json")]
    json: String,

    /// HTML report file (default: fasterp.html)
    #[arg(short = 'h', long, default_value = "fasterp.html")]
    html: String,

    /// Stats output format: compact (default), pretty, off, or jsonl
    #[arg(long, default_value = "compact", alias = "stats-format")]
    stats_format: String,

    /// Compression level for gzip output (1-9, default: 4)
    #[arg(short = 'z', long = "compression", alias = "compression_level", alias = "compression-level", value_parser = clap::value_parser!(u32).range(1..=9))]
    compression_level: Option<u32>,

    /// Use parallel compression for gzip output (disable with --parallel_compression false)
    #[arg(long, default_value = "true", action = clap::ArgAction::Set, alias = "parallel-compression")]
    parallel_compression: bool,

    /// Number of worker threads (default: 3)
    #[arg(short = 'w', long = "thread", alias = "threads", default_value = "3")]
    threads: Option<usize>,

    /// Batch size in bytes (default: 32 MiB)
    #[arg(long, default_value = "33554432", alias = "batch-bytes")]
    batch_bytes: usize,

    /// Maximum backlog of batches (default: threads+1)
    #[arg(long, alias = "max-backlog")]
    max_backlog: Option<usize>,

    /// Processing batch memory budget in MB (excludes profiling and adapter samples)
    #[arg(long, alias = "max-memory")]
    max_memory: Option<usize>,

    /// Skip k-mer counting for ceiling performance tests
    #[arg(long, alias = "no-kmer")]
    no_kmer: bool,

    /// Quality cutoff for sliding-window trimming (default: 20)
    /// Trim when mean quality in window falls below this value
    #[arg(short = 'M', long, default_value = "20", alias = "cut-mean-quality")]
    cut_mean_quality: u8,

    /// Sliding window size for quality trimming (default: 4)
    #[arg(short = 'W', long, default_value = "4", alias = "cut-window-size")]
    cut_window_size: usize,

    /// Enable quality trimming at 5' end (front)
    #[arg(short = '5', long, alias = "cut-front", alias = "cut_by_quality5")]
    cut_front: bool,

    /// Enable quality trimming at 3' end (tail)
    #[arg(short = '3', long, alias = "cut-tail", alias = "cut_by_quality3")]
    cut_tail: bool,

    /// Disable polyG/quality tail trimming (for backward compatibility)
    #[arg(long, alias = "disable-trim-tail")]
    disable_trim_tail: bool,

    /// Trim N bases from 5' (front) end
    #[arg(long, default_value = "0", alias = "trim-front")]
    trim_front: usize,

    /// Trim N bases from 3' (tail) end
    #[arg(long, default_value = "0", alias = "trim-tail")]
    trim_tail: usize,

    /// Trim N bases from front for read1 (fastp compatibility)
    #[arg(
        short = 'f',
        long = "trim_front1",
        default_value = "0",
        alias = "trim-front1"
    )]
    trim_front1: usize,

    /// Trim N bases from tail for read1 (fastp compatibility)
    #[arg(
        short = 't',
        long = "trim_tail1",
        default_value = "0",
        alias = "trim-tail1"
    )]
    trim_tail1: usize,

    /// Trim N bases from front for read2 (fastp compatibility)
    #[arg(short = 'F', long = "trim_front2", alias = "trim-front2")]
    trim_front2: Option<usize>,

    /// Trim N bases from tail for read2 (fastp compatibility)
    #[arg(short = 'T', long = "trim_tail2", alias = "trim-tail2")]
    trim_tail2: Option<usize>,

    /// Enable polyG tail trimming
    #[arg(short = 'g', long, alias = "trim-poly-g")]
    trim_poly_g: bool,

    /// Disable polyG tail trimming (for backward compatibility)
    #[arg(short = 'G', long, alias = "disable-trim-poly-g")]
    disable_trim_poly_g: bool,

    /// Enable generic polyX tail trimming (any homopolymer)
    #[arg(short = 'x', long, alias = "trim-poly-x")]
    trim_poly_x: bool,

    /// Minimum length for polyG/polyX detection (default: 10)
    #[arg(long, default_value = "10", alias = "poly-g-min-len")]
    poly_g_min_len: usize,

    /// Disable adapter trimming (adapter trimming is enabled by default)
    #[arg(
        short = 'A',
        long = "disable_adapter_trimming",
        alias = "disable-adapter-trimming"
    )]
    disable_adapter_trimming: bool,

    /// The adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped
    #[arg(short = 'a', long = "adapter_sequence", alias = "adapter-sequence")]
    adapter_sequence: Option<String>,

    /// Adapter sequence for read2
    #[arg(long = "adapter_sequence_r2", alias = "adapter-sequence-r2")]
    adapter_sequence_r2: Option<String>,

    /// Disable adapter auto-detection (enabled by default, use this to turn off)
    #[arg(
        long = "disable_adapter_detection",
        alias = "disable-adapter-detection"
    )]
    disable_adapter_detection: bool,

    /// Enable base correction using overlap analysis for paired-end data
    #[arg(short = 'c', long)]
    correction: bool,

    /// Minimum overlap length required for correction (default: 30)
    #[arg(long, default_value = "30", alias = "overlap-len-require")]
    overlap_len_require: usize,

    /// Maximum allowed differences in overlap region (default: 5)
    #[arg(long, default_value = "5", alias = "overlap-diff-limit")]
    overlap_diff_limit: usize,

    /// Maximum allowed difference percentage (0-100) (default: 20)
    #[arg(long, default_value = "20", alias = "overlap-diff-percent-limit")]
    overlap_diff_percent_limit: usize,

    /// Merge overlapping paired-end reads into single reads
    #[arg(short = 'm', long)]
    merge: bool,

    /// Output file for merged reads (required when --merge is enabled)
    #[arg(long, alias = "merged-out")]
    merged_out: Option<String>,

    /// Write unmerged reads to merged output file
    #[arg(long, alias = "include-unmerged")]
    include_unmerged: bool,

    /// Enable UMI preprocessing
    #[arg(short = 'U', long)]
    umi: bool,

    /// UMI location: `read1/read2/index1/index2/per_read/per_index` (default: read1)
    #[arg(long, default_value = "", alias = "umi-loc")]
    umi_loc: String,

    /// UMI length (required for extraction from reads)
    #[arg(long, default_value = "0", alias = "umi-len")]
    umi_len: usize,

    /// Prefix in read name for UMI (default: none)
    #[arg(long, default_value = "", alias = "umi-prefix")]
    umi_prefix: String,

    /// Skip N bases following UMI
    #[arg(long, default_value = "0", alias = "umi-skip")]
    umi_skip: usize,

    /// Enable deduplication to remove duplicate reads (fastp compatible)
    #[arg(short = 'D', long)]
    dedup: bool,

    /// Deduplication accuracy (1-6). Higher = more memory but fewer false positives (default: 3 for dedup, 1 otherwise)
    #[arg(long, default_value = "0", alias = "dup-calc-accuracy")]
    dup_calc_accuracy: u8,

    /// Don't evaluate duplication rate (saves time and memory)
    #[arg(long, alias = "dont-eval-duplication")]
    dont_eval_duplication: bool,

    /// Split output by limiting total number of files (2-999). Cannot be used with --split-by-lines
    #[arg(short = 's', long, conflicts_with = "split_by_lines")]
    split: Option<usize>,

    /// Split output by limiting lines per file (>=1000). Cannot be used with --split
    #[arg(short = 'S', long, conflicts_with = "split", alias = "split-by-lines")]
    split_by_lines: Option<usize>,

    /// Digits for file number padding (1-10, default: 4). 0 to disable padding
    #[arg(short = 'd', long, default_value = "4", alias = "split-prefix-digits")]
    split_prefix_digits: usize,

    /// Show help (--help; -h writes the HTML report, as in fastp).
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
    /// Print version.
    #[arg(short = 'v', long, action = clap::ArgAction::Version)]
    version: Option<bool>,
    /// Read FASTQ from standard input.
    #[arg(long)]
    stdin: bool,
    /// Write passing reads to standard output (interleaved for paired-end).
    #[arg(long)]
    stdout: bool,
    /// Convert Phred+64 input to Phred+33 output.
    #[arg(short = '6', long)]
    phred64: bool,
    /// Stop after this many reads or pairs (0 = all).
    #[arg(long, alias = "reads-to-process", default_value_t = 0)]
    reads_to_process: usize,
    /// Refuse to overwrite existing output or report files.
    #[arg(long, alias = "dont-overwrite")]
    dont_overwrite: bool,
    /// Fix MGI read identifiers ending in /1 or /2.
    #[arg(long, alias = "fix-mgi-id")]
    fix_mgi_id: bool,
    /// Print progress information.
    #[arg(short = 'V', long)]
    verbose: bool,
    /// Disable quality, average quality, and N-base filters.
    #[arg(short = 'Q', long, alias = "disable-quality-filtering")]
    disable_quality_filtering: bool,
    /// Discard reads longer than this (0 = disabled).
    #[arg(long, alias = "length-limit", default_value_t = 0)]
    length_limit: usize,
    /// Write failing reads with failure tags.
    #[arg(long, alias = "failed-out")]
    failed_out: Option<String>,
    /// Write perfectly matching paired overlap regions.
    #[arg(long, alias = "overlapped-out")]
    overlapped_out: Option<String>,
    /// Enable adapter auto-detection for paired-end input.
    #[arg(short = '2', long, alias = "detect-adapter-for-pe")]
    detect_adapter_for_pe: bool,
    /// Trim adapters from every sequence in a FASTA file.
    #[arg(long, alias = "adapter-fasta")]
    adapter_fasta: Option<String>,
    /// Allow a single gap in overlap analysis used for adapter trimming.
    #[arg(long, alias = "allow-gap-overlap-trimming")]
    allow_gap_overlap_trimming: bool,
    /// Adapter dimer length threshold; requires evidence of an adapter.
    #[arg(long, alias = "dimer-max-len", default_value_t = 2)]
    dimer_max_len: usize,
    /// Separate minimum length for polyX trimming.
    #[arg(long, alias = "poly-x-min-len", default_value_t = 10)]
    poly_x_min_len: usize,
    /// Cut from the first low quality window to the end.
    #[arg(
        short = 'r',
        long,
        alias = "cut-right",
        alias = "cut_by_quality_aggressive"
    )]
    cut_right: bool,
    #[arg(long, alias = "cut-front-window-size")]
    cut_front_window_size: Option<usize>,
    #[arg(long, alias = "cut-front-mean-quality")]
    cut_front_mean_quality: Option<u8>,
    #[arg(long, alias = "cut-tail-window-size")]
    cut_tail_window_size: Option<usize>,
    #[arg(long, alias = "cut-tail-mean-quality")]
    cut_tail_mean_quality: Option<u8>,
    #[arg(long, alias = "cut-right-window-size")]
    cut_right_window_size: Option<usize>,
    #[arg(long, alias = "cut-right-mean-quality")]
    cut_right_mean_quality: Option<u8>,
    /// Delimiter between the read name and UMI.
    #[arg(long, alias = "umi-delim", default_value = ":")]
    umi_delim: String,
    #[arg(long, alias = "filter-by-index1")]
    filter_by_index1: Option<String>,
    #[arg(long, alias = "filter-by-index2")]
    filter_by_index2: Option<String>,
    #[arg(long, alias = "filter-by-index-threshold", default_value_t = 0)]
    filter_by_index_threshold: usize,
    /// Enable overrepresented sequence analysis.
    #[arg(short = 'p', long, alias = "overrepresentation-analysis")]
    overrepresentation_analysis: bool,
    #[arg(
        short = 'P',
        long,
        alias = "overrepresentation-sampling",
        default_value_t = 20
    )]
    overrepresentation_sampling: usize,
    /// Report title.
    #[arg(
        short = 'R',
        long,
        alias = "report-title",
        default_value = "fasterp report"
    )]
    report_title: String,
    /// Deprecated fastp option; has no effect.
    #[arg(long, alias = "discard-unmerged", hide = true)]
    discard_unmerged: bool,
}

fn main() -> Result<()> {
    compat::run(Args::parse())
}
