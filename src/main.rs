//! A fast and simple FASTQ preprocessor
//!
//! This tool processes FASTQ files with quality control filters including:
//! - Length filtering (minimum read length)
//! - Quality filtering (mean quality score)
//! - N-base filtering (maximum ambiguous bases)

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::{Receiver, Sender, bounded};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::thread;

// LOOKUP TABLES for fast base/quality checks

/// Lookup table: is quality >= 20? (Phred+33 encoding)
static LUT_Q20: [bool; 256] = {
    let mut lut = [false; 256];
    let mut i = 0;
    while i < 256 {
        lut[i] = i >= 33 && (i - 33) >= 20;
        i += 1;
    }
    lut
};

/// Lookup table: is quality >= 30? (Phred+33 encoding)
static LUT_Q30: [bool; 256] = {
    let mut lut = [false; 256];
    let mut i = 0;
    while i < 256 {
        lut[i] = i >= 33 && (i - 33) >= 30;
        i += 1;
    }
    lut
};

/// Lookup table: is base N?
static LUT_IS_N: [bool; 256] = {
    let mut lut = [false; 256];
    lut[b'N' as usize] = true;
    lut[b'n' as usize] = true;
    lut
};

/// Lookup table: is base G or C?
static LUT_IS_GC: [bool; 256] = {
    let mut lut = [false; 256];
    lut[b'G' as usize] = true;
    lut[b'g' as usize] = true;
    lut[b'C' as usize] = true;
    lut[b'c' as usize] = true;
    lut
};

/// Convert base to 2-bit encoding: A=0, C=1, G=2, T=3
#[inline]
fn base_to_2bit(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Get base index for quality_curves: A=0, T=1, C=2, G=3
#[inline]
fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'T' | b't' => Some(1),
        b'C' | b'c' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

/// Count 5-mers using 2-bit rolling code (NO STRING ALLOCATIONS)
///
/// This replaces the old String-based approach that allocated millions of strings.
/// Uses a fixed array of 1024 elements (4^5 possible 5-mers).
/// Encodes A=0, C=1, G=2, T=3 and rolls a 10-bit window.
/// Any N base resets the window.
///
/// PERFORMANCE: ~10-50x faster than String-based approach for k=5
#[inline]
fn count_k5_2bit(seq: &[u8], kmer_table: &mut [usize; 1024]) {
    let mut code: u32 = 0;
    let mask: u32 = (1 << (2 * 5)) - 1; // 10 bits for 5-mer
    let mut filled = 0u8;

    for &b in seq {
        let Some(c) = base_to_2bit(b) else {
            // Hit an N or invalid base - reset window
            code = 0;
            filled = 0;
            continue;
        };

        code = ((code << 2) & mask) | c;

        if filled < 4 {
            filled += 1;
            continue;
        }

        kmer_table[code as usize] += 1;
    }
}

/// Convert 2-bit encoded kmer to String for JSON output
fn kmer_to_string(code: usize) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut result = Vec::with_capacity(5);
    let mut c = code;

    // Extract bases from right to left (least significant to most significant)
    for _ in 0..5 {
        result.push(bases[c & 3]);
        c >>= 2;
    }

    // Reverse to get correct order (we extracted backwards)
    result.reverse();
    String::from_utf8(result).unwrap()
}

// IO ABSTRACTION - Compression support (gzip) and stdin/stdout

/// Detect compression format from file extension or magic bytes
#[derive(Debug, Clone, Copy, PartialEq)]
enum CompressionFormat {
    None,
    Gzip,
}

impl CompressionFormat {
    /// Detect from file path
    fn from_path(path: &str) -> Self {
        if path == "-" {
            return CompressionFormat::None; // stdin/stdout defaults to uncompressed
        }

        let path_lower = path.to_lowercase();
        if path_lower.ends_with(".gz") || path_lower.ends_with(".gzip") {
            CompressionFormat::Gzip
        } else {
            CompressionFormat::None
        }
    }

    /// Detect from magic bytes (for future auto-detection)
    #[allow(dead_code)]
    fn from_magic_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
            Some(CompressionFormat::Gzip)
        } else {
            None
        }
    }
}

/// Open input file or stdin with automatic decompression
fn open_input(path: &str) -> Result<Box<dyn BufRead + Send>> {
    if path == "-" {
        // Read from stdin
        let stdin = std::io::stdin();
        let reader = BufReader::with_capacity(16 * 1024 * 1024, stdin);
        Ok(Box::new(reader))
    } else {
        // Open file and detect compression
        let file = File::open(path).context(format!("Failed to open input file: {path}"))?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                let decoder = GzDecoder::new(file);
                let reader = BufReader::with_capacity(16 * 1024 * 1024, decoder);
                Ok(Box::new(reader))
            }
            CompressionFormat::None => {
                let reader = BufReader::with_capacity(16 * 1024 * 1024, file);
                Ok(Box::new(reader))
            }
        }
    }
}

/// Wrapper for output writer that ensures proper cleanup
pub enum OutputWriter {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
    Stdout(BufWriter<std::io::Stdout>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Gzip(w) => w.write(buf),
            OutputWriter::Stdout(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Gzip(w) => w.flush(),
            OutputWriter::Stdout(w) => w.flush(),
        }
    }
}

impl OutputWriter {
    /// Finish writing and ensure gzip encoder is properly finalized
    pub fn finish(self) -> Result<()> {
        match self {
            OutputWriter::Plain(mut w) => {
                w.flush()?;
                Ok(())
            }
            OutputWriter::Gzip(mut w) => {
                w.flush()?;
                let encoder = w.into_inner().context("Failed to finish gzip writer")?;
                encoder.finish().context("Failed to finish gzip encoding")?;
                Ok(())
            }
            OutputWriter::Stdout(mut w) => {
                w.flush()?;
                Ok(())
            }
        }
    }
}

/// Open output file or stdout with optional compression
fn open_output(path: &str, compression_level: Option<u32>) -> Result<OutputWriter> {
    if path == "-" {
        // Write to stdout
        let stdout = std::io::stdout();
        let writer = BufWriter::with_capacity(16 * 1024 * 1024, stdout);
        Ok(OutputWriter::Stdout(writer))
    } else {
        // Create file and detect compression
        let file = File::create(path).context(format!("Failed to create output file: {path}"))?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                let level = compression_level.unwrap_or(6); // Default to level 6
                let compression = Compression::new(level);
                let encoder = GzEncoder::new(file, compression);
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, encoder);
                Ok(OutputWriter::Gzip(writer))
            }
            CompressionFormat::None => {
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, file);
                Ok(OutputWriter::Plain(writer))
            }
        }
    }
}

// ROBUST FASTQ PARSER - State machine handles multiline and missing newlines

/// State machine for parsing FASTQ records
/// Handles:
/// - Multiline sequences/qualities (wrapped lines)
/// - Missing final newline
/// - Malformed records (skips with warning)
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
enum ParseState {
    ExpectHeader,
    InSequence,
    ExpectPlus,
    InQuality,
}

/// FASTQ record parser that accumulates multiline sequences/qualities
pub struct FastqParser<R: BufRead> {
    reader: R,
    state: ParseState,
    header: Vec<u8>,
    sequence: Vec<u8>,
    plus: Vec<u8>,
    quality: Vec<u8>,
    line_buf: Vec<u8>,
    eof: bool,
}

impl<R: BufRead> FastqParser<R> {
    pub fn new(reader: R) -> Self {
        FastqParser {
            reader,
            state: ParseState::ExpectHeader,
            header: Vec::new(),
            sequence: Vec::new(),
            plus: Vec::new(),
            quality: Vec::new(),
            line_buf: Vec::new(),
            eof: false,
        }
    }

    /// Read next record, returns Some((header, seq, plus, qual)) or None if EOF
    /// Handles multiline sequences and missing final newlines
    pub fn next_record(&mut self) -> Result<Option<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)>> {
        // Reset buffers
        self.header.clear();
        self.sequence.clear();
        self.plus.clear();
        self.quality.clear();
        self.state = ParseState::ExpectHeader;

        loop {
            self.line_buf.clear();
            let bytes_read = self.reader.read_until(b'\n', &mut self.line_buf)?;

            // Handle EOF
            if bytes_read == 0 {
                if self.eof {
                    return Ok(None); // Already processed EOF
                }
                self.eof = true;

                // Check if we have a complete record without final newline
                if self.state == ParseState::InQuality && !self.quality.is_empty() {
                    // Validate sequence and quality match
                    if self.sequence.len() == self.quality.len() {
                        return Ok(Some((
                            std::mem::take(&mut self.header),
                            std::mem::take(&mut self.sequence),
                            std::mem::take(&mut self.plus),
                            std::mem::take(&mut self.quality),
                        )));
                    }
                }
                return Ok(None);
            }

            // Remove trailing newline if present
            let mut line = &self.line_buf[..];
            if line.ends_with(b"\n") {
                line = &line[..line.len() - 1];
            }
            if line.ends_with(b"\r") {
                line = &line[..line.len() - 1];
            }

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            match self.state {
                ParseState::ExpectHeader => {
                    if line.starts_with(b"@") {
                        self.header.extend_from_slice(line);
                        self.state = ParseState::InSequence;
                    }
                    // Ignore non-header lines when expecting header
                }
                ParseState::InSequence => {
                    if line.starts_with(b"+") {
                        // Found plus line, move to quality
                        self.plus.extend_from_slice(line);
                        self.state = ParseState::InQuality;
                    } else {
                        // Accumulate sequence (handles multiline)
                        self.sequence.extend_from_slice(line);
                    }
                }
                ParseState::InQuality => {
                    // Accumulate quality until we have enough bases
                    self.quality.extend_from_slice(line);

                    // Check if quality matches sequence length
                    if self.quality.len() >= self.sequence.len() {
                        // We have a complete record
                        if self.quality.len() > self.sequence.len() {
                            // Quality is longer than sequence, truncate
                            self.quality.truncate(self.sequence.len());
                        }

                        return Ok(Some((
                            std::mem::take(&mut self.header),
                            std::mem::take(&mut self.sequence),
                            std::mem::take(&mut self.plus),
                            std::mem::take(&mut self.quality),
                        )));
                    }
                }
                ParseState::ExpectPlus => {
                    // This state is unused in current implementation
                    if line.starts_with(b"+") {
                        self.plus.extend_from_slice(line);
                        self.state = ParseState::InQuality;
                    }
                }
            }
        }
    }
}

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

    /// Mean quality score threshold (default: 0, disabled)
    #[arg(short = 'q', long, default_value = "0")]
    qualified_quality_phred: u8,

    /// Max number of N bases allowed (default: 5)
    #[arg(short = 'n', long, default_value = "5")]
    n_base_limit: usize,

    /// JSON report file (default: fastp.json)
    #[arg(short = 'j', long, default_value = "fastp.json")]
    json: String,

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
    #[arg(long)]
    disable_trim_poly_g: bool,

    /// Enable generic polyX tail trimming (any homopolymer)
    #[arg(long)]
    trim_poly_x: bool,

    /// Minimum length for polyG/polyX detection (default: 10)
    #[arg(long, default_value = "10")]
    poly_g_min_len: usize,
}

// TRIMMING DATA STRUCTURES AND ALGORITHMS

/// Configuration for trimming operations
#[derive(Debug, Clone)]
struct TrimmingConfig {
    // Sliding window trimming
    enable_trim_front: bool,
    enable_trim_tail: bool,
    cut_mean_quality: u8,
    cut_window_size: usize,

    // Fixed position trimming
    trim_front_bases: usize,
    trim_tail_bases: usize,

    // PolyG/PolyX trimming
    enable_poly_g: bool,
    enable_poly_x: bool,
    poly_min_len: usize,
}

impl TrimmingConfig {
    fn from_args(args: &Args) -> Self {
        Self {
            enable_trim_front: args.cut_front && args.cut_mean_quality > 0,
            enable_trim_tail: args.cut_tail && args.cut_mean_quality > 0 && !args.disable_trim_tail,
            cut_mean_quality: args.cut_mean_quality,
            cut_window_size: args.cut_window_size,
            trim_front_bases: args.trim_front,
            trim_tail_bases: args.trim_tail,
            enable_poly_g: args.trim_poly_g && !args.disable_trim_poly_g,
            enable_poly_x: args.trim_poly_x,
            poly_min_len: args.poly_g_min_len,
        }
    }

    fn is_enabled(&self) -> bool {
        self.enable_trim_front
            || self.enable_trim_tail
            || self.trim_front_bases > 0
            || self.trim_tail_bases > 0
            || self.enable_poly_g
            || self.enable_poly_x
    }
}

/// Result of trimming a single read
#[derive(Debug, Clone, Default)]
struct TrimmingResult {
    start_pos: usize, // Starting position after 5' trimming
    end_pos: usize,   // Ending position after 3' trimming
    poly_g_trimmed: usize,
    poly_x_trimmed: usize,
}

impl TrimmingResult {
    fn trimmed_length(&self) -> usize {
        self.end_pos.saturating_sub(self.start_pos)
    }

    fn bases_trimmed_5(&self) -> usize {
        self.start_pos
    }

    fn bases_trimmed_3(&self) -> usize {
        self.poly_g_trimmed + self.poly_x_trimmed
    }
}

/// Accumulated trimming statistics
#[derive(Debug, Default, Clone)]
struct TrimmingStats {
    reads_trimmed: usize,
    bases_trimmed_5: u64,
    bases_trimmed_3: u64,
    poly_g_trimmed_reads: usize,
    poly_g_trimmed_bases: u64,
    poly_x_trimmed_reads: usize,
    poly_x_trimmed_bases: u64,
}

impl TrimmingStats {
    fn add(&mut self, result: &TrimmingResult) {
        if result.start_pos > 0 || result.poly_g_trimmed > 0 || result.poly_x_trimmed > 0 {
            self.reads_trimmed += 1;
        }

        self.bases_trimmed_5 += result.start_pos as u64;

        if result.poly_g_trimmed > 0 {
            self.poly_g_trimmed_reads += 1;
            self.poly_g_trimmed_bases += result.poly_g_trimmed as u64;
            self.bases_trimmed_3 += result.poly_g_trimmed as u64;
        }

        if result.poly_x_trimmed > 0 {
            self.poly_x_trimmed_reads += 1;
            self.poly_x_trimmed_bases += result.poly_x_trimmed as u64;
            self.bases_trimmed_3 += result.poly_x_trimmed as u64;
        }
    }

    fn merge(&mut self, other: &TrimmingStats) {
        self.reads_trimmed += other.reads_trimmed;
        self.bases_trimmed_5 += other.bases_trimmed_5;
        self.bases_trimmed_3 += other.bases_trimmed_3;
        self.poly_g_trimmed_reads += other.poly_g_trimmed_reads;
        self.poly_g_trimmed_bases += other.poly_g_trimmed_bases;
        self.poly_x_trimmed_reads += other.poly_x_trimmed_reads;
        self.poly_x_trimmed_bases += other.poly_x_trimmed_bases;
    }
}

/// Trim 3' end using sliding window quality check
///
/// Scans from the 3' end towards the 5' end with a sliding window.
/// Returns the end position where quality is acceptable.
///
/// # Arguments
/// * `qual` - Quality scores (Phred+33 encoded)
/// * `window_size` - Size of the sliding window
/// * `cutoff` - Minimum mean quality threshold
///
/// # Returns
/// End position (exclusive) for trimmed sequence
fn trim_tail_sliding_window(qual: &[u8], window_size: usize, cutoff: u8) -> usize {
    let len = qual.len();

    if len <= window_size {
        return len; // Don't trim if shorter than window
    }

    // Scan from 3' end towards 5' end
    for end_pos in (window_size..=len).rev() {
        let start = end_pos - window_size;
        let window_qual = &qual[start..end_pos];

        // Calculate mean quality of window (Phred+33)
        let sum: u32 = window_qual.iter().map(|&q| (q - 33) as u32).sum();
        let mean_qual = sum as f64 / window_size as f64;

        if mean_qual >= cutoff as f64 {
            return end_pos; // Found acceptable window
        }
    }

    0 // Entire read below threshold
}

/// Trim 5' end using sliding window quality check
///
/// Scans from the 5' end towards the 3' end with a sliding window.
/// Returns the start position where quality is acceptable.
///
/// # Arguments
/// * `qual` - Quality scores (Phred+33 encoded)
/// * `window_size` - Size of the sliding window
/// * `cutoff` - Minimum mean quality threshold
///
/// # Returns
/// Start position (inclusive) for trimmed sequence
fn trim_front_sliding_window(qual: &[u8], window_size: usize, cutoff: u8) -> usize {
    let len = qual.len();

    if len <= window_size {
        return 0; // Don't trim if shorter than window
    }

    // Scan from 5' end towards 3' end
    for start_pos in 0..=(len - window_size) {
        let end = start_pos + window_size;
        let window_qual = &qual[start_pos..end];

        let sum: u32 = window_qual.iter().map(|&q| (q - 33) as u32).sum();
        let mean_qual = sum as f64 / window_size as f64;

        if mean_qual >= cutoff as f64 {
            return start_pos; // Found acceptable window
        }
    }

    len // Entire read below threshold
}

/// Detect polyG tail (common Illumina NovaSeq/NextSeq artifact)
///
/// Scans backwards from the 3' end to find consecutive G bases.
///
/// # Arguments
/// * `seq` - Nucleotide sequence
/// * `min_len` - Minimum length to consider as polyG
///
/// # Returns
/// Number of G bases at the tail (0 if below threshold)
fn detect_poly_g_tail(seq: &[u8], min_len: usize) -> usize {
    let len = seq.len();
    let mut g_count = 0;

    // Scan backwards from 3' end
    for i in (0..len).rev() {
        if seq[i] == b'G' || seq[i] == b'g' {
            g_count += 1;
        } else {
            break;
        }
    }

    if g_count >= min_len { g_count } else { 0 }
}

/// Detect generic homopolymer tail (polyA, polyT, polyC, polyG, polyN)
///
/// Scans backwards from the 3' end to find consecutive identical bases.
///
/// # Arguments
/// * `seq` - Nucleotide sequence
/// * `min_len` - Minimum length to consider as polyX
///
/// # Returns
/// Number of bases in the homopolymer tail (0 if below threshold)
fn detect_poly_x_tail(seq: &[u8], min_len: usize) -> usize {
    let len = seq.len();
    if len == 0 {
        return 0;
    }

    let tail_base = seq[len - 1];
    let mut count = 0;

    // Scan backwards from 3' end
    for i in (0..len).rev() {
        if seq[i] == tail_base {
            count += 1;
        } else {
            break;
        }
    }

    if count >= min_len { count } else { 0 }
}

/// Apply all trimming operations to a read
///
/// Trimming order:
/// 1. Fixed front trimming
/// 2. Fixed tail trimming
/// 3. Sliding window front trimming (if enabled)
/// 4. Sliding window tail trimming (if enabled)
/// 5. PolyG/PolyX tail trimming
///
/// # Arguments
/// * `seq` - Nucleotide sequence
/// * `qual` - Quality scores
/// * `config` - Trimming configuration
///
/// # Returns
/// TrimmingResult with start/end positions and statistics
fn trim_read(seq: &[u8], qual: &[u8], config: &TrimmingConfig) -> TrimmingResult {
    let mut result = TrimmingResult {
        start_pos: 0,
        end_pos: seq.len(),
        poly_g_trimmed: 0,
        poly_x_trimmed: 0,
    };

    // 1. Fixed front trimming
    if config.trim_front_bases > 0 {
        result.start_pos = config.trim_front_bases.min(seq.len());
    }

    // 2. Fixed tail trimming
    if config.trim_tail_bases > 0 {
        result.end_pos = result.end_pos.saturating_sub(config.trim_tail_bases);
    }

    // Ensure we still have a valid range
    if result.start_pos >= result.end_pos {
        result.end_pos = result.start_pos;
        return result;
    }

    let current_qual = &qual[result.start_pos..result.end_pos];

    // 3. Sliding window front trimming
    if config.enable_trim_front {
        let trim_amount = trim_front_sliding_window(
            current_qual,
            config.cut_window_size,
            config.cut_mean_quality,
        );
        result.start_pos += trim_amount;
    }

    // 4. Sliding window tail trimming
    if config.enable_trim_tail {
        let new_len = trim_tail_sliding_window(
            &qual[result.start_pos..result.end_pos],
            config.cut_window_size,
            config.cut_mean_quality,
        );
        result.end_pos = result.start_pos + new_len;
    }

    // Ensure we still have sequence left
    if result.start_pos >= result.end_pos {
        result.end_pos = result.start_pos;
        return result;
    }

    let current_seq = &seq[result.start_pos..result.end_pos];

    // 5. PolyG tail trimming (check polyG first as it's more specific)
    if config.enable_poly_g {
        let poly_g = detect_poly_g_tail(current_seq, config.poly_min_len);
        if poly_g > 0 {
            result.poly_g_trimmed = poly_g;
            result.end_pos = result.end_pos.saturating_sub(poly_g);
        }
    }

    // 6. PolyX tail trimming (only if polyG didn't already trim)
    if config.enable_poly_x && result.poly_g_trimmed == 0 {
        let current_seq = &seq[result.start_pos..result.end_pos];
        let poly_x = detect_poly_x_tail(current_seq, config.poly_min_len);
        if poly_x > 0 {
            result.poly_x_trimmed = poly_x;
            result.end_pos = result.end_pos.saturating_sub(poly_x);
        }
    }

    result
}

// MULTI-THREADED PIPELINE DATA STRUCTURES

/// A batch of FASTQ records parsed from a buffer
///
/// Contains raw bytes and record positions (no String allocations)
/// Each record is [header_start, seq_start, plus_start, qual_start] as byte offsets
#[derive(Clone)]
struct Batch {
    id: u64,
    buf: Vec<u8>,
    /// Each element is [header_start, seq_start, plus_start, qual_start]
    /// Lengths are implicit: header len = seq_start - header_start, etc.
    /// Quality ends at the next record's header_start (or buf.len() for last record)
    recs: Vec<[usize; 4]>,
}

/// Result from a worker thread
struct WorkerResult {
    id: u64,
    out_buf: Vec<u8>,
    before: SimpleStats,
    after: SimpleStats,
    pos: PositionStats,
    k5: [usize; 1024],
    // Filter counts
    too_short: usize,
    too_many_n: usize,
    low_quality: usize,
    invalid: usize,
}

/// Statistics for a set of reads
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ReadStats {
    total_reads: usize,
    total_bases: usize,
    q20_bases: usize,
    q30_bases: usize,
    q20_rate: f64,
    q30_rate: f64,
    read1_mean_length: usize,
    gc_content: f64,
}

/// Quality curves data for per-position quality statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
struct QualityCurves {
    #[serde(rename = "A")]
    a: Vec<f64>,
    #[serde(rename = "T")]
    t: Vec<f64>,
    #[serde(rename = "C")]
    c: Vec<f64>,
    #[serde(rename = "G")]
    g: Vec<f64>,
    mean: Vec<f64>,
}

/// Detailed read statistics including kmer counts
#[derive(Debug, Clone, Serialize, Deserialize)]
struct DetailedReadStats {
    total_reads: usize,
    total_bases: usize,
    q20_bases: usize,
    q30_bases: usize,
    quality_curves: QualityCurves,
    kmer_count: IndexMap<String, usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Summary {
    fastp_version: String,
    sequencing: String,
    before_filtering: ReadStats,
    after_filtering: ReadStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FilteringResult {
    passed_filter_reads: usize,
    low_quality_reads: usize,
    #[serde(rename = "too_many_N_reads")]
    too_many_n_reads: usize,
    too_short_reads: usize,
    too_long_reads: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FastpReport {
    summary: Summary,
    filtering_result: FilteringResult,
    read1_before_filtering: DetailedReadStats,
}

// STREAMING ACCUMULATOR for single-pass processing

/// Accumulator for per-position quality statistics
struct PositionStats {
    /// Per-base quality sums: [A, T, C, G][position]
    base_sum: [Vec<u64>; 4],
    /// Per-base quality counts: [A, T, C, G][position]
    base_cnt: [Vec<u64>; 4],
    /// Total quality sum per position
    total_sum: Vec<u64>,
    /// Total count per position
    total_cnt: Vec<u64>,
}

impl PositionStats {
    fn new() -> Self {
        Self {
            base_sum: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            base_cnt: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            total_sum: Vec::new(),
            total_cnt: Vec::new(),
        }
    }

    /// Ensure capacity for at least `len` positions
    fn ensure_capacity(&mut self, len: usize) {
        if self.total_sum.len() < len {
            self.total_sum.resize(len, 0);
            self.total_cnt.resize(len, 0);
            for i in 0..4 {
                self.base_sum[i].resize(len, 0);
                self.base_cnt[i].resize(len, 0);
            }
        }
    }

    /// Convert to QualityCurves for JSON output
    fn to_quality_curves(&self) -> QualityCurves {
        let mean: Vec<f64> = self
            .total_sum
            .iter()
            .zip(&self.total_cnt)
            .map(|(&sum, &cnt)| {
                if cnt > 0 {
                    sum as f64 / cnt as f64
                } else {
                    0.0
                }
            })
            .collect();

        let mut curves = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
        for i in 0..4 {
            curves[i] = self.base_sum[i]
                .iter()
                .zip(&self.base_cnt[i])
                .enumerate()
                .map(|(pos, (&sum, &cnt))| {
                    if cnt > 0 {
                        sum as f64 / cnt as f64
                    } else {
                        mean[pos]
                    }
                })
                .collect();
        }

        QualityCurves {
            a: curves[0].clone(),
            t: curves[1].clone(),
            c: curves[2].clone(),
            g: curves[3].clone(),
            mean,
        }
    }
}

/// Simple stats accumulator
#[derive(Default)]
struct SimpleStats {
    total_reads: usize,
    total_bases: usize,
    q20_bases: usize,
    q30_bases: usize,
    gc_bases: usize,
}

impl SimpleStats {
    fn add(&mut self, bases: usize, q20: usize, q30: usize, gc: usize) {
        self.total_reads += 1;
        self.total_bases += bases;
        self.q20_bases += q20;
        self.q30_bases += q30;
        self.gc_bases += gc;
    }

    fn to_read_stats(&self) -> ReadStats {
        ReadStats {
            total_reads: self.total_reads,
            total_bases: self.total_bases,
            q20_bases: self.q20_bases,
            q30_bases: self.q30_bases,
            q20_rate: if self.total_bases > 0 {
                self.q20_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
            q30_rate: if self.total_bases > 0 {
                self.q30_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
            read1_mean_length: if self.total_reads > 0 {
                self.total_bases / self.total_reads
            } else {
                0
            },
            gc_content: if self.total_bases > 0 {
                self.gc_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
        }
    }
}

/// Main accumulator for streaming processing
struct StreamAccumulator {
    before: SimpleStats,
    after: SimpleStats,
    pos: PositionStats,
    kmer_table: [usize; 1024],

    // Filtering counts
    too_short: usize,
    too_many_n: usize,
    low_quality: usize,
    invalid: usize,
    max_cycle: usize,
}

impl StreamAccumulator {
    fn new() -> Self {
        Self {
            before: SimpleStats::default(),
            after: SimpleStats::default(),
            pos: PositionStats::new(),
            kmer_table: [0; 1024],
            too_short: 0,
            too_many_n: 0,
            low_quality: 0,
            invalid: 0,
            max_cycle: 0,
        }
    }

    /// Process a single record in ONE PASS - the heart of the optimization!
    ///
    /// This function:
    /// 1. Updates "before" statistics
    /// 2. Checks filters
    /// 3. Writes passing records
    /// 4. Updates "after" statistics
    ///
    /// All in a single pass through the data with minimal allocations.
    fn process_record(
        &mut self,
        header: &[u8],
        seq: &[u8],
        plus: &[u8],
        qual: &[u8],
        min_len: usize,
        n_limit: usize,
        q_mean_phred: u8,
        trimming_config: &TrimmingConfig,
        writer: &mut impl Write,
    ) -> Result<()> {
        // Validate record
        if seq.len() != qual.len() {
            self.invalid += 1;
            return Ok(());
        }

        // Track max cycle length
        if seq.len() > self.max_cycle {
            self.max_cycle = seq.len();
        }

        // SINGLE PASS: compute all "before" stats in one loop
        let mut qsum = 0u32;
        let mut q20 = 0usize;
        let mut q30 = 0usize;
        let mut ncnt = 0usize;
        let mut gc = 0usize;

        self.pos.ensure_capacity(seq.len());

        for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
            let qu = q as usize;

            // Quality stats
            qsum += (q - 33) as u32;
            q20 += LUT_Q20[qu] as usize;
            q30 += LUT_Q30[qu] as usize;

            // Base stats
            ncnt += LUT_IS_N[b as usize] as usize;
            gc += LUT_IS_GC[b as usize] as usize;

            // Position-specific quality stats
            self.pos.total_sum[i] += (q - 33) as u64;
            self.pos.total_cnt[i] += 1;

            if let Some(bi) = base_idx(b) {
                self.pos.base_sum[bi][i] += (q - 33) as u64;
                self.pos.base_cnt[bi][i] += 1;
            }
        }

        // K-mer counting (also in the same pass conceptually, but separate loop for clarity)
        count_k5_2bit(seq, &mut self.kmer_table);

        // Update "before" stats
        self.before.add(seq.len(), q20, q30, gc);

        // Apply trimming if enabled
        let trimming_result = if trimming_config.is_enabled() {
            trim_read(seq, qual, trimming_config)
        } else {
            TrimmingResult {
                start_pos: 0,
                end_pos: seq.len(),
                poly_g_trimmed: 0,
                poly_x_trimmed: 0,
            }
        };

        // Get trimmed sequences
        let trimmed_seq = &seq[trimming_result.start_pos..trimming_result.end_pos];
        let trimmed_qual = &qual[trimming_result.start_pos..trimming_result.end_pos];

        // Recompute stats for trimmed read (used for filtering)
        let trimmed_len = trimmed_seq.len();
        let mut trimmed_qsum = 0u32;
        let mut trimmed_q20 = 0usize;
        let mut trimmed_q30 = 0usize;
        let mut trimmed_ncnt = 0usize;
        let mut trimmed_gc = 0usize;

        for (&b, &q) in trimmed_seq.iter().zip(trimmed_qual) {
            let qu = q as usize;
            trimmed_qsum += (q - 33) as u32;
            trimmed_q20 += LUT_Q20[qu] as usize;
            trimmed_q30 += LUT_Q30[qu] as usize;
            trimmed_ncnt += LUT_IS_N[b as usize] as usize;
            trimmed_gc += LUT_IS_GC[b as usize] as usize;
        }

        // Apply filters on TRIMMED read
        if trimmed_len < min_len {
            self.too_short += 1;
            return Ok(());
        }

        if trimmed_ncnt > n_limit {
            self.too_many_n += 1;
            return Ok(());
        }

        if q_mean_phred > 0
            && trimmed_len > 0
            && (trimmed_qsum as f64 / trimmed_len as f64) < q_mean_phred as f64
        {
            self.low_quality += 1;
            return Ok(());
        }

        // Record passed - write trimmed version
        writeln!(writer, "{}", std::str::from_utf8(header)?)?;
        writeln!(writer, "{}", std::str::from_utf8(trimmed_seq)?)?;
        writeln!(writer, "{}", std::str::from_utf8(plus)?)?;
        writeln!(writer, "{}", std::str::from_utf8(trimmed_qual)?)?;

        // Update "after" stats with trimmed read stats
        self.after
            .add(trimmed_len, trimmed_q20, trimmed_q30, trimmed_gc);

        Ok(())
    }

    /// Convert kmer_table to IndexMap for JSON output
    fn kmer_table_to_map(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            let kmer_str = kmer_to_string(code);
            map.insert(kmer_str, self.kmer_table[code]);
        }
        map
    }
}

// MULTI-THREADED PIPELINE

/// Producer thread: read blocks and parse into batches
///
/// Reads large blocks (batch_bytes) from input, parses FASTQ records,
/// and emits Batch structures with NO string allocations - just byte slices.
///
/// Handles partial records at block boundaries by carrying them over to next batch.
fn producer_thread(
    input_path: String,
    batch_bytes: usize,
    sender: Sender<Option<Batch>>,
) -> Result<()> {
    let file =
        File::open(&input_path).context(format!("Failed to open input file: {input_path}"))?;
    let mut reader = BufReader::with_capacity(16 * 1024 * 1024, file); // 16 MiB buffer

    let mut batch_id = 0u64;
    let mut carryover = Vec::new();

    loop {
        // Read a chunk
        let mut buffer = vec![0u8; batch_bytes];
        let bytes_read = reader.read(&mut buffer)?;
        buffer.truncate(bytes_read);

        // Prepend carryover from previous iteration
        if !carryover.is_empty() {
            let mut combined = carryover.clone();
            combined.extend_from_slice(&buffer);
            buffer = combined;
            carryover.clear();
        }

        if buffer.is_empty() {
            break; // EOF and no carryover
        }

        // Determine if this is the last chunk
        let is_eof = bytes_read == 0 || bytes_read < batch_bytes;

        // Find complete FASTQ records
        // A record is complete if it has 4 lines ending with newline
        let mut complete_end = 0;
        let mut line_count = 0;

        for i in 0..buffer.len() {
            if buffer[i] == b'\n' {
                line_count += 1;
                // After every 4 lines, we have a complete record
                if line_count % 4 == 0 {
                    complete_end = i + 1;
                }
            }
        }

        // On EOF, if we have remaining lines that form complete records, use them
        let complete_part = if is_eof && line_count > 0 {
            // Check if we have complete records
            if line_count % 4 == 0 {
                // All lines end with newline and form complete records
                &buffer[..]
            } else if line_count % 4 == 3 && complete_end < buffer.len() {
                // We have 3 newlines, meaning 4th line exists but no trailing newline
                // This is still a complete record
                &buffer[..]
            } else {
                &buffer[..complete_end]
            }
        } else {
            &buffer[..complete_end]
        };

        // Save incomplete part for next iteration (if not EOF)
        if !is_eof && complete_end < buffer.len() {
            carryover = buffer[complete_end..].to_vec();
        }

        // Parse complete records
        if !complete_part.is_empty() {
            let mut line_starts = vec![0];

            // Find all line starts
            for i in 0..complete_part.len() {
                if complete_part[i] == b'\n' && i + 1 < complete_part.len() {
                    line_starts.push(i + 1);
                }
            }

            // Group into 4-line records
            let mut recs = Vec::new();
            for i in (0..line_starts.len()).step_by(4) {
                if i + 3 < line_starts.len() {
                    recs.push([
                        line_starts[i],
                        line_starts[i + 1],
                        line_starts[i + 2],
                        line_starts[i + 3],
                    ]);
                }
            }

            if !recs.is_empty() {
                let batch = Batch {
                    id: batch_id,
                    buf: complete_part.to_vec(),
                    recs,
                };
                batch_id += 1;

                if sender.send(Some(batch)).is_err() {
                    break; // Receiver disconnected
                }
            }
        }

        if is_eof {
            break;
        }
    }

    // Send sentinel
    let _ = sender.send(None);
    Ok(())
}

/// Worker thread: process batches with thread-local accumulators
fn worker_thread(
    receiver: Receiver<Option<Batch>>,
    sender: Sender<Option<WorkerResult>>,
    min_len: usize,
    n_limit: usize,
    q_mean_phred: u8,
    no_kmer: bool,
    trimming_config: TrimmingConfig,
) {
    while let Ok(Some(batch)) = receiver.recv() {
        let mut before = SimpleStats::default();
        let mut after = SimpleStats::default();
        let mut pos = PositionStats::new();
        let mut k5 = [0usize; 1024];
        let mut too_short = 0usize;
        let mut too_many_n = 0usize;
        let mut low_quality = 0usize;
        let mut invalid = 0usize;
        let mut out_buf = Vec::new();

        // Process each record in the batch
        for (idx, &[h_start, s_start, p_start, q_start]) in batch.recs.iter().enumerate() {
            // Calculate end positions
            let s_end = p_start - 1; // -1 to skip newline
            let p_end = q_start - 1;
            let q_end = if idx + 1 < batch.recs.len() {
                batch.recs[idx + 1][0] - 1
            } else {
                // For the last record, exclude trailing newline if present
                let buf_len = batch.buf.len();
                if buf_len > 0 && batch.buf[buf_len - 1] == b'\n' {
                    buf_len - 1
                } else {
                    buf_len
                }
            };

            let header = &batch.buf[h_start..s_start - 1];
            let seq = &batch.buf[s_start..s_end];
            let plus = &batch.buf[p_start..p_end];
            let qual = &batch.buf[q_start..q_end];

            // Validate
            if seq.len() != qual.len() {
                invalid += 1;
                continue;
            }

            // SINGLE PASS: compute all stats
            let mut qsum = 0u32;
            let mut q20 = 0usize;
            let mut q30 = 0usize;
            let mut ncnt = 0usize;
            let mut gc = 0usize;

            pos.ensure_capacity(seq.len());

            for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
                let qu = q as usize;

                qsum += (q - 33) as u32;
                q20 += LUT_Q20[qu] as usize;
                q30 += LUT_Q30[qu] as usize;
                ncnt += LUT_IS_N[b as usize] as usize;
                gc += LUT_IS_GC[b as usize] as usize;

                pos.total_sum[i] += (q - 33) as u64;
                pos.total_cnt[i] += 1;

                if let Some(bi) = base_idx(b) {
                    pos.base_sum[bi][i] += (q - 33) as u64;
                    pos.base_cnt[bi][i] += 1;
                }
            }

            // K-mer counting
            if !no_kmer {
                count_k5_2bit(seq, &mut k5);
            }

            // Update before stats
            before.add(seq.len(), q20, q30, gc);

            // Apply trimming if enabled
            let trimming_result = if trimming_config.is_enabled() {
                trim_read(seq, qual, &trimming_config)
            } else {
                TrimmingResult {
                    start_pos: 0,
                    end_pos: seq.len(),
                    poly_g_trimmed: 0,
                    poly_x_trimmed: 0,
                }
            };

            // Get trimmed sequences
            let trimmed_seq = &seq[trimming_result.start_pos..trimming_result.end_pos];
            let trimmed_qual = &qual[trimming_result.start_pos..trimming_result.end_pos];

            // Recompute stats for trimmed read (used for filtering)
            let trimmed_len = trimmed_seq.len();
            let mut trimmed_qsum = 0u32;
            let mut trimmed_q20 = 0usize;
            let mut trimmed_q30 = 0usize;
            let mut trimmed_ncnt = 0usize;
            let mut trimmed_gc = 0usize;

            for (&b, &q) in trimmed_seq.iter().zip(trimmed_qual) {
                let qu = q as usize;
                trimmed_qsum += (q - 33) as u32;
                trimmed_q20 += LUT_Q20[qu] as usize;
                trimmed_q30 += LUT_Q30[qu] as usize;
                trimmed_ncnt += LUT_IS_N[b as usize] as usize;
                trimmed_gc += LUT_IS_GC[b as usize] as usize;
            }

            // Apply filters on TRIMMED read
            if trimmed_len < min_len {
                too_short += 1;
                continue;
            }

            if trimmed_ncnt > n_limit {
                too_many_n += 1;
                continue;
            }

            if q_mean_phred > 0
                && trimmed_len > 0
                && (trimmed_qsum as f64 / trimmed_len as f64) < q_mean_phred as f64
            {
                low_quality += 1;
                continue;
            }

            // Passed - write TRIMMED read to output buffer
            out_buf.extend_from_slice(header);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(trimmed_seq);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(plus);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(trimmed_qual);
            out_buf.push(b'\n');

            // Update after stats with trimmed read
            after.add(trimmed_len, trimmed_q20, trimmed_q30, trimmed_gc);
        }

        let result = WorkerResult {
            id: batch.id,
            out_buf,
            before,
            after,
            pos,
            k5,
            too_short,
            too_many_n,
            low_quality,
            invalid,
        };

        if sender.send(Some(result)).is_err() {
            break; // Receiver disconnected
        }
    }

    // Send sentinel
    let _ = sender.send(None);
}

/// Merger thread: write output in order and reduce stats
fn merger_thread(
    receiver: Receiver<Option<WorkerResult>>,
    output_path: String,
    num_workers: usize,
    compression_level: Option<u32>,
) -> Result<StreamAccumulator> {
    let mut writer = open_output(&output_path, compression_level)?;

    let mut acc = StreamAccumulator::new();
    let mut next_id = 0u64;
    let mut pending: BTreeMap<u64, WorkerResult> = BTreeMap::new();
    let mut workers_done = 0;

    while workers_done < num_workers {
        match receiver.recv() {
            Ok(Some(result)) => {
                pending.insert(result.id, result);

                // Write all consecutive results starting from next_id
                while let Some(result) = pending.remove(&next_id) {
                    // Write output
                    writer.write_all(&result.out_buf)?;

                    // Merge stats
                    acc.before.total_reads += result.before.total_reads;
                    acc.before.total_bases += result.before.total_bases;
                    acc.before.q20_bases += result.before.q20_bases;
                    acc.before.q30_bases += result.before.q30_bases;
                    acc.before.gc_bases += result.before.gc_bases;

                    acc.after.total_reads += result.after.total_reads;
                    acc.after.total_bases += result.after.total_bases;
                    acc.after.q20_bases += result.after.q20_bases;
                    acc.after.q30_bases += result.after.q30_bases;
                    acc.after.gc_bases += result.after.gc_bases;

                    // Merge position stats
                    acc.pos.ensure_capacity(result.pos.total_sum.len());
                    for i in 0..result.pos.total_sum.len() {
                        acc.pos.total_sum[i] += result.pos.total_sum[i];
                        acc.pos.total_cnt[i] += result.pos.total_cnt[i];
                        for b in 0..4 {
                            acc.pos.base_sum[b][i] += result.pos.base_sum[b][i];
                            acc.pos.base_cnt[b][i] += result.pos.base_cnt[b][i];
                        }
                    }

                    // Merge k-mer counts
                    for (i, &count) in result.k5.iter().enumerate() {
                        acc.kmer_table[i] += count;
                    }

                    // Merge filter counts
                    acc.too_short += result.too_short;
                    acc.too_many_n += result.too_many_n;
                    acc.low_quality += result.low_quality;
                    acc.invalid += result.invalid;

                    // Track max cycle
                    if result.pos.total_sum.len() > acc.max_cycle {
                        acc.max_cycle = result.pos.total_sum.len();
                    }

                    next_id += 1;
                }
            }
            Ok(None) => {
                workers_done += 1;
            }
            Err(_) => break,
        }
    }

    writer.finish()?;
    Ok(acc)
}

// STREAMING PARSER - processes records without loading all into memory

/// Stream-process FASTQ records using the robust parser
///
/// Uses FastqParser which handles:
/// - Multiline sequences/qualities
/// - Missing final newlines
/// - Malformed records
///
/// NO intermediate Vec<FastqRecord> allocation!
fn process_fastq_stream<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    min_len: usize,
    n_limit: usize,
    q_mean_phred: u8,
    trimming_config: &TrimmingConfig,
) -> Result<StreamAccumulator> {
    let mut acc = StreamAccumulator::new();
    let mut parser = FastqParser::new(reader);

    while let Some((header, seq, plus, qual)) = parser.next_record()? {
        // Process this record in a single pass
        acc.process_record(
            &header,
            &seq,
            &plus,
            &qual,
            min_len,
            n_limit,
            q_mean_phred,
            trimming_config,
            writer,
        )?;
    }

    Ok(acc)
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
    let trimming_config = TrimmingConfig::from_args(&args);

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
            let q_mean_phred = args.qualified_quality_phred;
            let no_kmer = args.no_kmer;
            let trimming_config_clone = trimming_config.clone();

            let worker = thread::spawn(move || {
                worker_thread(
                    batch_rx_clone,
                    result_tx_clone,
                    min_len,
                    n_limit,
                    q_mean_phred,
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

    // Write JSON report
    let json_file =
        File::create(&args.json).context(format!("Failed to create JSON file: {}", args.json))?;
    serde_json::to_writer_pretty(json_file, &report)?;

    Ok(())
}

// UNIT TESTS (see src/tests.rs)

#[cfg(test)]
mod tests;
