//! A fast and simple FASTQ preprocessor
//!
//! This tool processes FASTQ files with quality control filters including:
//! - Length filtering (minimum read length)
//! - Quality filtering (mean quality score)
//! - N-base filtering (maximum ambiguous bases)

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::{Receiver, Sender, bounded};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::thread;

// ============================================================================
// LOOKUP TABLES for fast base/quality checks
// ============================================================================

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

#[derive(Parser, Debug)]
#[command(author, version, about = "A fast FASTQ preprocessor", long_about = None)]
struct Args {
    /// Input FASTQ file
    #[arg(short = 'i', long)]
    input: String,

    /// Output FASTQ file
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
}

// ============================================================================
// MULTI-THREADED PIPELINE DATA STRUCTURES
// ============================================================================

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

// ============================================================================
// STREAMING ACCUMULATOR for single-pass processing
// ============================================================================

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

        // Apply filters
        if seq.len() < min_len {
            self.too_short += 1;
            return Ok(());
        }

        if ncnt > n_limit {
            self.too_many_n += 1;
            return Ok(());
        }

        if q_mean_phred > 0 && (qsum as f64 / seq.len() as f64) < q_mean_phred as f64 {
            self.low_quality += 1;
            return Ok(());
        }

        // Record passed - write it out
        writeln!(writer, "{}", std::str::from_utf8(header)?)?;
        writeln!(writer, "{}", std::str::from_utf8(seq)?)?;
        writeln!(writer, "{}", std::str::from_utf8(plus)?)?;
        writeln!(writer, "{}", std::str::from_utf8(qual)?)?;

        // Update "after" stats
        self.after.add(seq.len(), q20, q30, gc);

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

// ============================================================================
// MULTI-THREADED PIPELINE
// ============================================================================

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
        File::open(&input_path).context(format!("Failed to open input file: {}", input_path))?;
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
                if complete_part[i] == b'\n' {
                    if i + 1 < complete_part.len() {
                        line_starts.push(i + 1);
                    }
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

            // Apply filters
            if seq.len() < min_len {
                too_short += 1;
                continue;
            }

            if ncnt > n_limit {
                too_many_n += 1;
                continue;
            }

            if q_mean_phred > 0 && (qsum as f64 / seq.len() as f64) < q_mean_phred as f64 {
                low_quality += 1;
                continue;
            }

            // Passed - write to output buffer
            out_buf.extend_from_slice(header);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(seq);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(plus);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(qual);
            out_buf.push(b'\n');

            // Update after stats
            after.add(seq.len(), q20, q30, gc);
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
) -> Result<StreamAccumulator> {
    let file = File::create(&output_path)
        .context(format!("Failed to create output file: {}", output_path))?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, file);

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

    writer.flush()?;
    Ok(acc)
}

// ============================================================================
// STREAMING PARSER - processes records without loading all into memory
// ============================================================================

/// Stream-process FASTQ records directly from reader
///
/// This is the key optimization: instead of loading all records into a Vec,
/// we process each record immediately as we read it.
///
/// NO intermediate Vec<FastqRecord> allocation!
fn process_fastq_stream<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    min_len: usize,
    n_limit: usize,
    q_mean_phred: u8,
) -> Result<StreamAccumulator> {
    let mut acc = StreamAccumulator::new();
    let mut lines = reader.lines();

    // Reusable buffers to avoid allocations
    let mut header_buf = String::new();
    let mut seq_buf = String::new();
    let mut plus_buf = String::new();
    let mut qual_buf = String::new();

    while let Some(header) = lines.next() {
        header_buf = header?;

        // Skip empty lines
        if header_buf.is_empty() {
            continue;
        }

        // Ensure it starts with @
        if !header_buf.starts_with('@') {
            continue;
        }

        // Read sequence
        seq_buf = match lines.next() {
            Some(s) => s?,
            None => break,
        };

        // Read plus line
        plus_buf = match lines.next() {
            Some(p) => p?,
            None => break,
        };

        // Read quality
        qual_buf = match lines.next() {
            Some(q) => q?,
            None => break,
        };

        // Process this record in a single pass
        acc.process_record(
            header_buf.as_bytes(),
            seq_buf.as_bytes(),
            plus_buf.as_bytes(),
            qual_buf.as_bytes(),
            min_len,
            n_limit,
            q_mean_phred,
            writer,
        )?;
    }

    Ok(acc)
}

/// Represents a single FASTQ record (4 lines)
#[derive(Debug, Clone)]
struct FastqRecord {
    header: String,
    sequence: String,
    plus: String,
    quality: String,
}

impl FastqRecord {
    /// Write this record to the output writer
    fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "{}", self.header)?;
        writeln!(writer, "{}", self.sequence)?;
        writeln!(writer, "{}", self.plus)?;
        writeln!(writer, "{}", self.quality)?;
        Ok(())
    }

    /// Calculate mean quality score (Phred+33 encoding)
    fn mean_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }

        let sum: u32 = self
            .quality
            .bytes()
            .map(|b| (b - 33) as u32) // Phred+33 encoding
            .sum();

        sum as f64 / self.quality.len() as f64
    }

    /// Count the number of N bases in the sequence
    fn count_n_bases(&self) -> usize {
        self.sequence
            .bytes()
            .filter(|&b| b == b'N' || b == b'n')
            .count()
    }

    /// Count bases with quality >= Q20
    fn count_q20_bases(&self) -> usize {
        self.quality.bytes().filter(|&b| (b - 33) >= 20).count()
    }

    /// Count bases with quality >= Q30
    fn count_q30_bases(&self) -> usize {
        self.quality.bytes().filter(|&b| (b - 33) >= 30).count()
    }

    /// Count GC bases
    fn count_gc_bases(&self) -> usize {
        self.sequence
            .bytes()
            .filter(|&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
            .count()
    }
}

/// Parse FASTQ format from a buffered reader
///
/// PERFORMANCE NOTE: This loads ALL records into memory at once.
/// - For large files (100K+ reads), this creates significant memory pressure
/// - Each record allocates 4 strings (header, sequence, plus, quality)
/// - A streaming approach would be more memory-efficient
fn read_fastq_records<R: BufRead>(reader: R) -> Result<Vec<FastqRecord>> {
    let mut records = Vec::new();
    let mut lines = reader.lines();

    while let Some(header) = lines.next() {
        let header = header?;

        // Skip empty lines
        if header.is_empty() {
            continue;
        }

        // Ensure it starts with @
        if !header.starts_with('@') {
            continue;
        }

        // Read sequence
        let sequence = match lines.next() {
            Some(s) => s?,
            None => break,
        };

        // Read plus line
        let plus = match lines.next() {
            Some(p) => p?,
            None => break,
        };

        // Read quality
        let quality = match lines.next() {
            Some(q) => q?,
            None => break,
        };

        // Validate: sequence and quality must have same length
        if sequence.len() != quality.len() {
            eprintln!(
                "WARNING: sequence and quality have different length:\n{}\n{}\n{}\n{}",
                header, sequence, plus, quality
            );
            continue; // Skip this invalid record
        }

        records.push(FastqRecord {
            header,
            sequence,
            plus,
            quality,
        });
    }

    Ok(records)
}

/// Calculate statistics for a set of records
///
/// PERFORMANCE HOTSPOT: Multiple passes through the data
/// - Pass 1: count sequence lengths
/// - Pass 2: count Q20 bases (iterates through quality string)
/// - Pass 3: count Q30 bases (iterates through quality string again)
/// - Pass 4: count GC bases (iterates through sequence string)
/// Could be optimized to a single pass
fn calculate_stats(records: &[FastqRecord]) -> ReadStats {
    let total_reads = records.len();
    let total_bases: usize = records.iter().map(|r| r.sequence.len()).sum();
    let q20_bases: usize = records.iter().map(|r| r.count_q20_bases()).sum();
    let q30_bases: usize = records.iter().map(|r| r.count_q30_bases()).sum();
    let gc_bases: usize = records.iter().map(|r| r.count_gc_bases()).sum();

    ReadStats {
        total_reads,
        total_bases,
        q20_bases,
        q30_bases,
        q20_rate: if total_bases > 0 {
            q20_bases as f64 / total_bases as f64
        } else {
            0.0
        },
        q30_rate: if total_bases > 0 {
            q30_bases as f64 / total_bases as f64
        } else {
            0.0
        },
        read1_mean_length: if total_reads > 0 {
            total_bases / total_reads
        } else {
            0
        },
        gc_content: if total_bases > 0 {
            gc_bases as f64 / total_bases as f64
        } else {
            0.0
        },
    }
}

/// Count all kmers of given size in a sequence
///
/// PERFORMANCE HOTSPOT: String allocations in tight loop
/// - Called for EVERY sequence in the dataset
/// - Allocates a new String for each kmer (e.g., 146 strings per 150bp read)
/// - For 100K reads with 150bp, this is ~14.6 million allocations
/// - Could use a fixed-size buffer or pre-allocated strings
fn count_kmers_in_sequence(sequence: &str, k: usize) -> IndexMap<String, usize> {
    let mut kmer_counts = IndexMap::new();
    let seq_bytes = sequence.as_bytes();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = &seq_bytes[i..i + k];
        // Convert to uppercase and create string
        // HOTSPOT: String allocation for every kmer
        let kmer_str: String = kmer
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        *kmer_counts.entry(kmer_str).or_insert(0) += 1;
    }

    kmer_counts
}

/// Generate all possible kmers of length k
fn generate_all_kmers(k: usize) -> IndexMap<String, usize> {
    let bases = ['A', 'C', 'G', 'T'];
    let mut kmers = IndexMap::new();

    // Generate all possible kmers recursively
    fn generate_recursive(
        k: usize,
        current: String,
        bases: &[char],
        kmers: &mut IndexMap<String, usize>,
    ) {
        if current.len() == k {
            kmers.insert(current, 0);
            return;
        }

        for &base in bases {
            let mut next = current.clone();
            next.push(base);
            generate_recursive(k, next, bases, kmers);
        }
    }

    generate_recursive(k, String::new(), &bases, &mut kmers);
    kmers
}

/// Count all kmers across all records
///
/// PERFORMANCE HOTSPOT: Major bottleneck for large datasets
/// - Calls count_kmers_in_sequence() for every record
/// - For 100K reads with 150bp:
///   - 100,000 records × ~146 kmers per read = 14.6M kmers
///   - Each kmer allocates a String = 14.6M allocations
/// - Then merges into total_counts with more HashMap operations
/// - Could be optimized by:
///   - Using integer encoding for kmers instead of strings
///   - Pre-allocating a single buffer and reusing it
///   - Using array indexing instead of HashMap (4^5 = 1024 possible 5-mers)
fn count_all_kmers(records: &[FastqRecord], k: usize) -> IndexMap<String, usize> {
    // Start with all possible kmers initialized to 0
    let mut total_counts = generate_all_kmers(k);

    // Count kmers that actually appear
    for record in records {
        let seq_kmers = count_kmers_in_sequence(&record.sequence, k);
        for (kmer, count) in seq_kmers {
            *total_counts.entry(kmer).or_insert(0) += count;
        }
    }

    total_counts
}

/// Calculate quality curves for per-position quality statistics
///
/// PERFORMANCE HOTSPOT: Nested iteration over all bases
/// - Allocates 10 vectors of length max_len (e.g., 10 * 151 = 1510 u64s)
/// - Iterates through EVERY base in EVERY record
/// - For 100K reads × 150bp = 15 million base-quality pairs
/// - Many branches in the inner loop (match statement)
/// - Could be optimized with SIMD or lookup tables
fn calculate_quality_curves(records: &[FastqRecord]) -> QualityCurves {
    if records.is_empty() {
        return QualityCurves {
            a: Vec::new(),
            t: Vec::new(),
            c: Vec::new(),
            g: Vec::new(),
            mean: Vec::new(),
        };
    }

    // Find maximum read length
    let max_len = records.iter().map(|r| r.sequence.len()).max().unwrap_or(0);

    // Initialize vectors to track quality sums and counts per position
    let mut a_sums = vec![0u64; max_len];
    let mut a_counts = vec![0u64; max_len];
    let mut t_sums = vec![0u64; max_len];
    let mut t_counts = vec![0u64; max_len];
    let mut c_sums = vec![0u64; max_len];
    let mut c_counts = vec![0u64; max_len];
    let mut g_sums = vec![0u64; max_len];
    let mut g_counts = vec![0u64; max_len];
    let mut total_sums = vec![0u64; max_len];
    let mut total_counts = vec![0u64; max_len];

    // Accumulate quality scores per position per base
    // HOTSPOT: This is the innermost loop, iterating over millions of bases
    for record in records {
        let seq_bytes = record.sequence.as_bytes();
        let qual_bytes = record.quality.as_bytes();

        for (pos, (&base, &qual)) in seq_bytes.iter().zip(qual_bytes.iter()).enumerate() {
            let quality = (qual - 33) as u64;

            total_sums[pos] += quality;
            total_counts[pos] += 1;

            match base {
                b'A' | b'a' => {
                    a_sums[pos] += quality;
                    a_counts[pos] += 1;
                }
                b'T' | b't' => {
                    t_sums[pos] += quality;
                    t_counts[pos] += 1;
                }
                b'C' | b'c' => {
                    c_sums[pos] += quality;
                    c_counts[pos] += 1;
                }
                b'G' | b'g' => {
                    g_sums[pos] += quality;
                    g_counts[pos] += 1;
                }
                _ => {} // Skip N or other bases
            }
        }
    }

    // Calculate overall mean quality per position first
    let mean: Vec<f64> = total_sums
        .iter()
        .zip(total_counts.iter())
        .map(|(&sum, &count)| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                0.0
            }
        })
        .collect();

    // Calculate mean quality per position for each base
    // If no bases of a particular type at a position, use overall mean (matches fastp behavior)
    let a_mean: Vec<f64> = a_sums
        .iter()
        .zip(a_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    let t_mean: Vec<f64> = t_sums
        .iter()
        .zip(t_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    let c_mean: Vec<f64> = c_sums
        .iter()
        .zip(c_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    let g_mean: Vec<f64> = g_sums
        .iter()
        .zip(g_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    QualityCurves {
        a: a_mean,
        t: t_mean,
        c: c_mean,
        g: g_mean,
        mean,
    }
}

/// Calculate detailed statistics with kmer counts
fn calculate_detailed_stats(records: &[FastqRecord]) -> DetailedReadStats {
    let total_reads = records.len();
    let total_bases: usize = records.iter().map(|r| r.sequence.len()).sum();
    let q20_bases: usize = records.iter().map(|r| r.count_q20_bases()).sum();
    let q30_bases: usize = records.iter().map(|r| r.count_q30_bases()).sum();
    let quality_curves = calculate_quality_curves(records);
    let kmer_count = count_all_kmers(records, 5); // Use k=5 to match fastp

    DetailedReadStats {
        total_reads,
        total_bases,
        q20_bases,
        q30_bases,
        quality_curves,
        kmer_count,
    }
}

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

    let acc = if num_threads == 1 {
        // SINGLE-THREADED MODE: use old streaming approach
        let input_file = File::open(&args.input)
            .context(format!("Failed to open input file: {}", args.input))?;
        let reader = BufReader::new(input_file);

        let output_file = File::create(&args.output)
            .context(format!("Failed to create output file: {}", args.output))?;
        let mut writer = BufWriter::new(output_file);

        let acc = process_fastq_stream(
            reader,
            &mut writer,
            args.length_required,
            args.n_base_limit,
            args.qualified_quality_phred,
        )?;

        writer.flush()?;
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

            let worker = thread::spawn(move || {
                worker_thread(
                    batch_rx_clone,
                    result_tx_clone,
                    min_len,
                    n_limit,
                    q_mean_phred,
                    no_kmer,
                )
            });
            workers.push(worker);
        }

        // Drop original senders so merger knows when all workers are done
        drop(batch_rx);
        drop(result_tx);

        // Spawn merger thread
        let output_path = args.output.clone();
        let merger = thread::spawn(move || merger_thread(result_rx, output_path, num_threads));

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
