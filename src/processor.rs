//! Single-threaded FASTQ processing
//!
//! This module provides the streaming single-threaded processing path:
//! - ParseState: FASTQ parsing state machine
//! - FastqParser: Handles multiline sequences and missing newlines
//! - StreamAccumulator: Single-threaded stats accumulator
//! - process_fastq_stream: Main single-threaded entry point

use anyhow::Result;
use indexmap::IndexMap;
use std::io::{BufRead, Write};

use crate::kmer::*;
use crate::simd;
use crate::stats::*;
use crate::trimming::*;
use crate::util;

/// Calculate sequence complexity as percentage of bases different from next base
///
/// Complexity is defined as: (count of positions where base[i] != base[i+1]) / (length - 1) * 100
/// This matches fastp's low complexity filter algorithm.
///
/// Returns complexity percentage (0-100)
#[inline]
fn calculate_complexity(seq: &[u8]) -> usize {
    if seq.len() <= 1 {
        return 100; // Single base or empty is considered max complexity
    }

    let mut different_count = 0;
    for i in 0..seq.len() - 1 {
        if seq[i] != seq[i + 1] {
            different_count += 1;
        }
    }

    // Calculate percentage: (different_count / (len - 1)) * 100
    // Use integer math to avoid floating point
    (different_count * 100) / (seq.len() - 1)
}

/// State machine for parsing FASTQ records
/// Handles:
/// - Multiline sequences/qualities (wrapped lines)
/// - Missing final newline
/// - Malformed records (skips with warning)
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub(crate) enum ParseState {
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
        // Pre-allocate for typical FASTQ record sizes
        const TYPICAL_READ_LEN: usize = 150;
        const TYPICAL_HEADER_LEN: usize = 64;
        const TYPICAL_LINE_LEN: usize = 256;

        FastqParser {
            reader,
            state: ParseState::ExpectHeader,
            header: Vec::with_capacity(TYPICAL_HEADER_LEN),
            sequence: Vec::with_capacity(TYPICAL_READ_LEN),
            plus: Vec::with_capacity(4), // Usually just "+"
            quality: Vec::with_capacity(TYPICAL_READ_LEN),
            line_buf: Vec::with_capacity(TYPICAL_LINE_LEN),
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

/// Main accumulator for streaming processing
pub(crate) struct StreamAccumulator {
    pub before: SimpleStats,
    pub after: SimpleStats,
    pub pos: PositionStats,
    pub kmer_table: [usize; 1024],

    // Filtering counts
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
    pub max_cycle: usize,
}

/// Paired-end accumulator for streaming processing
pub(crate) struct PairedEndAccumulator {
    pub before_r1: SimpleStats,
    pub before_r2: SimpleStats,
    pub after_r1: SimpleStats,
    pub after_r2: SimpleStats,
    pub pos_r1: PositionStats,
    pub pos_r2: PositionStats,
    pub kmer_table_r1: [usize; 1024],
    pub kmer_table_r2: [usize; 1024],

    // Filtering counts
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
    pub max_cycle_r1: usize,
    pub max_cycle_r2: usize,
}

impl StreamAccumulator {
    pub(crate) fn new() -> Self {
        Self {
            before: SimpleStats::default(),
            after: SimpleStats::default(),
            pos: PositionStats::new(),
            kmer_table: [0; 1024],
            too_short: 0,
            too_many_n: 0,
            low_quality: 0,
            low_complexity: 0,
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
    pub(crate) fn process_record(
        &mut self,
        header: &[u8],
        seq: &[u8],
        plus: &[u8],
        qual: &[u8],
        min_len: usize,
        n_limit: usize,
        qualified_quality_phred: u8,
        unqualified_percent_limit: usize,
        average_qual: u8,
        low_complexity_filter: bool,
        complexity_threshold: usize,
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

        // Compute stats - use SIMD when available, otherwise single-pass
        let (_qsum, q20, q30, _ncnt, gc) = if simd::is_simd_available() {
            // SIMD path: compute basic stats fast, then position-specific
            let stats = simd::compute_stats(seq, qual, 0); // 0 = don't count unqualified for before stats

            self.pos.ensure_capacity(seq.len());

            // SAFETY: We've validated seq.len() == qual.len() above, and ensured capacity
            // Extract raw pointers ONCE before loop to avoid repeated Vec::as_mut_slice calls
            unsafe {
                let seq_ptr = seq.as_ptr();
                let qual_ptr = qual.as_ptr();
                let len = seq.len();

                // Get raw pointers to all position stat arrays (ONCE, not per iteration!)
                let total_sum_ptr = self.pos.total_sum.as_mut_ptr();
                let total_cnt_ptr = self.pos.total_cnt.as_mut_ptr();
                let base_sum_ptrs = [
                    self.pos.base_sum[0].as_mut_ptr(),
                    self.pos.base_sum[1].as_mut_ptr(),
                    self.pos.base_sum[2].as_mut_ptr(),
                    self.pos.base_sum[3].as_mut_ptr(),
                ];
                let base_cnt_ptrs = [
                    self.pos.base_cnt[0].as_mut_ptr(),
                    self.pos.base_cnt[1].as_mut_ptr(),
                    self.pos.base_cnt[2].as_mut_ptr(),
                    self.pos.base_cnt[3].as_mut_ptr(),
                ];

                util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                    let qval = (q - 33) as u64;
                    // Direct pointer arithmetic - no Vec::as_mut_slice overhead!
                    *total_sum_ptr.add(i) += qval;
                    *total_cnt_ptr.add(i) += 1;

                    if let Some(bi) = base_idx(b) {
                        *base_sum_ptrs[bi].add(i) += qval;
                        *base_cnt_ptrs[bi].add(i) += 1;
                    }
                });
            }

            (stats.qsum, stats.q20, stats.q30, stats.ncnt, stats.gc)
        } else {
            // Non-SIMD path: single pass for both basic and position stats
            let mut qsum = 0u32;
            let mut q20 = 0usize;
            let mut q30 = 0usize;
            let mut ncnt = 0usize;
            let mut gc = 0usize;

            self.pos.ensure_capacity(seq.len());

            // SAFETY: We've validated seq.len() == qual.len() above, and ensured capacity
            // Extract raw pointers ONCE before loop to avoid repeated Vec::as_mut_slice calls
            unsafe {
                let seq_ptr = seq.as_ptr();
                let qual_ptr = qual.as_ptr();
                let len = seq.len();

                // Get raw pointers to all position stat arrays (ONCE, not per iteration!)
                let total_sum_ptr = self.pos.total_sum.as_mut_ptr();
                let total_cnt_ptr = self.pos.total_cnt.as_mut_ptr();
                let base_sum_ptrs = [
                    self.pos.base_sum[0].as_mut_ptr(),
                    self.pos.base_sum[1].as_mut_ptr(),
                    self.pos.base_sum[2].as_mut_ptr(),
                    self.pos.base_sum[3].as_mut_ptr(),
                ];
                let base_cnt_ptrs = [
                    self.pos.base_cnt[0].as_mut_ptr(),
                    self.pos.base_cnt[1].as_mut_ptr(),
                    self.pos.base_cnt[2].as_mut_ptr(),
                    self.pos.base_cnt[3].as_mut_ptr(),
                ];

                util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                    let qval = (q - 33) as u32;
                    qsum += qval;
                    if q >= 53 {
                        q20 += 1;
                    }
                    if q >= 63 {
                        q30 += 1;
                    }
                    if b == b'N' || b == b'n' {
                        ncnt += 1;
                    }
                    if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                        gc += 1;
                    }

                    // Direct pointer arithmetic - no Vec::as_mut_slice overhead!
                    *total_sum_ptr.add(i) += qval as u64;
                    *total_cnt_ptr.add(i) += 1;

                    if let Some(bi) = base_idx(b) {
                        *base_sum_ptrs[bi].add(i) += qval as u64;
                        *base_cnt_ptrs[bi].add(i) += 1;
                    }
                });
            }

            (qsum, q20, q30, ncnt, gc)
        };

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

        // Early length check - skip expensive stats computation for reads that are too short
        let trimmed_len = trimmed_seq.len();
        if trimmed_len < min_len {
            self.too_short += 1;
            return Ok(());
        }

        // Recompute stats for trimmed read (used for filtering) - SIMD accelerated
        let trimmed_stats = simd::compute_stats(trimmed_seq, trimmed_qual, qualified_quality_phred);
        let trimmed_qsum = trimmed_stats.qsum;
        let trimmed_q20 = trimmed_stats.q20;
        let trimmed_q30 = trimmed_stats.q30;
        let trimmed_ncnt = trimmed_stats.ncnt;
        let trimmed_gc = trimmed_stats.gc;
        let unqualified_count = trimmed_stats.unqualified;

        // Apply remaining filters on TRIMMED read

        if trimmed_ncnt > n_limit {
            self.too_many_n += 1;
            return Ok(());
        }

        // Check unqualified percent (fastp -q/-u logic)
        // unqualified_count already computed by SIMD above
        if qualified_quality_phred > 0 && trimmed_len > 0 {
            // Avoid division to prevent rounding issues: check if 100*unqualified > limit*len
            if 100 * unqualified_count > unqualified_percent_limit * trimmed_len {
                self.low_quality += 1;
                return Ok(());
            }
        }

        // Check average quality (fastp -e logic)
        if average_qual > 0 && trimmed_len > 0 {
            let mean_qual = trimmed_qsum as f64 / trimmed_len as f64;
            if mean_qual < average_qual as f64 {
                self.low_quality += 1;
                return Ok(());
            }
        }

        // Check low complexity (fastp -y/-Y logic)
        if low_complexity_filter && trimmed_len > 0 {
            let complexity = calculate_complexity(trimmed_seq);
            if complexity < complexity_threshold {
                self.low_complexity += 1;
                return Ok(());
            }
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
    /// Uses static string cache - only 1024 allocations total (reused across calls)
    pub(crate) fn kmer_table_to_map(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // kmer_to_str returns &'static str from cache, we convert to String for the map
            // This still allocates, but only once per unique kmer string (cached in kmer_to_str)
            map.insert(
                crate::kmer::kmer_to_str(code).to_string(),
                self.kmer_table[code],
            );
        }
        map
    }
}

impl PairedEndAccumulator {
    pub(crate) fn new() -> Self {
        Self {
            before_r1: SimpleStats::default(),
            before_r2: SimpleStats::default(),
            after_r1: SimpleStats::default(),
            after_r2: SimpleStats::default(),
            pos_r1: PositionStats::new(),
            pos_r2: PositionStats::new(),
            kmer_table_r1: [0; 1024],
            kmer_table_r2: [0; 1024],
            too_short: 0,
            too_many_n: 0,
            low_quality: 0,
            low_complexity: 0,
            invalid: 0,
            max_cycle_r1: 0,
            max_cycle_r2: 0,
        }
    }

    /// Convert kmer tables to IndexMaps
    pub(crate) fn kmer_table_to_map_r1(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            map.insert(
                crate::kmer::kmer_to_str(code).to_string(),
                self.kmer_table_r1[code],
            );
        }
        map
    }

    pub(crate) fn kmer_table_to_map_r2(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            map.insert(
                crate::kmer::kmer_to_str(code).to_string(),
                self.kmer_table_r2[code],
            );
        }
        map
    }
}

/// Single-threaded streaming FASTQ processing
///
/// Processes FASTQ records one-by-one with zero intermediate allocations.
///
/// NO intermediate Vec<FastqRecord> allocation!
pub(crate) fn process_fastq_stream<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    min_len: usize,
    n_limit: usize,
    qualified_quality_phred: u8,
    unqualified_percent_limit: usize,
    average_qual: u8,
    low_complexity_filter: bool,
    complexity_threshold: usize,
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
            qualified_quality_phred,
            unqualified_percent_limit,
            average_qual,
            low_complexity_filter,
            complexity_threshold,
            trimming_config,
            writer,
        )?;
    }

    Ok(acc)
}

/// Paired-end streaming FASTQ processing
///
/// Processes two FASTQ files simultaneously, maintaining read pair synchronization.
/// Read pairs pass/fail together - if either fails filtering, both are discarded.
pub(crate) fn process_paired_fastq_stream<R1: BufRead, R2: BufRead, W1: Write, W2: Write>(
    reader1: R1,
    reader2: R2,
    writer1: &mut W1,
    writer2: &mut W2,
    min_len: usize,
    n_limit: usize,
    qualified_quality_phred: u8,
    unqualified_percent_limit: usize,
    average_qual: u8,
    low_complexity_filter: bool,
    complexity_threshold: usize,
    trimming_config_r1: &TrimmingConfig,
    trimming_config_r2: &TrimmingConfig,
    overlap_config: Option<&crate::overlap::OverlapConfig>,
) -> Result<PairedEndAccumulator> {
    let mut acc = PairedEndAccumulator::new();
    let mut parser1 = FastqParser::new(reader1);
    let mut parser2 = FastqParser::new(reader2);

    loop {
        let record1 = parser1.next_record()?;
        let record2 = parser2.next_record()?;

        // Check for EOF - both files must end at the same time
        match (record1, record2) {
            (None, None) => break, // Both files ended - normal termination
            (Some(_), None) | (None, Some(_)) => {
                anyhow::bail!("Read1 and Read2 files have different numbers of records");
            }
            (Some((header1, seq1, plus1, qual1)), Some((header2, seq2, plus2, qual2))) => {
                // Validate both records
                if seq1.len() != qual1.len() || seq2.len() != qual2.len() {
                    acc.invalid += 1;
                    continue;
                }

                // Track max cycle lengths
                if seq1.len() > acc.max_cycle_r1 {
                    acc.max_cycle_r1 = seq1.len();
                }
                if seq2.len() > acc.max_cycle_r2 {
                    acc.max_cycle_r2 = seq2.len();
                }

                // Process read1 statistics (before filtering)
                let stats1 =
                    process_read_stats(&seq1, &qual1, &mut acc.pos_r1, &mut acc.kmer_table_r1)?;
                acc.before_r1
                    .add(seq1.len(), stats1.q20, stats1.q30, stats1.gc);

                // Process read2 statistics (before filtering)
                let stats2 =
                    process_read_stats(&seq2, &qual2, &mut acc.pos_r2, &mut acc.kmer_table_r2)?;
                acc.before_r2
                    .add(seq2.len(), stats2.q20, stats2.q30, stats2.gc);

                // Apply trimming to both reads
                let trim_result1 = if trimming_config_r1.is_enabled() {
                    trim_read(&seq1, &qual1, trimming_config_r1)
                } else {
                    TrimmingResult {
                        start_pos: 0,
                        end_pos: seq1.len(),
                        poly_g_trimmed: 0,
                        poly_x_trimmed: 0,
                    }
                };

                let trim_result2 = if trimming_config_r2.is_enabled() {
                    // Use adapter_seq_r2 for read2 if available
                    let adapter_r2 = trimming_config_r2.adapter_config.adapter_seq_r2.as_deref();
                    trim_read_with_adapter(&seq2, &qual2, trimming_config_r2, adapter_r2)
                } else {
                    TrimmingResult {
                        start_pos: 0,
                        end_pos: seq2.len(),
                        poly_g_trimmed: 0,
                        poly_x_trimmed: 0,
                    }
                };

                let trimmed_seq1 = &seq1[trim_result1.start_pos..trim_result1.end_pos];
                let trimmed_qual1 = &qual1[trim_result1.start_pos..trim_result1.end_pos];
                let trimmed_seq2 = &seq2[trim_result2.start_pos..trim_result2.end_pos];
                let trimmed_qual2 = &qual2[trim_result2.start_pos..trim_result2.end_pos];

                // Apply base correction using overlap analysis (if enabled)
                let (final_seq1, final_qual1, final_seq2, final_qual2);
                let mut seq1_corrected;
                let mut qual1_corrected;
                let mut seq2_corrected;
                let mut qual2_corrected;

                if let Some(config) = overlap_config {
                    // Create mutable copies for correction
                    seq1_corrected = trimmed_seq1.to_vec();
                    qual1_corrected = trimmed_qual1.to_vec();
                    seq2_corrected = trimmed_seq2.to_vec();
                    qual2_corrected = trimmed_qual2.to_vec();

                    // Detect overlap and correct mismatches
                    if let Some(overlap) =
                        crate::overlap::detect_overlap(&seq1_corrected, &seq2_corrected, config)
                    {
                        let _correction_stats = crate::overlap::correct_by_overlap(
                            &mut seq1_corrected,
                            &mut qual1_corrected,
                            &mut seq2_corrected,
                            &mut qual2_corrected,
                            &overlap,
                        );
                        // TODO: Track correction statistics in accumulator
                    }

                    final_seq1 = &seq1_corrected[..];
                    final_qual1 = &qual1_corrected[..];
                    final_seq2 = &seq2_corrected[..];
                    final_qual2 = &qual2_corrected[..];
                } else {
                    // No correction - use trimmed sequences directly
                    final_seq1 = trimmed_seq1;
                    final_qual1 = trimmed_qual1;
                    final_seq2 = trimmed_seq2;
                    final_qual2 = trimmed_qual2;
                }

                // Check if either read is too short
                if final_seq1.len() < min_len || final_seq2.len() < min_len {
                    acc.too_short += 1;
                    continue;
                }

                // Recompute stats for final reads (after trimming and correction)
                let trimmed_stats1 =
                    simd::compute_stats(final_seq1, final_qual1, qualified_quality_phred);
                let trimmed_stats2 =
                    simd::compute_stats(final_seq2, final_qual2, qualified_quality_phred);

                // Check N-base filter for both reads
                if trimmed_stats1.ncnt > n_limit || trimmed_stats2.ncnt > n_limit {
                    acc.too_many_n += 1;
                    continue;
                }

                // Check unqualified percent for both reads
                let mut fail_quality = false;
                if qualified_quality_phred > 0 {
                    if final_seq1.len() > 0
                        && 100 * trimmed_stats1.unqualified
                            > unqualified_percent_limit * final_seq1.len()
                    {
                        fail_quality = true;
                    }
                    if final_seq2.len() > 0
                        && 100 * trimmed_stats2.unqualified
                            > unqualified_percent_limit * final_seq2.len()
                    {
                        fail_quality = true;
                    }
                }

                // Check average quality for both reads
                if average_qual > 0 {
                    if final_seq1.len() > 0 {
                        let mean_qual1 = trimmed_stats1.qsum as f64 / final_seq1.len() as f64;
                        if mean_qual1 < average_qual as f64 {
                            fail_quality = true;
                        }
                    }
                    if final_seq2.len() > 0 {
                        let mean_qual2 = trimmed_stats2.qsum as f64 / final_seq2.len() as f64;
                        if mean_qual2 < average_qual as f64 {
                            fail_quality = true;
                        }
                    }
                }

                if fail_quality {
                    acc.low_quality += 1;
                    continue;
                }

                // Check low complexity for both reads
                if low_complexity_filter {
                    let mut fail_complexity = false;
                    if final_seq1.len() > 0 {
                        let complexity1 = calculate_complexity(final_seq1);
                        if complexity1 < complexity_threshold {
                            fail_complexity = true;
                        }
                    }
                    if final_seq2.len() > 0 {
                        let complexity2 = calculate_complexity(final_seq2);
                        if complexity2 < complexity_threshold {
                            fail_complexity = true;
                        }
                    }

                    if fail_complexity {
                        acc.low_complexity += 1;
                        continue;
                    }
                }

                // Both reads passed - write them
                writeln!(writer1, "{}", std::str::from_utf8(&header1)?)?;
                writeln!(writer1, "{}", std::str::from_utf8(final_seq1)?)?;
                writeln!(writer1, "{}", std::str::from_utf8(&plus1)?)?;
                writeln!(writer1, "{}", std::str::from_utf8(final_qual1)?)?;

                writeln!(writer2, "{}", std::str::from_utf8(&header2)?)?;
                writeln!(writer2, "{}", std::str::from_utf8(final_seq2)?)?;
                writeln!(writer2, "{}", std::str::from_utf8(&plus2)?)?;
                writeln!(writer2, "{}", std::str::from_utf8(final_qual2)?)?;

                // Update "after" stats
                acc.after_r1.add(
                    final_seq1.len(),
                    trimmed_stats1.q20,
                    trimmed_stats1.q30,
                    trimmed_stats1.gc,
                );
                acc.after_r2.add(
                    final_seq2.len(),
                    trimmed_stats2.q20,
                    trimmed_stats2.q30,
                    trimmed_stats2.gc,
                );
            }
        }
    }

    Ok(acc)
}

/// Helper function to process read statistics
struct ReadStats {
    q20: usize,
    q30: usize,
    gc: usize,
}

fn process_read_stats(
    seq: &[u8],
    qual: &[u8],
    pos: &mut PositionStats,
    kmer_table: &mut [usize; 1024],
) -> Result<ReadStats> {
    // Compute basic stats using SIMD if available
    let stats = if simd::is_simd_available() {
        let s = simd::compute_stats(seq, qual, 0);

        pos.ensure_capacity(seq.len());

        // SAFETY: We've validated seq.len() == qual.len(), and ensured capacity
        unsafe {
            let seq_ptr = seq.as_ptr();
            let qual_ptr = qual.as_ptr();
            let len = seq.len();

            let total_sum_ptr = pos.total_sum.as_mut_ptr();
            let total_cnt_ptr = pos.total_cnt.as_mut_ptr();
            let base_sum_ptrs = [
                pos.base_sum[0].as_mut_ptr(),
                pos.base_sum[1].as_mut_ptr(),
                pos.base_sum[2].as_mut_ptr(),
                pos.base_sum[3].as_mut_ptr(),
            ];
            let base_cnt_ptrs = [
                pos.base_cnt[0].as_mut_ptr(),
                pos.base_cnt[1].as_mut_ptr(),
                pos.base_cnt[2].as_mut_ptr(),
                pos.base_cnt[3].as_mut_ptr(),
            ];

            util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                let qval = (q - 33) as u64;
                *total_sum_ptr.add(i) += qval;
                *total_cnt_ptr.add(i) += 1;

                if let Some(bi) = base_idx(b) {
                    *base_sum_ptrs[bi].add(i) += qval;
                    *base_cnt_ptrs[bi].add(i) += 1;
                }
            });
        }

        ReadStats {
            q20: s.q20,
            q30: s.q30,
            gc: s.gc,
        }
    } else {
        // Non-SIMD path
        let mut q20 = 0usize;
        let mut q30 = 0usize;
        let mut gc = 0usize;

        pos.ensure_capacity(seq.len());

        unsafe {
            let seq_ptr = seq.as_ptr();
            let qual_ptr = qual.as_ptr();
            let len = seq.len();

            let total_sum_ptr = pos.total_sum.as_mut_ptr();
            let total_cnt_ptr = pos.total_cnt.as_mut_ptr();
            let base_sum_ptrs = [
                pos.base_sum[0].as_mut_ptr(),
                pos.base_sum[1].as_mut_ptr(),
                pos.base_sum[2].as_mut_ptr(),
                pos.base_sum[3].as_mut_ptr(),
            ];
            let base_cnt_ptrs = [
                pos.base_cnt[0].as_mut_ptr(),
                pos.base_cnt[1].as_mut_ptr(),
                pos.base_cnt[2].as_mut_ptr(),
                pos.base_cnt[3].as_mut_ptr(),
            ];

            util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                let qval = (q - 33) as u32;
                if q >= 53 {
                    q20 += 1;
                }
                if q >= 63 {
                    q30 += 1;
                }
                if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                    gc += 1;
                }

                *total_sum_ptr.add(i) += qval as u64;
                *total_cnt_ptr.add(i) += 1;

                if let Some(bi) = base_idx(b) {
                    *base_sum_ptrs[bi].add(i) += qval as u64;
                    *base_cnt_ptrs[bi].add(i) += 1;
                }
            });
        }

        ReadStats { q20, q30, gc }
    };

    // K-mer counting
    count_k5_2bit(seq, kmer_table);

    Ok(stats)
}
