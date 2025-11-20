//! Single-threaded FASTQ processing
//!
//! This module provides the streaming single-threaded processing path:
//! - `ParseState`: FASTQ parsing state machine
//! - `FastqParser`: Handles multiline sequences and missing newlines
//! - `StreamAccumulator`: Single-threaded stats accumulator
//! - `process_fastq_stream`: Main single-threaded entry point

use anyhow::Result;
use indexmap::IndexMap;
use std::io::{BufRead, Write};

use crate::kmer::{base_idx, count_k5_2bit};
use crate::simd;
use crate::stats::{PositionStats, SimpleStats};
use crate::trimming::{TrimmingConfig, TrimmingResult, trim_read, trim_read_with_adapter};
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
    #[allow(clippy::type_complexity)]
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
    pub pos_after: PositionStats,
    pub kmer_table: [usize; 1024],
    pub kmer_table_after: [usize; 1024],

    // Filtering counts
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
    pub max_cycle: usize,

    // Adapter trimming stats
    pub adapter_trimmed_reads: usize,
    pub adapter_trimmed_bases: usize,
}

/// Paired-end accumulator for streaming processing
pub(crate) struct PairedEndAccumulator {
    pub before_r1: SimpleStats,
    pub before_r2: SimpleStats,
    pub after_r1: SimpleStats,
    pub after_r2: SimpleStats,
    pub pos_r1: PositionStats,
    pub pos_r2: PositionStats,
    pub pos_r1_after: PositionStats,
    pub pos_r2_after: PositionStats,
    pub kmer_table_r1: [usize; 1024],
    pub kmer_table_r2: [usize; 1024],
    pub kmer_table_r1_after: [usize; 1024],
    pub kmer_table_r2_after: [usize; 1024],

    // Filtering counts
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
    pub duplicated: usize,
    pub max_cycle_r1: usize,
    pub max_cycle_r2: usize,

    // Adapter trimming stats
    pub adapter_trimmed_reads: usize,
    pub adapter_trimmed_bases: usize,

    // Insert size stats (paired-end only)
    pub insert_size_histogram: Vec<usize>,
    pub insert_size_unknown: usize,
}

impl StreamAccumulator {
    pub(crate) fn new() -> Self {
        Self {
            before: SimpleStats::default(),
            after: SimpleStats::default(),
            pos: PositionStats::new(),
            pos_after: PositionStats::new(),
            kmer_table: [0; 1024],
            kmer_table_after: [0; 1024],
            too_short: 0,
            too_many_n: 0,
            low_quality: 0,
            low_complexity: 0,
            invalid: 0,
            max_cycle: 0,
            adapter_trimmed_reads: 0,
            adapter_trimmed_bases: 0,
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
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_lines)]
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
        let simd_available = simd::is_simd_available();
        let (_qsum, q20, q30, q40, _ncnt, gc) = if simd_available {
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

                let qual_hist_ptr = self.pos.qual_hist.as_mut_ptr();
                util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                    let quality_score = u64::from(q - 33);
                    // Direct pointer arithmetic - no Vec::as_mut_slice overhead!
                    *total_sum_ptr.add(i) += quality_score;
                    *total_cnt_ptr.add(i) += 1;

                    if let Some(bi) = base_idx(b) {
                        *base_sum_ptrs[bi].add(i) += quality_score;
                        *base_cnt_ptrs[bi].add(i) += 1;
                    }

                    // Track quality histogram
                    if quality_score < 94 {
                        *qual_hist_ptr.add(quality_score as usize) += 1;
                    }
                });
            }

            (
                stats.qsum, stats.q20, stats.q30, stats.q40, stats.ncnt, stats.gc,
            )
        } else {
            // Non-SIMD path: single pass for both basic and position stats
            let mut qsum = 0u32;
            let mut q20 = 0usize;
            let mut q30 = 0usize;
            let mut q40 = 0usize;
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

                let qual_hist_ptr = self.pos.qual_hist.as_mut_ptr();
                util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                    let quality_u32 = u32::from(q - 33);
                    qsum += quality_u32;
                    if q >= 53 {
                        q20 += 1;
                    }
                    if q >= 63 {
                        q30 += 1;
                    }
                    if q >= 73 {
                        q40 += 1;
                    }
                    if b == b'N' || b == b'n' {
                        ncnt += 1;
                    }
                    if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                        gc += 1;
                    }

                    // Direct pointer arithmetic - no Vec::as_mut_slice overhead!
                    *total_sum_ptr.add(i) += u64::from(quality_u32);
                    *total_cnt_ptr.add(i) += 1;

                    if let Some(bi) = base_idx(b) {
                        *base_sum_ptrs[bi].add(i) += u64::from(quality_u32);
                        *base_cnt_ptrs[bi].add(i) += 1;
                    }

                    // Track quality histogram
                    if (quality_u32 as usize) < 94 {
                        *qual_hist_ptr.add(quality_u32 as usize) += 1;
                    }
                });
            }

            (qsum, q20, q30, q40, ncnt, gc)
        };

        // K-mer counting (also in the same pass conceptually, but separate loop for clarity)
        count_k5_2bit(seq, &mut self.kmer_table);

        // Update "before" stats
        self.before.add(seq.len(), q20, q30, q40, gc);

        // Apply trimming if enabled
        let trimming_result = if trimming_config.is_enabled() {
            trim_read(seq, qual, trimming_config)
        } else {
            TrimmingResult {
                start_pos: 0,
                end_pos: seq.len(),
                poly_g_trimmed: 0,
                poly_x_trimmed: 0,
                adapter_trimmed: false,
                adapter_bases_trimmed: 0,
            }
        };

        // Track adapter trimming stats
        if trimming_result.adapter_trimmed {
            self.adapter_trimmed_reads += 1;
            self.adapter_trimmed_bases += trimming_result.adapter_bases_trimmed;

            // Debug output for investigating base count difference
            if std::env::var("FASTERP_DEBUG_ADAPTER").is_ok() {
                eprintln!(
                    "ADAPTER_TRIM: read={} seq_len={} start={} end={} bases_trimmed={}",
                    self.before.total_reads,
                    seq.len(),
                    trimming_result.start_pos,
                    trimming_result.end_pos,
                    trimming_result.adapter_bases_trimmed
                );
            }
        }

        // Get trimmed sequences
        let trimmed_seq = &seq[trimming_result.start_pos..trimming_result.end_pos];
        let trimmed_qual = &qual[trimming_result.start_pos..trimming_result.end_pos];
        let trimmed_len = trimmed_seq.len();

        // Recompute stats for trimmed read (used for filtering) - SIMD accelerated
        let trimmed_stats = simd::compute_stats(trimmed_seq, trimmed_qual, qualified_quality_phred);
        let trimmed_qsum = trimmed_stats.qsum;
        let trimmed_q20 = trimmed_stats.q20;
        let trimmed_q30 = trimmed_stats.q30;
        let trimmed_q40 = trimmed_stats.q40;
        let trimmed_ncnt = trimmed_stats.ncnt;
        let trimmed_gc = trimmed_stats.gc;
        let unqualified_count = trimmed_stats.unqualified;

        // Apply filters in fastp order: quality FIRST, then N-bases, then length

        // Check unqualified percent (fastp -q/-u logic)
        // Match fastp's exact formula: lowQualNum > (unqualifiedPercentLimit * rlen / 100.0)
        if qualified_quality_phred > 0 && trimmed_len > 0 {
            let threshold = (unqualified_percent_limit * trimmed_len) as f64 / 100.0;
            if (unqualified_count as f64) > threshold {
                self.low_quality += 1;
                return Ok(());
            }
        }

        // Check average quality (fastp -e logic)
        // Use integer division to match fastp exactly
        if average_qual > 0 && trimmed_len > 0 {
            let mean_qual = trimmed_qsum / trimmed_len as u32;
            if mean_qual < u32::from(average_qual) {
                self.low_quality += 1;
                return Ok(());
            }
        }

        // Check N-base filter (fastp checks this AFTER quality but BEFORE length)
        if trimmed_ncnt > n_limit {
            self.too_many_n += 1;
            return Ok(());
        }

        // Length check (fastp -l logic) - AFTER quality checks to match fastp behavior
        if trimmed_len < min_len {
            self.too_short += 1;
            return Ok(());
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
        self.after.add(
            trimmed_len,
            trimmed_q20,
            trimmed_q30,
            trimmed_q40,
            trimmed_gc,
        );

        // Track position stats for reads that passed filters (after filtering)
        track_position_stats(trimmed_seq, trimmed_qual, &mut self.pos_after);

        // K-mer counting for after-filtering data
        count_k5_2bit(trimmed_seq, &mut self.kmer_table_after);

        Ok(())
    }

    /// Convert `kmer_table` to `IndexMap` for JSON output
    /// Uses static string cache - only 1024 allocations total (reused across calls)
    pub(crate) fn kmer_table_to_map(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // Only include kmers that were actually seen (count > 0)
            if self.kmer_table[code] > 0 {
                map.insert(
                    crate::kmer::kmer_to_str(code).to_string(),
                    self.kmer_table[code],
                );
            }
        }
        map
    }

    /// Convert `kmer_table_after` to `IndexMap` for JSON output
    pub(crate) fn kmer_table_after_to_map(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // Only include kmers that were actually seen (count > 0)
            if self.kmer_table_after[code] > 0 {
                map.insert(
                    crate::kmer::kmer_to_str(code).to_string(),
                    self.kmer_table_after[code],
                );
            }
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
            pos_r1_after: PositionStats::new(),
            pos_r2_after: PositionStats::new(),
            kmer_table_r1: [0; 1024],
            kmer_table_r2: [0; 1024],
            kmer_table_r1_after: [0; 1024],
            kmer_table_r2_after: [0; 1024],
            too_short: 0,
            too_many_n: 0,
            low_quality: 0,
            low_complexity: 0,
            invalid: 0,
            duplicated: 0,
            max_cycle_r1: 0,
            max_cycle_r2: 0,
            adapter_trimmed_reads: 0,
            adapter_trimmed_bases: 0,
            insert_size_histogram: vec![0; 512], // Track insert sizes up to 512bp
            insert_size_unknown: 0,
        }
    }

    /// Convert kmer tables to `IndexMaps`
    pub(crate) fn kmer_table_to_map_r1(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // Only include kmers that were actually seen (count > 0)
            if self.kmer_table_r1[code] > 0 {
                map.insert(
                    crate::kmer::kmer_to_str(code).to_string(),
                    self.kmer_table_r1[code],
                );
            }
        }
        map
    }

    pub(crate) fn kmer_table_to_map_r2(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // Only include kmers that were actually seen (count > 0)
            if self.kmer_table_r2[code] > 0 {
                map.insert(
                    crate::kmer::kmer_to_str(code).to_string(),
                    self.kmer_table_r2[code],
                );
            }
        }
        map
    }

    pub(crate) fn kmer_table_to_map_r1_after(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // Only include kmers that were actually seen (count > 0)
            if self.kmer_table_r1_after[code] > 0 {
                map.insert(
                    crate::kmer::kmer_to_str(code).to_string(),
                    self.kmer_table_r1_after[code],
                );
            }
        }
        map
    }

    pub(crate) fn kmer_table_to_map_r2_after(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            // Only include kmers that were actually seen (count > 0)
            if self.kmer_table_r2_after[code] > 0 {
                map.insert(
                    crate::kmer::kmer_to_str(code).to_string(),
                    self.kmer_table_r2_after[code],
                );
            }
        }
        map
    }

    /// Calculate the peak (mode) insert size from the histogram
    pub(crate) fn calculate_insert_size_peak(&self) -> usize {
        if self.insert_size_histogram.is_empty() {
            return 0;
        }

        let mut max_count = 0;
        let mut peak = 0;
        for (size, &count) in self.insert_size_histogram.iter().enumerate() {
            if count > max_count {
                max_count = count;
                peak = size;
            }
        }
        peak
    }
}

/// Single-threaded streaming FASTQ processing
///
/// Processes FASTQ records one-by-one with zero intermediate allocations.
///
/// NO intermediate Vec<FastqRecord> allocation!
#[allow(clippy::too_many_arguments)]
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
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)]
pub(crate) fn process_paired_fastq_stream<
    R1: BufRead,
    R2: BufRead,
    W1: Write,
    W2: Write,
    WM: Write,
>(
    reader1: R1,
    reader2: R2,
    writer1: &mut W1,
    writer2: &mut W2,
    mut merged_writer: Option<&mut WM>,
    merge_enabled: bool,
    include_unmerged: bool,
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
    umi_config: Option<&crate::umi::UmiConfig>,
    dedup_config: Option<&crate::dedup::DedupConfig>,
    dup_detector: Option<&std::sync::Arc<crate::bloom::DuplicateDetector>>,
) -> Result<PairedEndAccumulator> {
    let mut acc = PairedEndAccumulator::new();
    let mut parser1 = FastqParser::new(reader1);
    let mut parser2 = FastqParser::new(reader2);

    // Initialize dedup tracker if deduplication is enabled
    let mut dedup_tracker = if let Some(cfg) = dedup_config {
        if cfg.enabled {
            Some(crate::dedup::DedupTracker::new(cfg.accuracy))
        } else {
            None
        }
    } else {
        None
    };

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

                // Track duplicates using Bloom filter (before any processing)
                // This is for statistics only, not filtering
                if let Some(detector) = dup_detector {
                    detector.check_pair(&seq1, &seq2);
                }

                // Extract UMI if enabled (happens BEFORE any other processing)
                let (
                    final_header1,
                    final_header2,
                    final_seq1,
                    final_qual1,
                    final_seq2,
                    final_qual2,
                );
                let header1_with_umi;
                let header2_with_umi;
                let seq1_after_umi;
                let qual1_after_umi;
                let seq2_after_umi;
                let qual2_after_umi;

                if let Some(umi_cfg) = umi_config {
                    match umi_cfg.location {
                        crate::umi::UmiLocation::Read1 => {
                            // Extract UMI from read1
                            if let Some(extraction) = crate::umi::extract_umi(&seq1, umi_cfg) {
                                // Add UMI to both headers
                                header1_with_umi = crate::umi::add_umi_to_header(
                                    &header1,
                                    &extraction.umi_seq,
                                    &umi_cfg.prefix,
                                );
                                header2_with_umi = crate::umi::add_umi_to_header(
                                    &header2,
                                    &extraction.umi_seq,
                                    &umi_cfg.prefix,
                                );

                                // Trim UMI from read1 sequence and quality
                                seq1_after_umi = seq1[extraction.end_pos..].to_vec();
                                qual1_after_umi = qual1[extraction.end_pos..].to_vec();
                                seq2_after_umi = seq2.clone();
                                qual2_after_umi = qual2.clone();

                                final_header1 = &header1_with_umi[..];
                                final_header2 = &header2_with_umi[..];
                                final_seq1 = &seq1_after_umi[..];
                                final_qual1 = &qual1_after_umi[..];
                                final_seq2 = &seq2_after_umi[..];
                                final_qual2 = &qual2_after_umi[..];
                            } else {
                                // UMI extraction failed - use original (convert to owned for consistency)
                                header1_with_umi = header1.clone();
                                header2_with_umi = header2.clone();
                                seq1_after_umi = seq1.clone();
                                qual1_after_umi = qual1.clone();
                                seq2_after_umi = seq2.clone();
                                qual2_after_umi = qual2.clone();

                                final_header1 = &header1_with_umi[..];
                                final_header2 = &header2_with_umi[..];
                                final_seq1 = &seq1_after_umi[..];
                                final_qual1 = &qual1_after_umi[..];
                                final_seq2 = &seq2_after_umi[..];
                                final_qual2 = &qual2_after_umi[..];
                            }
                        }
                        crate::umi::UmiLocation::Read2 => {
                            // Extract UMI from read2
                            if let Some(extraction) = crate::umi::extract_umi(&seq2, umi_cfg) {
                                // Add UMI to both headers
                                header1_with_umi = crate::umi::add_umi_to_header(
                                    &header1,
                                    &extraction.umi_seq,
                                    &umi_cfg.prefix,
                                );
                                header2_with_umi = crate::umi::add_umi_to_header(
                                    &header2,
                                    &extraction.umi_seq,
                                    &umi_cfg.prefix,
                                );

                                // Trim UMI from read2 sequence and quality
                                seq1_after_umi = seq1.clone();
                                qual1_after_umi = qual1.clone();
                                seq2_after_umi = seq2[extraction.end_pos..].to_vec();
                                qual2_after_umi = qual2[extraction.end_pos..].to_vec();

                                final_header1 = &header1_with_umi[..];
                                final_header2 = &header2_with_umi[..];
                                final_seq1 = &seq1_after_umi[..];
                                final_qual1 = &qual1_after_umi[..];
                                final_seq2 = &seq2_after_umi[..];
                                final_qual2 = &qual2_after_umi[..];
                            } else {
                                // UMI extraction failed - use original (convert to owned for consistency)
                                header1_with_umi = header1.clone();
                                header2_with_umi = header2.clone();
                                seq1_after_umi = seq1.clone();
                                qual1_after_umi = qual1.clone();
                                seq2_after_umi = seq2.clone();
                                qual2_after_umi = qual2.clone();

                                final_header1 = &header1_with_umi[..];
                                final_header2 = &header2_with_umi[..];
                                final_seq1 = &seq1_after_umi[..];
                                final_qual1 = &qual1_after_umi[..];
                                final_seq2 = &seq2_after_umi[..];
                                final_qual2 = &qual2_after_umi[..];
                            }
                        }
                        _ => {
                            // Other UMI locations not yet supported - use original
                            header1_with_umi = header1.clone();
                            header2_with_umi = header2.clone();
                            seq1_after_umi = seq1.clone();
                            qual1_after_umi = qual1.clone();
                            seq2_after_umi = seq2.clone();
                            qual2_after_umi = qual2.clone();

                            final_header1 = &header1_with_umi[..];
                            final_header2 = &header2_with_umi[..];
                            final_seq1 = &seq1_after_umi[..];
                            final_qual1 = &qual1_after_umi[..];
                            final_seq2 = &seq2_after_umi[..];
                            final_qual2 = &qual2_after_umi[..];
                        }
                    }
                } else {
                    // No UMI processing - use original (convert to owned for consistency)
                    header1_with_umi = header1.clone();
                    header2_with_umi = header2.clone();
                    seq1_after_umi = seq1.clone();
                    qual1_after_umi = qual1.clone();
                    seq2_after_umi = seq2.clone();
                    qual2_after_umi = qual2.clone();

                    final_header1 = &header1_with_umi[..];
                    final_header2 = &header2_with_umi[..];
                    final_seq1 = &seq1_after_umi[..];
                    final_qual1 = &qual1_after_umi[..];
                    final_seq2 = &seq2_after_umi[..];
                    final_qual2 = &qual2_after_umi[..];
                }

                // Track max cycle lengths (use final sequences after UMI removal)
                if final_seq1.len() > acc.max_cycle_r1 {
                    acc.max_cycle_r1 = final_seq1.len();
                }
                if final_seq2.len() > acc.max_cycle_r2 {
                    acc.max_cycle_r2 = final_seq2.len();
                }

                // Process read1 statistics (before filtering, after UMI removal)
                let stats1 = process_read_stats(
                    final_seq1,
                    final_qual1,
                    &mut acc.pos_r1,
                    &mut acc.kmer_table_r1,
                );
                acc.before_r1.add(
                    final_seq1.len(),
                    stats1.q20,
                    stats1.q30,
                    stats1.q40,
                    stats1.gc,
                );

                // Process read2 statistics (before filtering, after UMI removal)
                let stats2 = process_read_stats(
                    final_seq2,
                    final_qual2,
                    &mut acc.pos_r2,
                    &mut acc.kmer_table_r2,
                );
                acc.before_r2.add(
                    final_seq2.len(),
                    stats2.q20,
                    stats2.q30,
                    stats2.q40,
                    stats2.gc,
                );

                // Always detect overlap for paired-end adapter trimming (fastp does this)
                // Overlap-based trimming works even when adapters aren't detected via k-mer analysis
                let overlap_result = {
                    let overlap_detect_config = crate::overlap::OverlapConfig {
                        min_overlap_len: 30, // fastp default (options.cpp:overlapRequire = 30)
                        max_diff: 5,
                        max_diff_percent: 20,
                    };
                    crate::overlap::detect_overlap(final_seq1, final_seq2, &overlap_detect_config)
                };

                // Track insert size from overlap detection
                if let Some(ref overlap) = overlap_result {
                    if overlap.overlapped {
                        // Insert size calculation depends on offset (fastp logic):
                        // - offset > 0: reads don't fully overlap, use traditional formula
                        // - offset <= 0: reads fully overlap with adapters, insert size = overlap_len
                        let insert_size = if overlap.offset > 0 {
                            final_seq1.len() + final_seq2.len() - overlap.overlap_len
                        } else {
                            overlap.overlap_len
                        };
                        if insert_size < acc.insert_size_histogram.len() {
                            acc.insert_size_histogram[insert_size] += 1;
                        }
                    } else {
                        acc.insert_size_unknown += 1;
                    }
                } else {
                    acc.insert_size_unknown += 1;
                }

                // Try overlap-based adapter trimming first
                // Always use overlap-based trimming when overlap is detected (fastp does this)
                let overlap_trimmed = if let Some(ref overlap) = overlap_result {
                    crate::overlap::trim_by_overlap_analysis(
                        final_seq1.len(),
                        final_seq2.len(),
                        overlap,
                    )
                } else {
                    None
                };

                // Apply trimming to both reads (using sequences after UMI removal)
                let mut trim_result1;
                let mut trim_result2;

                if let Some((overlap_trim_len1, overlap_trim_len2)) = overlap_trimmed {
                    // Overlap-based adapter trimming succeeded
                    // Still need to apply other trimming (poly-G, quality, etc.)
                    // but skip adapter trimming since we already did it with overlap

                    // Temporarily disable adapter trimming for other trims
                    let mut temp_config1 = trimming_config_r1.clone();
                    let mut temp_config2 = trimming_config_r2.clone();
                    temp_config1.adapter_config.adapter_seq = None;
                    temp_config2.adapter_config.adapter_seq = None;
                    temp_config2.adapter_config.adapter_seq_r2 = None;

                    trim_result1 = if trimming_config_r1.is_enabled() {
                        trim_read(final_seq1, final_qual1, &temp_config1)
                    } else {
                        TrimmingResult {
                            start_pos: 0,
                            end_pos: final_seq1.len(),
                            poly_g_trimmed: 0,
                            poly_x_trimmed: 0,
                            adapter_trimmed: false,
                            adapter_bases_trimmed: 0,
                        }
                    };

                    trim_result2 = if trimming_config_r2.is_enabled() {
                        trim_read(final_seq2, final_qual2, &temp_config2)
                    } else {
                        TrimmingResult {
                            start_pos: 0,
                            end_pos: final_seq2.len(),
                            poly_g_trimmed: 0,
                            poly_x_trimmed: 0,
                            adapter_trimmed: false,
                            adapter_bases_trimmed: 0,
                        }
                    };

                    // Apply overlap-based adapter trim lengths
                    // Track the end position before overlap trim to calculate adapter bases correctly
                    let end_before_overlap1 = trim_result1.end_pos;
                    let end_before_overlap2 = trim_result2.end_pos;

                    trim_result1.end_pos = std::cmp::min(trim_result1.end_pos, overlap_trim_len1);
                    trim_result2.end_pos = std::cmp::min(trim_result2.end_pos, overlap_trim_len2);

                    // Mark as adapter trimmed (only count bases trimmed beyond other trims)
                    if trim_result1.end_pos < end_before_overlap1 {
                        trim_result1.adapter_trimmed = true;
                        trim_result1.adapter_bases_trimmed =
                            end_before_overlap1 - trim_result1.end_pos;
                    }
                    if trim_result2.end_pos < end_before_overlap2 {
                        trim_result2.adapter_trimmed = true;
                        trim_result2.adapter_bases_trimmed =
                            end_before_overlap2 - trim_result2.end_pos;
                    }
                } else {
                    // No overlap-based trimming, use sequence-based adapter trimming
                    trim_result1 = if trimming_config_r1.is_enabled() {
                        trim_read(final_seq1, final_qual1, trimming_config_r1)
                    } else {
                        TrimmingResult {
                            start_pos: 0,
                            end_pos: final_seq1.len(),
                            poly_g_trimmed: 0,
                            poly_x_trimmed: 0,
                            adapter_trimmed: false,
                            adapter_bases_trimmed: 0,
                        }
                    };

                    trim_result2 = if trimming_config_r2.is_enabled() {
                        // Use adapter_seq_r2 for read2 if available
                        let adapter_r2 =
                            trimming_config_r2.adapter_config.adapter_seq_r2.as_deref();
                        trim_read_with_adapter(
                            final_seq2,
                            final_qual2,
                            trimming_config_r2,
                            adapter_r2,
                        )
                    } else {
                        TrimmingResult {
                            start_pos: 0,
                            end_pos: final_seq2.len(),
                            poly_g_trimmed: 0,
                            poly_x_trimmed: 0,
                            adapter_trimmed: false,
                            adapter_bases_trimmed: 0,
                        }
                    };
                }

                // Track adapter trimming stats for paired-end
                if trim_result1.adapter_trimmed {
                    acc.adapter_trimmed_reads += 1;
                    acc.adapter_trimmed_bases += trim_result1.adapter_bases_trimmed;
                }
                if trim_result2.adapter_trimmed {
                    acc.adapter_trimmed_reads += 1;
                    acc.adapter_trimmed_bases += trim_result2.adapter_bases_trimmed;
                }

                let trimmed_seq1 = &final_seq1[trim_result1.start_pos..trim_result1.end_pos];
                let trimmed_qual1 = &final_qual1[trim_result1.start_pos..trim_result1.end_pos];
                let trimmed_seq2 = &final_seq2[trim_result2.start_pos..trim_result2.end_pos];
                let trimmed_qual2 = &final_qual2[trim_result2.start_pos..trim_result2.end_pos];

                // Apply base correction using overlap analysis (if enabled)
                let (corrected_seq1, corrected_qual1, corrected_seq2, corrected_qual2);
                let mut seq1_corrected_buf;
                let mut qual1_corrected_buf;
                let mut seq2_corrected_buf;
                let mut qual2_corrected_buf;

                if let Some(config) = overlap_config {
                    // Create mutable copies for correction
                    seq1_corrected_buf = trimmed_seq1.to_vec();
                    qual1_corrected_buf = trimmed_qual1.to_vec();
                    seq2_corrected_buf = trimmed_seq2.to_vec();
                    qual2_corrected_buf = trimmed_qual2.to_vec();

                    // Detect overlap on trimmed sequences for correction
                    // (trimmed sequences may have different overlap than original)
                    if let Some(overlap) = crate::overlap::detect_overlap(
                        &seq1_corrected_buf,
                        &seq2_corrected_buf,
                        config,
                    ) {
                        let _correction_stats = crate::overlap::correct_by_overlap(
                            &mut seq1_corrected_buf,
                            &mut qual1_corrected_buf,
                            &mut seq2_corrected_buf,
                            &mut qual2_corrected_buf,
                            &overlap,
                        );
                        // TODO: Track correction statistics in accumulator
                    }

                    corrected_seq1 = &seq1_corrected_buf[..];
                    corrected_qual1 = &qual1_corrected_buf[..];
                    corrected_seq2 = &seq2_corrected_buf[..];
                    corrected_qual2 = &qual2_corrected_buf[..];
                } else {
                    // No correction - use trimmed sequences directly
                    corrected_seq1 = trimmed_seq1;
                    corrected_qual1 = trimmed_qual1;
                    corrected_seq2 = trimmed_seq2;
                    corrected_qual2 = trimmed_qual2;
                }

                // Check if either read is too short
                if corrected_seq1.len() < min_len || corrected_seq2.len() < min_len {
                    acc.too_short += 1;
                    continue;
                }

                // Recompute stats for final reads (after UMI removal, trimming and correction)
                let trimmed_stats1 =
                    simd::compute_stats(corrected_seq1, corrected_qual1, qualified_quality_phred);
                let trimmed_stats2 =
                    simd::compute_stats(corrected_seq2, corrected_qual2, qualified_quality_phred);

                // Mimic fastp's per-read filtering with max() logic
                // In fastp, all quality and N-base checks are inside qualfilter.enabled block
                // Filter R1
                let mut result1_fail_quality = false;
                let mut result1_fail_n = false;

                // R1: Check quality filters (fastp checks quality, then N-bases in else-if chain)
                // Use fastp's exact comparison logic with floating point to match behavior
                if !corrected_seq1.is_empty() {
                    let rlen1 = corrected_seq1.len();
                    let unqual_threshold1 = (unqualified_percent_limit * rlen1) as f64 / 100.0;
                    if (trimmed_stats1.unqualified as f64) > unqual_threshold1 {
                        result1_fail_quality = true;
                    } else if average_qual > 0 {
                        let mean_qual1 = f64::from(trimmed_stats1.qsum) / rlen1 as f64;
                        if mean_qual1 < f64::from(average_qual) {
                            result1_fail_quality = true;
                        }
                    } else if trimmed_stats1.ncnt > n_limit {
                        // Only check N-bases if quality checks passed (fastp's else if logic)
                        result1_fail_n = true;
                    }
                }

                // Filter R2
                let mut result2_fail_quality = false;
                let mut result2_fail_n = false;

                // R2: Check quality filters (fastp checks quality, then N-bases in else-if chain)
                // Use fastp's exact comparison logic with floating point to match behavior
                if !corrected_seq2.is_empty() {
                    let rlen2 = corrected_seq2.len();
                    let unqual_threshold2 = (unqualified_percent_limit * rlen2) as f64 / 100.0;
                    if (trimmed_stats2.unqualified as f64) > unqual_threshold2 {
                        result2_fail_quality = true;
                    } else if average_qual > 0 {
                        let mean_qual2 = f64::from(trimmed_stats2.qsum) / rlen2 as f64;
                        if mean_qual2 < f64::from(average_qual) {
                            result2_fail_quality = true;
                        }
                    } else if trimmed_stats2.ncnt > n_limit {
                        // Only check N-bases if quality checks passed (fastp's else if logic)
                        result2_fail_n = true;
                    }
                }

                // Use max() logic like fastp: FAIL_QUALITY (20) > FAIL_N_BASE (12)
                // If either read fails quality, count as low_quality
                if result1_fail_quality || result2_fail_quality {
                    acc.low_quality += 1;
                    continue;
                }
                // Otherwise, if either read fails N-base, count as too_many_n
                if result1_fail_n || result2_fail_n {
                    acc.too_many_n += 1;
                    continue;
                }

                // Check low complexity for both reads
                if low_complexity_filter {
                    let mut fail_complexity = false;
                    if !corrected_seq1.is_empty() {
                        let complexity1 = calculate_complexity(corrected_seq1);
                        if complexity1 < complexity_threshold {
                            fail_complexity = true;
                        }
                    }
                    if !corrected_seq2.is_empty() {
                        let complexity2 = calculate_complexity(corrected_seq2);
                        if complexity2 < complexity_threshold {
                            fail_complexity = true;
                        }
                    }

                    if fail_complexity {
                        acc.low_complexity += 1;
                        continue;
                    }
                }

                // Check for duplicates if deduplication is enabled
                if let Some(ref mut tracker) = dedup_tracker {
                    if tracker.is_duplicate_pe(corrected_seq1, corrected_seq2) {
                        // This is a duplicate - skip it
                        acc.duplicated += 1;
                        continue;
                    }
                }

                // Try to merge if merge mode enabled
                if merge_enabled {
                    if let Some(ref mut mwriter) = merged_writer {
                        // Detect overlap for merging
                        if let Some(overlap) = crate::overlap::detect_overlap(
                            corrected_seq1,
                            corrected_seq2,
                            &crate::overlap::OverlapConfig::default(),
                        ) {
                            if overlap.overlapped {
                                // Merge the reads
                                let (merged_header, merged_seq, merged_qual) =
                                    crate::overlap::merge_reads(
                                        corrected_seq1,
                                        corrected_qual1,
                                        final_header1,
                                        corrected_seq2,
                                        corrected_qual2,
                                        &overlap,
                                    );

                                // Write merged read
                                writeln!(mwriter, "{}", std::str::from_utf8(&merged_header)?)?;
                                writeln!(mwriter, "{}", std::str::from_utf8(&merged_seq)?)?;
                                writeln!(mwriter, "+")?;
                                writeln!(mwriter, "{}", std::str::from_utf8(&merged_qual)?)?;

                                // Update "after" stats for merged read
                                let merged_stats = simd::compute_stats(
                                    &merged_seq,
                                    &merged_qual,
                                    qualified_quality_phred,
                                );
                                acc.after_r1.add(
                                    merged_seq.len(),
                                    merged_stats.q20,
                                    merged_stats.q30,
                                    merged_stats.q40,
                                    merged_stats.gc,
                                );

                                // Track position stats for merged read
                                track_position_stats(
                                    &merged_seq,
                                    &merged_qual,
                                    &mut acc.pos_r1_after,
                                );

                                // Track kmers for merged read
                                count_k5_2bit(&merged_seq, &mut acc.kmer_table_r1_after);

                                continue; // Skip normal R1/R2 output
                            }
                        }

                        // Merge failed or no overlap - write unmerged if flag set
                        if include_unmerged {
                            writeln!(mwriter, "{}", std::str::from_utf8(final_header1)?)?;
                            writeln!(mwriter, "{}", std::str::from_utf8(corrected_seq1)?)?;
                            writeln!(mwriter, "{}", std::str::from_utf8(&plus1)?)?;
                            writeln!(mwriter, "{}", std::str::from_utf8(corrected_qual1)?)?;

                            writeln!(mwriter, "{}", std::str::from_utf8(final_header2)?)?;
                            writeln!(mwriter, "{}", std::str::from_utf8(corrected_seq2)?)?;
                            writeln!(mwriter, "{}", std::str::from_utf8(&plus2)?)?;
                            writeln!(mwriter, "{}", std::str::from_utf8(corrected_qual2)?)?;

                            // Update "after" stats
                            acc.after_r1.add(
                                corrected_seq1.len(),
                                trimmed_stats1.q20,
                                trimmed_stats1.q30,
                                trimmed_stats1.q40,
                                trimmed_stats1.gc,
                            );
                            acc.after_r2.add(
                                corrected_seq2.len(),
                                trimmed_stats2.q20,
                                trimmed_stats2.q30,
                                trimmed_stats2.q40,
                                trimmed_stats2.gc,
                            );

                            // Track position stats for unmerged reads
                            track_position_stats(
                                corrected_seq1,
                                corrected_qual1,
                                &mut acc.pos_r1_after,
                            );
                            track_position_stats(
                                corrected_seq2,
                                corrected_qual2,
                                &mut acc.pos_r2_after,
                            );

                            // Track kmers for unmerged reads
                            count_k5_2bit(corrected_seq1, &mut acc.kmer_table_r1_after);
                            count_k5_2bit(corrected_seq2, &mut acc.kmer_table_r2_after);

                            continue; // Skip normal R1/R2 output
                        }
                    }
                }

                // Normal output (merge disabled, or merge failed without include_unmerged)
                writeln!(writer1, "{}", std::str::from_utf8(final_header1)?)?;
                writeln!(writer1, "{}", std::str::from_utf8(corrected_seq1)?)?;
                writeln!(writer1, "{}", std::str::from_utf8(&plus1)?)?;
                writeln!(writer1, "{}", std::str::from_utf8(corrected_qual1)?)?;

                writeln!(writer2, "{}", std::str::from_utf8(final_header2)?)?;
                writeln!(writer2, "{}", std::str::from_utf8(corrected_seq2)?)?;
                writeln!(writer2, "{}", std::str::from_utf8(&plus2)?)?;
                writeln!(writer2, "{}", std::str::from_utf8(corrected_qual2)?)?;

                // Update "after" stats
                acc.after_r1.add(
                    corrected_seq1.len(),
                    trimmed_stats1.q20,
                    trimmed_stats1.q30,
                    trimmed_stats1.q40,
                    trimmed_stats1.gc,
                );
                acc.after_r2.add(
                    corrected_seq2.len(),
                    trimmed_stats2.q20,
                    trimmed_stats2.q30,
                    trimmed_stats2.q40,
                    trimmed_stats2.gc,
                );

                // Track position stats for reads that passed filters (after filtering)
                track_position_stats(corrected_seq1, corrected_qual1, &mut acc.pos_r1_after);
                track_position_stats(corrected_seq2, corrected_qual2, &mut acc.pos_r2_after);

                // Track kmers for reads that passed filters
                count_k5_2bit(corrected_seq1, &mut acc.kmer_table_r1_after);
                count_k5_2bit(corrected_seq2, &mut acc.kmer_table_r2_after);
            }
        }
    }

    Ok(acc)
}

// EXPERT OPTIMIZATION:
// Map ASCII bases to 0..3 index instantly.
// 255 acts as the "None" sentinel to avoid an Option enum (which adds branching).
const BASE_LUT: [u8; 256] = {
    let mut lut = [255u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'T' as usize] = 1;
    lut[b'C' as usize] = 2;
    lut[b'G' as usize] = 3;
    // Handle lower case too if necessary, though fastp is usually upper.
    lut[b'a' as usize] = 0;
    lut[b't' as usize] = 1;
    lut[b'c' as usize] = 2;
    lut[b'g' as usize] = 3;
    lut
};

/// High-performance stats tracker
/// Uses raw pointers and manual unrolling to maximize instruction throughput.
pub(crate) fn track_position_stats(seq: &[u8], qual: &[u8], pos: &mut PositionStats) {
    let len = seq.len();
    pos.ensure_capacity(len);

    unsafe {
        let seq_ptr = seq.as_ptr();
        let qual_ptr = qual.as_ptr();

        let total_sum_ptr = pos.total_sum.as_mut_ptr();
        let total_cnt_ptr = pos.total_cnt.as_mut_ptr();
        let qual_hist_ptr = pos.qual_hist.as_mut_ptr();

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

        for i in 0..len {
            let b = *seq_ptr.add(i);
            let q = *qual_ptr.add(i);
            let qual_val = (q.wrapping_sub(33)) as usize;
            let qual_u64 = qual_val as u64;

            *total_sum_ptr.add(i) += qual_u64;
            *total_cnt_ptr.add(i) += 1;

            if qual_val < 94 {
                *qual_hist_ptr.add(qual_val) += 1;
            }

            let bi = *BASE_LUT.get_unchecked(b as usize);
            if bi != 255 {
                let bi_idx = bi as usize;
                *base_sum_ptrs.get_unchecked(bi_idx).add(i) += qual_u64;
                *base_cnt_ptrs.get_unchecked(bi_idx).add(i) += 1;
            }
        }
    }
}

/// Helper function to process read statistics
struct ReadStats {
    q20: usize,
    q30: usize,
    q40: usize,
    gc: usize,
}

fn process_read_stats(
    seq: &[u8],
    qual: &[u8],
    pos: &mut PositionStats,
    kmer_table: &mut [usize; 1024],
) -> ReadStats {
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

            let qual_hist_ptr = pos.qual_hist.as_mut_ptr();
            util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                let quality_score = u64::from(q - 33);
                *total_sum_ptr.add(i) += quality_score;
                *total_cnt_ptr.add(i) += 1;

                if let Some(bi) = base_idx(b) {
                    *base_sum_ptrs[bi].add(i) += quality_score;
                    *base_cnt_ptrs[bi].add(i) += 1;
                }

                // Track quality histogram
                if quality_score < 94 {
                    *qual_hist_ptr.add(quality_score as usize) += 1;
                }
            });
        }

        ReadStats {
            q20: s.q20,
            q30: s.q30,
            q40: s.q40,
            gc: s.gc,
        }
    } else {
        // Non-SIMD path
        let mut q20 = 0usize;
        let mut q30 = 0usize;
        let mut q40 = 0usize;
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

            let qual_hist_ptr = pos.qual_hist.as_mut_ptr();
            util::loop_seq_qual_indexed(seq_ptr, qual_ptr, len, |i, b, q| {
                let quality_u32 = u32::from(q - 33);
                if q >= 53 {
                    q20 += 1;
                }
                if q >= 63 {
                    q30 += 1;
                }
                if q >= 73 {
                    q40 += 1;
                }
                if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                    gc += 1;
                }

                *total_sum_ptr.add(i) += u64::from(quality_u32);
                *total_cnt_ptr.add(i) += 1;

                if let Some(bi) = base_idx(b) {
                    *base_sum_ptrs[bi].add(i) += u64::from(quality_u32);
                    *base_cnt_ptrs[bi].add(i) += 1;
                }

                // Track quality histogram
                if (quality_u32 as usize) < 94 {
                    *qual_hist_ptr.add(quality_u32 as usize) += 1;
                }
            });
        }

        ReadStats { q20, q30, q40, gc }
    };

    // K-mer counting
    count_k5_2bit(seq, kmer_table);

    stats
}
