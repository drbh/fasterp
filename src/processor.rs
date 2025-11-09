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
    pub invalid: usize,
    pub max_cycle: usize,
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
        let (qsum, q20, q30, ncnt, gc) = if simd::is_simd_available() {
            // SIMD path: compute basic stats fast, then position-specific
            let stats = simd::compute_stats(seq, qual);

            self.pos.ensure_capacity(seq.len());
            for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
                self.pos.total_sum[i] += (q - 33) as u64;
                self.pos.total_cnt[i] += 1;

                if let Some(bi) = base_idx(b) {
                    self.pos.base_sum[bi][i] += (q - 33) as u64;
                    self.pos.base_cnt[bi][i] += 1;
                }
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
            for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
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

                self.pos.total_sum[i] += qval as u64;
                self.pos.total_cnt[i] += 1;

                if let Some(bi) = base_idx(b) {
                    self.pos.base_sum[bi][i] += qval as u64;
                    self.pos.base_cnt[bi][i] += 1;
                }
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

        // Recompute stats for trimmed read (used for filtering) - SIMD accelerated
        let trimmed_len = trimmed_seq.len();
        let trimmed_stats = simd::compute_stats(trimmed_seq, trimmed_qual);
        let trimmed_qsum = trimmed_stats.qsum;
        let trimmed_q20 = trimmed_stats.q20;
        let trimmed_q30 = trimmed_stats.q30;
        let trimmed_ncnt = trimmed_stats.ncnt;
        let trimmed_gc = trimmed_stats.gc;

        // Apply filters on TRIMMED read
        if trimmed_len < min_len {
            self.too_short += 1;
            return Ok(());
        }

        if trimmed_ncnt > n_limit {
            self.too_many_n += 1;
            return Ok(());
        }

        // Check unqualified percent (fastp -q/-u logic)
        if qualified_quality_phred > 0 && trimmed_len > 0 {
            let qual_threshold = qualified_quality_phred + 33; // Convert to ASCII
            let unqualified_count = trimmed_qual.iter().filter(|&&q| q < qual_threshold).count();

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
    pub(crate) fn kmer_table_to_map(&self) -> IndexMap<String, usize> {
        let mut map = IndexMap::new();
        for code in 0..1024 {
            let kmer_str = kmer_to_string(code);
            map.insert(kmer_str, self.kmer_table[code]);
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
            trimming_config,
            writer,
        )?;
    }

    Ok(acc)
}
