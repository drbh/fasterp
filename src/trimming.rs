//! Quality trimming algorithms
//!
//! This module provides various trimming strategies for FASTQ reads:
//! - Fixed position trimming (front and tail)
//! - Sliding window quality-based trimming
//! - PolyG/PolyX tail detection and removal

/// Configuration for trimming operations
#[derive(Debug, Clone)]
pub(crate) struct TrimmingConfig {
    // Sliding window trimming
    pub enable_trim_front: bool,
    pub enable_trim_tail: bool,
    pub cut_mean_quality: u8,
    pub cut_window_size: usize,

    // Fixed position trimming
    pub trim_front_bases: usize,
    pub trim_tail_bases: usize,

    // Length limits
    pub max_len: usize, // Maximum read length (trim tail if longer, 0 = disabled)

    // PolyG/PolyX trimming
    pub enable_poly_g: bool,
    pub enable_poly_x: bool,
    pub poly_min_len: usize,
}

impl TrimmingConfig {
    pub(crate) fn is_enabled(&self) -> bool {
        self.enable_trim_front
            || self.enable_trim_tail
            || self.trim_front_bases > 0
            || self.trim_tail_bases > 0
            || self.max_len > 0
            || self.enable_poly_g
            || self.enable_poly_x
    }
}

/// Result of trimming a single read
#[derive(Debug, Clone, Default)]
pub(crate) struct TrimmingResult {
    pub start_pos: usize, // Starting position after 5' trimming
    pub end_pos: usize,   // Ending position after 3' trimming
    pub poly_g_trimmed: usize,
    pub poly_x_trimmed: usize,
}

impl TrimmingResult {
    #[allow(dead_code)]
    pub(crate) fn trimmed_length(&self) -> usize {
        self.end_pos.saturating_sub(self.start_pos)
    }

    #[allow(dead_code)]
    pub(crate) fn bases_trimmed_5(&self) -> usize {
        self.start_pos
    }

    #[allow(dead_code)]
    pub(crate) fn bases_trimmed_3(&self) -> usize {
        self.poly_g_trimmed + self.poly_x_trimmed
    }
}

/// Accumulated trimming statistics
#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
pub(crate) struct TrimmingStats {
    pub reads_trimmed: usize,
    pub bases_trimmed_5: u64,
    pub bases_trimmed_3: u64,
    pub poly_g_trimmed_reads: usize,
    pub poly_g_trimmed_bases: u64,
    pub poly_x_trimmed_reads: usize,
    pub poly_x_trimmed_bases: u64,
}

#[allow(dead_code)]
impl TrimmingStats {
    pub(crate) fn add(&mut self, result: &TrimmingResult) {
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

    pub(crate) fn merge(&mut self, other: &TrimmingStats) {
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
pub(crate) fn trim_tail_sliding_window(qual: &[u8], window_size: usize, cutoff: u8) -> usize {
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
pub(crate) fn trim_front_sliding_window(qual: &[u8], window_size: usize, cutoff: u8) -> usize {
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
pub(crate) fn detect_poly_g_tail(seq: &[u8], min_len: usize) -> usize {
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
pub(crate) fn detect_poly_x_tail(seq: &[u8], min_len: usize) -> usize {
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
pub(crate) fn trim_read(seq: &[u8], qual: &[u8], config: &TrimmingConfig) -> TrimmingResult {
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

    // 3. Maximum length trimming (trim from tail if read is too long)
    if config.max_len > 0 {
        let current_len = result.end_pos - result.start_pos;
        if current_len > config.max_len {
            result.end_pos = result.start_pos + config.max_len;
        }
    }

    // Ensure we still have a valid range
    if result.start_pos >= result.end_pos {
        result.end_pos = result.start_pos;
        return result;
    }

    let current_qual = &qual[result.start_pos..result.end_pos];

    // 4. Sliding window front trimming
    if config.enable_trim_front {
        let trim_amount = trim_front_sliding_window(
            current_qual,
            config.cut_window_size,
            config.cut_mean_quality,
        );
        result.start_pos += trim_amount;
    }

    // 5. Sliding window tail trimming
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

    // 6. PolyG tail trimming (check polyG first as it's more specific)
    if config.enable_poly_g {
        let poly_g = detect_poly_g_tail(current_seq, config.poly_min_len);
        if poly_g > 0 {
            result.poly_g_trimmed = poly_g;
            result.end_pos = result.end_pos.saturating_sub(poly_g);
        }
    }

    // 7. PolyX tail trimming (only if polyG didn't already trim)
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
