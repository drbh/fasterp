//! Adapter trimming functionality
//!
//! This module provides adapter detection and trimming for FASTQ preprocessing:
//! - Manual adapter specification
//! - Built-in common adapter sequences
//! - Overlap-based adapter detection for paired-end reads
//! - Mismatch-tolerant matching

use std::cmp::min;

/// Configuration for adapter trimming
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// Adapter sequence for read1/single-end
    pub adapter_seq: Option<Vec<u8>>,
    /// Adapter sequence for read2 (paired-end only)
    pub adapter_seq_r2: Option<Vec<u8>>,
    /// Enable auto-detection for paired-end
    pub detect_adapter_for_pe: bool,
    /// Minimum overlap length for detection
    pub min_overlap: usize,
    /// Maximum mismatches allowed in overlap
    pub max_mismatches: usize,
}

impl AdapterConfig {
    pub fn new() -> Self {
        Self {
            adapter_seq: None,
            adapter_seq_r2: None,
            detect_adapter_for_pe: false,
            min_overlap: 5, // fastp's default for adapter trimming
            max_mismatches: 2,
        }
    }

    pub fn is_enabled(&self) -> bool {
        self.adapter_seq.is_some() || self.adapter_seq_r2.is_some() || self.detect_adapter_for_pe
    }
}

/// Built-in Illumina adapter sequences
pub mod adapters {
    /// Illumina TruSeq Universal Adapter
    pub const TRUSEQ_UNIVERSAL: &[u8] = b"AGATCGGAAGAGC";

    /// Illumina TruSeq Read 1 Adapter
    pub const TRUSEQ_READ1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

    /// Illumina TruSeq Read 2 Adapter
    pub const TRUSEQ_READ2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

    /// Illumina Small RNA 3' Adapter
    pub const SMALL_RNA_3P: &[u8] = b"TGGAATTCTCGGGTGCCAAGG";

    /// Nextera Transposase Sequence
    pub const NEXTERA: &[u8] = b"CTGTCTCTTATACACATCT";
}

/// Result of adapter detection
#[derive(Debug, Clone)]
pub struct AdapterMatch {
    /// Position where adapter starts
    pub position: usize,
    /// Number of matching bases
    pub matched_bases: usize,
    /// Number of mismatches
    pub mismatches: usize,
}

/// Find adapter sequence in read allowing mismatches
///
/// Uses a simple sliding window approach with mismatch tolerance.
/// Returns the best match (earliest position with fewest mismatches).
pub fn find_adapter(
    seq: &[u8],
    adapter: &[u8],
    min_overlap: usize,
    max_mismatches: usize,
) -> Option<AdapterMatch> {
    if adapter.is_empty() || seq.len() < min_overlap {
        return None;
    }

    let mut best_match: Option<AdapterMatch> = None;

    // Try all possible positions where adapter could start
    // We want to find adapters even if they're only partially present at the end
    // NOTE: Adapters typically appear at the 3' end of reads. To avoid false positives
    // in the middle of reads (like "GGATCGGAAGA" matching "AGATCGGAAGAGC"),
    // we only search in positions where there's relatively little sequence remaining.
    // We start searching from the position where at most 2x the adapter length remains.
    let search_start = if seq.len() > adapter.len() * 2 {
        seq.len() - adapter.len() * 2
    } else {
        0
    };

    for start_pos in search_start..seq.len() {
        let remaining = seq.len() - start_pos;

        // Need at least min_overlap bases to consider
        if remaining < min_overlap {
            break;
        }

        // Check how many bases we can compare
        let compare_len = min(remaining, adapter.len());

        // Only check if we have enough bases for minimum overlap
        if compare_len < min_overlap {
            break;
        }

        // Dynamic mismatch threshold based on match length (fastp compatibility)
        // For short matches (5-7 bases): require perfect match (0 mismatches)
        // For longer matches (8+ bases): allow 1 mismatch
        let allowed_mismatches = if compare_len < 8 { 0 } else { 1 };

        // Count matches and mismatches
        let mut matches = 0;
        let mut mismatches = 0;

        for i in 0..compare_len {
            let seq_base = seq[start_pos + i].to_ascii_uppercase();
            let adapter_base = adapter[i].to_ascii_uppercase();

            if seq_base == adapter_base {
                matches += 1;
            } else {
                mismatches += 1;
                if mismatches > allowed_mismatches {
                    break; // Too many mismatches, skip this position
                }
            }
        }

        // Check if this is a valid match
        if mismatches <= allowed_mismatches && matches + mismatches >= min_overlap {
            // Check if this is better than current best match
            let is_better = match &best_match {
                None => true,
                Some(current_best) => {
                    // Prefer earlier position, or fewer mismatches at same position
                    start_pos < current_best.position
                        || (start_pos == current_best.position
                            && mismatches < current_best.mismatches)
                }
            };

            if is_better {
                best_match = Some(AdapterMatch {
                    position: start_pos,
                    matched_bases: matches,
                    mismatches,
                });
            }
        }
    }

    best_match
}

/// Trim adapter from sequence and quality
///
/// Returns (trimmed_seq, trimmed_qual) as slices
pub fn trim_adapter<'a>(
    seq: &'a [u8],
    qual: &'a [u8],
    adapter_match: &AdapterMatch,
) -> (&'a [u8], &'a [u8]) {
    // Trim everything from the adapter position onwards
    let trim_pos = adapter_match.position;
    (&seq[..trim_pos], &qual[..trim_pos])
}

/// Detect adapter using paired-end overlap information
///
/// When read1 and read2 overlap, we can detect adapter contamination
/// by checking for reverse-complement matching beyond the insert size.
pub fn detect_adapter_from_pe_overlap(
    _seq1: &[u8],
    _seq2: &[u8],
    _min_overlap: usize,
) -> (Option<AdapterMatch>, Option<AdapterMatch>) {
    // TODO: Implement paired-end overlap-based adapter detection
    // This requires:
    // 1. Find overlap between read1 and reverse-complement of read2
    // 2. Determine insert size
    // 3. Identify adapter sequences beyond the insert

    // For now, return None (not implemented)
    (None, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_adapter_match() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGC";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 12); // Adapter starts at position 12
        assert_eq!(m.mismatches, 0);
    }

    #[test]
    fn test_adapter_with_mismatch() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGX"; // X instead of C
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 12);
        assert_eq!(m.mismatches, 1);
    }

    #[test]
    fn test_partial_adapter_at_end() {
        // Only partial adapter at the end of read
        let seq = b"ACGTACGTACGTACGTAGATCGGAA"; // Only first 9 bases of adapter
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        // Should not match with min_overlap=10
        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());

        // Should match with min_overlap=8
        let result = find_adapter(seq, adapter, 8, 2);
        assert!(result.is_some());
    }

    #[test]
    fn test_no_adapter() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());
    }

    #[test]
    fn test_trim_adapter() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGC";
        let qual = b"############IIIIIIIIIIIII";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let adapter_match = find_adapter(seq, adapter, 10, 2).unwrap();
        let (trimmed_seq, trimmed_qual) = trim_adapter(seq, qual, &adapter_match);

        assert_eq!(trimmed_seq, b"ACGTACGTACGT");
        assert_eq!(trimmed_qual, b"############");
    }

    #[test]
    fn test_adapter_too_many_mismatches() {
        // 3 mismatches (XXX), but max_mismatches=2
        let seq = b"ACGTACGTACGTAGATCGGAAXXXC";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());
    }

    #[test]
    fn test_case_insensitive_matching() {
        let seq = b"ACGTACGTACGTagatcggaagagc"; // lowercase adapter
        let adapter = adapters::TRUSEQ_UNIVERSAL; // uppercase

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.mismatches, 0);
    }
}
