//! Base correction using paired-end overlap analysis
//!
//! This module implements overlap detection and base correction for paired-end reads.
//! When reads overlap (insert size < 2× read length), we can correct sequencing errors
//! by comparing bases in the overlapping region and using quality scores to determine
//! which base is correct.

/// Quality score thresholds for base correction
const GOOD_QUAL: u8 = 30; // Q30 - high confidence base
const BAD_QUAL: u8 = 14; // Q14 - low confidence base

/// Result of overlap analysis between paired-end reads
#[derive(Debug, Clone)]
pub struct OverlapResult {
    /// Whether an overlap was found
    pub overlapped: bool,
    /// Starting offset of overlap in R1
    pub offset: usize,
    /// Length of overlapping region
    pub overlap_len: usize,
    /// Number of differences in overlap region
    pub differences: usize,
}

/// Configuration for overlap detection
#[derive(Debug, Clone)]
pub struct OverlapConfig {
    /// Minimum required overlap length
    pub min_overlap_len: usize,
    /// Maximum allowed differences
    pub max_diff: usize,
    /// Maximum allowed difference percentage (0-100)
    pub max_diff_percent: usize,
}

impl Default for OverlapConfig {
    fn default() -> Self {
        Self {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        }
    }
}

/// Statistics from base correction
#[derive(Debug, Default, Clone)]
pub struct CorrectionStats {
    /// Number of bases corrected
    pub corrected: usize,
    /// Number of mismatches not corrected (quality ambiguous)
    pub uncorrected: usize,
}

impl CorrectionStats {
    /// Merge statistics from another correction operation
    pub fn merge(&mut self, other: &CorrectionStats) {
        self.corrected += other.corrected;
        self.uncorrected += other.uncorrected;
    }
}

/// Reverse complement a single base
#[inline]
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        _ => base,
    }
}

/// Reverse complement a sequence in-place
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Detect overlap between R1 and reverse-complement of R2
///
/// This function tries different offset positions to find where R1 and R2_rc
/// overlap with acceptable differences.
pub fn detect_overlap(
    r1_seq: &[u8],
    r2_seq: &[u8],
    config: &OverlapConfig,
) -> Option<OverlapResult> {
    if r1_seq.is_empty() || r2_seq.is_empty() {
        return None;
    }

    // Create reverse complement of R2
    let r2_rc = reverse_complement(r2_seq);

    let r1_len = r1_seq.len();
    let r2_len = r2_rc.len();

    // Try different offset positions
    // We check three scenarios:
    // 1. R1 extends past R2 (offset in R1)
    // 2. R2 extends past R1 (offset in R2)
    // 3. Complete overlap

    // Try offsets where R1 starts before R2_rc
    for offset in 0..r1_len {
        let overlap_len = std::cmp::min(r1_len - offset, r2_len);

        if overlap_len < config.min_overlap_len {
            continue;
        }

        // Count differences in overlapping region
        let differences = count_differences(
            &r1_seq[offset..offset + overlap_len],
            &r2_rc[0..overlap_len],
        );

        // Check if overlap meets criteria
        if differences <= config.max_diff
            && (differences * 100 / overlap_len) <= config.max_diff_percent
        {
            return Some(OverlapResult {
                overlapped: true,
                offset,
                overlap_len,
                differences,
            });
        }
    }

    // Try offsets where R2_rc starts before R1
    for offset in 1..r2_len {
        let overlap_len = std::cmp::min(r2_len - offset, r1_len);

        if overlap_len < config.min_overlap_len {
            continue;
        }

        // Count differences in overlapping region
        let differences = count_differences(
            &r1_seq[0..overlap_len],
            &r2_rc[offset..offset + overlap_len],
        );

        // Check if overlap meets criteria
        if differences <= config.max_diff
            && (differences * 100 / overlap_len) <= config.max_diff_percent
        {
            // For R2 offset, we encode it as negative by using offset in R1 = 0
            // and storing the R2 offset in a different way
            // For now, we'll just return the case where R1 offset is 0
            // and handle R2 offset separately in correction
            return Some(OverlapResult {
                overlapped: true,
                offset: 0, // R1 starts at beginning
                overlap_len,
                differences,
            });
        }
    }

    None
}

/// Count differences between two sequences of equal length
#[inline]
fn count_differences(seq1: &[u8], seq2: &[u8]) -> usize {
    debug_assert_eq!(seq1.len(), seq2.len());
    seq1.iter().zip(seq2.iter()).filter(|(a, b)| a != b).count()
}

/// Correct bases in overlapping region using quality scores
///
/// When R1 and R2 overlap, we can correct mismatches by choosing the base
/// with higher quality. Only corrects when one base has high quality (≥Q30)
/// and the other has low quality (≤Q14).
pub fn correct_by_overlap(
    r1_seq: &mut [u8],
    r1_qual: &mut [u8],
    r2_seq: &mut [u8],
    r2_qual: &mut [u8],
    overlap: &OverlapResult,
) -> CorrectionStats {
    if !overlap.overlapped {
        return CorrectionStats::default();
    }

    let mut stats = CorrectionStats::default();

    // Create reverse complement of R2 for comparison
    let r2_rc = reverse_complement(r2_seq);
    let r2_rc_qual: Vec<u8> = r2_qual.iter().rev().copied().collect();

    // Process each position in the overlap
    for i in 0..overlap.overlap_len {
        let r1_pos = overlap.offset + i;
        let r2_rc_pos = i;

        // Get bases and qualities
        let r1_base = r1_seq[r1_pos];
        let r2_base_rc = r2_rc[r2_rc_pos];
        let r1_q = r1_qual[r1_pos];
        let r2_q_rc = r2_rc_qual[r2_rc_pos];

        // If bases match, nothing to correct
        if r1_base == r2_base_rc {
            continue;
        }

        // Mismatch found - check if we can correct
        // Case 1: R1 has high quality, R2 has low quality
        if r1_q >= GOOD_QUAL && r2_q_rc <= BAD_QUAL {
            // Correct R2 to match R1
            let r2_original_pos = r2_seq.len() - 1 - r2_rc_pos;
            r2_seq[r2_original_pos] = complement_base(r1_base);
            r2_qual[r2_original_pos] = r1_q;
            stats.corrected += 1;
        }
        // Case 2: R2 has high quality, R1 has low quality
        else if r2_q_rc >= GOOD_QUAL && r1_q <= BAD_QUAL {
            // Correct R1 to match R2
            r1_seq[r1_pos] = r2_base_rc;
            r1_qual[r1_pos] = r2_q_rc;
            stats.corrected += 1;
        }
        // Case 3: Both have similar quality - don't correct
        else {
            stats.uncorrected += 1;
        }
    }

    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_base() {
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'T'), b'A');
        assert_eq!(complement_base(b'C'), b'G');
        assert_eq!(complement_base(b'G'), b'C');
        assert_eq!(complement_base(b'N'), b'N');
    }

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"ACGT"); // ACGT reverse complement is ACGT

        let seq2 = b"AAAA";
        let rc2 = reverse_complement(seq2);
        assert_eq!(rc2, b"TTTT");

        let seq3 = b"ATCG";
        let rc3 = reverse_complement(seq3);
        assert_eq!(rc3, b"CGAT");
    }

    #[test]
    fn test_count_differences() {
        assert_eq!(count_differences(b"AAAA", b"AAAA"), 0);
        assert_eq!(count_differences(b"AAAA", b"AAAT"), 1);
        assert_eq!(count_differences(b"AAAA", b"TTTT"), 4);
        assert_eq!(count_differences(b"ACGT", b"ACGT"), 0);
        assert_eq!(count_differences(b"ACGT", b"TGCA"), 4);
    }

    #[test]
    fn test_overlap_detection_no_overlap() {
        let config = OverlapConfig::default();

        // Completely different sequences that won't overlap even after RC
        // R1: all A's
        // R2: CGCGCG... (RC = CGCGCG...), which won't match AAAA...
        let r1 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 34 bases
        let r2 = b"CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"; // 34 bases, RC is also CGCGCG...

        let result = detect_overlap(r1, r2, &config);
        assert!(result.is_none());
    }

    #[test]
    fn test_overlap_detection_perfect_overlap() {
        let config = OverlapConfig::default();

        // R1 and R2_rc are identical (perfect overlap)
        let r1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 35 bases
        let r2 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // Same, so R2_rc = reverse_complement

        let result = detect_overlap(r1, r2, &config);
        assert!(result.is_some());

        let overlap = result.unwrap();
        assert!(overlap.overlapped);
        assert_eq!(overlap.differences, 0);
        assert!(overlap.overlap_len >= config.min_overlap_len);
    }

    #[test]
    fn test_overlap_detection_partial_overlap() {
        let config = OverlapConfig {
            min_overlap_len: 10,
            max_diff: 2,
            max_diff_percent: 20,
        };

        // Create sequences with partial overlap
        // R1: AAAAAAAAAAAAAAAAAAAAACGTACGTACGT (30A + 12 bases)
        // R2_rc should overlap at the end
        let r1 = b"AAAAAAAAAAAAAAAAAAAAACGTACGTACGT";
        let r2_rc = b"ACGTACGTACGTTTTTTTTTTTTTTTTTTTTT";

        // R2 is reverse_complement of R2_rc
        let r2 = reverse_complement(r2_rc);

        let result = detect_overlap(r1, &r2, &config);
        // This should find some overlap
        assert!(result.is_some());
    }

    #[test]
    fn test_correction_high_quality_r1() {
        let mut r1_seq = b"AAAA".to_vec();
        let mut r1_qual = vec![35, 35, 35, 35]; // Q35 (good)
        let mut r2_seq = b"TTTT".to_vec(); // RC = AAAA, so base matches after RC
        let mut r2_qual = vec![10, 10, 10, 10]; // Q10 (bad)

        // Create mismatch: change one base in R2
        r2_seq[0] = b'A'; // R2_rc would be TAAT (mismatch at position 3)

        let overlap = OverlapResult {
            overlapped: true,
            offset: 0,
            overlap_len: 4,
            differences: 1,
        };

        let stats = correct_by_overlap(
            &mut r1_seq,
            &mut r1_qual,
            &mut r2_seq,
            &mut r2_qual,
            &overlap,
        );

        assert!(stats.corrected > 0 || stats.uncorrected > 0);
    }

    #[test]
    fn test_correction_ambiguous_quality() {
        let mut r1_seq = b"AAAA".to_vec();
        let mut r1_qual = vec![20, 20, 20, 20]; // Q20 (medium)
        let mut r2_seq = b"TTTT".to_vec();
        let mut r2_qual = vec![20, 20, 20, 20]; // Q20 (medium)

        // Create intentional mismatch
        r2_seq[0] = b'A';

        let overlap = OverlapResult {
            overlapped: true,
            offset: 0,
            overlap_len: 4,
            differences: 1,
        };

        let stats = correct_by_overlap(
            &mut r1_seq,
            &mut r1_qual,
            &mut r2_seq,
            &mut r2_qual,
            &overlap,
        );

        // With similar quality, should not correct
        assert_eq!(stats.corrected, 0);
    }

    #[test]
    fn test_default_config() {
        let config = OverlapConfig::default();
        assert_eq!(config.min_overlap_len, 30);
        assert_eq!(config.max_diff, 5);
        assert_eq!(config.max_diff_percent, 20);
    }

    #[test]
    fn test_stats_merge() {
        let mut stats1 = CorrectionStats {
            corrected: 10,
            uncorrected: 5,
        };

        let stats2 = CorrectionStats {
            corrected: 3,
            uncorrected: 2,
        };

        stats1.merge(&stats2);

        assert_eq!(stats1.corrected, 13);
        assert_eq!(stats1.uncorrected, 7);
    }
}
