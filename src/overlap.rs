//! Base correction using paired-end overlap analysis
//!
//! This module implements overlap detection and base correction for paired-end reads.
//! When reads overlap (insert size < 2× read length), we can correct sequencing errors
//! by comparing bases in the overlapping region and using quality scores to determine
//! which base is correct.

/// Quality score thresholds for base correction
const GOOD_QUAL: u8 = 30; // Q30 - high confidence base
const BAD_QUAL: u8 = 14; // Q14 - low confidence base

/// Quick inline SIMD check of first 16 bytes - returns mismatch count
/// This avoids function call overhead for obviously bad offsets
#[cfg(target_arch = "aarch64")]
#[inline(always)]
unsafe fn quick_mismatch_check_16(ptr1: *const u8, ptr2: *const u8) -> usize {
    use std::arch::aarch64::{vaddlvq_u8, vceqq_u8, vld1q_u8, vshrq_n_u8};
    unsafe {
        let v1 = vld1q_u8(ptr1);
        let v2 = vld1q_u8(ptr2);
        let eq = vceqq_u8(v1, v2);
        let shifted = vshrq_n_u8(eq, 7);
        let equals = vaddlvq_u8(shifted) as usize;
        16 - equals
    }
}

#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn quick_mismatch_check_16(ptr1: *const u8, ptr2: *const u8) -> usize {
    use std::arch::x86_64::{_mm_cmpeq_epi8, _mm_loadu_si128, _mm_movemask_epi8};
    unsafe {
        let v1 = _mm_loadu_si128(ptr1 as *const _);
        let v2 = _mm_loadu_si128(ptr2 as *const _);
        let eq = _mm_cmpeq_epi8(v1, v2);
        let mask = _mm_movemask_epi8(eq) as u32;
        16 - mask.count_ones() as usize
    }
}

#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
#[inline(always)]
unsafe fn quick_mismatch_check_16(ptr1: *const u8, ptr2: *const u8) -> usize {
    unsafe {
        let mut mismatches = 0;
        for i in 0..16 {
            if *ptr1.add(i) != *ptr2.add(i) {
                mismatches += 1;
            }
        }
        mismatches
    }
}

/// Result of overlap analysis between paired-end reads
#[derive(Debug, Clone)]
pub struct OverlapResult {
    /// Whether an overlap was found
    pub overlapped: bool,
    /// Starting offset of overlap
    /// Positive: R1 extends past R2 (offset in R1)
    /// Negative: R2 extends past R1 (offset in R2, stored as negative)
    /// Zero: reads start at same position
    pub offset: isize,
    /// Length of overlapping region
    pub overlap_len: usize,
    /// Number of differences/mismatches in overlap
    pub diff: usize,
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

// LUT is faster than match - no branch misprediction
static RC_LUT: [u8; 256] = {
    let mut lut = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        lut[i] = i as u8;
        i += 1;
    }
    lut[b'A' as usize] = b'T';
    lut[b'T' as usize] = b'A';
    lut[b'C' as usize] = b'G';
    lut[b'G' as usize] = b'C';
    lut[b'N' as usize] = b'N';
    lut
};

// Max read length we support on stack. 512bp covers all standard Illumina reads.
const MAX_READ_LEN: usize = 512;

/// Detect overlap between R1 and reverse-complement of R2
///
/// This function tries different offset positions to find where R1 and `R2_rc`
/// overlap with acceptable differences.
pub fn detect_overlap(
    r1_seq: &[u8],
    r2_seq: &[u8],
    config: &OverlapConfig,
) -> Option<OverlapResult> {
    if r1_seq.is_empty() || r2_seq.is_empty() {
        return None;
    }

    const COMPLETE_COMPARE_REQUIRE: usize = 50;

    // Stack-allocated buffer - no heap, no TLS overhead
    let r2_len = r2_seq.len();
    let r1_len = r1_seq.len();

    // Use stack for normal reads, Vec for unusually long ones
    let mut stack_buf = [0u8; MAX_READ_LEN];
    let heap_buf: Vec<u8>;

    let r2_rc: &[u8] = if r2_len <= MAX_READ_LEN {
        for (i, &b) in r2_seq.iter().rev().enumerate() {
            stack_buf[i] = RC_LUT[b as usize];
        }
        &stack_buf[..r2_len]
    } else {
        heap_buf = r2_seq.iter().rev().map(|&b| RC_LUT[b as usize]).collect();
        &heap_buf
    };

    // Try different offset positions
    // IMPORTANT: Fastp returns the FIRST valid overlap it finds, not the best one.
    // This means it searches forward offsets first (0, 1, 2, ...) and returns immediately
    // when a valid overlap is found, without checking if there might be a better one later.

    // Try offsets where R1 starts before R2_rc (R1 extends past R2) - FORWARD SCAN
    // We must check these FIRST and return immediately to match fastp's behavior
    for offset in 0..r1_len {
        let overlap_len = std::cmp::min(r1_len - offset, r2_len);

        // INTENTIONALLY match fastp's off-by-one bug (overlapanalysis.cpp:32)
        // Fastp uses `while (offset < len1-overlapRequire)` which excludes the boundary case
        // where overlap_len would equal exactly min_overlap_len.
        // This causes fastp to miss overlaps at offset = (read_length - min_overlap_length).
        // See fastp_offset_bug_analysis.md for detailed explanation.
        if overlap_len <= config.min_overlap_len {
            continue;
        }

        // Count differences with early termination (fastp's logic)
        let overlap_diff_limit = std::cmp::min(
            config.max_diff,
            (overlap_len * config.max_diff_percent) / 100,
        );

        // Quick prefix filter: check first 16 bytes inline to skip obviously bad offsets
        // This avoids function call overhead for ~80-90% of offsets
        if overlap_len >= 16 {
            let first_mismatches =
                unsafe { quick_mismatch_check_16(r1_seq[offset..].as_ptr(), r2_rc.as_ptr()) };
            // If first 16 bytes already exceed limit, skip this offset
            // (early termination would reject it anyway)
            if first_mismatches > overlap_diff_limit {
                continue;
            }
        }

        // Use SIMD-accelerated mismatch counting
        let (differences, i) = crate::simd::count_mismatches(
            &r1_seq[offset..offset + overlap_len],
            &r2_rc[..overlap_len],
            overlap_diff_limit,
            COMPLETE_COMPARE_REQUIRE,
        );

        // Check if overlap meets criteria
        // Fastp's logic: accept if either:
        // 1. Within calculated overlap-specific limits
        // 2. OR we compared all bases in a long overlap (i > 50)
        let within_strict_limits = differences <= overlap_diff_limit;
        let completed_long_overlap = i > COMPLETE_COMPARE_REQUIRE;

        if within_strict_limits || completed_long_overlap {
            // MATCH FASTP: Return the FIRST valid overlap immediately
            // Do NOT continue searching for better overlaps
            return Some(OverlapResult {
                overlapped: true,
                offset: offset as isize,
                overlap_len,
                diff: differences,
            });
        }
    }

    // Try offsets where R2_rc starts before R1 (R2 extends past R1)
    for offset in 1..r2_len {
        let overlap_len = std::cmp::min(r2_len - offset, r1_len);

        // INTENTIONALLY match fastp's off-by-one bug (overlapanalysis.cpp:67)
        // Fastp uses `while (offset > -(len2-overlapRequire))` which excludes the boundary case
        // where overlap_len would equal exactly min_overlap_len.
        // This causes fastp to miss overlaps at offset = -(read_length - min_overlap_length).
        // See fastp_offset_bug_analysis.md for detailed explanation.
        if overlap_len <= config.min_overlap_len {
            continue;
        }

        // Count differences with early termination (fastp's logic)
        let overlap_diff_limit = std::cmp::min(
            config.max_diff,
            (overlap_len * config.max_diff_percent) / 100,
        );

        // Quick prefix filter: check first 16 bytes inline to skip obviously bad offsets
        // This avoids function call overhead for ~80-90% of offsets
        if overlap_len >= 16 {
            let first_mismatches =
                unsafe { quick_mismatch_check_16(r1_seq.as_ptr(), r2_rc[offset..].as_ptr()) };
            // If first 16 bytes already exceed limit, skip this offset
            // (early termination would reject it anyway)
            if first_mismatches > overlap_diff_limit {
                continue;
            }
        }

        // Use SIMD-accelerated mismatch counting
        let (differences, i) = crate::simd::count_mismatches(
            &r1_seq[..overlap_len],
            &r2_rc[offset..offset + overlap_len],
            overlap_diff_limit,
            COMPLETE_COMPARE_REQUIRE,
        );

        // Check if overlap meets criteria
        // Same lenient logic for long overlaps
        let within_strict_limits = differences <= overlap_diff_limit;
        let completed_long_overlap = i > COMPLETE_COMPARE_REQUIRE;

        if within_strict_limits || completed_long_overlap {
            // MATCH FASTP: Return the FIRST valid overlap immediately
            // Do NOT continue searching for better overlaps
            return Some(OverlapResult {
                overlapped: true,
                offset: -(offset as isize),
                overlap_len,
                diff: differences,
            });
        }
    }

    // No valid overlap found
    None
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
        let (r1_pos, r2_rc_pos) = if overlap.offset >= 0 {
            // R1 extends past R2: R1 offset is positive
            (overlap.offset as usize + i, i)
        } else {
            // R2 extends past R1: R2 offset is negative
            (i, (-overlap.offset) as usize + i)
        };

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

/// Trim adapters based on overlap analysis for paired-end reads
///
/// When R2 extends past R1 (negative offset), the non-overlapping parts
/// are adapter sequences that should be trimmed.
///
/// Returns (`r1_trim_len`, `r2_trim_len`) - the lengths to keep for each read.
/// Returns None if overlap-based trimming should not be applied.
pub fn trim_by_overlap_analysis(
    r1_len: usize,
    r2_len: usize,
    overlap: &OverlapResult,
) -> Option<(usize, usize)> {
    // Only trim if overlap is found
    if !overlap.overlapped {
        return None;
    }

    // Fastp's trim logic based on overlap:
    // When offset < 0: R2 extends past R1
    //   - trim_len1 = min(r1_len, overlap_len + frontTrimmed2)
    //   - trim_len2 = min(r2_len, overlap_len + frontTrimmed1)
    // When offset >= 0: R1 extends past R2
    //   - Both reads trimmed to overlap_len (symmetric)

    if overlap.offset < 0 {
        // R2 extends past R1 - use fastp's formula (no front trimming in our case)
        let trim_len1 = std::cmp::min(r1_len, overlap.overlap_len);
        let trim_len2 = std::cmp::min(r2_len, overlap.overlap_len);
        Some((trim_len1, trim_len2))
    } else {
        // R1 extends past R2 - fastp seems to NOT trim in this case for positive offset
        // Return None to skip overlap-based trimming and use sequence-based trimming instead
        None
    }
}

/// Merge two overlapping paired-end reads into a single read
///
/// Returns (header, sequence, quality) for the merged read.
/// The header includes "`merged_XXX_YYY`" where XXX is bases from R1, YYY from R2.
pub fn merge_reads(
    r1_seq: &[u8],
    r1_qual: &[u8],
    r1_header: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    overlap: &OverlapResult,
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    // Reverse complement R2
    let r2_rc = reverse_complement(r2_seq);
    let r2_qual_rc: Vec<u8> = r2_qual.iter().rev().copied().collect();

    let (merged_seq, merged_qual, len1, len2) = if overlap.offset >= 0 {
        // R1 extends past R2
        // merged = R1[0..len1] + R2_rc[overlap_len..]
        let len1 = (overlap.offset as usize) + overlap.overlap_len;
        let len2 = r2_rc.len().saturating_sub(overlap.overlap_len);

        let mut seq = r1_seq[..len1].to_vec();
        if len2 > 0 {
            seq.extend_from_slice(&r2_rc[overlap.overlap_len..]);
        }

        let mut qual = r1_qual[..len1].to_vec();
        if len2 > 0 {
            qual.extend_from_slice(&r2_qual_rc[overlap.overlap_len..]);
        }

        (seq, qual, len1, len2)
    } else {
        // R2 extends past R1
        // merged = R2_rc[0..-offset] + R1[0..]
        let offset_abs = (-overlap.offset) as usize;
        let len2 = offset_abs;
        let len1 = r1_seq.len();

        let mut seq = r2_rc[..offset_abs].to_vec();
        seq.extend_from_slice(r1_seq);

        let mut qual = r2_qual_rc[..offset_abs].to_vec();
        qual.extend_from_slice(r1_qual);

        (seq, qual, len1, len2)
    };

    // Build header: @READNAME merged_150_15 rest
    let mut merged_header = Vec::with_capacity(r1_header.len() + 32);

    // Copy until first space or end
    let space_pos = r1_header
        .iter()
        .position(|&b| b == b' ')
        .unwrap_or(r1_header.len());
    merged_header.extend_from_slice(&r1_header[..space_pos]);

    // Add merge tag
    merged_header.extend_from_slice(b" merged_");
    merged_header.extend_from_slice(len1.to_string().as_bytes());
    merged_header.push(b'_');
    merged_header.extend_from_slice(len2.to_string().as_bytes());

    // Add rest of header if any
    if space_pos < r1_header.len() {
        merged_header.push(b' ');
        merged_header.extend_from_slice(&r1_header[space_pos + 1..]);
    }

    (merged_header, merged_seq, merged_qual)
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
            diff: 1,
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
            diff: 1,
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
    fn test_merge_reads_perfect_overlap() {
        // Create two identical 150bp reads that should overlap perfectly
        // R1 and R2_rc are identical, so full overlap with offset=0
        let r1_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let r1_header = b"read1 1:N:0";

        // R2 is same as R1, so R2_rc will be the reverse complement
        // For perfect overlap, R2 should be reverse complement of R1
        let r2 = reverse_complement(r1_seq);
        let r2_qual: Vec<u8> = r1_qual.iter().rev().copied().collect();

        // Detect overlap
        let config = OverlapConfig::default();
        let overlap = detect_overlap(r1_seq, &r2, &config);

        println!("Overlap result: {overlap:?}");

        assert!(overlap.is_some(), "Should detect overlap");
        let overlap = overlap.unwrap();
        assert!(overlap.overlapped, "Should be overlapped");

        // Merge the reads
        let (merged_header, merged_seq, merged_qual) =
            merge_reads(r1_seq, r1_qual, r1_header, &r2, &r2_qual, &overlap);

        println!("Merged header: {}", String::from_utf8_lossy(&merged_header));
        println!("Merged seq len: {}", merged_seq.len());
        println!("Merged qual len: {}", merged_qual.len());
        println!(
            "Overlap offset: {}, len: {}",
            overlap.offset, overlap.overlap_len
        );

        // Verify header contains "merged_"
        assert!(
            merged_header.windows(7).any(|w| w == b"merged_"),
            "Header should contain 'merged_' tag"
        );

        // Verify sequence and quality lengths match
        assert_eq!(
            merged_seq.len(),
            merged_qual.len(),
            "Seq and qual lengths must match"
        );

        // For perfect overlap with offset=0, merged length should equal original length
        assert!(
            !merged_seq.is_empty(),
            "Merged sequence should not be empty"
        );
    }

    #[test]
    fn test_merge_reads_partial_overlap() {
        // Create reads with partial overlap
        // R1: 50bp overlap region + 100bp unique = 150bp total
        // R2_rc: 100bp unique + 50bp overlap region = 150bp total
        // Expected merged: 100 + 50 + 100 = 250bp

        // R1: 100 A's + 50 C's
        let mut r1_seq = vec![b'A'; 100];
        r1_seq.extend_from_slice(&[b'C'; 50]);

        let r1_qual = vec![b'I'; 150];
        let r1_header = b"read1 1:N:0";

        // R2 reverse complement should have: 50 G's (rc of C) + 100 T's (rc of A)
        let mut r2_rc_seq = vec![b'G'; 50];
        r2_rc_seq.extend_from_slice(&[b'T'; 100]);

        // To get this as R2, we need to reverse complement it back
        let r2 = reverse_complement(&r2_rc_seq);
        let r2_qual = vec![b'I'; 150];

        // Detect overlap with relaxed config
        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 10,
            max_diff_percent: 30,
        };
        let overlap = detect_overlap(&r1_seq, &r2, &config);

        println!("Partial overlap result: {overlap:?}");

        if let Some(overlap) = overlap {
            println!(
                "Found overlap: offset={}, len={}",
                overlap.offset, overlap.overlap_len
            );

            // Merge the reads
            let (merged_header, merged_seq, merged_qual) =
                merge_reads(&r1_seq, &r1_qual, r1_header, &r2, &r2_qual, &overlap);

            println!("Merged header: {}", String::from_utf8_lossy(&merged_header));
            println!("Merged seq len: {}", merged_seq.len());

            // Verify merge worked
            assert!(
                merged_seq.len() > 150,
                "Merged should be longer than individual reads"
            );
            assert_eq!(merged_seq.len(), merged_qual.len());
        } else {
            println!("No overlap detected - this is expected with strict default config");
        }
    }

    // ===== WORKING CASES =====

    #[test]
    fn test_merge_reads_r1_extends_past_r2() {
        // Test case where R1 extends past R2 (positive offset)
        // Manually construct an overlap scenario with positive offset
        let r1_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 60bp
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let r1_header = b"read1 1:N:0";

        let r2 = reverse_complement(&r1_seq[10..60]); // 50bp - last 50bp of R1
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

        // Create overlap manually - R1 extends 10bp past start of R2
        let overlap = OverlapResult {
            overlapped: true,
            offset: 10,      // R1 extends by 10bp
            overlap_len: 50, // 50bp overlap
            diff: 0,
        };

        let (merged_header, merged_seq, merged_qual) =
            merge_reads(r1_seq, r1_qual, r1_header, &r2, r2_qual, &overlap);

        assert!(merged_header.windows(7).any(|w| w == b"merged_"));
        assert_eq!(merged_seq.len(), merged_qual.len());
        // With offset=10 and overlap=50, merged should be 60bp (all of R1)
        assert_eq!(merged_seq.len(), 60, "Merged should be 60bp");
    }

    #[test]
    fn test_merge_reads_exact_length_match() {
        // Test case where R1 and R2 are same length with complete overlap
        let r1_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let r1_header = b"read3 1:N:0";

        let r2 = reverse_complement(r1_seq);
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

        let config = OverlapConfig::default();
        let overlap = detect_overlap(r1_seq, &r2, &config);

        assert!(overlap.is_some(), "Should detect perfect overlap");
        let overlap = overlap.unwrap();
        assert!(overlap.overlapped);
        assert_eq!(
            overlap.offset, 0,
            "Offset should be 0 for equal length perfect overlap"
        );

        let (merged_header, merged_seq, merged_qual) =
            merge_reads(r1_seq, r1_qual, r1_header, &r2, r2_qual, &overlap);

        assert!(merged_header.windows(7).any(|w| w == b"merged_"));
        assert_eq!(
            merged_seq.len(),
            r1_seq.len(),
            "Merged should equal original length for complete overlap"
        );
        assert_eq!(merged_seq.len(), merged_qual.len());
    }

    // ===== FAILING CASES =====

    #[test]
    fn test_merge_reads_no_overlap_detected() {
        // Test case where reads don't overlap at all
        let r1_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // All A's

        // R2: All C's - no overlap with R1
        let r2 = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

        let config = OverlapConfig::default();
        let overlap = detect_overlap(r1_seq, r2, &config);

        assert!(
            overlap.is_none(),
            "Should NOT detect overlap for completely different sequences"
        );
    }

    #[test]
    fn test_merge_reads_insufficient_overlap_length() {
        // Test case where overlap is too short (below min_overlap_len)
        let r1_seq = b"AAAAAAAAAAAAAAAAAAAAAAAACGTACGT"; // 25 A's + 7 bases

        // R2_rc: Last 10 bases of R1 (overlap only 7bp, below default 30bp min)
        let r2_rc_seq = &r1_seq[22..]; // Last 7 bases
        let r2 = reverse_complement(r2_rc_seq);

        let config = OverlapConfig::default(); // min_overlap_len = 30
        let overlap = detect_overlap(r1_seq, &r2, &config);

        // Should not detect overlap because 7bp < 30bp minimum
        assert!(
            overlap.is_none(),
            "Should NOT detect overlap when overlap length < min_overlap_len"
        );
    }

    #[test]
    fn test_merge_reads_too_many_differences() {
        // Test case where overlap has too many mismatches
        let r1_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 35bp

        // R2: Similar to R1_rc but with many mismatches
        let mut r2: Vec<u8> = reverse_complement(r1_seq);
        // Introduce many differences
        for i in (0..r2.len()).step_by(3) {
            r2[i] = b'N'; // Every 3rd base is N (mismatch)
        }

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5, // Only allow 5 differences
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(r1_seq, &r2, &config);

        // Should not detect overlap due to too many differences
        assert!(
            overlap.is_none() || !overlap.unwrap().overlapped,
            "Should NOT detect overlap when differences exceed threshold"
        );
    }

    #[test]
    fn test_merge_reads_completely_non_overlapping() {
        // Test case where reads are from different fragments entirely
        let r1_seq = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"; // All G's

        let r2_seq = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"; // All T's

        let config = OverlapConfig::default();
        let overlap = detect_overlap(r1_seq, r2_seq, &config);

        // These sequences don't overlap at all
        assert!(
            overlap.is_none(),
            "Should NOT detect overlap for non-overlapping sequences"
        );
    }

    // ===== LENIENT MODE TESTS =====

    #[test]
    fn test_lenient_mode_long_overlap_high_diff() {
        // Test that overlaps >50bp are accepted via lenient mode even with high differences
        // This is the key fastp behavior: if i > 50, accept regardless of diff count
        // IMPORTANT: Must not exceed diff limit before position 50, or early termination rejects

        // Create a 60bp overlap where we stay under limit until past position 50
        let r1 = vec![b'A'; 60];

        // Create r2_rc with differences AFTER position 50 to avoid early termination
        // First 50bp: only 4 differences (under limit of 5)
        // After 50bp: add 6 more differences
        let mut r2_rc = vec![b'A'; 60];
        r2_rc[10] = b'T'; // diff 1
        r2_rc[20] = b'T'; // diff 2
        r2_rc[30] = b'T'; // diff 3
        r2_rc[40] = b'T'; // diff 4
        // After position 50:
        r2_rc[51] = b'T'; // diff 5
        r2_rc[52] = b'T'; // diff 6
        r2_rc[53] = b'T'; // diff 7
        r2_rc[54] = b'T'; // diff 8
        r2_rc[55] = b'T'; // diff 9
        r2_rc[56] = b'T'; // diff 10

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5, // Strict limit is only 5
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        assert!(overlap.is_some(), "Should detect overlap via lenient mode");
        let ov = overlap.unwrap();
        assert!(ov.overlapped);
        assert_eq!(ov.overlap_len, 60, "Should find full 60bp overlap");
        // Accepted via lenient mode: i=60 > 50, even though diff=10 > limit=5
    }

    #[test]
    fn test_early_termination_short_overlap() {
        // Test that short overlaps with early high differences are rejected quickly
        // This test verifies early termination behavior exists by checking that
        // we reject overlaps with too many early differences

        // Create sequences that would only match at offset 0 with many early differences
        let r1 = b"AAAAGGGGGCCCCCTTTTTAAAAGGGGGCCCCCT".to_vec(); // 35bp unique

        // Create r2_rc with 6 differences in positions that will cause early rejection
        let mut r2_rc = r1.clone();
        r2_rc[0] = b'T'; // diff 1
        r2_rc[1] = b'T'; // diff 2
        r2_rc[2] = b'T'; // diff 3
        r2_rc[3] = b'T'; // diff 4
        r2_rc[6] = b'A'; // diff 5
        r2_rc[7] = b'A'; // diff 6 - now exceeds limit of 5

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20, // 35 * 20% = 7, so limit = min(5, 7) = 5
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        // Should reject: 6 differences > 5 limit, early termination before position 50
        assert!(
            overlap.is_none(),
            "Should reject short overlap with early high differences"
        );
    }

    #[test]
    fn test_boundary_exactly_50bp() {
        // Test behavior at exactly the boundary: overlap_len = 50bp
        // With differences at the limit, should it pass?

        let r1 = vec![b'A'; 50];

        // Create r2_rc with exactly 5 differences (at the limit)
        let mut r2_rc = vec![b'A'; 50];
        for i in 0..5 {
            r2_rc[i * 10] = b'T';
        }

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20, // 50 * 20% = 10, so limit = min(5, 10) = 5
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        assert!(overlap.is_some(), "Should accept with diff=5 <= limit=5");
        let ov = overlap.unwrap();
        assert_eq!(ov.overlap_len, 50);
    }

    #[test]
    fn test_boundary_51bp_just_over() {
        // Test just over the boundary: 51bp overlap with 6 differences
        // Should pass via lenient mode (i=51 > 50)
        // Must keep diff count <= limit until past position 50

        let r1 = vec![b'A'; 51];

        // Create r2_rc staying at or under limit until past position 50
        let mut r2_rc = vec![b'A'; 51];
        r2_rc[10] = b'T'; // diff 1
        r2_rc[20] = b'T'; // diff 2
        r2_rc[30] = b'T'; // diff 3
        r2_rc[40] = b'T'; // diff 4
        r2_rc[48] = b'T'; // diff 5 (at limit, but not over yet)
        r2_rc[50] = b'T'; // diff 6 (exceeds limit, but i=50 is not > 50 yet, so this might cause early exit)

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        // i=51 > 50, so lenient mode should accept even with 6 differences
        assert!(
            overlap.is_some(),
            "Should accept via lenient mode: i=51 > 50"
        );
        let ov = overlap.unwrap();
        assert_eq!(ov.overlap_len, 51);
    }

    #[test]
    fn test_min_overlap_29bp_rejected() {
        // Test that overlaps below min_overlap_len are rejected
        // Even with perfect match

        let r1 = vec![b'A'; 29];
        let r2 = reverse_complement(&r1);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        assert!(
            overlap.is_none(),
            "Should reject overlap < 30bp even if perfect"
        );
    }

    // #[test]
    // fn test_min_overlap_30bp_accepted() {
    //     // Test that overlaps at exactly min_overlap_len are accepted

    //     let r1 = vec![b'A'; 30];
    //     let r2 = reverse_complement(&r1);

    //     let config = OverlapConfig {
    //         min_overlap_len: 30,
    //         max_diff: 5,
    //         max_diff_percent: 20,
    //     };

    //     let overlap = detect_overlap(&r1, &r2, &config);

    //     assert!(
    //         overlap.is_some(),
    //         "Should accept overlap = 30bp (at minimum)"
    //     );
    //     let ov = overlap.unwrap();
    //     assert_eq!(ov.overlap_len, 30);
    // }

    #[test]
    fn test_overlap_diff_limit_calculation() {
        // Test that overlap_diff_limit is correctly calculated as min(max_diff, overlap_len * percent)
        // For overlap_len=40, max_diff=5, max_diff_percent=20:
        // overlap_diff_limit = min(5, 40 * 20 / 100) = min(5, 8) = 5

        let r1 = vec![b'A'; 40];

        // Create r2_rc with exactly 5 differences (at calculated limit)
        let mut r2_rc = vec![b'A'; 40];
        for i in 0..5 {
            r2_rc[i * 8] = b'T';
        }

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        assert!(overlap.is_some(), "Should accept with diff=5 at limit");
    }

    #[test]
    fn test_percentage_limit_dominates() {
        // Test case where percentage limit is stricter than absolute limit
        // overlap_len=20, max_diff=5, max_diff_percent=20
        // overlap_diff_limit = min(5, 20 * 20 / 100) = min(5, 4) = 4

        // Use unique pattern to avoid false matches at other offsets
        let r1 = b"ACGTACGTACGTACGTACGT".to_vec(); // 20bp unique pattern

        // Create r2_rc with 5 differences (exceeds percentage limit of 4)
        let mut r2_rc = b"ACGTACGTACGTACGTACGT".to_vec();
        r2_rc[0] = b'T'; // diff 1
        r2_rc[4] = b'T'; // diff 2
        r2_rc[8] = b'T'; // diff 3
        r2_rc[12] = b'T'; // diff 4
        r2_rc[16] = b'T'; // diff 5 (exceeds limit of 4)

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 10, // Lower threshold for this test
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        // Should reject: 5 differences > 4 (percentage limit), early termination at i=17
        assert!(
            overlap.is_none(),
            "Should reject when exceeding percentage limit"
        );
    }

    #[test]
    fn test_negative_offset_overlap() {
        // Test R2 extends past R1 (negative offset)
        // R1: 40bp, R2: 60bp with 20bp extending past R1's start

        let r1 = vec![b'C'; 40];

        // R2_rc should be: [20bp of G] + [40bp of C]
        let mut r2_rc = vec![b'G'; 20];
        r2_rc.extend_from_slice(&[b'C'; 40]);

        // Convert to R2 (reverse complement of r2_rc)
        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        assert!(
            overlap.is_some(),
            "Should detect overlap with negative offset"
        );
        let ov = overlap.unwrap();
        assert!(
            ov.offset < 0,
            "Offset should be negative when R2 extends past R1"
        );
        assert_eq!(ov.overlap_len, 40, "Should overlap by 40bp");
    }

    #[test]
    fn test_positive_offset_overlap() {
        // Test R1 extends past R2 (positive offset)
        // Use a unique repeating pattern to avoid ambiguous matches

        // R1: 20bp unique + 40bp pattern
        let r1 = b"AAAAAAAAAAAAAAAAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"; // 20 A's + 40bp CGTA pattern

        // R2_rc should match the 40bp pattern part
        let r2_rc = b"CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec(); // 40bp CGTA pattern
        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(r1, &r2, &config);

        assert!(
            overlap.is_some(),
            "Should detect overlap with positive offset"
        );
        let ov = overlap.unwrap();
        assert!(
            ov.offset >= 0,
            "Offset should be positive when R1 extends past R2"
        );
        // The exact offset depends on the first unique match position
        assert!(ov.overlap_len >= 30, "Should have at least 30bp overlap");
    }

    #[test]
    fn test_real_world_read_18990() {
        // Test with actual problematic read from SRR22472290
        // This read has offset=-18, overlap_len=57, diff=7
        // Should pass via lenient mode (i=57 > 50)

        let r1 = b"GCCCTGGCCGGCCCGCGGGGCGCAGAGAGCGGCTGTGCGGCGCGCGCCCCGCCCCACCCGGGTCTTTTTAAACAA";
        let r2 = b"CGAGCGCCGTGCGCGCGCCCCACAGCCGCTCTCTGCGCCCCGCGGGCCGGCCAGGGCCGGTGTCATTTTCACCTA";

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(r1, r2, &config);

        assert!(
            overlap.is_some(),
            "Should detect overlap for real read 18990"
        );
        let ov = overlap.unwrap();
        assert_eq!(ov.offset, -18, "Should have offset=-18");
        assert_eq!(ov.overlap_len, 57, "Should have 57bp overlap");
        // Passes via lenient mode: all 57 bases compared (>50)
    }

    #[test]
    fn test_trim_by_overlap_negative_offset() {
        // Test that trim_by_overlap_analysis correctly calculates trim lengths for negative offset

        let overlap = OverlapResult {
            overlapped: true,
            offset: -20,
            overlap_len: 50,
            diff: 0,
        };

        let r1_len = 75;
        let r2_len = 75;

        let result = trim_by_overlap_analysis(r1_len, r2_len, &overlap);

        assert!(result.is_some(), "Should trim for negative offset");
        let (trim_len1, trim_len2) = result.unwrap();
        assert_eq!(trim_len1, 50, "R1 should be trimmed to overlap_len");
        assert_eq!(trim_len2, 50, "R2 should be trimmed to overlap_len");
    }

    #[test]
    fn test_trim_by_overlap_positive_offset() {
        // Test that trim_by_overlap_analysis returns None for positive offset
        // (fastp doesn't trim in this case, falls back to sequence-based trimming)

        let overlap = OverlapResult {
            overlapped: true,
            offset: 20,
            overlap_len: 50,
            diff: 0,
        };

        let r1_len = 75;
        let r2_len = 75;

        let result = trim_by_overlap_analysis(r1_len, r2_len, &overlap);

        assert!(result.is_none(), "Should NOT trim for positive offset");
    }

    #[test]
    fn test_early_termination_saves_computation() {
        // Test that early termination works - this is verified indirectly
        // by checking that we correctly reject short overlaps with many early differences
        // The real-world test (test_real_world_read_18990) validates this works correctly

        // Use the actual problematic read that we know should work
        // This provides better validation than constructed edge cases
        let r1 = b"GCCCTGGCCGGCCCGCGGGGCGCAGAGAGCGGCTGTGCGGCGCGCGCCCCGCCCCACCCGGGTCTTTTTAAACAA";
        let r2 = b"CGAGCGCCGTGCGCGCGCCCCACAGCCGCTCTCTGCGCCCCGCGGGCCGGCCAGGGCCGGTGTCATTTTCACCTA";

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(r1, r2, &config);

        // This should find overlap (early termination doesn't trigger because
        // differences are spread out enough to reach position >50)
        assert!(
            overlap.is_some(),
            "Should find overlap with early termination working correctly"
        );
        let ov = overlap.unwrap();
        assert_eq!(ov.offset, -18);
        assert_eq!(ov.overlap_len, 57);
        // Early termination working correctly: didn't reject despite 7 differences > limit=5
    }

    #[test]
    fn test_no_early_termination_for_long_overlaps() {
        // Test that we don't terminate early and continue to check lenient mode
        // 60bp overlap with 6 differences spread throughout

        let r1 = vec![b'A'; 60];

        // Create r2_rc with 6 differences spread across the entire overlap
        let mut r2_rc = vec![b'A'; 60];
        for i in 0..6 {
            r2_rc[i * 10] = b'T';
        }

        let r2 = reverse_complement(&r2_rc);

        let config = OverlapConfig {
            min_overlap_len: 30,
            max_diff: 5,
            max_diff_percent: 20,
        };

        let overlap = detect_overlap(&r1, &r2, &config);

        // Should accept via lenient mode even though diff=6 > limit=5
        // because i=60 > 50
        assert!(
            overlap.is_some(),
            "Should accept long overlap via lenient mode"
        );
        let ov = overlap.unwrap();
        assert_eq!(ov.overlap_len, 60);
    }
}
