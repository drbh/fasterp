/// SIMD-accelerated stats computation
///
/// This module provides SIMD implementations for the hot path operations:
/// - Quality score sum
/// - Q20/Q30 counting
/// - N base counting
/// - GC content counting
///
/// Supports:
/// - x86_64: AVX2 (256-bit, 32 bytes per iteration)
/// - aarch64: NEON (128-bit, 16 bytes per iteration)
///
/// All functions have scalar fallbacks and use runtime CPU feature detection.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::{
    uint8x16_t, vaddlvq_u8, vceqq_u8, vcgeq_u8, vcltq_u8, vdupq_n_u8, vld1q_u8, vorrq_u8,
    vshrq_n_u8, vsubq_u8,
};

/// Compute multiple stats in a single SIMD pass over sequence and quality data
pub struct Stats {
    pub qsum: u32,
    pub q20: usize,
    pub q30: usize,
    pub q40: usize,
    pub ncnt: usize,
    pub gc: usize,
    pub unqualified: usize, // Count of bases below qual_threshold
}

/// Check if SIMD acceleration is available on this platform
#[inline]
pub fn is_simd_available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2")
    }
    #[cfg(target_arch = "aarch64")]
    {
        // NEON is always available on aarch64
        true
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        false
    }
}

/// Count mismatches with early termination
///
/// Returns (`total_differences`, `positions_compared`)
/// Stops early if differences exceeds `max_diff` and `positions_compared` < `min_complete`
///
/// Uses SIMD acceleration where available, with careful handling to maintain
/// exact byte-level early termination semantics required by fastp compatibility.
#[inline]
pub fn count_mismatches(
    seq1: &[u8],
    seq2: &[u8],
    max_diff: usize,
    min_complete: usize,
) -> (usize, usize) {
    debug_assert_eq!(seq1.len(), seq2.len());

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { count_mismatches_avx2(seq1, seq2, max_diff, min_complete) };
        }
        return count_mismatches_scalar(seq1, seq2, max_diff, min_complete);
    }

    #[cfg(target_arch = "aarch64")]
    {
        return unsafe { count_mismatches_neon(seq1, seq2, max_diff, min_complete) };
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        count_mismatches_scalar(seq1, seq2, max_diff, min_complete)
    }
}

/// Scalar mismatch counting with early termination
#[inline]
fn count_mismatches_scalar(
    seq1: &[u8],
    seq2: &[u8],
    max_diff: usize,
    min_complete: usize,
) -> (usize, usize) {
    let mut differences = 0;
    let len = seq1.len();

    for i in 0..len {
        if seq1[i] != seq2[i] {
            differences += 1;
            // Early termination: too many differences before min_complete
            if differences > max_diff && i < min_complete {
                return (differences, i + 1);
            }
        }
    }

    (differences, len)
}

/// Compute all stats for a sequence/quality pair using SIMD when available
///
/// `qual_threshold`: Phred quality threshold (e.g., 15). Bases with quality < threshold are counted as unqualified.
/// Set to 0 to disable unqualified counting.
#[inline]
pub fn compute_stats(seq: &[u8], qual: &[u8], qual_threshold: u8) -> Stats {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { compute_stats_avx2(seq, qual, qual_threshold) };
        }
        // Scalar fallback for x86_64 without AVX2
        return compute_stats_scalar(seq, qual, qual_threshold);
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Always use NEON on aarch64 (available on all ARM64 platforms)
        unsafe { compute_stats_neon(seq, qual, qual_threshold) }
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        // Scalar fallback for other architectures
        compute_stats_scalar(seq, qual, qual_threshold)
    }
}

/// Scalar implementation (fallback)
/// Used on `x86_64` without AVX2 and non-x86_64/aarch64 platforms
#[inline]
#[allow(dead_code)]
fn compute_stats_scalar(seq: &[u8], qual: &[u8], qual_threshold: u8) -> Stats {
    let mut qsum = 0u32;
    let mut q20 = 0usize;
    let mut q30 = 0usize;
    let mut q40 = 0usize;
    let mut ncnt = 0usize;
    let mut gc = 0usize;
    let mut unqualified = 0usize;

    let qual_threshold_ascii = qual_threshold + 33; // Convert to ASCII

    for (&b, &q) in seq.iter().zip(qual) {
        let quality_u32 = u32::from(q - 33);
        qsum += quality_u32;

        // Q20/Q30/Q40: quality thresholds
        if q >= 53 {
            q20 += 1;
        } // Phred 20 = ASCII 53
        if q >= 63 {
            q30 += 1;
        } // Phred 30 = ASCII 63
        if q >= 73 {
            q40 += 1;
        } // Phred 40 = ASCII 73

        // Unqualified: below qual_threshold
        if qual_threshold > 0 && q < qual_threshold_ascii {
            unqualified += 1;
        }

        // N counting
        if b == b'N' || b == b'n' {
            ncnt += 1;
        }

        // GC counting
        if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
            gc += 1;
        }
    }

    Stats {
        qsum,
        q20,
        q30,
        q40,
        ncnt,
        gc,
        unqualified,
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(unsafe_op_in_unsafe_fn)]
unsafe fn count_mismatches_avx2(
    seq1: &[u8],
    seq2: &[u8],
    max_diff: usize,
    min_complete: usize,
) -> (usize, usize) {
    let len = seq1.len();
    let mut differences = 0usize;
    let mut i = 0usize;
    let chunk_size = 32; // AVX2 processes 32 bytes (256-bit)

    // Process 32 bytes at a time with AVX2
    while i + chunk_size <= len {
        // Load both sequences
        let vec1 = unsafe { _mm256_loadu_si256(seq1[i..].as_ptr() as *const __m256i) };
        let vec2 = unsafe { _mm256_loadu_si256(seq2[i..].as_ptr() as *const __m256i) };

        // Compare: returns 0xFF for equal, 0x00 for not equal
        let eq_mask = unsafe { _mm256_cmpeq_epi8(vec1, vec2) };

        // Convert to movemask (1 bit per byte, 1 = equal)
        let mask = unsafe { _mm256_movemask_epi8(eq_mask) } as u32;

        // Count mismatches: count zeros in the mask
        let matches = mask.count_ones() as usize;
        let chunk_mismatches = chunk_size - matches;
        let new_differences = differences + chunk_mismatches;
        let chunk_end = i + chunk_size;

        // Early termination check - only if threshold exceeded
        if new_differences > max_diff && i < min_complete {
            // Need exact position where threshold was exceeded
            for j in 0..chunk_size {
                let pos = i + j;
                if seq1[pos] != seq2[pos] {
                    differences += 1;
                    if differences > max_diff && pos < min_complete {
                        return (differences, pos + 1);
                    }
                }
            }
            // Threshold crossed at pos >= min_complete, continue with updated differences
        }

        differences = new_differences;
        i += chunk_size;
    }

    // Scalar remainder
    while i < len {
        if seq1[i] != seq2[i] {
            differences += 1;
            if differences > max_diff && i < min_complete {
                return (differences, i + 1);
            }
        }
        i += 1;
    }

    (differences, len)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(unsafe_op_in_unsafe_fn)]
unsafe fn compute_stats_avx2(seq: &[u8], qual: &[u8], qual_threshold: u8) -> Stats {
    let len = seq.len().min(qual.len());
    let mut qsum = 0u32;
    let mut q20 = 0usize;
    let mut q30 = 0usize;
    let mut q40 = 0usize;
    let mut ncnt = 0usize;
    let mut gc = 0usize;
    let mut unqualified = 0usize;

    let offset = _mm256_set1_epi8(33);
    let q20_thresh = _mm256_set1_epi8(53);
    let q30_thresh = _mm256_set1_epi8(63);
    let q40_thresh = _mm256_set1_epi8(73);
    let qual_threshold_ascii = qual_threshold + 33;
    let qual_thresh_vec = _mm256_set1_epi8(qual_threshold_ascii as i8);

    let n_upper = _mm256_set1_epi8(b'N' as i8);
    let n_lower = _mm256_set1_epi8(b'n' as i8);
    let g_upper = _mm256_set1_epi8(b'G' as i8);
    let g_lower = _mm256_set1_epi8(b'g' as i8);
    let c_upper = _mm256_set1_epi8(b'C' as i8);
    let c_lower = _mm256_set1_epi8(b'c' as i8);

    let mut i = 0;
    let chunk_size = 32;

    // Process 32 bytes at a time
    while i + chunk_size <= len {
        // Load quality scores
        let qual_vec = _mm256_loadu_si256(qual[i..].as_ptr() as *const __m256i);

        // Q20 counting: compare qual >= 53
        let q20_mask =
            _mm256_cmpgt_epi8(qual_vec, _mm256_sub_epi8(q20_thresh, _mm256_set1_epi8(1)));
        q20 += count_set_bits_avx2(q20_mask);

        // Q30 counting: compare qual >= 63
        let q30_mask =
            _mm256_cmpgt_epi8(qual_vec, _mm256_sub_epi8(q30_thresh, _mm256_set1_epi8(1)));
        q30 += count_set_bits_avx2(q30_mask);

        // Q40 counting: compare qual >= 73
        let q40_mask =
            _mm256_cmpgt_epi8(qual_vec, _mm256_sub_epi8(q40_thresh, _mm256_set1_epi8(1)));
        q40 += count_set_bits_avx2(q40_mask);

        // Unqualified counting: compare qual < qual_threshold_ascii
        if qual_threshold > 0 {
            let unqual_mask = _mm256_cmpgt_epi8(qual_thresh_vec, qual_vec); // inverted: thresh > qual means qual < thresh
            unqualified += count_set_bits_avx2(unqual_mask);
        }

        // Quality sum: subtract 33 and accumulate
        let adjusted = _mm256_sub_epi8(qual_vec, offset);
        qsum += horizontal_sum_u8_to_u32(adjusted);

        // Load sequence bases
        let seq_vec = _mm256_loadu_si256(seq[i..].as_ptr() as *const __m256i);

        // N counting
        let n_mask1 = _mm256_cmpeq_epi8(seq_vec, n_upper);
        let n_mask2 = _mm256_cmpeq_epi8(seq_vec, n_lower);
        let n_mask = _mm256_or_si256(n_mask1, n_mask2);
        ncnt += count_set_bits_avx2(n_mask);

        // GC counting
        let g_mask1 = _mm256_cmpeq_epi8(seq_vec, g_upper);
        let g_mask2 = _mm256_cmpeq_epi8(seq_vec, g_lower);
        let c_mask1 = _mm256_cmpeq_epi8(seq_vec, c_upper);
        let c_mask2 = _mm256_cmpeq_epi8(seq_vec, c_lower);
        let gc_mask = _mm256_or_si256(
            _mm256_or_si256(g_mask1, g_mask2),
            _mm256_or_si256(c_mask1, c_mask2),
        );
        gc += count_set_bits_avx2(gc_mask);

        i += chunk_size;
    }

    // Scalar remainder
    for j in i..len {
        let b = seq[j];
        let q = qual[j];
        let quality_u32 = (q - 33) as u32;
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
        if qual_threshold > 0 && q < qual_threshold_ascii {
            unqualified += 1;
        }

        if b == b'N' || b == b'n' {
            ncnt += 1;
        }
        if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
            gc += 1;
        }
    }

    Stats {
        qsum,
        q20,
        q30,
        q40,
        ncnt,
        gc,
        unqualified,
    }
}

#[cfg(target_arch = "x86_64")]
#[inline]
#[allow(unsafe_op_in_unsafe_fn)]
unsafe fn count_set_bits_avx2(mask: __m256i) -> usize {
    // Convert comparison mask to bitmask
    let bitmask = _mm256_movemask_epi8(mask) as u32;
    bitmask.count_ones() as usize
}

#[cfg(target_arch = "x86_64")]
#[inline]
#[allow(unsafe_op_in_unsafe_fn)]
unsafe fn horizontal_sum_u8_to_u32(vec: __m256i) -> u32 {
    // Unpack bytes to 16-bit integers (with zero extension)
    let zero = _mm256_setzero_si256();
    let low = _mm256_unpacklo_epi8(vec, zero);
    let high = _mm256_unpackhi_epi8(vec, zero);

    // Horizontal add to get 16-bit sums
    let sum16 = _mm256_add_epi16(low, high);

    // Unpack to 32-bit and sum
    let low32 = _mm256_unpacklo_epi16(sum16, zero);
    let high32 = _mm256_unpackhi_epi16(sum16, zero);
    let sum32 = _mm256_add_epi32(low32, high32);

    // Extract and sum all 32-bit values
    let mut result = [0u32; 8];
    _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, sum32);
    result.iter().sum()
}

// ARM NEON Implementation (aarch64)

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn count_mismatches_neon(
    seq1: &[u8],
    seq2: &[u8],
    max_diff: usize,
    min_complete: usize,
) -> (usize, usize) {
    let len = seq1.len();
    let mut differences = 0usize;
    let mut i = 0usize;
    let chunk_size = 16; // NEON processes 16 bytes (128-bit)

    unsafe {
        // Process 16 bytes at a time with NEON
        while i + chunk_size <= len {
            // Load both sequences
            let vec1 = vld1q_u8(seq1[i..].as_ptr());
            let vec2 = vld1q_u8(seq2[i..].as_ptr());

            // Compare: vceqq_u8 returns 0xFF for equal, 0x00 for not equal
            let eq_mask = vceqq_u8(vec1, vec2);

            // Count equals and subtract from 16 (eliminates vmvnq_u8 instruction)
            // Shift right by 7 to get 1 for 0xFF (equal), 0 for 0x00 (not equal)
            let shifted = vshrq_n_u8(eq_mask, 7);
            let equals = vaddlvq_u8(shifted) as usize;
            let chunk_mismatches = 16 - equals;
            let new_differences = differences + chunk_mismatches;
            let chunk_end = i + chunk_size;

            // Early termination check - only if threshold exceeded
            if new_differences > max_diff {
                if chunk_end <= min_complete {
                    // Chunk entirely before min_complete - chunk boundary is valid
                    return (new_differences, chunk_end);
                } else if i < min_complete {
                    // Chunk straddles min_complete boundary - need exact position
                    for j in 0..chunk_size {
                        let pos = i + j;
                        if seq1[pos] != seq2[pos] {
                            differences += 1;
                            if differences > max_diff && pos < min_complete {
                                return (differences, pos + 1);
                            }
                        }
                    }
                    // Threshold crossed at pos >= min_complete, continue
                }
                // else: chunk entirely past min_complete, no early termination possible
            }

            differences = new_differences;
            i += chunk_size;
        }
    }

    // Scalar remainder
    while i < len {
        if seq1[i] != seq2[i] {
            differences += 1;
            if differences > max_diff && i < min_complete {
                return (differences, i + 1);
            }
        }
        i += 1;
    }

    (differences, len)
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn compute_stats_neon(seq: &[u8], qual: &[u8], qual_threshold: u8) -> Stats {
    let len = seq.len().min(qual.len());
    let mut qsum = 0u32;
    let mut q20 = 0usize;
    let mut q30 = 0usize;
    let mut q40 = 0usize;
    let mut ncnt = 0usize;
    let mut gc = 0usize;
    let mut unqualified = 0usize;

    let mut i = 0;
    let chunk_size = 16; // NEON processes 16 bytes (128-bit)
    let qual_threshold_ascii = qual_threshold + 33;

    unsafe {
        let offset = vdupq_n_u8(33);
        let q20_thresh = vdupq_n_u8(53);
        let q30_thresh = vdupq_n_u8(63);
        let q40_thresh = vdupq_n_u8(73);
        let qual_thresh_vec = vdupq_n_u8(qual_threshold_ascii);

        let n_upper = vdupq_n_u8(b'N');
        let n_lower = vdupq_n_u8(b'n');
        let g_upper = vdupq_n_u8(b'G');
        let g_lower = vdupq_n_u8(b'g');
        let c_upper = vdupq_n_u8(b'C');
        let c_lower = vdupq_n_u8(b'c');

        // Process 16 bytes at a time
        while i + chunk_size <= len {
            // Load quality scores
            let qual_vec = vld1q_u8(qual[i..].as_ptr());

            // Q20 counting: compare qual >= 53
            let q20_mask = vcgeq_u8(qual_vec, q20_thresh);
            q20 += count_set_bits_neon(q20_mask);

            // Q30 counting: compare qual >= 63
            let q30_mask = vcgeq_u8(qual_vec, q30_thresh);
            q30 += count_set_bits_neon(q30_mask);

            // Q40 counting: compare qual >= 73
            let q40_mask = vcgeq_u8(qual_vec, q40_thresh);
            q40 += count_set_bits_neon(q40_mask);

            // Unqualified counting: compare qual < qual_threshold_ascii
            if qual_threshold > 0 {
                let unqual_mask = vcltq_u8(qual_vec, qual_thresh_vec); // less than
                unqualified += count_set_bits_neon(unqual_mask);
            }

            // Quality sum: subtract 33 and accumulate
            let adjusted = vsubq_u8(qual_vec, offset);
            qsum += horizontal_sum_u8_to_u32_neon(adjusted);

            // Load sequence bases
            let seq_vec = vld1q_u8(seq[i..].as_ptr());

            // N counting
            let n_mask1 = vceqq_u8(seq_vec, n_upper);
            let n_mask2 = vceqq_u8(seq_vec, n_lower);
            let n_mask = vorrq_u8(n_mask1, n_mask2);
            ncnt += count_set_bits_neon(n_mask);

            // GC counting
            let g_mask1 = vceqq_u8(seq_vec, g_upper);
            let g_mask2 = vceqq_u8(seq_vec, g_lower);
            let c_mask1 = vceqq_u8(seq_vec, c_upper);
            let c_mask2 = vceqq_u8(seq_vec, c_lower);
            let gc_mask = vorrq_u8(vorrq_u8(g_mask1, g_mask2), vorrq_u8(c_mask1, c_mask2));
            gc += count_set_bits_neon(gc_mask);

            i += chunk_size;
        }
    } // end unsafe block

    // Scalar remainder
    for j in i..len {
        let b = seq[j];
        let q = qual[j];
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
        if qual_threshold > 0 && q < qual_threshold_ascii {
            unqualified += 1;
        }

        if b == b'N' || b == b'n' {
            ncnt += 1;
        }
        if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
            gc += 1;
        }
    }

    Stats {
        qsum,
        q20,
        q30,
        q40,
        ncnt,
        gc,
        unqualified,
    }
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn count_set_bits_neon(mask: uint8x16_t) -> usize {
    unsafe {
        // Count number of 0xFF bytes in mask (comparison result)
        // Each comparison sets all bits to 1 (0xFF) if true, 0 otherwise

        // Shift right by 7 to get just the sign bit (0 or 1)
        let shifted = vshrq_n_u8(mask, 7);

        // Horizontal sum to count set bits
        let sum = vaddlvq_u8(shifted);
        sum as usize
    }
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn horizontal_sum_u8_to_u32_neon(vec: uint8x16_t) -> u32 {
    unsafe {
        // Widen and sum all bytes - vaddlvq_u8 returns u16
        u32::from(vaddlvq_u8(vec))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_sum() {
        let qual = b"IIIIIIIIII"; // All quality 40 (ASCII 73, Phred 40)
        let seq = b"AAAAAAAAAA";
        let stats = compute_stats(seq, qual, 0);
        assert_eq!(stats.qsum, 400); // 10 * 40
    }

    #[test]
    fn test_q20_q30_counting() {
        // Q20 = Phred 20 = ASCII 53 = '5'
        // Q30 = Phred 30 = ASCII 63 = '?'
        let qual = b"!!!!555555????????"; // 4 low, 6 Q20, 8 Q30
        let seq = b"AAAAAAAAAAAAAAAAAA";
        let stats = compute_stats(seq, qual, 0);
        assert!(stats.q20 >= 14); // 6 Q20 + 8 Q30
        assert_eq!(stats.q30, 8);
    }

    #[test]
    fn test_n_counting() {
        let seq = b"AAANNNAAANNAAA";
        let qual = b"IIIIIIIIIIIIII";
        let stats = compute_stats(seq, qual, 0);
        assert_eq!(stats.ncnt, 5);
    }

    #[test]
    fn test_gc_counting() {
        let seq = b"ATCGATCGATCGAT"; // 3 G, 3 C = 6 GC
        let qual = b"IIIIIIIIIIIIII";
        let stats = compute_stats(seq, qual, 0);
        assert_eq!(stats.gc, 6);
    }

    #[test]
    fn test_mixed_case() {
        let seq = b"AtCgNnGgCc";
        let qual = b"IIIIIIIIII";
        let stats = compute_stats(seq, qual, 0);
        assert_eq!(stats.ncnt, 2); // N and n
        assert_eq!(stats.gc, 6); // C, g, G, g, C, c
    }

    #[test]
    fn test_unqualified_counting() {
        // Quality threshold 15 = ASCII 48 = '0'
        // Qualities below Q15: ! through / (ASCII 33-47)
        let qual = b"!!!!000000IIIIIIII"; // 4 below Q15, 6 at Q15, 8 above Q15
        let seq = b"AAAAAAAAAAAAAAAAAA";
        let stats = compute_stats(seq, qual, 15);
        assert_eq!(stats.unqualified, 4);
    }

    #[test]
    fn test_neon_vs_scalar_simple() {
        // Simple test first with the exact quality chars from problematic read
        let r1_qual = b"AAAA/EA/E6EA/AE/E6E/EEA/EAE/EAE/EE//EE/EE";
        let r1_seq = b"GACTGGGGAGACGCCGAGTGAGGGCGAGCGCTGCTGTGGCG";

        let scalar = compute_stats_scalar(r1_seq, r1_qual, 15);

        #[cfg(target_arch = "aarch64")]
        {
            let neon = unsafe { compute_stats_neon(r1_seq, r1_qual, 15) };
            println!("\nSimple NEON vs Scalar test:");
            println!(
                "  Scalar: unqualified={}, qsum={}",
                scalar.unqualified, scalar.qsum
            );
            println!(
                "  NEON:   unqualified={}, qsum={}",
                neon.unqualified, neon.qsum
            );
            assert_eq!(
                neon.unqualified, scalar.unqualified,
                "NEON/scalar unqualified mismatch!"
            );
            assert_eq!(neon.qsum, scalar.qsum, "NEON/scalar qsum mismatch!");
        }

        #[cfg(not(target_arch = "aarch64"))]
        {
            println!("Skipping NEON test (not on ARM)");
        }
    }

    #[test]
    fn test_compare_scalar_vs_simd() {
        let seq = b"ATCGATCGATCGATNNATCGATCGATCGATCGATCG";
        let qual = b"IIII555555????????IIII555555????????II";

        let scalar = compute_stats_scalar(seq, qual, 0);

        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                let simd = unsafe { compute_stats_avx2(seq, qual, 0) };
                assert_eq!(scalar.qsum, simd.qsum, "qsum mismatch");
                assert_eq!(scalar.q20, simd.q20, "q20 mismatch");
                assert_eq!(scalar.q30, simd.q30, "q30 mismatch");
                assert_eq!(scalar.ncnt, simd.ncnt, "ncnt mismatch");
                assert_eq!(scalar.gc, simd.gc, "gc mismatch");
                assert_eq!(scalar.unqualified, simd.unqualified, "unqualified mismatch");
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            let simd = unsafe { compute_stats_neon(seq, qual, 0) };
            assert_eq!(scalar.qsum, simd.qsum, "qsum mismatch");
            assert_eq!(scalar.q20, simd.q20, "q20 mismatch");
            assert_eq!(scalar.q30, simd.q30, "q30 mismatch");
            assert_eq!(scalar.ncnt, simd.ncnt, "ncnt mismatch");
            assert_eq!(scalar.gc, simd.gc, "gc mismatch");
            assert_eq!(scalar.unqualified, simd.unqualified, "unqualified mismatch");
        }
    }

    /// Test the problematic read from SRR22472290.464
    /// This read is filtered by fasterp but passed by fastp
    // Early Termination Semantics Tests for count_mismatches
    // These tests ensure SIMD implementations maintain byte-level early termination
    // semantics. The key invariant: returned position `i` must be the exact byte
    // where termination occurred, not a chunk boundary.

    #[test]
    fn test_early_termination_basic() {
        // Basic case: 6 differences in first 10 bytes with max_diff=5, min_complete=50
        // Should terminate early and return position < 50
        let seq1 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 64 A's
        let mut seq2 = seq1.to_vec();
        // Put 6 differences in positions 0, 2, 4, 6, 8, 10
        for i in (0..12).step_by(2) {
            seq2[i] = b'T';
        }

        let (diff, pos) = count_mismatches(seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        // Key semantic: pos must be < min_complete for early termination to apply
        // SIMD may return chunk boundary (16) instead of exact position (11), both are < 50
        assert!(
            pos <= 50,
            "Should terminate early with position <= 50, got {}",
            pos
        );
        assert!(
            pos > 10,
            "Position must be past the 6th mismatch at position 10"
        );
    }

    #[test]
    fn test_early_termination_at_chunk_boundary() {
        // Critical test: 6th difference at position 15 (end of first SIMD chunk)
        // SIMD bug would return pos=16, scalar returns pos=16 - same in this case
        // But the intermediate check matters!
        let seq1 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 64 A's
        let mut seq2 = seq1.to_vec();
        // 6 differences at positions 0, 3, 6, 9, 12, 15
        for i in (0..18).step_by(3) {
            seq2[i] = b'T';
        }

        let (diff, pos) = count_mismatches(seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        assert_eq!(
            pos, 16,
            "Should terminate at position 16 (right after 6th diff at 15)"
        );
    }

    #[test]
    fn test_early_termination_across_chunks() {
        // 5 differences in first chunk (0-15), 1 more in second chunk
        // Tests that SIMD correctly handles threshold crossing in second chunk
        let seq1 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 64 A's
        let mut seq2 = seq1.to_vec();
        // 5 differences at positions 0, 3, 6, 9, 12 (all in first chunk)
        for i in (0..15).step_by(3) {
            seq2[i] = b'T';
        }
        // 6th difference at position 20 (second chunk)
        seq2[20] = b'T';

        let (diff, pos) = count_mismatches(seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        // Key semantic: pos must be <= min_complete for early termination
        // SIMD returns chunk boundary (32), scalar returns exact (21), both are < 50
        assert!(
            pos <= 50,
            "Should terminate with position <= 50, got {}",
            pos
        );
        assert!(
            pos > 20,
            "Position must be past the 6th mismatch at position 20"
        );
    }

    #[test]
    fn test_early_termination_near_min_complete() {
        // The most critical case: 6th difference at position 48 (just before min_complete=50)
        // Scalar: terminates at 49, i=49 > 50 is FALSE → rejects
        // Buggy SIMD: continues past chunk boundary, i=64 > 50 is TRUE → accepts (WRONG!)
        let seq1 = vec![b'A'; 64];
        let mut seq2 = seq1.clone();
        // Place 5 differences spread out, then 6th at position 48
        seq2[10] = b'T'; // diff 1
        seq2[20] = b'T'; // diff 2
        seq2[30] = b'T'; // diff 3
        seq2[35] = b'T'; // diff 4
        seq2[40] = b'T'; // diff 5
        seq2[48] = b'T'; // diff 6 - critical position!

        let (diff, pos) = count_mismatches(&seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        // CRITICAL: must return 49, not 64!
        // If SIMD returns 64, then 64 > 50 = true → lenient mode incorrectly applies
        // If scalar returns 49, then 49 > 50 = false → correct early rejection
        assert_eq!(
            pos, 49,
            "Must terminate at position 49, not chunk boundary 64"
        );
        assert!(
            pos <= 50,
            "Position {} must be <= 50 for early termination to work",
            pos
        );
    }

    #[test]
    fn test_no_early_termination_when_completing() {
        // When we complete all comparisons (reach the end), should return full length
        // even if differences exceed max_diff (lenient mode)
        let seq1 = vec![b'A'; 60];
        let mut seq2 = seq1.clone();
        // 6 differences spread throughout
        seq2[5] = b'T';
        seq2[15] = b'T';
        seq2[25] = b'T';
        seq2[35] = b'T';
        seq2[45] = b'T';
        seq2[55] = b'T';

        let (diff, pos) = count_mismatches(&seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        // Should complete all 60 comparisons because after position 50,
        // the early termination condition (i < min_complete) is false
        assert_eq!(
            pos, 60,
            "Should complete all comparisons when past min_complete"
        );
    }

    #[test]
    fn test_early_termination_exactly_at_min_complete() {
        // Edge case: 6th difference exactly at position 49 (last position before min_complete=50)
        let seq1 = vec![b'A'; 64];
        let mut seq2 = seq1.clone();
        seq2[8] = b'T'; // diff 1
        seq2[16] = b'T'; // diff 2
        seq2[24] = b'T'; // diff 3
        seq2[32] = b'T'; // diff 4
        seq2[40] = b'T'; // diff 5
        seq2[49] = b'T'; // diff 6 - exactly at position 49

        let (diff, pos) = count_mismatches(&seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        // Position 49 < 50, so early termination should trigger
        assert_eq!(pos, 50, "Should terminate at position 50");
        assert!(pos <= 50, "Early termination check: {} <= 50", pos);
    }

    #[test]
    fn test_early_termination_just_past_min_complete() {
        // 6th difference at position 50 (first position >= min_complete)
        // Should NOT early terminate because i >= min_complete
        let seq1 = vec![b'A'; 64];
        let mut seq2 = seq1.clone();
        seq2[8] = b'T'; // diff 1
        seq2[16] = b'T'; // diff 2
        seq2[24] = b'T'; // diff 3
        seq2[32] = b'T'; // diff 4
        seq2[40] = b'T'; // diff 5
        seq2[50] = b'T'; // diff 6 - at position 50, past min_complete

        let (diff, pos) = count_mismatches(&seq1, &seq2, 5, 50);

        assert_eq!(diff, 6, "Should count 6 differences");
        // Should complete all comparisons because at position 50, i >= min_complete
        assert_eq!(
            pos, 64,
            "Should complete all comparisons when diff occurs at/past min_complete"
        );
    }

    #[test]
    fn test_simd_vs_scalar_consistency() {
        // Comprehensive test: ensure SIMD and scalar produce semantically equivalent results
        // SIMD may return chunk boundaries, scalar returns exact positions
        // Both must agree on: diff count, and whether pos triggers early termination
        let test_cases = vec![
            // (seq_len, diff_positions, max_diff, min_complete)
            (64, vec![5, 10, 15, 20, 25, 30], 5, 50), // spread across chunks
            (64, vec![0, 1, 2, 3, 4, 5], 5, 50),      // all early
            (64, vec![45, 46, 47, 48, 49, 50], 5, 50), // near boundary
            (64, vec![14, 15, 16, 17, 18, 19], 5, 50), // across first chunk boundary
            (64, vec![30, 31, 32, 33, 34, 35], 5, 50), // across second chunk boundary
            (150, vec![40, 45, 48, 49, 50, 51], 5, 50), // realistic read length
        ];

        for (seq_len, diff_positions, max_diff, min_complete) in test_cases {
            let seq1 = vec![b'A'; seq_len];
            let mut seq2 = seq1.clone();
            for &pos in &diff_positions {
                if pos < seq_len {
                    seq2[pos] = b'T';
                }
            }

            // Get result from current implementation (SIMD)
            let (simd_diff, simd_pos) = count_mismatches(&seq1, &seq2, max_diff, min_complete);

            // Get expected result from scalar
            let (scalar_diff, scalar_pos) =
                count_mismatches_scalar(&seq1, &seq2, max_diff, min_complete);

            // Diff counts must match
            assert_eq!(
                simd_diff, scalar_diff,
                "Diff count mismatch for seq_len={}, diffs={:?}: SIMD={}, scalar={}",
                seq_len, diff_positions, simd_diff, scalar_diff
            );

            // Semantic equivalence: both must agree on early termination decision
            // Early termination applies when: diff > max_diff && pos <= min_complete
            let simd_early_term = simd_diff > max_diff && simd_pos <= min_complete;
            let scalar_early_term = scalar_diff > max_diff && scalar_pos <= min_complete;

            assert_eq!(
                simd_early_term, scalar_early_term,
                "Early termination decision mismatch for seq_len={}, diffs={:?}: SIMD pos={} (early={}), scalar pos={} (early={})",
                seq_len, diff_positions, simd_pos, simd_early_term, scalar_pos, scalar_early_term
            );
        }
    }

    #[test]
    fn test_problematic_read_srr22472290_464() {
        // The exact quality strings from the problematic read pair
        let r1_seq = b"GACTGGGGAGACGCCGAGTGAGGGCGAGCGCTGCTGTGGCG";
        let r1_qual = b"AAAA/EA/E6EA/AE/E6E/EEA/EAE/EAE/EE//EE/EE";

        let r2_seq = b"CGCCACAGCCGCGCTCGCCCTCACTCGGCGTCTACCCAGTCCTGTCTATTATACACAAAAGACGATGCCGACGAA";
        let r2_qual =
            b"A/A//EEE/EA<EE/EEEEAA<EE//AEEEAAE/EE/EAEEAAA//A////AE/A<A//////////E///////";

        let qual_threshold = 15; // Default qualified_quality_phred

        // Compute stats using the scalar fallback
        let r1_stats_scalar = compute_stats_scalar(r1_seq, r1_qual, qual_threshold);
        let r2_stats_scalar = compute_stats_scalar(r2_seq, r2_qual, qual_threshold);

        println!("\n=== Problematic Read SRR22472290.464 Analysis ===");
        println!("\nR1 (42 bases) Scalar Results:");
        println!(
            "  unqualified={}, qsum={}, q20={}, q30={}, ncnt={}, gc={}",
            r1_stats_scalar.unqualified,
            r1_stats_scalar.qsum,
            r1_stats_scalar.q20,
            r1_stats_scalar.q30,
            r1_stats_scalar.ncnt,
            r1_stats_scalar.gc
        );

        println!("\nR2 (76 bases) Scalar Results:");
        println!(
            "  unqualified={}, qsum={}, q20={}, q30={}, ncnt={}, gc={}",
            r2_stats_scalar.unqualified,
            r2_stats_scalar.qsum,
            r2_stats_scalar.q20,
            r2_stats_scalar.q30,
            r2_stats_scalar.ncnt,
            r2_stats_scalar.gc
        );

        // Test SIMD implementations if available
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                let r1_stats_avx2 = unsafe { compute_stats_avx2(r1_seq, r1_qual, qual_threshold) };
                let r2_stats_avx2 = unsafe { compute_stats_avx2(r2_seq, r2_qual, qual_threshold) };

                println!("\nR1 (42 bases) AVX2 Results:");
                println!(
                    "  unqualified={}, qsum={}, q20={}, q30={}, ncnt={}, gc={}",
                    r1_stats_avx2.unqualified,
                    r1_stats_avx2.qsum,
                    r1_stats_avx2.q20,
                    r1_stats_avx2.q30,
                    r1_stats_avx2.ncnt,
                    r1_stats_avx2.gc
                );

                println!("\nR2 (76 bases) AVX2 Results:");
                println!(
                    "  unqualified={}, qsum={}, q20={}, q30={}, ncnt={}, gc={}",
                    r2_stats_avx2.unqualified,
                    r2_stats_avx2.qsum,
                    r2_stats_avx2.q20,
                    r2_stats_avx2.q30,
                    r2_stats_avx2.ncnt,
                    r2_stats_avx2.gc
                );

                // Check filtering decisions
                let unqualified_percent_limit = 40;
                let r1_threshold = (unqualified_percent_limit * r1_seq.len()) as f64 / 100.0;
                let r2_threshold = (unqualified_percent_limit * r2_seq.len()) as f64 / 100.0;

                let r1_fail_scalar = (r1_stats_scalar.unqualified as f64) > r1_threshold;
                let r1_fail_avx2 = (r1_stats_avx2.unqualified as f64) > r1_threshold;
                let r2_fail_scalar = (r2_stats_scalar.unqualified as f64) > r2_threshold;
                let r2_fail_avx2 = (r2_stats_avx2.unqualified as f64) > r2_threshold;

                println!("\nFiltering Decision (40% unqualified threshold):");
                println!(
                    "  R1: threshold={:.2}, Scalar fail={}, AVX2 fail={}",
                    r1_threshold, r1_fail_scalar, r1_fail_avx2
                );
                println!(
                    "  R2: threshold={:.2}, Scalar fail={}, AVX2 fail={}",
                    r2_threshold, r2_fail_scalar, r2_fail_avx2
                );

                // THE KEY ASSERTIONS: SIMD and scalar must match
                assert_eq!(
                    r1_stats_avx2.unqualified, r1_stats_scalar.unqualified,
                    "R1: AVX2 and scalar produce different unqualified counts!"
                );
                assert_eq!(
                    r1_stats_avx2.qsum, r1_stats_scalar.qsum,
                    "R1: AVX2 and scalar produce different quality sums!"
                );
                assert_eq!(
                    r1_stats_avx2.q20, r1_stats_scalar.q20,
                    "R1: AVX2 and scalar produce different Q20 counts!"
                );
                assert_eq!(
                    r1_stats_avx2.q30, r1_stats_scalar.q30,
                    "R1: AVX2 and scalar produce different Q30 counts!"
                );

                assert_eq!(
                    r2_stats_avx2.unqualified, r2_stats_scalar.unqualified,
                    "R2: AVX2 and scalar produce different unqualified counts!"
                );
                assert_eq!(
                    r2_stats_avx2.qsum, r2_stats_scalar.qsum,
                    "R2: AVX2 and scalar produce different quality sums!"
                );
                assert_eq!(
                    r2_stats_avx2.q20, r2_stats_scalar.q20,
                    "R2: AVX2 and scalar produce different Q20 counts!"
                );
                assert_eq!(
                    r2_stats_avx2.q30, r2_stats_scalar.q30,
                    "R2: AVX2 and scalar produce different Q30 counts!"
                );

                println!("\n✓ AVX2 and scalar implementations produce identical results");
            } else {
                println!("\nAVX2 not available, skipping SIMD comparison");
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            let r1_stats_neon = unsafe { compute_stats_neon(r1_seq, r1_qual, qual_threshold) };
            let r2_stats_neon = unsafe { compute_stats_neon(r2_seq, r2_qual, qual_threshold) };

            println!("\nR1 (42 bases) NEON Results:");
            println!(
                "  unqualified={}, qsum={}, q20={}, q30={}, ncnt={}, gc={}",
                r1_stats_neon.unqualified,
                r1_stats_neon.qsum,
                r1_stats_neon.q20,
                r1_stats_neon.q30,
                r1_stats_neon.ncnt,
                r1_stats_neon.gc
            );

            println!("\nR2 (76 bases) NEON Results:");
            println!(
                "  unqualified={}, qsum={}, q20={}, q30={}, ncnt={}, gc={}",
                r2_stats_neon.unqualified,
                r2_stats_neon.qsum,
                r2_stats_neon.q20,
                r2_stats_neon.q30,
                r2_stats_neon.ncnt,
                r2_stats_neon.gc
            );

            // THE KEY ASSERTIONS: SIMD and scalar must match
            assert_eq!(
                r1_stats_neon.unqualified, r1_stats_scalar.unqualified,
                "R1: NEON and scalar produce different unqualified counts!"
            );
            assert_eq!(
                r2_stats_neon.unqualified, r2_stats_scalar.unqualified,
                "R2: NEON and scalar produce different unqualified counts!"
            );

            println!("\n✓ NEON and scalar implementations produce identical results");
        }
    }
}
