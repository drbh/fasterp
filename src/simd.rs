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
use std::arch::aarch64::*;

/// Compute multiple stats in a single SIMD pass over sequence and quality data
pub struct Stats {
    pub qsum: u32,
    pub q20: usize,
    pub q30: usize,
    pub ncnt: usize,
    pub gc: usize,
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

/// Compute all stats for a sequence/quality pair using SIMD when available
#[inline]
pub fn compute_stats(seq: &[u8], qual: &[u8]) -> Stats {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { compute_stats_avx2(seq, qual) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        return unsafe { compute_stats_neon(seq, qual) };
    }

    // Scalar fallback
    compute_stats_scalar(seq, qual)
}

/// Scalar implementation (fallback)
#[inline]
fn compute_stats_scalar(seq: &[u8], qual: &[u8]) -> Stats {
    let mut qsum = 0u32;
    let mut q20 = 0usize;
    let mut q30 = 0usize;
    let mut ncnt = 0usize;
    let mut gc = 0usize;

    for (&b, &q) in seq.iter().zip(qual) {
        let qval = (q - 33) as u32;
        qsum += qval;

        // Q20/Q30: quality thresholds
        if q >= 53 {
            q20 += 1;
        } // Phred 20 = ASCII 53
        if q >= 63 {
            q30 += 1;
        } // Phred 30 = ASCII 63

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
        ncnt,
        gc,
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compute_stats_avx2(seq: &[u8], qual: &[u8]) -> Stats {
    let len = seq.len().min(qual.len());
    let mut qsum = 0u32;
    let mut q20 = 0usize;
    let mut q30 = 0usize;
    let mut ncnt = 0usize;
    let mut gc = 0usize;

    let offset = _mm256_set1_epi8(33);
    let q20_thresh = _mm256_set1_epi8(53);
    let q30_thresh = _mm256_set1_epi8(63);

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
    }

    Stats {
        qsum,
        q20,
        q30,
        ncnt,
        gc,
    }
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn count_set_bits_avx2(mask: __m256i) -> usize {
    // Convert comparison mask to bitmask
    let bitmask = _mm256_movemask_epi8(mask) as u32;
    bitmask.count_ones() as usize
}

#[cfg(target_arch = "x86_64")]
#[inline]
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

// ============================================================================
// ARM NEON Implementation (aarch64)
// ============================================================================

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn compute_stats_neon(seq: &[u8], qual: &[u8]) -> Stats {
    let len = seq.len().min(qual.len());
    let mut qsum = 0u32;
    let mut q20 = 0usize;
    let mut q30 = 0usize;
    let mut ncnt = 0usize;
    let mut gc = 0usize;

    let mut i = 0;
    let chunk_size = 16; // NEON processes 16 bytes (128-bit)

    unsafe {
        let offset = vdupq_n_u8(33);
        let q20_thresh = vdupq_n_u8(53);
        let q30_thresh = vdupq_n_u8(63);

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
    }

    Stats {
        qsum,
        q20,
        q30,
        ncnt,
        gc,
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
        vaddlvq_u8(vec) as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_sum() {
        let qual = b"IIIIIIIIII"; // All quality 40 (ASCII 73, Phred 40)
        let seq = b"AAAAAAAAAA";
        let stats = compute_stats(seq, qual);
        assert_eq!(stats.qsum, 400); // 10 * 40
    }

    #[test]
    fn test_q20_q30_counting() {
        // Q20 = Phred 20 = ASCII 53 = '5'
        // Q30 = Phred 30 = ASCII 63 = '?'
        let qual = b"!!!!555555????????"; // 4 low, 6 Q20, 8 Q30
        let seq = b"AAAAAAAAAAAAAAAAAA";
        let stats = compute_stats(seq, qual);
        assert!(stats.q20 >= 14); // 6 Q20 + 8 Q30
        assert_eq!(stats.q30, 8);
    }

    #[test]
    fn test_n_counting() {
        let seq = b"AAANNNAAANNAAA";
        let qual = b"IIIIIIIIIIIIII";
        let stats = compute_stats(seq, qual);
        assert_eq!(stats.ncnt, 5);
    }

    #[test]
    fn test_gc_counting() {
        let seq = b"ATCGATCGATCGAT"; // 3 G, 3 C = 6 GC
        let qual = b"IIIIIIIIIIIIII";
        let stats = compute_stats(seq, qual);
        assert_eq!(stats.gc, 6);
    }

    #[test]
    fn test_mixed_case() {
        let seq = b"AtCgNnGgCc";
        let qual = b"IIIIIIIIII";
        let stats = compute_stats(seq, qual);
        assert_eq!(stats.ncnt, 2); // N and n
        assert_eq!(stats.gc, 6); // C, g, G, g, C, c
    }

    #[test]
    fn test_compare_scalar_vs_simd() {
        let seq = b"ATCGATCGATCGATNNATCGATCGATCGATCGATCG";
        let qual = b"IIII555555????????IIII555555????????II";

        let scalar = compute_stats_scalar(seq, qual);

        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                let simd = unsafe { compute_stats_avx2(seq, qual) };
                assert_eq!(scalar.qsum, simd.qsum, "qsum mismatch");
                assert_eq!(scalar.q20, simd.q20, "q20 mismatch");
                assert_eq!(scalar.q30, simd.q30, "q30 mismatch");
                assert_eq!(scalar.ncnt, simd.ncnt, "ncnt mismatch");
                assert_eq!(scalar.gc, simd.gc, "gc mismatch");
            }
        }
    }
}
