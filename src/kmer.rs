//! Kmer counting and base encoding
//!
//! This module provides efficient kmer counting using 2-bit encoding
//! to avoid string allocations. Also includes lookup tables for fast
//! base and quality checks.

// LOOKUP TABLES for fast base/quality checks

/// Lookup table: is quality >= 20? (Phred+33 encoding)
#[allow(dead_code)]
pub(crate) static LUT_Q20: [bool; 256] = {
    let mut lut = [false; 256];
    let mut i = 0;
    while i < 256 {
        lut[i] = i >= 33 && (i - 33) >= 20;
        i += 1;
    }
    lut
};

/// Lookup table: is quality >= 30? (Phred+33 encoding)
#[allow(dead_code)]
pub(crate) static LUT_Q30: [bool; 256] = {
    let mut lut = [false; 256];
    let mut i = 0;
    while i < 256 {
        lut[i] = i >= 33 && (i - 33) >= 30;
        i += 1;
    }
    lut
};

/// Lookup table: is base N?
#[allow(dead_code)]
pub(crate) static LUT_IS_N: [bool; 256] = {
    let mut lut = [false; 256];
    lut[b'N' as usize] = true;
    lut[b'n' as usize] = true;
    lut
};

/// Lookup table: is base G or C?
#[allow(dead_code)]
pub(crate) static LUT_IS_GC: [bool; 256] = {
    let mut lut = [false; 256];
    lut[b'G' as usize] = true;
    lut[b'g' as usize] = true;
    lut[b'C' as usize] = true;
    lut[b'c' as usize] = true;
    lut
};

/// Lookup table: base to 2-bit encoding (A=0, C=1, G=2, T=3)
/// Invalid bases map to 255 (sentinel value)
const LUT_BASE_TO_2BIT: [u8; 256] = {
    let mut lut = [255u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut
};

/// Lookup table: base index for `quality_curves` (A=0, T=1, C=2, G=3)
/// Invalid bases map to 255 (sentinel value)
const LUT_BASE_IDX: [u8; 256] = {
    let mut lut = [255u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'T' as usize] = 1;
    lut[b't' as usize] = 1;
    lut[b'C' as usize] = 2;
    lut[b'c' as usize] = 2;
    lut[b'G' as usize] = 3;
    lut[b'g' as usize] = 3;
    lut
};

/// Convert base to 2-bit encoding: A=0, C=1, G=2, T=3
/// Uses lookup table for O(1) performance
#[allow(clippy::inline_always)]
#[inline(always)]
pub(crate) fn base_to_2bit(b: u8) -> Option<u32> {
    let code = LUT_BASE_TO_2BIT[b as usize];
    if code == 255 {
        None
    } else {
        Some(u32::from(code))
    }
}

/// Get base index for `quality_curves`: A=0, T=1, C=2, G=3
/// Uses lookup table for O(1) performance
#[allow(clippy::inline_always)]
#[inline(always)]
pub(crate) fn base_idx(b: u8) -> Option<usize> {
    let idx = LUT_BASE_IDX[b as usize];
    if idx == 255 { None } else { Some(idx as usize) }
}

/// Count 5-mers using 2-bit rolling code
///
/// This replaces the old String-based approach that allocated millions of strings.
/// Uses a fixed array of 1024 elements (4^5 possible 5-mers).
/// Encodes A=0, C=1, G=2, T=3 and rolls a 10-bit window.
/// Any N base resets the window.
#[inline]
pub(crate) fn count_k5_2bit(seq: &[u8], kmer_table: &mut [usize; 1024]) {
    if seq.len() < 5 {
        return;
    }

    let mask: u32 = (1 << (2 * 5)) - 1; // 10 bits for 5-mer
    let mut code: u32 = 0;
    let mut valid_count: usize = 0;

    // SAFETY: All accesses are bounds-checked by the loop structure
    // and LUT_BASE_TO_2BIT is 256 elements (all u8 values valid)
    unsafe {
        for i in 0..seq.len() {
            let b = *seq.get_unchecked(i);
            let c = *LUT_BASE_TO_2BIT.get_unchecked(b as usize);

            if c == 255 {
                // Hit an N or invalid base - reset window
                code = 0;
                valid_count = 0;
                continue;
            }

            code = ((code << 2) & mask) | u32::from(c);
            valid_count += 1;

            if valid_count >= 5 {
                *kmer_table.get_unchecked_mut(code as usize) += 1;
            }
        }
    }
}

/// Static cache of all 1024 possible 5-mer strings
/// Initialized once on first use, then reused
static KMER_CACHE: std::sync::OnceLock<Vec<String>> = std::sync::OnceLock::new();

/// Convert 2-bit encoded kmer to &str for JSON output
/// Uses static cache - initialized once, then zero allocations
#[allow(clippy::inline_always)]
#[inline(always)]
pub(crate) fn kmer_to_str(code: usize) -> &'static str {
    let cache = KMER_CACHE.get_or_init(|| {
        let mut kmers = Vec::with_capacity(1024);
        for c in 0..1024 {
            kmers.push(kmer_to_string_impl(c));
        }
        kmers
    });

    // SAFETY: The cache is never mutated after initialization
    // and lives for the entire program lifetime
    unsafe { std::mem::transmute::<&str, &'static str>(cache[code].as_str()) }
}

/// Internal implementation of `kmer_to_string`
fn kmer_to_string_impl(code: usize) -> String {
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
