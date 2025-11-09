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

/// Convert base to 2-bit encoding: A=0, C=1, G=2, T=3
#[inline]
pub(crate) fn base_to_2bit(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Get base index for quality_curves: A=0, T=1, C=2, G=3
#[inline]
pub(crate) fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'T' | b't' => Some(1),
        b'C' | b'c' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

/// Count 5-mers using 2-bit rolling code
///
/// This replaces the old String-based approach that allocated millions of strings.
/// Uses a fixed array of 1024 elements (4^5 possible 5-mers).
/// Encodes A=0, C=1, G=2, T=3 and rolls a 10-bit window.
/// Any N base resets the window.
#[inline]
pub(crate) fn count_k5_2bit(seq: &[u8], kmer_table: &mut [usize; 1024]) {
    let mut code: u32 = 0;
    let mask: u32 = (1 << (2 * 5)) - 1; // 10 bits for 5-mer
    let mut filled = 0u8;

    for &b in seq {
        let Some(c) = base_to_2bit(b) else {
            // Hit an N or invalid base - reset window
            code = 0;
            filled = 0;
            continue;
        };

        code = ((code << 2) & mask) | c;

        if filled < 4 {
            filled += 1;
            continue;
        }

        kmer_table[code as usize] += 1;
    }
}

/// Convert 2-bit encoded kmer to String for JSON output
pub(crate) fn kmer_to_string(code: usize) -> String {
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
