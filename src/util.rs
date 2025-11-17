//! Unsafe fast helpers for hot paths
//!
//! All unsafe operations are encapsulated in small, auditable functions.
//! Callers must validate lengths/capacity before calling these functions.

/// Fast iteration over seq, qual, and position index without bounds checks
///
/// # Safety
/// - `seq` and `qual` must point to valid memory of at least `len` bytes
/// - `len` must be accurate
/// - Caller must ensure seq and qual are aligned and valid for the entire length
#[allow(clippy::inline_always)]
#[inline(always)]
pub(crate) unsafe fn loop_seq_qual_indexed<F: FnMut(usize, u8, u8)>(
    seq: *const u8,
    qual: *const u8,
    len: usize,
    mut f: F,
) {
    // SAFETY: Caller ensures seq and qual are valid for len bytes
    unsafe {
        let mut i = 0;

        // Process 4 bytes at a time (unrolled for better performance)
        while i + 4 <= len {
            f(i, *seq.add(i), *qual.add(i));
            f(i + 1, *seq.add(i + 1), *qual.add(i + 1));
            f(i + 2, *seq.add(i + 2), *qual.add(i + 2));
            f(i + 3, *seq.add(i + 3), *qual.add(i + 3));
            i += 4;
        }

        // Process remaining bytes
        while i < len {
            f(i, *seq.add(i), *qual.add(i));
            i += 1;
        }
    }
}
