//! Locate all line endings in a block with one SIMD comparison per byte.
#[cfg(any(
    all(target_arch = "aarch64", target_endian = "little"),
    target_arch = "x86_64"
))]
pub(super) struct Offsets<'a> {
    bytes: &'a [u8],
    position: usize,
    matches: u64,
}

#[cfg(any(
    all(target_arch = "aarch64", target_endian = "little"),
    target_arch = "x86_64"
))]
pub(super) fn offsets(bytes: &[u8]) -> Offsets<'_> {
    Offsets {
        bytes,
        position: 0,
        matches: 0,
    }
}

#[cfg(any(
    all(target_arch = "aarch64", target_endian = "little"),
    target_arch = "x86_64"
))]
impl Iterator for Offsets<'_> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        while self.matches == 0 {
            if let Some((chunk, rest)) = self.bytes.split_first_chunk::<64>() {
                self.matches = newline_mask(chunk);
                self.bytes = rest;
                self.position += chunk.len();
            } else {
                let index = memchr::memchr(b'\n', self.bytes)?;
                self.bytes = &self.bytes[index + 1..];
                self.position += index + 1;
                return Some(self.position - 1);
            }
        }
        let index = self.matches.trailing_zeros() as usize;
        self.matches &= self.matches - 1;
        Some(self.position - 64 + index)
    }
}

#[cfg(not(any(
    all(target_arch = "aarch64", target_endian = "little"),
    target_arch = "x86_64"
)))]
pub(super) fn offsets(bytes: &[u8]) -> memchr::Memchr<'_> {
    memchr::memchr_iter(b'\n', bytes)
}

#[cfg(all(target_arch = "aarch64", target_endian = "little"))]
#[inline]
fn newline_mask(bytes: &[u8; 64]) -> u64 {
    use std::arch::aarch64::*;
    // Four in-bounds 16-byte loads. Pairwise addition collects eight
    // disjoint bit weights into each byte of the resulting mask.
    unsafe {
        let weights = vcreate_u8(0x8040_2010_0804_0201);
        let weights = vcombine_u8(weights, weights);
        let mask = |offset| {
            let input = vld1q_u8(bytes.as_ptr().add(offset));
            vandq_u8(vceqq_u8(input, vdupq_n_u8(b'\n')), weights)
        };
        let pairs = vpaddq_u8(vpaddq_u8(mask(0), mask(16)), vpaddq_u8(mask(32), mask(48)));
        let octets = vpadd_u8(vget_low_u8(pairs), vget_high_u8(pairs));
        vget_lane_u64(vreinterpret_u64_u8(octets), 0)
    }
}

#[cfg(target_arch = "x86_64")]
#[inline]
fn newline_mask(bytes: &[u8; 64]) -> u64 {
    use std::arch::x86_64::*;
    // SSE2 is baseline on x86-64. Each load covers a full 16-byte chunk.
    unsafe {
        let newline = _mm_set1_epi8(b'\n' as i8);
        let mut mask = 0;
        for (i, chunk) in bytes.chunks_exact(16).enumerate() {
            let data = _mm_loadu_si128(chunk.as_ptr().cast());
            mask |= (_mm_movemask_epi8(_mm_cmpeq_epi8(data, newline)) as u64) << (i * 16);
        }
        mask
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_memchr_at_every_lane_and_tail_length() {
        let mut bytes = [b'X'; 200];
        for period in [1, 2, 7, 16, 39, 63, 64, 79, 127] {
            for (i, b) in bytes.iter_mut().enumerate() {
                *b = if i % period == 0 { b'\n' } else { b'X' };
            }
            for end in 0..=bytes.len() {
                for start in 0..4.min(end + 1) {
                    let data = &bytes[start..end];
                    assert!(offsets(data).eq(memchr::memchr_iter(b'\n', data)));
                }
            }
        }
    }
}
