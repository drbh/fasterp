//! Two-hash duplication evaluation for at most 512 bases. Split each prime
//! into three bytes so one dot product counts 16 bases per coefficient byte.
use std::arch::aarch64::*;

pub(super) struct PairHash {
    prime_bytes: [[[u8; 512]; 3]; 2],
    position_bias: [[u64; 2]; 513],
}

impl PairHash {
    pub(super) fn new(primes: &[u64]) -> Box<Self> {
        let mut table = Box::new(Self {
            prime_bytes: [[[0; 512]; 3]; 2],
            position_bias: [[0; 2]; 513],
        });
        for pos in 0..512 {
            for hash in 0..2 {
                let prime = primes[pos * 2 + hash];
                assert!(prime < 1 << 24);
                for byte in 0..3 {
                    table.prime_bytes[hash][byte][pos] = (prime >> (8 * byte)) as u8;
                }
                table.position_bias[pos + 1][hash] =
                    table.position_bias[pos][hash] + pos as u64 * prime;
            }
        }
        table
    }

    pub(super) fn positions(&self, r1: &[u8], r2: Option<&[u8]>, primes: &[u64]) -> [u64; 8] {
        if std::arch::is_aarch64_feature_detected!("dotprod") {
            unsafe { self.positions_dotprod(r1, r2, primes) }
        } else {
            self.positions_neon(r1, r2, primes)
        }
    }

    pub(super) fn positions_neon(&self, r1: &[u8], r2: Option<&[u8]>, primes: &[u64]) -> [u64; 8] {
        // NEON is available on the supported AArch64 targets.
        unsafe { self.accumulate::<false>(r1, r2, primes) }
    }

    #[target_feature(enable = "dotprod")]
    unsafe fn positions_dotprod(&self, r1: &[u8], r2: Option<&[u8]>, primes: &[u64]) -> [u64; 8] {
        unsafe { self.accumulate::<true>(r1, r2, primes) }
    }

    #[inline(always)]
    unsafe fn accumulate<const DOT: bool>(
        &self,
        r1: &[u8],
        r2: Option<&[u8]>,
        primes: &[u64],
    ) -> [u64; 8] {
        // TBL returns zero for out-of-range indices. Adding the default code
        // afterwards maps every byte except uppercase A/C/G/T to exactly 13.
        const DIFFERENCES: [u8; 32] = {
            let mut v = [0; 32];
            v[0] = 7u8.wrapping_sub(13);
            v[2] = 74 - 13;
            v[6] = 31 - 13;
            v[19] = 222 - 13;
            v
        };
        let len = r1.len() + r2.map_or(0, <[u8]>::len);
        let [first, second] = self.position_bias[len];
        let mut result = [first, second, 0, 0, 0, 0, 0, 0];
        unsafe {
            let table = vld1q_u8_x2(DIFFERENCES.as_ptr());
            let mut sums = [[vdupq_n_u32(0); 3]; 2];
            let mut offset = 0;
            for read in std::iter::once(r1).chain(r2) {
                let mut chunks = read.chunks_exact(16);
                for chunk in &mut chunks {
                    let indices = vsubq_u8(vld1q_u8(chunk.as_ptr()), vdupq_n_u8(b'A'));
                    let codes = vaddq_u8(vqtbl2q_u8(table, indices), vdupq_n_u8(13));
                    for (hash, parts) in sums.iter_mut().enumerate() {
                        for (byte, sum) in parts.iter_mut().enumerate() {
                            // len <= 512 and this chunk contains all 16 bytes.
                            let weights =
                                vld1q_u8(self.prime_bytes[hash][byte].as_ptr().add(offset));
                            if DOT {
                                // The stable Rust toolchain has no UDOT intrinsic.
                                std::arch::asm!(
                                    "udot {sum:v}.4s, {codes:v}.16b, {weights:v}.16b",
                                    sum = inout(vreg) *sum, codes = in(vreg) codes, weights = in(vreg) weights,
                                    options(pure, nomem, nostack)
                                );
                            } else {
                                let lo = vmull_u8(vget_low_u8(codes), vget_low_u8(weights));
                                let hi = vmull_high_u8(codes, weights);
                                *sum = vpadalq_u16(vpadalq_u16(*sum, lo), hi);
                            }
                        }
                    }
                    offset += 16;
                }
                for &base in chunks.remainder() {
                    let code = match base {
                        b'A' => 7,
                        b'C' => 74,
                        b'G' => 31,
                        b'T' => 222,
                        _ => 13,
                    };
                    for i in 0..2 {
                        result[i] += primes[offset * 2 + i] * code;
                    }
                    offset += 1;
                }
            }
            // Even the sum across lanes fits in u32 (512 * 222 * 255).
            for (hash, parts) in sums.into_iter().enumerate() {
                for (byte, sum) in parts.into_iter().enumerate() {
                    result[hash] += u64::from(vaddvq_u32(sum)) << (8 * byte);
                }
            }
        }
        result
    }
}
