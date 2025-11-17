//! Bloom filter based duplicate detection matching fastp's implementation
//!
//! Uses multiple hash functions with prime number based rolling hash
//! to detect exact duplicate reads/pairs with low memory footprint.

use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};

const PRIME_ARRAY_LEN: usize = 1 << 9; // 512 primes per buffer

/// Duplicate detector using Bloom filter
/// Memory layout and algorithm matches fastp for exact parity
pub struct DuplicateDetector {
    /// Bit arrays for Bloom filter (one per buffer)
    /// Each byte stores 8 bits, accessed atomically
    bit_arrays: Vec<Vec<AtomicU8>>,
    /// Buffer size in bits (for modulo operations)
    buf_len_in_bits: u64,
    /// Number of hash buffers (determines number of hash functions)
    buf_num: usize,
    /// Prime numbers for hashing (`buf_num` * `PRIME_ARRAY_LEN`)
    prime_arrays: Vec<u64>,
    /// Mask for prime array indexing
    offset_mask: usize,
    /// Total reads processed
    total_reads: AtomicUsize,
    /// Duplicate reads found
    duplicate_reads: AtomicUsize,
}

impl DuplicateDetector {
    /// Create new duplicate detector with specified accuracy level
    /// Levels 1-6 use 1GB, 2GB, 4GB, 8GB, 16GB, 24GB respectively
    pub fn new(accuracy_level: u8) -> Self {
        // Base: 512MB per buffer
        let mut buf_len_in_bytes: usize = 1 << 29; // 536870912 bytes = 512MB
        let mut buf_num: usize = 2;

        // Scale memory based on accuracy level (matching fastp)
        match accuracy_level {
            2 => {
                buf_len_in_bytes *= 2;
            } // 1GB × 2 = 2GB
            3 => {
                buf_len_in_bytes *= 2;
                buf_num *= 2;
            } // 1GB × 4 = 4GB
            4 => {
                buf_len_in_bytes *= 4;
                buf_num *= 2;
            } // 2GB × 4 = 8GB
            5 => {
                buf_len_in_bytes *= 8;
                buf_num *= 2;
            } // 4GB × 4 = 16GB
            6 => {
                buf_len_in_bytes *= 8;
                buf_num = 3;
            } // 4GB × 6 = 24GB
            _ => {} // Level 1 and default: 512MB × 2 = 1GB
        }

        let buf_len_in_bits = (buf_len_in_bytes as u64) << 3; // Convert bytes to bits
        let offset_mask = PRIME_ARRAY_LEN * buf_num - 1;

        // Allocate bit arrays (one per buffer)
        let mut bit_arrays = Vec::with_capacity(buf_num);
        for _ in 0..buf_num {
            let mut buf = Vec::with_capacity(buf_len_in_bytes);
            for _ in 0..buf_len_in_bytes {
                buf.push(AtomicU8::new(0));
            }
            bit_arrays.push(buf);
        }

        // Generate prime numbers
        let prime_arrays = Self::init_prime_arrays(buf_num);

        Self {
            bit_arrays,
            buf_len_in_bits,
            buf_num,
            prime_arrays,
            offset_mask,
            total_reads: AtomicUsize::new(0),
            duplicate_reads: AtomicUsize::new(0),
        }
    }

    /// Generate prime numbers for hashing (matches fastp's initPrimeArrays)
    fn init_prime_arrays(buf_num: usize) -> Vec<u64> {
        let total_primes = buf_num * PRIME_ARRAY_LEN;
        let mut primes = Vec::with_capacity(total_primes);
        let mut number: u64 = 10000;

        while primes.len() < total_primes {
            number += 1;
            if Self::is_prime(number) {
                primes.push(number);
                number += 10000; // Jump ahead for next prime
            }
        }

        primes
    }

    /// Simple primality test
    fn is_prime(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        let sqrt_n = (n as f64).sqrt() as u64;
        for i in 2..=sqrt_n {
            if n % i == 0 {
                return false;
            }
        }
        true
    }

    /// Hash sequence to integer vector (matches fastp's seq2intvector)
    /// For paired-end, call twice: once for R1 (`pos_offset=0`) and once for R2 (`pos_offset=R1.len()`)
    fn seq2intvector(&self, data: &[u8], pos_offset: usize, output: &mut [u64]) {
        for (p, &base_char) in data.iter().enumerate() {
            let base: u64 = match base_char {
                b'A' => 7,
                b'T' => 222,
                b'C' => 74,
                b'G' => 31,
                _ => 13, // N or any other character
            };

            for (i, output_item) in output.iter_mut().take(self.buf_num).enumerate() {
                let offset = ((p + pos_offset) * self.buf_num + i) & self.offset_mask;
                *output_item += self.prime_arrays[offset] * (base + (p + pos_offset) as u64);
            }
        }
    }

    /// Apply Bloom filter: check if all bits are set, and set them if not
    /// Returns true if this is a duplicate (all bits were already set)
    fn apply_bloom_filter(&self, positions: &[u64]) -> bool {
        let mut is_dup = true;

        for (i, &position) in positions.iter().enumerate().take(self.buf_num) {
            let pos = position % self.buf_len_in_bits;
            let byte_pos = (pos >> 3) as usize; // Divide by 8 to get byte position
            let bit_offset = (pos & 0x07) as u8; // Remainder is bit position within byte
            let bit_mask = 1u8 << bit_offset;

            // Atomically fetch old value and set the bit
            let old_value = self.bit_arrays[i][byte_pos].fetch_or(bit_mask, Ordering::Relaxed);

            // Check if bit was already set
            is_dup = is_dup && (old_value & bit_mask) != 0;
        }

        is_dup
    }

    /// Check if a paired-end read is duplicate
    /// Hashes R1 and R2 concatenated together
    pub fn check_pair(&self, seq1: &[u8], seq2: &[u8]) -> bool {
        let mut positions = vec![0u64; self.buf_num];

        // Hash R1 starting at position 0
        self.seq2intvector(seq1, 0, &mut positions);
        // Hash R2 starting at position R1.len()
        self.seq2intvector(seq2, seq1.len(), &mut positions);

        let is_dup = self.apply_bloom_filter(&positions);

        self.total_reads.fetch_add(1, Ordering::Relaxed);
        if is_dup {
            self.duplicate_reads.fetch_add(1, Ordering::Relaxed);
        }

        is_dup
    }

    /// Get duplication rate
    pub fn get_dup_rate(&self) -> f64 {
        let total = self.total_reads.load(Ordering::Relaxed);
        if total == 0 {
            return 0.0;
        }
        let dups = self.duplicate_reads.load(Ordering::Relaxed);
        dups as f64 / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paired_duplicate_detection() {
        let detector = DuplicateDetector::new(1);

        let r1 = b"ATCGATCGATCGATCGATCG";
        let r2 = b"GCTAGCTAGCTAGCTAGCTA";
        let r1_diff = b"AAAAAAAAAAAAAAAAAAAA";

        // First pair should not be duplicate
        assert!(!detector.check_pair(r1, r2));

        // Same pair should be duplicate
        assert!(detector.check_pair(r1, r2));

        // Different R1 with same R2 should not be duplicate
        assert!(!detector.check_pair(r1_diff, r2));
    }

    #[test]
    fn test_prime_generation() {
        let primes = DuplicateDetector::init_prime_arrays(2);
        assert_eq!(primes.len(), 2 * PRIME_ARRAY_LEN);

        // Check first few primes are actually prime
        for &p in &primes[..10] {
            assert!(DuplicateDetector::is_prime(p));
        }
    }
}
