//! Bloom filter based duplicate detection matching fastp's implementation
//!
//! Uses multiple hash functions with prime number based rolling hash
//! to detect exact duplicate reads/pairs with low memory footprint.

use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};

const PRIME_ARRAY_LEN: usize = 1 << 9; // 512 primes per buffer

// DATA ORIENTED OPTIMIZATION:
// A static Lookup Table (LUT) maps DNA bases to their Fastp integer values.
// This allows O(1) conversion without 'if/else' or 'match' branching,
// preventing pipeline stalls on the CPU.
static BASE_LUT: [u64; 256] = {
    let mut lut = [13u64; 256]; // Default to 13 (N/other)
    lut[b'A' as usize] = 7;
    lut[b'T' as usize] = 222;
    lut[b'C' as usize] = 74;
    lut[b'G' as usize] = 31;
    lut
};

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

    /// Optimized hash accumulation with linear memory access
    ///
    /// Changes from original seq2intvector:
    /// 1. Uses LUT for branchless base mapping (no match/if)
    /// 2. Linear memory access pattern for CPU prefetching
    /// 3. Uses unsafe get_unchecked (bounds are mathematically guaranteed by mask)
    #[inline(always)]
    fn accumulate_hash(&self, seq: &[u8], start_offset: usize, accumulators: &mut [u64]) {
        let buf_num = self.buf_num;
        let mask = self.offset_mask;

        // Calculate the starting index in the prime array.
        // We track this linearly to avoid multiplication in the loop.
        let mut current_prime_base_idx = start_offset * buf_num;

        for &byte in seq {
            // 1. LUT Lookup (No Branching)
            // SAFETY: byte is u8, so it is always a valid index for [u64; 256]
            let base = unsafe { *BASE_LUT.get_unchecked(byte as usize) };

            // Calculate position component once per base
            let pos = current_prime_base_idx / buf_num;
            let val_component = base + pos as u64;

            // 2. Update Accumulators
            // We iterate buf_num times (2-6). LLVM will unroll this small loop.
            for i in 0..buf_num {
                // Handle the circular buffer wrap-around using bitwise AND
                let prime_idx = (current_prime_base_idx + i) & mask;

                // SAFETY: prime_arrays length matches mask+1, and mask guarantees bounds
                let prime = unsafe { *self.prime_arrays.get_unchecked(prime_idx) };

                // Update the accumulator
                accumulators[i] += prime * val_component;
            }

            // Stride forward by buf_num (Linear Memory Access)
            current_prime_base_idx += buf_num;
        }
    }

    /// Check if a paired-end read is duplicate
    ///
    /// Optimized to use:
    /// - Stack allocation instead of heap (Vec -> fixed array)
    /// - Linear memory access pattern
    /// - Inlined Bloom filter logic
    pub fn check_pair(&self, seq1: &[u8], seq2: &[u8]) -> bool {
        // OPTIMIZATION: Stack Allocation
        // Instead of vec![0u64; self.buf_num] (heap), we use a fixed array.
        // Max buf_num is 6 (accuracy level 6). 16 provides safety margin.
        let mut positions = [0u64; 16];

        // Slice to actual size needed
        let active_accs = &mut positions[..self.buf_num];

        // Hash R1 starting at position 0
        self.accumulate_hash(seq1, 0, active_accs);

        // Hash R2 starting at position seq1.len()
        self.accumulate_hash(seq2, seq1.len(), active_accs);

        // Apply Bloom filter (inlined for reduced overhead)
        let mut is_dup = true;

        for (i, &position) in active_accs.iter().enumerate() {
            let pos = position % self.buf_len_in_bits;
            let byte_pos = (pos >> 3) as usize;
            let bit_mask = 1u8 << (pos & 0x07);

            // SAFETY: bit_arrays are initialized to buf_len_in_bytes.
            // The modulo arithmetic guarantees we are in bounds.
            let bit_array = unsafe { self.bit_arrays.get_unchecked(i) };
            let byte_ptr = unsafe { bit_array.get_unchecked(byte_pos) };

            let old_value = byte_ptr.fetch_or(bit_mask, Ordering::Relaxed);

            // Match fastp: overwrite (not AND) - only last buffer's result matters
            is_dup = (old_value & bit_mask) != 0;
        }

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
