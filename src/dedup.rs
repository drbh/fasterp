//! Sequence-based deduplication module
//!
//! This module implements fastp-compatible deduplication using exact sequence matching.
//! Reads are considered duplicates if ALL base pairs are identical (hash-based comparison).
//!
//! Note: This is INDEPENDENT from UMI processing. fastp's --dedup uses sequence matching,
//! not UMI-aware deduplication.

use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Configuration for deduplication
#[derive(Debug, Clone)]
pub struct DedupConfig {
    /// Enable deduplication
    pub enabled: bool,
    /// Accuracy level (higher = more memory, lower false positive rate)
    /// Default: 5 (matches fastp default)
    pub accuracy: u8,
}

impl Default for DedupConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            accuracy: 5,
        }
    }
}

/// Deduplication tracker using hash-based exact matching
///
/// For single-end: tracks sequence hashes
/// For paired-end: tracks combined (seq1, seq2) hashes
pub struct DedupTracker {
    /// Set of seen sequence hashes
    seen_hashes: HashSet<u64>,
    /// Number of duplicates found
    pub duplicates: usize,
    /// Accuracy level (affects hash precision)
    accuracy: u8,
}

impl DedupTracker {
    /// Create new deduplication tracker
    pub fn new(accuracy: u8) -> Self {
        Self {
            seen_hashes: HashSet::new(),
            duplicates: 0,
            accuracy,
        }
    }

    /// Check if single-end read is a duplicate
    /// Returns true if this is a duplicate (already seen)
    pub fn is_duplicate_se(&mut self, seq: &[u8]) -> bool {
        let hash = self.hash_sequence(seq);
        if self.seen_hashes.contains(&hash) {
            self.duplicates += 1;
            true
        } else {
            self.seen_hashes.insert(hash);
            false
        }
    }

    /// Check if paired-end read pair is a duplicate
    /// Returns true if this pair is a duplicate (already seen)
    pub fn is_duplicate_pe(&mut self, seq1: &[u8], seq2: &[u8]) -> bool {
        let hash = self.hash_pair(seq1, seq2);
        if self.seen_hashes.contains(&hash) {
            self.duplicates += 1;
            true
        } else {
            self.seen_hashes.insert(hash);
            false
        }
    }

    /// Hash a single sequence
    fn hash_sequence(&self, seq: &[u8]) -> u64 {
        let mut hasher = DefaultHasher::new();

        // Use only part of sequence based on accuracy level to reduce memory
        // Higher accuracy = use more of the sequence
        let sample_size = self.calculate_sample_size(seq.len());
        let sample_seq = if seq.len() <= sample_size {
            seq
        } else {
            &seq[..sample_size]
        };

        sample_seq.hash(&mut hasher);
        hasher.finish()
    }

    /// Hash a pair of sequences (for paired-end)
    fn hash_pair(&self, seq1: &[u8], seq2: &[u8]) -> u64 {
        let mut hasher = DefaultHasher::new();

        let sample_size1 = self.calculate_sample_size(seq1.len());
        let sample_size2 = self.calculate_sample_size(seq2.len());

        let sample_seq1 = if seq1.len() <= sample_size1 {
            seq1
        } else {
            &seq1[..sample_size1]
        };

        let sample_seq2 = if seq2.len() <= sample_size2 {
            seq2
        } else {
            &seq2[..sample_size2]
        };

        sample_seq1.hash(&mut hasher);
        sample_seq2.hash(&mut hasher);
        hasher.finish()
    }

    /// Calculate how much of the sequence to use based on accuracy level
    /// Higher accuracy = larger sample = more memory but fewer false positives
    fn calculate_sample_size(&self, seq_len: usize) -> usize {
        // accuracy 1-6 maps to using different fractions of the sequence
        // fastp default is 5
        let fraction = match self.accuracy {
            1 => 16, // 1/16 of sequence
            2 => 8,  // 1/8 of sequence
            3 => 4,  // 1/4 of sequence
            4 => 2,  // 1/2 of sequence
            5 => 1,  // full sequence (fastp default)
            _ => 1,  // full sequence for accuracy >= 6
        };

        if fraction == 1 {
            seq_len
        } else {
            std::cmp::max(seq_len / fraction, 32) // minimum 32 bases
        }
    }

    /// Get number of unique reads seen
    pub fn unique_count(&self) -> usize {
        self.seen_hashes.len()
    }

    /// Clear the tracker (for memory management in streaming)
    pub fn clear(&mut self) {
        self.seen_hashes.clear();
        // Keep duplicates count
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dedup_single_end_exact_duplicate() {
        let mut tracker = DedupTracker::new(5);

        let seq = b"ACGTACGTACGTACGT";

        // First occurrence - not a duplicate
        assert!(!tracker.is_duplicate_se(seq));

        // Second occurrence - is a duplicate
        assert!(tracker.is_duplicate_se(seq));

        assert_eq!(tracker.duplicates, 1);
        assert_eq!(tracker.unique_count(), 1);
    }

    #[test]
    fn test_dedup_single_end_different_sequences() {
        let mut tracker = DedupTracker::new(5);

        let seq1 = b"ACGTACGTACGTACGT";
        let seq2 = b"GGGGCCCCAAAATTTT";

        assert!(!tracker.is_duplicate_se(seq1));
        assert!(!tracker.is_duplicate_se(seq2));

        assert_eq!(tracker.duplicates, 0);
        assert_eq!(tracker.unique_count(), 2);
    }

    #[test]
    fn test_dedup_paired_end_exact_duplicate() {
        let mut tracker = DedupTracker::new(5);

        let seq1 = b"ACGTACGTACGTACGT";
        let seq2 = b"TTTTAAAACCCCGGGG";

        // First occurrence - not a duplicate
        assert!(!tracker.is_duplicate_pe(seq1, seq2));

        // Second occurrence - is a duplicate
        assert!(tracker.is_duplicate_pe(seq1, seq2));

        assert_eq!(tracker.duplicates, 1);
        assert_eq!(tracker.unique_count(), 1);
    }

    #[test]
    fn test_dedup_paired_end_different_r1() {
        let mut tracker = DedupTracker::new(5);

        let seq1a = b"ACGTACGTACGTACGT";
        let seq1b = b"GGGGACGTACGTACGT"; // Different R1
        let seq2 = b"TTTTAAAACCCCGGGG";

        assert!(!tracker.is_duplicate_pe(seq1a, seq2));
        assert!(!tracker.is_duplicate_pe(seq1b, seq2)); // Different pair

        assert_eq!(tracker.duplicates, 0);
        assert_eq!(tracker.unique_count(), 2);
    }

    #[test]
    fn test_dedup_paired_end_different_r2() {
        let mut tracker = DedupTracker::new(5);

        let seq1 = b"ACGTACGTACGTACGT";
        let seq2a = b"TTTTAAAACCCCGGGG";
        let seq2b = b"TTTTAAAAGGGGCCCC"; // Different R2

        assert!(!tracker.is_duplicate_pe(seq1, seq2a));
        assert!(!tracker.is_duplicate_pe(seq1, seq2b)); // Different pair

        assert_eq!(tracker.duplicates, 0);
        assert_eq!(tracker.unique_count(), 2);
    }

    #[test]
    fn test_dedup_accuracy_levels() {
        // Test that different accuracy levels work
        for accuracy in 1..=6 {
            let mut tracker = DedupTracker::new(accuracy);
            let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";

            assert!(!tracker.is_duplicate_se(seq));
            assert!(tracker.is_duplicate_se(seq));
            assert_eq!(tracker.duplicates, 1);
        }
    }

    #[test]
    fn test_default_config() {
        let config = DedupConfig::default();
        assert!(!config.enabled);
        assert_eq!(config.accuracy, 5);
    }
}
