//! Matching semantics updated against OpenGene/fastp v1.4.0; see LICENSE.fastp.
//! Adapter trimming functionality
//!
//! This module provides adapter detection and trimming for FASTQ preprocessing:
//! - Manual adapter specification
//! - Built-in common adapter sequences
//! - Overlap-based adapter detection for paired-end reads
//! - Auto-detection using PE overlap analysis
//! - Mismatch-tolerant matching

use std::collections::BTreeMap;

/// Configuration for adapter trimming
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// Adapter sequence for read1/single-end
    pub adapter_seq: Option<Vec<u8>>,
    /// Adapter sequence for read2 (paired-end only)
    pub adapter_seq_r2: Option<Vec<u8>>,
    /// Enable auto-detection for paired-end
    pub detect_adapter_for_pe: bool,
    /// Minimum overlap length for detection
    pub min_overlap: usize,
    /// Maximum mismatches allowed in overlap
    pub max_mismatches: usize,
}

impl AdapterConfig {
    pub fn new() -> Self {
        Self {
            adapter_seq: None,
            adapter_seq_r2: None,
            detect_adapter_for_pe: false, // fastp does NOT auto-detect for PE reads, only uses overlap-based trimming
            min_overlap: 4,               // fastp's default for adapter trimming
            max_mismatches: 2,
        }
    }

    pub fn is_enabled(&self) -> bool {
        self.adapter_seq.is_some() || self.adapter_seq_r2.is_some() || self.detect_adapter_for_pe
    }
}

/// Result of adapter auto-detection
#[derive(Debug, Clone)]
pub struct AdapterDetectionResult {
    /// Detected adapter for read 1
    pub adapter_r1: Option<Vec<u8>>,
    /// Detected adapter for read 2
    pub adapter_r2: Option<Vec<u8>>,
    /// Confidence score (fraction of reads supporting this adapter)
    pub confidence: f64,
    /// Number of reads analyzed
    pub reads_analyzed: usize,
}

/// Detect adapters from paired-end reads using overlap analysis
///
/// This function samples read pairs and uses overlap detection to identify
/// adapter sequences. When R1 and R2 overlap, any bases beyond the overlap
/// are considered adapter contamination.
///
/// # Arguments
/// * `r1_seqs` - Iterator of R1 sequences
/// * `r2_seqs` - Iterator of R2 sequences
/// * `sample_size` - Maximum number of read pairs to analyze
/// * `min_frequency` - Minimum fraction of reads that must support an adapter (default: 0.01 = 1%)
///
/// # Returns
/// `AdapterDetectionResult` with detected adapters and confidence scores
pub fn detect_adapters_from_pe_reads<'a>(
    r1_seqs: impl Iterator<Item = &'a [u8]>,
    r2_seqs: impl Iterator<Item = &'a [u8]>,
    sample_size: usize,
    min_frequency: f64,
) -> AdapterDetectionResult {
    use crate::overlap::{detect_overlap, reverse_complement, OverlapConfig};

    let overlap_config = OverlapConfig {
        min_overlap_len: 10, // Lower threshold for detection
        max_diff: 5,
        max_diff_percent: 20,
    };

    // Use BTreeMap instead of HashMap for deterministic iteration order
    let mut adapter_r1_counts: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    let mut adapter_r2_counts: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    let mut total_reads = 0;
    let mut reads_with_adapters = 0;

    // Sample read pairs and extract adapters
    for (r1_seq, r2_seq) in r1_seqs.zip(r2_seqs).take(sample_size) {
        total_reads += 1;

        // Detect overlap between R1 and R2
        if let Some(overlap_result) = detect_overlap(r1_seq, r2_seq, &overlap_config) {
            let r2_rc = reverse_complement(r2_seq);

            // Extract adapters based on offset direction
            if overlap_result.offset >= 0 {
                // R1 extends past R2: adapter is at end of R1
                // SAFETY: offset is guaranteed to be non-negative by the if condition
                let r1_adapter_start = usize::try_from(overlap_result.offset).unwrap_or(0)
                    + overlap_result.overlap_len;
                if r1_adapter_start < r1_seq.len() {
                    let adapter = r1_seq[r1_adapter_start..].to_vec();
                    if !adapter.is_empty() && adapter.len() >= 3 {
                        *adapter_r1_counts.entry(adapter).or_insert(0) += 1;
                        reads_with_adapters += 1;
                    }
                }
            } else {
                // R2 extends past R1: adapters are at ends of both reads
                // Adapter in R1 is at the end (beyond overlap)
                if overlap_result.overlap_len < r1_seq.len() {
                    let adapter = r1_seq[overlap_result.overlap_len..].to_vec();
                    if !adapter.is_empty() && adapter.len() >= 3 {
                        *adapter_r1_counts.entry(adapter).or_insert(0) += 1;
                        reads_with_adapters += 1;
                    }
                }

                // Adapter in R2_rc is at the end (beyond overlap)
                // SAFETY: offset is negative here, so negating makes it positive
                let r2_offset = usize::try_from(-overlap_result.offset).unwrap_or(0);
                let r2_adapter_start = r2_offset + overlap_result.overlap_len;
                if r2_adapter_start < r2_rc.len() {
                    let adapter_rc = r2_rc[r2_adapter_start..].to_vec();
                    if !adapter_rc.is_empty() && adapter_rc.len() >= 3 {
                        // Store the original R2 adapter (reverse complement of what we extracted)
                        let adapter = reverse_complement(&adapter_rc);
                        *adapter_r2_counts.entry(adapter).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    // Find the most common adapters
    // Threshold calculation: total_reads * min_frequency
    // Using multiplication then floor division to avoid floating point
    let min_freq_percent = (min_frequency * 1000.0) as u64;
    let threshold = ((total_reads as u64 * min_freq_percent) / 1000) as usize;

    let adapter_r1 = find_consensus_adapter(&adapter_r1_counts, threshold);
    let adapter_r2 = find_consensus_adapter(&adapter_r2_counts, threshold);

    let confidence = if total_reads > 0 {
        f64::from(reads_with_adapters) / total_reads as f64
    } else {
        0.0
    };

    AdapterDetectionResult {
        adapter_r1,
        adapter_r2,
        confidence,
        reads_analyzed: total_reads,
    }
}

/// Find consensus adapter from frequency counts
fn find_consensus_adapter(
    adapter_counts: &BTreeMap<Vec<u8>, usize>,
    threshold: usize,
) -> Option<Vec<u8>> {
    if adapter_counts.is_empty() {
        return None;
    }

    // Find the most frequent adapter
    let mut best_adapter: Option<Vec<u8>> = None;
    let mut best_count = threshold;

    for (adapter, &count) in adapter_counts {
        if count > best_count {
            best_count = count;
            best_adapter = Some(adapter.clone());
        }
    }

    // Return the most frequently detected adapter without extending it
    // This matches fastp's behavior which uses the detected fragments directly
    best_adapter
}

//
// ========== OPTIMIZED ADAPTER DETECTION (2-BIT FLAT COUNTER) ==========
//

mod optimized {
    #[cfg(feature = "native")]
    use rayon::prelude::*;

    /// 2-bit encoding for DNA bases: A=00, C=01, G=10, T=11
    #[inline]
    const fn base_to_2bit(base: u8) -> Option<u8> {
        match base {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None, // N or other
        }
    }

    /// Rolling k-mer encoder for K=10 (20 bits in u32)
    struct RollingKmerEncoder {
        kmer_bits: u32,
        kmer_len: usize,
        k: usize,
        mask: u32,
    }

    impl RollingKmerEncoder {
        #[inline]
        fn new(k: usize) -> Self {
            assert!(k <= 16, "K must be <= 16 for u32 encoding");
            let mask = (1u32 << (k * 2)) - 1;
            Self {
                kmer_bits: 0,
                kmer_len: 0,
                k,
                mask,
            }
        }

        #[inline]
        fn push(&mut self, base: u8) -> Option<u32> {
            if let Some(bits) = base_to_2bit(base) {
                self.kmer_bits = ((self.kmer_bits << 2) | u32::from(bits)) & self.mask;
                self.kmer_len += 1;
                if self.kmer_len >= self.k {
                    Some(self.kmer_bits)
                } else {
                    None
                }
            } else {
                self.kmer_len = 0;
                self.kmer_bits = 0;
                None
            }
        }

        #[inline]
        fn reset(&mut self) {
            self.kmer_len = 0;
            self.kmer_bits = 0;
        }
    }

    /// Fast filters using both byte slice and pre-computed 2-bit value
    #[inline(always)]
    fn should_skip_kmer(kmer_bytes: &[u8], kmer_bits: u32) -> bool {
        const K: usize = 10;

        // Filter 1: Low complexity (runs of 5+ same base)
        let mut run_len = 1;
        let mut prev_bits = kmer_bits & 0b11;
        for i in 1..K {
            let curr_bits = (kmer_bits >> (i * 2)) & 0b11;
            if curr_bits == prev_bits {
                run_len += 1;
                if run_len >= 5 {
                    return true;
                }
            } else {
                run_len = 1;
            }
            prev_bits = curr_bits;
        }

        // Filter 2: GC-rich (GC count >= K-2)
        let gc_count = kmer_bytes
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
            .count();
        if gc_count >= K - 2 {
            return true;
        }

        // Filter 3: Diversity (at least 3 adjacent different bases)
        let mut diff_count = 0;
        for i in 0..K - 1 {
            let curr = (kmer_bits >> (i * 2)) & 0b11;
            let next = (kmer_bits >> ((i + 1) * 2)) & 0b11;
            if curr != next {
                diff_count += 1;
            }
        }
        if diff_count < 3 {
            return true;
        }

        false
    }

    /// Count K=10 k-mers from read tails using 2-bit encoding (single-threaded)
    fn count_k10_2bit_tail(seqs: &[&[u8]], start: usize) -> Vec<u32> {
        const K: usize = 10;
        const ARRAY_SIZE: usize = 1 << (K * 2); // 4^10 = 1,048,576

        let mut counts = vec![0u32; ARRAY_SIZE];
        let mut encoder = RollingKmerEncoder::new(K);

        for seq in seqs {
            if seq.len() < start + K {
                continue;
            }

            encoder.reset();

            for (pos, &base) in seq[start..].iter().enumerate() {
                if let Some(kmer_idx) = encoder.push(base) {
                    let kmer_slice = &seq[(start + pos - K + 1)..=(start + pos)];

                    if should_skip_kmer(kmer_slice, kmer_idx) {
                        continue;
                    }

                    counts[kmer_idx as usize] = counts[kmer_idx as usize].saturating_add(1);
                }
            }
        }

        counts
    }

    /// Parallel k-mer counting using rayon (native) or serial (WASM)
    #[cfg(feature = "native")]
    pub fn count_k10_2bit_parallel(seqs: &[&[u8]], start: usize) -> Vec<u32> {
        const K: usize = 10;
        const ARRAY_SIZE: usize = 1 << (K * 2);
        const CHUNK_SIZE: usize = 10_000;

        let thread_counts: Vec<Vec<u32>> = seqs
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| count_k10_2bit_tail(chunk, start))
            .collect();

        let mut final_counts = vec![0u32; ARRAY_SIZE];
        for thread_count in thread_counts {
            for (i, count) in thread_count.iter().enumerate() {
                final_counts[i] = final_counts[i].saturating_add(*count);
            }
        }

        final_counts
    }

    /// Serial k-mer counting for WASM
    #[cfg(not(feature = "native"))]
    pub fn count_k10_2bit_parallel(seqs: &[&[u8]], start: usize) -> Vec<u32> {
        // Fall back to serial processing on WASM
        count_k10_2bit_tail(seqs, start)
    }

    /// Find top N k-mer candidates by count
    pub fn find_top_kmers(counts: &[u32], top_n: usize, threshold: u32) -> Vec<(u32, u32)> {
        let mut candidates: Vec<(u32, u32)> = counts
            .iter()
            .enumerate()
            .filter(|&(_, count)| *count >= threshold)
            .map(|(idx, count)| (idx as u32, *count))
            .collect();

        // Sort by count DESC, then by k-mer index ASC for deterministic tiebreaking
        candidates.sort_by(|a, b| match b.1.cmp(&a.1) {
            std::cmp::Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        });
        candidates.truncate(top_n);
        candidates
    }

    /// Known adapter with pre-encoded windows
    pub struct KnownAdapter {
        pub full_seq: &'static [u8],
        pub encoded_windows: Vec<u32>,
    }

    /// Encode all K=10 windows in a sequence
    fn encode_all_windows(seq: &[u8], k: usize) -> Vec<u32> {
        let mut encoder = RollingKmerEncoder::new(k);
        let mut windows = Vec::new();

        for &base in seq {
            if let Some(kmer_bits) = encoder.push(base) {
                windows.push(kmer_bits);
            }
        }

        windows
    }

    /// Pre-encoded known adapters (lazy initialization)
    pub static KNOWN_ADAPTERS: std::sync::LazyLock<Vec<KnownAdapter>> =
        std::sync::LazyLock::new(|| {
            let adapter_map = crate::adapter::adapters::get_known_adapters();
            let mut result = Vec::new();

            // Sort keys for deterministic iteration order
            let mut keys: Vec<_> = adapter_map.keys().collect();
            keys.sort();

            for seq_str in keys {
                let seq = seq_str.as_bytes();
                let windows = encode_all_windows(seq, 10);

                result.push(KnownAdapter {
                    full_seq: seq,
                    encoded_windows: windows,
                });
            }

            result
        });

    /// Match candidate k-mer against known adapters
    pub fn match_known_adapter(kmer_bits: u32) -> Option<&'static KnownAdapter> {
        KNOWN_ADAPTERS
            .iter()
            .find(|&adapter| adapter.encoded_windows.contains(&kmer_bits))
    }
}

/// Detect adapters from single-end reads using optimized 2-bit flat counter
///
/// This function uses a 2-bit encoded flat array for k-mer counting, which is much
/// faster than the HashMap-based approach. It parallelizes counting with rayon.
///
/// # Arguments
/// * `seqs` - Iterator of sequences
/// * `sample_size` - Maximum number of reads to analyze
/// * `min_frequency` - Minimum fraction of reads that must contain a k-mer
///
/// # Returns
/// `AdapterDetectionResult` with detected adapter
pub fn detect_adapter_from_se_reads<'a>(
    seqs: impl Iterator<Item = &'a [u8]>,
    sample_size: usize,
    _min_frequency: f64,
) -> AdapterDetectionResult {
    use optimized::{count_k10_2bit_parallel, find_top_kmers, match_known_adapter};

    const K: usize = 10;
    const MIN_READS_FOR_DETECTION: usize = 10000;
    const FOLD_THRESHOLD: usize = 20;
    const START_POS: usize = 20;

    // Collect sequences into a vector for parallel processing
    let seq_vec: Vec<&[u8]> = seqs.take(sample_size).collect();
    let total_reads = seq_vec.len();

    // Fastp requires at least 10,000 reads for detection
    if total_reads < MIN_READS_FOR_DETECTION {
        return AdapterDetectionResult {
            adapter_r1: None,
            adapter_r2: None,
            confidence: 0.0,
            reads_analyzed: total_reads,
        };
    }

    // Count k-mers using parallel 2-bit flat counter
    let counts = count_k10_2bit_parallel(&seq_vec, START_POS);

    // Calculate total k-mers counted
    let total_kmers: u64 = counts.iter().map(|&c| u64::from(c)).sum();

    if total_kmers == 0 {
        return AdapterDetectionResult {
            adapter_r1: None,
            adapter_r2: None,
            confidence: 0.0,
            reads_analyzed: total_reads,
        };
    }

    // Calculate fold threshold like fastp
    let size = 1usize << (K * 2); // 4^10 = 1048576 possible k-mers
    let fold_threshold = ((total_kmers as usize * FOLD_THRESHOLD) / size).max(10) as u32;

    // Find top k-mers that meet the threshold
    let top_kmers = find_top_kmers(&counts, 100, fold_threshold);

    if top_kmers.is_empty() {
        return AdapterDetectionResult {
            adapter_r1: None,
            adapter_r2: None,
            confidence: 0.0,
            reads_analyzed: total_reads,
        };
    }

    // Try to match top k-mers against known adapters
    for (kmer_bits, count) in &top_kmers {
        if let Some(adapter) = match_known_adapter(*kmer_bits) {
            return AdapterDetectionResult {
                adapter_r1: Some(adapter.full_seq.to_vec()),
                adapter_r2: None,
                confidence: f64::from(*count) / total_reads as f64,
                reads_analyzed: total_reads,
            };
        }
    }

    // No adapter detected
    AdapterDetectionResult {
        adapter_r1: None,
        adapter_r2: None,
        confidence: 0.0,
        reads_analyzed: total_reads,
    }
}

/// Built-in Illumina adapter sequences
pub mod adapters {
    use std::collections::HashMap;
    use std::sync::OnceLock;

    /// Add Illumina `TruSeq` and basic adapters
    fn add_truseq_adapters(map: &mut HashMap<&'static str, &'static str>) {
        map.insert(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            ">Illumina TruSeq Adapter Read 1",
        );
        map.insert(
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            ">Illumina TruSeq Adapter Read 2",
        );
        map.insert(
            "GATCGTCGGACTGTAGAACTCTGAACGTGTAGA",
            ">Illumina Small RNA Adapter Read 2",
        );
        map.insert("AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA", ">Illumina DpnII expression PCR Primer 2 | >Illumina NlaIII expression PCR Primer 2 | >Illumina Small RNA PCR Primer 2 | >Illumina DpnII Gex PCR Primer 2 | >Illumina NlaIII Gex PCR Primer 2");
        map.insert(
            "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
            ">Illumina RNA PCR Primer",
        );
        map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1 | >Illumina Single End PCR Primer 1 | >Illumina Paried End PCR Primer 1 | >Illumina Multiplexing PCR Primer 1.01 | >TruSeq Universal Adapter | >TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1");
        map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">pcr_dimer");
        map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT", ">PCR_Primers");
        map.insert("ACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">Illumina Single End Sequencing Primer | >Illumina Paired End Adapter 1 | >Illumina Paried End Sequencing Primer 1 | >Illumina Multiplexing Adapter 2 | >Illumina Multiplexing Read1 Sequencing Primer");
        map.insert(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            ">PE2_rc | >TruSeq3_IndexedAdapter | >PE2_rc | >TruSeq3_IndexedAdapter",
        );
        map.insert(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
            ">Reverse_adapter",
        );
        map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", ">TruSeq2_PE_r");
        map.insert(
            "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
            ">PCR_Primer2_rc",
        );
        map.insert(
            "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA",
            ">PhiX_read1_adapter",
        );
        map.insert(
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
            ">PE1_rc | >TruSeq3_UniversalAdapter | >PE1_rc | >TruSeq3_UniversalAdapter",
        );
        map.insert(
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
            ">PCR_Primer1_rc",
        );
        map.insert(
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA",
            ">PhiX_read2_adapter",
        );
        map.insert("AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", ">TruSeq2_SE");
    }

    /// Add Illumina RNA PCR Primer Index adapters
    #[allow(clippy::too_many_lines)]
    fn add_rna_pcr_index_adapters(map: &mut HashMap<&'static str, &'static str>) {
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 35",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 10",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 10",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 2",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 2",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 38",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 27",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 25",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 31",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 44",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 37",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 6",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 6",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 5",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 5",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 23",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 30",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 21",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 42",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 33",
        );
        map.insert("CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT", ">PrefixPE/2 | >PCR_Primer2 | >Illumina Paired End PCR Primer 2 | >PrefixPE/2 | >PCR_Primer2");
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 22",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 1",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 1",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 17",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 9",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 9",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 47",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 28",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 45",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 7",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 7",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 34",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 3",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 3",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 18",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 24",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 26",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 43",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 14",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 16",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 20",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 11",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 11",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 39",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 41",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 12",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 12",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 29",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 8",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 8",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 46",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 40",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 15",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 32",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 48",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC",
            ">Illumina PCR Primer Index 4",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 4",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 36",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 13",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
            ">RNA PCR Primer, Index 19",
        );
        map.insert(
            "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
            ">Illumina Single End Adapter 2 | >Illumina Single End PCR Primer 2",
        );
    }

    /// Get the comprehensive known adapters database
    /// Returns a `HashMap` mapping adapter sequences to their descriptions
    #[allow(clippy::too_many_lines)]
    pub fn get_known_adapters() -> &'static HashMap<&'static str, &'static str> {
        static KNOWN_ADAPTERS: OnceLock<HashMap<&'static str, &'static str>> = OnceLock::new();
        KNOWN_ADAPTERS.get_or_init(|| {
            let mut map = HashMap::new();
            add_truseq_adapters(&mut map);
            add_rna_pcr_index_adapters(&mut map);
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", ">Illumina TruSeq Adapter Read 2");
            map.insert("GATCGTCGGACTGTAGAACTCTGAACGTGTAGA", ">Illumina Small RNA Adapter Read 2");
            map.insert("AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA", ">Illumina DpnII expression PCR Primer 2 | >Illumina NlaIII expression PCR Primer 2 | >Illumina Small RNA PCR Primer 2 | >Illumina DpnII Gex PCR Primer 2 | >Illumina NlaIII Gex PCR Primer 2");
            map.insert("AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA", ">Illumina RNA PCR Primer");
            map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1 | >Illumina Single End PCR Primer 1 | >Illumina Paried End PCR Primer 1 | >Illumina Multiplexing PCR Primer 1.01 | >TruSeq Universal Adapter | >TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1");
            map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">pcr_dimer");
            map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT", ">PCR_Primers");
            map.insert("ACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">Illumina Single End Sequencing Primer | >Illumina Paired End Adapter 1 | >Illumina Paried End Sequencing Primer 1 | >Illumina Multiplexing Adapter 2 | >Illumina Multiplexing Read1 Sequencing Primer");
            map.insert("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", ">PE2_rc | >TruSeq3_IndexedAdapter | >PE2_rc | >TruSeq3_IndexedAdapter");
            map.insert("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG", ">Reverse_adapter");
            map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", ">TruSeq2_PE_r");
            map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">PCR_Primer2_rc");
            map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA", ">PhiX_read1_adapter");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA", ">PE1_rc | >TruSeq3_UniversalAdapter | >PE1_rc | >TruSeq3_UniversalAdapter");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">PCR_Primer1_rc");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA", ">PhiX_read2_adapter");
            map.insert("AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", ">TruSeq2_SE");
            map.insert("CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 35");
            map.insert("CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC", ">Illumina PCR Primer Index 10");
            map.insert("CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 10");
            map.insert("CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC", ">Illumina PCR Primer Index 2");
            map.insert("CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 2");
            map.insert("CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 38");
            map.insert("CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 27");
            map.insert("CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 25");
            map.insert("CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 31");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 44");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 37");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC", ">Illumina PCR Primer Index 6");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 6");
            map.insert("CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC", ">Illumina PCR Primer Index 5");
            map.insert("CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 5");
            map.insert("CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 23");
            map.insert("CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 30");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 21");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 42");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 33");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT", ">PrefixPE/2 | >PCR_Primer2 | >Illumina Paired End PCR Primer 2 | >PrefixPE/2 | >PCR_Primer2");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 22");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC", ">Illumina PCR Primer Index 1");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 1");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 17");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC", ">Illumina PCR Primer Index 9");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 9");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 47");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 28");
            map.insert("CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 45");
            map.insert("CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC", ">Illumina PCR Primer Index 7");
            map.insert("CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 7");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 34");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC", ">Illumina PCR Primer Index 3");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 3");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 18");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 24");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 26");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 43");
            map.insert("CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 14");
            map.insert("CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 16");
            map.insert("CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 20");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC", ">Illumina PCR Primer Index 11");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 11");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 39");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 41");
            map.insert("CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC", ">Illumina PCR Primer Index 12");
            map.insert("CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 12");
            map.insert("CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 29");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC", ">Illumina PCR Primer Index 8");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 8");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 46");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 40");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 15");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 32");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 48");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC", ">Illumina PCR Primer Index 4");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 4");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 36");
            map.insert("CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 13");
            map.insert("CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 19");
            map.insert("CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT", ">Illumina Single End Adapter 2 | >Illumina Single End PCR Primer 2");
            map.insert("CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT", ">ABI Solid3 Adapter B");
            map.insert("CCGACAGGTTCAGAGTTCTACAGTCCGACATG", ">Illumina NlaIII expression Sequencing Primer | >Illumina NlaIII Gex Sequencing Primer");
            map.insert("CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | >I7_Primer_Nextera_XT_Index_Kit_v2_N711 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | >I7_Primer_Nextera_XT_Index_Kit_v2_N711");
            map.insert("CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N716");
            map.insert("CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N724");
            map.insert("CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | >I7_Primer_Nextera_XT_Index_Kit_v2_N703 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | >I7_Primer_Nextera_XT_Index_Kit_v2_N703");
            map.insert("CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N715");
            map.insert("CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N722");
            map.insert("CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N708");
            map.insert("CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N726");
            map.insert("CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | >I7_Primer_Nextera_XT_Index_Kit_v2_N710 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | >I7_Primer_Nextera_XT_Index_Kit_v2_N710");
            map.insert("CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N727");
            map.insert("CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N720");
            map.insert("CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | >I7_Primer_Nextera_XT_Index_Kit_v2_N702 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | >I7_Primer_Nextera_XT_Index_Kit_v2_N702");
            map.insert("CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | >I7_Primer_Nextera_XT_Index_Kit_v2_N707 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | >I7_Primer_Nextera_XT_Index_Kit_v2_N707");
            map.insert("CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N719");
            map.insert("CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N709");
            map.insert("CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N714");
            map.insert("CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | >I7_Primer_Nextera_XT_Index_Kit_v2_N705 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | >I7_Primer_Nextera_XT_Index_Kit_v2_N705");
            map.insert("CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N718");
            map.insert("CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | >I7_Primer_Nextera_XT_Index_Kit_v2_N712 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | >I7_Primer_Nextera_XT_Index_Kit_v2_N712");
            map.insert("CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | >I7_Primer_Nextera_XT_Index_Kit_v2_N701 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | >I7_Primer_Nextera_XT_Index_Kit_v2_N701");
            map.insert("CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N721");
            map.insert("CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N723");
            map.insert("CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | >I7_Primer_Nextera_XT_Index_Kit_v2_N706 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | >I7_Primer_Nextera_XT_Index_Kit_v2_N706");
            map.insert("CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | >I7_Primer_Nextera_XT_Index_Kit_v2_N704 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | >I7_Primer_Nextera_XT_Index_Kit_v2_N704");
            map.insert("CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N729");
            map.insert("CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N728");
            map.insert("CGACAGGTTCAGAGTTCTACAGTCCGACGATC", ">Illumina DpnII expression Sequencing Primer | >Illumina Small RNA Sequencing Primer | >Illumina DpnII Gex Sequencing Primer");
            map.insert("CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT", ">Illumina Paired End Sequencing Primer 2");
            map.insert("CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT", ">Clontech Universal Primer Mix Long");
            map.insert("CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Adapter_Nextera_No_Barcode");
            map.insert("CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Adapter_Nextera");
            map.insert("CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG", ">ABI Solid3 Adapter A");
            map.insert("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", ">I7_Nextera_Transposase_1 | >Trans2_rc | >I7_Nextera_Transposase_1 | >Trans2_rc");
            map.insert("CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC", ">I7_Nextera_Transposase_2");
            map.insert("CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC", ">I5_Nextera_Transposase_2");
            map.insert("CTGTCTCTTATACACATCTGACGCTGCCGACGA", ">I5_Nextera_Transposase_1 | >Trans1_rc | >I5_Nextera_Transposase_1 | >Trans1_rc");
            map.insert("GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S516");
            map.insert("GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | >I5_Primer_Nextera_XT_Index_Kit_v2_S503 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | >I5_Primer_Nextera_XT_Index_Kit_v2_S503");
            map.insert("GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S515");
            map.insert("GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | >I5_Primer_Nextera_XT_Index_Kit_v2_S508 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | >I5_Primer_Nextera_XT_Index_Kit_v2_S508");
            map.insert("GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | >I5_Primer_Nextera_XT_Index_Kit_v2_S502 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | >I5_Primer_Nextera_XT_Index_Kit_v2_S502");
            map.insert("GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S520");
            map.insert("GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S510");
            map.insert("GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S511");
            map.insert("GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S513");
            map.insert("GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | >I5_Primer_Nextera_XT_Index_Kit_v2_S505 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | >I5_Primer_Nextera_XT_Index_Kit_v2_S505");
            map.insert("GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S518");
            map.insert("GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]501");
            map.insert("GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S521");
            map.insert("GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | >I5_Primer_Nextera_XT_Index_Kit_v2_S507 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | >I5_Primer_Nextera_XT_Index_Kit_v2_S507");
            map.insert("GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | >I5_Primer_Nextera_XT_Index_Kit_v2_S506 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | >I5_Primer_Nextera_XT_Index_Kit_v2_S506");
            map.insert("GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S522");
            map.insert("GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]504");
            map.insert("GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | >I5_Primer_Nextera_XT_Index_Kit_v2_S517 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | >I5_Primer_Nextera_XT_Index_Kit_v2_S517");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", ">Nextera_LMP_Read1_External_Adapter | >Illumina Multiplexing Index Sequencing Primer");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_5 | >TruSeq Adapter, Index 5");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_25");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 25");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_8 | >TruSeq Adapter, Index 8");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_13");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 13");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_14");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 14");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_1_6 | >TruSeq Adapter, Index 1");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_15");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 15");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_27");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 27");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_7 | >TruSeq Adapter, Index 7");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 23");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_16");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 16");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_2 | >TruSeq Adapter, Index 2");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_22");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 22");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_12 | >TruSeq Adapter, Index 12");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_23");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_9 | >TruSeq Adapter, Index 9");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_6 | >TruSeq Adapter, Index 6");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_11 | >TruSeq Adapter, Index 11");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_18_7");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 18");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_19");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 19");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_20");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 20");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_21");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 21");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_10 | >TruSeq Adapter, Index 10");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_4 | >TruSeq Adapter, Index 4");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_3 | >TruSeq Adapter, Index 3");
            map.insert("GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", ">Illumina Paired End Adapter 2");
            map.insert("GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", ">Nextera_LMP_Read2_External_Adapter");
            map.insert("GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", ">Illumina Single End Adapter 1");
            map.insert("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", ">Trans2");
            map.insert("GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", ">PrefixPE/2 | >PE2 | >Illumina Multiplexing PCR Primer 2.01 | >Illumina Multiplexing Read2 Sequencing Primer | >PrefixPE/2 | >PE2");
            map.insert("TACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">PrefixPE/1 | >PE1 | >PrefixPE/1 | >PE1");
            map.insert("TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT", ">RNA_PCR_Primer_(RP1)_part_#_15013198");
            map.insert("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", ">Trans1");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_5_(RPI5)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_25_(RPI25)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_8_(RPI8)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_13_(RPI13)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_14_(RPI14)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_1_(RPI1)_2,9");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_26_(RPI26)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_15_(RPI15)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_27_(RPI27)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_28_(RPI28)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_29_(RPI29)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_30_(RPI30)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_31_(RPI31)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_32_(RPI32)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_7_(RPI7)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_33_(RPI33)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_34_(RPI34)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_35_(RPI35)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_36_(RPI36)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_16_(RPI16)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_2_(RPI2)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_37_(RPI37)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_22_(RPI22)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_38_(RPI38)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_39_(RPI39)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_40_(RPI40)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_12_(RPI12)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_41_(RPI41)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_23_(RPI23)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_9_(RPI9)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_6_(RPI6)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_11_(RPI11)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_24_(RPI24)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_17_(RPI17)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_18_(RPI18)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_19_(RPI19)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_20_(RPI20)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_21_(RPI21)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_42_(RPI42)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_43_(RPI43)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_10_(RPI10)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_44_(RPI44)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_45_(RPI45)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_46_(RPI46)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_47_(RPI47)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_48_(RPI48)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_4_(RPI4)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_3_(RPI3)");
            map.insert("TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC", ">FlowCell1");
            map.insert("TTTTTTTTTTCAAGCAGAAGACGGCATACGA", ">FlowCell2");
            map.insert("AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA", ">MGI/BGI adapter (forward)");
            map.insert("AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG", ">MGI/BGI adapter (reverse)");
            map.insert("AACTGTAGGCACCATCAAT", ">QIASeq miRNA adapter");
            map
        })
    }
}

/// Type of adapter match found
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatchType {
    /// Exact match with allowed mismatches
    Exact,
    /// Match with single insertion in read sequence
    Insertion,
    /// Match with single deletion from adapter sequence
    Deletion,
}

/// Result of adapter detection
#[derive(Debug, Clone)]
pub struct AdapterMatch {
    /// Position where adapter starts (0 if originally negative)
    pub position: usize,
    /// Number of mismatches
    pub mismatches: usize,
    /// Type of match found
    pub match_type: MatchType,
    /// Overlap length between adapter and read (for correct base counting)
    pub overlap_len: usize,
    /// Whether this was an A-tailing match (original position was negative)
    pub is_atailing: bool,
}

/// Find adapter sequence in read allowing mismatches
///
/// Uses a three-stage matching approach (exact, insertion, deletion).
/// Returns the first accepted match in upstream's stage and position order.
pub fn find_adapter(
    seq: &[u8],
    adapter: &[u8],
    min_overlap: usize,
    _max_mismatches: usize,
) -> Option<AdapterMatch> {
    if adapter.len() < min_overlap || seq.len() <= min_overlap {
        return None;
    }
    let start = if adapter.len() >= 16 {
        -4
    } else if adapter.len() >= 12 {
        -3
    } else if adapter.len() >= 8 {
        -2
    } else {
        0
    };
    for pos in start..(seq.len() - min_overlap) as isize {
        let len = ((seq.len() as isize - pos) as usize).min(adapter.len());
        let skip = (-pos).max(0) as usize;
        let begin = (pos + skip as isize) as usize;
        let (diff, _) = crate::simd::count_mismatches(
            &seq[begin..begin + len - skip],
            &adapter[skip..len],
            len / 8,
            len - skip,
        );
        if diff <= len / 8 {
            return Some(AdapterMatch {
                position: pos.max(0) as usize,
                mismatches: diff,
                match_type: MatchType::Exact,
                overlap_len: if pos < 0 {
                    (adapter.len() as isize + pos) as usize
                } else {
                    len
                },
                is_atailing: pos < 0,
            });
        }
    }
    for insertion in [true, false] {
        let bound = seq
            .len()
            .saturating_sub(min_overlap + usize::from(insertion));
        let mut previous_len = None;
        for pos in 0..bound {
            let len = if insertion {
                (seq.len() - pos - 1).min(adapter.len())
            } else {
                (seq.len() - pos).min(adapter.len().saturating_sub(1))
            };
            if len < 8 {
                continue;
            }
            // fastp 1.4.0 passes the original read/adapter pointers to its gap
            // matcher, without adding `pos`. Until the overlap length changes,
            // these calls would repeat the same failed comparison.
            if previous_len == Some(len) {
                continue;
            }
            previous_len = Some(len);
            let diff = if insertion {
                one_gap_match(seq, adapter, len, len / 8 - 1, true)
            } else {
                one_gap_match(adapter, seq, len, len / 8 - 1, true)
            };
            if let Some(diff) = diff {
                return Some(AdapterMatch {
                    position: pos,
                    mismatches: diff,
                    match_type: if insertion {
                        MatchType::Insertion
                    } else {
                        MatchType::Deletion
                    },
                    overlap_len: len,
                    is_atailing: false,
                });
            }
        }
    }
    None
}

/// Hamming distance with one internal insertion in `ins`, matching fastp's
/// Matcher. The gap itself is not counted as a mismatching base.
pub(crate) fn one_gap_diff(ins: &[u8], normal: &[u8], len: usize, limit: usize) -> Option<usize> {
    one_gap_match(ins, normal, len, limit, false)
}

fn one_gap_match(
    ins: &[u8],
    normal: &[u8],
    len: usize,
    limit: usize,
    first: bool,
) -> Option<usize> {
    if len < 2 || ins.len() <= len || normal.len() < len {
        return None;
    }
    // Keep fastp's prefix/tail rejection rule, including for the minimum-diff
    // overlap mode. Reject impossible first positions before counting a suffix.
    let tail = usize::from(ins[len] != normal[len - 1]);
    if usize::from(ins[0] != normal[0]) + tail > limit {
        return None;
    }
    let (mut suffix, _) =
        crate::simd::count_mismatches(&ins[1..=len], &normal[..len], usize::MAX, 0);
    let mut prefix = 0;
    let mut best = usize::MAX;
    for i in 1..len {
        prefix += usize::from(ins[i - 1] != normal[i - 1]);
        suffix -= usize::from(ins[i] != normal[i - 1]);
        if prefix + tail > limit {
            return None;
        }
        let diff = prefix + suffix;
        if first && diff <= limit {
            return Some(diff);
        }
        best = best.min(diff);
    }
    (best <= limit).then_some(best)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Keep the original scalar matcher as an independent test reference.
    fn reference_find_adapter(
        seq: &[u8],
        adapter: &[u8],
        min_overlap: usize,
        _max_mismatches: usize,
    ) -> Option<AdapterMatch> {
        if adapter.len() < min_overlap || seq.len() <= min_overlap {
            return None;
        }
        let start = if adapter.len() >= 16 {
            -4
        } else if adapter.len() >= 12 {
            -3
        } else if adapter.len() >= 8 {
            -2
        } else {
            0
        };
        for pos in start..(seq.len() - min_overlap) as isize {
            let len = ((seq.len() as isize - pos) as usize).min(adapter.len());
            let skip = (-pos).max(0) as usize;
            let begin = (pos + skip as isize) as usize;
            let diff = seq[begin..begin + len - skip]
                .iter()
                .zip(&adapter[skip..len])
                .filter(|(a, b)| a != b)
                .count();
            if diff <= len / 8 {
                return Some(AdapterMatch {
                    position: pos.max(0) as usize,
                    mismatches: diff,
                    match_type: MatchType::Exact,
                    overlap_len: if pos < 0 {
                        (adapter.len() as isize + pos) as usize
                    } else {
                        len
                    },
                    is_atailing: pos < 0,
                });
            }
        }
        for insertion in [true, false] {
            let bound = seq
                .len()
                .saturating_sub(min_overlap + usize::from(insertion));
            for pos in 0..bound {
                let len = if insertion {
                    (seq.len() - pos - 1).min(adapter.len())
                } else {
                    (seq.len() - pos).min(adapter.len().saturating_sub(1))
                };
                if len < 8 {
                    continue;
                }
                let diff = if insertion {
                    reference_one_gap_match(seq, adapter, len, len / 8 - 1, true)
                } else {
                    reference_one_gap_match(adapter, seq, len, len / 8 - 1, true)
                };
                if let Some(diff) = diff {
                    return Some(AdapterMatch {
                        position: pos,
                        mismatches: diff,
                        match_type: if insertion {
                            MatchType::Insertion
                        } else {
                            MatchType::Deletion
                        },
                        overlap_len: len,
                        is_atailing: false,
                    });
                }
            }
        }
        None
    }

    fn reference_one_gap_match(
        ins: &[u8],
        normal: &[u8],
        len: usize,
        limit: usize,
        first: bool,
    ) -> Option<usize> {
        if len < 2 || ins.len() <= len || normal.len() < len {
            return None;
        }
        let mut suffix = vec![0usize; len + 1];
        for i in (0..len).rev() {
            suffix[i] = suffix[i + 1] + usize::from(ins[i + 1] != normal[i]);
        }
        let mut prefix = 0;
        let mut best = usize::MAX;
        for i in 1..len {
            prefix += usize::from(ins[i - 1] != normal[i - 1]);
            if prefix + suffix[len - 1] > limit {
                return None;
            }
            let diff = prefix + suffix[i];
            if first && diff <= limit {
                return Some(diff);
            }
            best = best.min(diff);
        }
        (best <= limit).then_some(best)
    }

    fn match_fields(value: Option<AdapterMatch>) -> Option<(usize, usize, MatchType, usize, bool)> {
        value.map(|m| {
            (
                m.position,
                m.mismatches,
                m.match_type,
                m.overlap_len,
                m.is_atailing,
            )
        })
    }

    #[test]
    fn gap_match_matches_scalar_oracle() {
        // Exhaust every binary sequence pair, gap position, mismatch budget and
        // selection mode at short lengths, including the late rejection rule.
        for len in 2..=6 {
            for bits in 0..(1usize << (2 * len + 1)) {
                let ins: Vec<u8> = (0..=len).map(|i| b'A' + ((bits >> i) & 1) as u8).collect();
                let normal: Vec<u8> = (0..len)
                    .map(|i| b'A' + ((bits >> (len + 1 + i)) & 1) as u8)
                    .collect();
                for limit in 0..=len {
                    for first in [true, false] {
                        assert_eq!(
                            one_gap_match(&ins, &normal, len, limit, first),
                            reference_one_gap_match(&ins, &normal, len, limit, first),
                            "len={len} bits={bits} limit={limit} first={first}"
                        );
                    }
                }
            }
        }
        // Longer inputs exercise the SIMD paths and realistic adapter lengths.
        let mut state = 0x3b9d_81a2_f605_c74eu64;
        let mut random = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as usize
        };
        for _ in 0..4000 {
            let len = 8 + random() % 249;
            let normal: Vec<u8> = (0..len).map(|_| b"ACGTN"[random() % 5]).collect();
            let mut ins = normal.clone();
            ins.insert(1 + random() % (len - 1), b"ACGTN"[random() % 5]);
            for _ in 0..random() % (len / 4 + 1) {
                ins[random() % (len + 1)] = b"ACGTN"[random() % 5];
            }
            for first in [true, false] {
                let limit = random() % (len / 4 + 1);
                assert_eq!(
                    one_gap_match(&ins, &normal, len, limit, first),
                    reference_one_gap_match(&ins, &normal, len, limit, first)
                );
            }
        }
    }

    #[test]
    fn adapter_search_matches_scalar_oracle() {
        let mut state = 0xf3c6_705a_92b1_e84du64;
        let mut random = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as usize
        };
        for case in 0..6000 {
            let len = 8 + random() % 193;
            let adapter: Vec<u8> = (0..len).map(|_| b"ACGTN"[random() % 5]).collect();
            let mut seq: Vec<u8> = (0..random() % 200)
                .map(|_| b"ACGTN"[random() % 5])
                .collect();
            match case % 5 {
                0 => seq.extend_from_slice(&adapter),
                1 => {
                    seq = adapter.clone();
                    seq.insert(1 + random() % (len - 1), b'N');
                }
                2 => {
                    seq = adapter.clone();
                    seq.remove(1 + random() % (len - 1));
                }
                3 => seq = adapter[random() % 5..].to_vec(),
                _ => (),
            }
            for _ in 0..random() % 8 {
                if !seq.is_empty() {
                    let pos = random() % seq.len();
                    seq[pos] = b"ACGTN"[random() % 5];
                }
            }
            let min_overlap = 1 + random() % 15;
            assert_eq!(
                match_fields(find_adapter(&seq, &adapter, min_overlap, 0)),
                match_fields(reference_find_adapter(&seq, &adapter, min_overlap, 0)),
                "case={case} min_overlap={min_overlap}"
            );
        }
    }

    #[test]
    fn test_exact_adapter_match() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGC";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 12); // Adapter starts at position 12
        assert_eq!(m.mismatches, 0);
    }

    #[test]
    fn test_adapter_with_mismatch() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGX"; // X instead of C
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 12);
        assert_eq!(m.mismatches, 1);
    }

    #[test]
    fn test_partial_adapter_at_end() {
        // Only partial adapter at the end of read
        let seq = b"ACGTACGTACGTACGTAGATCGGAA"; // Only first 9 bases of adapter
        let adapter = b"AGATCGGAAGAGC";

        // Should not match with min_overlap=10
        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());

        // Should match with min_overlap=8
        let result = find_adapter(seq, adapter, 8, 2);
        assert!(result.is_some());
    }

    #[test]
    fn test_no_adapter() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());
    }

    #[test]
    fn test_adapter_too_many_mismatches() {
        // 3 mismatches (XXX), but max_mismatches=2
        let seq = b"ACGTACGTACGTAGATCGGAAXXXC";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());
    }

    #[test]
    fn test_matching_preserves_fastp_case_sensitivity() {
        let seq = b"ACGTACGTACGTagatcggaagagc"; // lowercase adapter
        let adapter = b"AGATCGGAAGAGC"; // uppercase

        let result = find_adapter(seq, adapter, 10, 2);
        // fastp compares original sequence bytes against the adapter.
        assert!(result.is_none());
    }

    #[test]
    fn test_insertion_in_read() {
        // Sequence has extra base (X) that needs to be skipped
        let seq = b"AGATCGGXAAGAGCTTTTTTTT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 5, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 0);
        assert_eq!(m.match_type, MatchType::Insertion);
    }

    #[test]
    fn test_deletion_from_adapter() {
        // seq.2033 case: AGATCGGAGCTC vs AGATCGGAAGAGC (missing 'A' at position 8)
        // This has too many mismatches (diff_limit = cmplen/8 - 1 = 0 for cmplen=12)
        // Fastp does NOT match this at position 0 - verified by testing
        let seq = b"AGATCGGAGCTCACGGATCAGGTGAAT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 5, 2);

        // Fastp finds a match at position 18, not position 0
        // The beginning doesn't match due to strict deletion matching criteria
        if let Some(m) = result {
            assert_ne!(
                m.position, 0,
                "Should not match at position 0 with deletion logic"
            );
        }
    }

    #[test]
    fn test_exact_match_preferred_over_indel() {
        // When both exact and indel match exist, prefer exact
        let seq = b"AGATCGGAAGAGCTTTTTTTT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 5, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.match_type, MatchType::Exact);
    }

    // #[test]
    // fn test_earlier_position_preferred() {
    //     // Earlier position should be preferred over later position
    //     let seq = b"AAAAAAGATCGGAAGAGCXXXXAGATCGGAAGAGC";
    //     let adapter = b"AGATCGGAAGAGC";

    //     let result = find_adapter(seq, adapter, 5, 2);
    //     assert!(result.is_some());

    //     let m = result.unwrap();
    //     assert_eq!(m.position, 6); // First occurrence
    // }

    #[test]
    fn test_adapter_detection_no_overlap() {
        // Reads without overlap should not detect adapters
        let r1_seqs = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".as_slice(),
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".as_slice(),
        ];
        let r2_seqs = vec![
            b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".as_slice(),
            b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".as_slice(),
        ];

        let result =
            detect_adapters_from_pe_reads(r1_seqs.into_iter(), r2_seqs.into_iter(), 100, 0.01);

        // Should not detect adapters when reads don't overlap
        assert!(result.adapter_r1.is_none() || result.confidence < 0.01);
    }

    #[test]
    fn test_adapter_detection_with_overlap() {
        use crate::overlap::reverse_complement;

        // Create reads that overlap with adapter contamination
        // R1: 30bp insert + 10bp adapter
        // R2: RC of same 30bp insert (no adapter for simplicity)
        let insert = b"ACGTACGTACGTACGTACGTACGTACGTAC"; // 30bp
        let adapter = b"AGATCGGAAG"; // 10bp adapter

        let mut r1 = insert.to_vec();
        r1.extend_from_slice(adapter);

        let r2 = reverse_complement(insert);

        let r1_seqs = vec![r1.as_slice(); 100]; // 100 identical reads
        let r2_seqs = vec![r2.as_slice(); 100];

        let result =
            detect_adapters_from_pe_reads(r1_seqs.into_iter(), r2_seqs.into_iter(), 100, 0.01);

        // Should detect adapter
        assert!(result.adapter_r1.is_some());
        assert!(result.confidence > 0.5); // High confidence
        assert_eq!(result.reads_analyzed, 100);
    }

    #[test]
    fn test_consensus_adapter_selection() {
        let mut adapter_counts = BTreeMap::new();
        adapter_counts.insert(b"AGATCGG".to_vec(), 5);
        adapter_counts.insert(b"AGATCGGAAGAGC".to_vec(), 100); // Most frequent
        adapter_counts.insert(b"TTTTTTT".to_vec(), 3);

        let result = find_consensus_adapter(&adapter_counts, 10);

        assert!(result.is_some());
        let adapter = result.unwrap();
        // Should pick the most frequent one
        assert!(adapter.len() >= 7); // At least AGATCGG
    }

    #[test]
    fn test_consensus_adapter_threshold() {
        let mut adapter_counts = BTreeMap::new();
        adapter_counts.insert(b"AGATCGG".to_vec(), 5); // Below threshold of 10

        let result = find_consensus_adapter(&adapter_counts, 10);

        // Should return None when no adapter meets threshold
        assert!(result.is_none());
    }
}
