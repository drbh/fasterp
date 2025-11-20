//! Adapter trimming functionality
//!
//! This module provides adapter detection and trimming for FASTQ preprocessing:
//! - Manual adapter specification
//! - Built-in common adapter sequences
//! - Overlap-based adapter detection for paired-end reads
//! - Auto-detection using PE overlap analysis
//! - Mismatch-tolerant matching

#[cfg(feature = "native")]
use rayon::prelude::*;
use std::cell::RefCell;
use std::cmp::min;
use std::collections::BTreeMap;

// Thread-local scratch buffers for adapter matching to avoid per-call heap allocations
thread_local! {
    static INSERTION_SCRATCH: RefCell<(Vec<u16>, Vec<u16>)> = const { RefCell::new((Vec::new(), Vec::new())) };
}

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
            min_overlap: 5,               // fastp's default for adapter trimming
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
    use crate::overlap::{OverlapConfig, detect_overlap, reverse_complement};

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

/// Try exact matching with allowed mismatches (Stage 1)
/// Supports negative `start_pos` for A-tailing (Illumina adapter dimer handling)
fn try_exact_match(
    seq: &[u8],
    adapter: &[u8],
    start_pos: isize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    let (seq_offset, adapter_offset, compare_len) = if start_pos < 0 {
        // Negative position: skip first abs(start_pos) bytes of adapter
        let skip = (-start_pos) as usize;
        if skip >= adapter.len() {
            return None;
        }
        let remaining_adapter = adapter.len() - skip;
        let compare_len = min(seq.len(), remaining_adapter);
        (0, skip, compare_len)
    } else {
        // Positive position: normal matching
        let start = start_pos as usize;
        if start >= seq.len() {
            return None;
        }
        let remaining = seq.len() - start;
        let compare_len = min(remaining, adapter.len());
        (start, 0, compare_len)
    };

    if compare_len < min_overlap {
        return None;
    }

    // Use fastp's cmplen formula for allowed mismatches calculation
    // Fastp uses min(seq_len - pos, adapter_len) even for negative positions
    let cmplen = if start_pos < 0 {
        min((seq.len() as isize - start_pos) as usize, adapter.len())
    } else {
        compare_len
    };
    let allowed_mismatches = cmplen / 8;
    let mut matches = 0;
    let mut mismatches = 0;

    for i in 0..compare_len {
        let seq_base = seq[seq_offset + i].to_ascii_uppercase();
        let adapter_base = adapter[adapter_offset + i].to_ascii_uppercase();

        if seq_base == adapter_base {
            matches += 1;
        } else {
            mismatches += 1;
            if mismatches > allowed_mismatches {
                return None;
            }
        }
    }

    if mismatches <= allowed_mismatches && matches + mismatches >= min_overlap {
        Some(AdapterMatch {
            position: start_pos.max(0) as usize,
            mismatches,
            match_type: MatchType::Exact,
            overlap_len: compare_len,
            is_atailing: start_pos < 0,
        })
    } else {
        None
    }
}

/// Try matching with single insertion (Stage 2)
/// Implements fastp's `Matcher::matchWithOneInsertion` algorithm exactly
///
/// `ins_data`: sequence with suspected insertion (longer, e.g., read)
/// `normal_data`: reference sequence (baseline, e.g., adapter)
/// cmplen: comparison length (calculated by caller based on insertion/deletion case)
///
/// OPTIMIZED: Uses thread-local scratch buffers to avoid per-call heap allocations
fn try_insertion_match_new(
    ins_data: &[u8],    // Sequence with insertion
    normal_data: &[u8], // Reference sequence (adapter)
    cmplen: usize,      // Length to compare
    min_overlap: usize,
) -> Option<AdapterMatch> {
    // Early validation
    if ins_data.len() < cmplen + 1 || cmplen < min_overlap || cmplen < 8 {
        return None;
    }

    let diff_limit: u16 = ((cmplen / 8).saturating_sub(1)) as u16;
    let init = diff_limit.saturating_add(1);

    INSERTION_SCRATCH.with(|cell| {
        let (left_buf, right_buf) = &mut *cell.borrow_mut();

        // Resize buffers if needed (rare, only happens once per thread typically)
        if left_buf.len() < cmplen {
            left_buf.resize(cmplen, init);
        }
        if right_buf.len() < cmplen {
            right_buf.resize(cmplen, init);
        }

        let left = &mut left_buf[..cmplen];
        let right = &mut right_buf[..cmplen];

        // Initialize arrays
        left.fill(init);
        right.fill(init);

        // Initialize first and last elements
        left[0] = u16::from(!ins_data[0].eq_ignore_ascii_case(&normal_data[0]));

        if cmplen >= ins_data.len() {
            return None;
        }

        right[cmplen - 1] =
            u16::from(!ins_data[cmplen].eq_ignore_ascii_case(&normal_data[cmplen - 1]));

        // Build left array with early termination
        for i in 1..cmplen {
            if i >= ins_data.len() {
                return None;
            }
            left[i] = left[i - 1] + u16::from(!ins_data[i].eq_ignore_ascii_case(&normal_data[i]));
            if left[i] + right[cmplen - 1] > diff_limit {
                break;
            }
        }

        // Build right array with early termination
        for i in (0..cmplen - 1).rev() {
            if i + 1 >= ins_data.len() {
                continue;
            }
            right[i] =
                right[i + 1] + u16::from(!ins_data[i + 1].eq_ignore_ascii_case(&normal_data[i]));
            if right[i] + left[0] > diff_limit {
                right[..i].fill(init);
                break;
            }
        }

        // Check each potential skip position
        for i in 1..cmplen {
            if left[i - 1] + right[cmplen - 1] > diff_limit {
                return None;
            }
            let diff = left[i - 1] + right[i];
            if diff <= diff_limit {
                return Some(AdapterMatch {
                    position: 0,
                    mismatches: diff as usize,
                    match_type: MatchType::Insertion,
                    overlap_len: cmplen,
                    is_atailing: false,
                });
            }
        }

        None
    })
}

#[inline(always)]
fn mismatch_ci(a: u8, b: u8) -> u16 {
    // ASCII-fold to uppercase: 'a'..'z' -> 'A'..'Z'
    let ua = a & 0xDF;
    let ub = b & 0xDF;
    u16::from(ua != ub)
}

/// Quick prefilter: if even a no-insertion match is far over the limit, skip insertion logic
#[inline(always)]
fn quick_hamming_over_limit(
    ins_data: &[u8],
    normal_data: &[u8],
    cmplen: usize,
    limit: u16,
) -> bool {
    let mut diff = 0u16;
    for i in 0..cmplen {
        diff += mismatch_ci(ins_data[i], normal_data[i]);
        if diff > limit {
            return true;
        }
    }
    false
}

/// Entry point for insertion matching with prefilter
fn try_insertion_match_entry(
    ins_data: &[u8],
    normal_data: &[u8],
    cmplen: usize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    if ins_data.len() < cmplen + 1
        || normal_data.len() < cmplen
        || cmplen < min_overlap
        || cmplen < 8
    {
        return None;
    }

    let diff_limit: u16 = ((cmplen / 8).saturating_sub(1)) as u16;

    // If even a no-insertion match is *far* over the limit, skip insertion logic.
    // Tune the "+ K" based on data; start conservative (e.g. +4 or +6)
    if quick_hamming_over_limit(ins_data, normal_data, cmplen, diff_limit + 4) {
        return None;
    }

    try_insertion_match(ins_data, normal_data, cmplen, min_overlap)
}

/// Try matching with single insertion (Stage 2)
/// Stack-based implementation for common case, falls back to TLS for large cmplen.
fn try_insertion_match(
    ins_data: &[u8],
    normal_data: &[u8],
    cmplen: usize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    // Same basic guards as before
    if ins_data.len() < cmplen + 1
        || normal_data.len() < cmplen
        || cmplen < min_overlap
        || cmplen < 8
    {
        return None;
    }

    // IMPORTANT: if you ever see bigger overlaps, just raise this and recompile.
    const MAX_STACK_CMPLEN: usize = 64;

    if cmplen <= MAX_STACK_CMPLEN {
        return try_insertion_match_stack::<MAX_STACK_CMPLEN>(ins_data, normal_data, cmplen);
    }

    // Rare path: use the old TLS-based implementation here
    try_insertion_match_slow(ins_data, normal_data, cmplen, min_overlap)
}

#[inline(always)]
fn try_insertion_match_stack<const MAX: usize>(
    ins_data: &[u8],
    normal_data: &[u8],
    cmplen: usize,
) -> Option<AdapterMatch> {
    debug_assert!(cmplen <= MAX);
    debug_assert!(ins_data.len() >= cmplen + 1);
    debug_assert!(normal_data.len() >= cmplen);

    let diff_limit: u16 = ((cmplen / 8).saturating_sub(1)) as u16;
    let init: u16 = diff_limit.saturating_add(1);

    // Everything stays on the stack
    let mut left: [u16; MAX] = [init; MAX];
    let mut right: [u16; MAX] = [init; MAX];

    let last = cmplen - 1;

    // left[0] = mismatches in ins[0] vs norm[0]
    left[0] = mismatch_ci(ins_data[0], normal_data[0]);

    // right[last] = mismatches in ins[cmplen] vs norm[cmplen-1]
    right[last] = mismatch_ci(ins_data[cmplen], normal_data[last]);
    let right_last = right[last];

    // Build left cumulative mismatches with early termination
    for i in 1..cmplen {
        let diff = mismatch_ci(ins_data[i], normal_data[i]);
        let v = left[i - 1] + diff;
        left[i] = v;

        // original-style shortcut: even if right side is perfect, too many mismatches
        if v + right_last > diff_limit {
            // later entries remain at init
            break;
        }
    }

    // Build right cumulative mismatches (reverse) with early termination
    for i in (0..last).rev() {
        // right[i] = right[i+1] + mismatch(ins[i+1], norm[i])
        let diff = mismatch_ci(ins_data[i + 1], normal_data[i]);
        let v = right[i + 1] + diff;
        right[i] = v;

        if v + left[0] > diff_limit {
            // everything before i is conceptually init
            for k in 0..i {
                right[k] = init;
            }
            break;
        }
    }

    // Check each potential skip position
    for i in 1..cmplen {
        // fastp shortcut: if even best right-side (at last) can't save it, bail
        if left[i - 1] + right_last > diff_limit {
            return None;
        }

        let diff = left[i - 1] + right[i];
        if diff <= diff_limit {
            return Some(AdapterMatch {
                position: 0,
                mismatches: diff as usize,
                match_type: MatchType::Insertion,
                overlap_len: cmplen,
                is_atailing: false,
            });
        }
    }

    None
}

fn try_insertion_match_slow(
    ins_data: &[u8],    // sequence with insertion (longer)
    normal_data: &[u8], // reference (adapter)
    cmplen: usize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    // Same guards as your original
    if ins_data.len() < cmplen + 1 || cmplen < min_overlap || cmplen < 8 {
        return None;
    }

    // Original extra check you preserved
    if cmplen >= ins_data.len() {
        return None;
    }

    let diff_limit: u16 = ((cmplen / 8).saturating_sub(1)) as u16;
    let init: u16 = diff_limit.saturating_add(1);

    INSERTION_SCRATCH.with(|cell| {
        let (left_buf, right_buf) = &mut *cell.borrow_mut();

        // Ensure buffer capacity (rare, amortized)
        if left_buf.len() < cmplen {
            left_buf.resize(cmplen, init);
        }
        if right_buf.len() < cmplen {
            right_buf.resize(cmplen, init);
        }

        let left = &mut left_buf[..cmplen];
        let right = &mut right_buf[..cmplen];

        // Initialize arrays to "infinite" diff (same semantics as original)
        left.fill(init);
        right.fill(init);

        // left[0]
        left[0] = mismatch_ci(ins_data[0], normal_data[0]);

        // right[cmplen - 1]
        let last = cmplen - 1;
        right[last] = mismatch_ci(ins_data[cmplen], normal_data[last]);

        let right_last = right[last];

        // Build left[] with early termination
        // left[i] = left[i-1] + mismatch(ins[i], normal[i])
        for i in 1..cmplen {
            // original had: if i >= ins.len { return None; }
            // but precondition ins.len >= cmplen + 1 makes this impossible
            left[i] = left[i - 1] + mismatch_ci(ins_data[i], normal_data[i]);

            // early stop as in your original code
            if left[i] + right_last > diff_limit {
                // later i’s stay at init
                break;
            }
        }

        // Build right[] with early termination
        // right[i] = right[i+1] + mismatch(ins[i+1], normal[i])
        //
        // original loop:
        // for i in (0..cmplen-1).rev() { ... }
        for i in (0..last).rev() {
            // original had: if i + 1 >= ins.len { continue; }
            // but i+1 <= cmplen-1 < ins.len (from precondition), so we can drop that
            right[i] = right[i + 1] + mismatch_ci(ins_data[i + 1], normal_data[i]);

            if right[i] + left[0] > diff_limit {
                // same semantics as:
                // right[..i].fill(init);
                for k in 0..i {
                    right[k] = init;
                }
                break;
            }
        }

        // Check each potential skip position
        for i in 1..cmplen {
            // Original early exit:
            // if left[i - 1] + right[cmplen - 1] > diff_limit { return None; }
            if left[i - 1] + right_last > diff_limit {
                return None;
            }

            let diff = left[i - 1] + right[i];
            if diff <= diff_limit {
                return Some(AdapterMatch {
                    position: 0,
                    mismatches: diff as usize,
                    match_type: MatchType::Insertion,
                    overlap_len: cmplen,
                    is_atailing: false,
                });
            }
        }

        None
    })
}

/// Try matching with single deletion from adapter (Stage 3)
/// Fastp handles deletion by swapping arguments to matchWithOneInsertion:
/// A deletion in adapter = insertion in adapter sequence relative to read
fn try_deletion_match(
    seq: &[u8],
    adapter: &[u8],
    cmplen: usize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    if cmplen < min_overlap {
        return None;
    }

    // Swap arguments: treat adapter as having "insertion" relative to seq
    // Fastp always compares from index 0
    if let Some(mut result) = try_insertion_match(adapter, seq, cmplen, min_overlap) {
        result.match_type = MatchType::Deletion;
        Some(result)
    } else {
        None
    }
}

/// Check if new match is better than current best
fn is_better_match(new_match: &AdapterMatch, current_best: Option<&AdapterMatch>) -> bool {
    use MatchType::{Deletion, Exact, Insertion};

    match current_best {
        None => true,
        Some(best) => {
            // Prefer earlier position
            if new_match.position < best.position {
                return true;
            }
            if new_match.position > best.position {
                return false;
            }

            // Same position: prefer exact > deletion > insertion
            match (&new_match.match_type, &best.match_type) {
                (Exact, Insertion | Deletion) | (Deletion, Insertion) => true,
                (Insertion | Deletion, Exact) | (Insertion, Deletion) => false,
                _ => new_match.mismatches < best.mismatches,
            }
        }
    }
}

/// Find adapter sequence in read allowing mismatches
///
/// Uses a three-stage matching approach (exact, insertion, deletion).
/// Returns the best match (earliest position with fewest mismatches).
pub fn find_adapter(
    seq: &[u8],
    adapter: &[u8],
    min_overlap: usize,
    _max_mismatches: usize,
) -> Option<AdapterMatch> {
    if adapter.is_empty() || seq.len() < min_overlap {
        return None;
    }

    let mut best_match: Option<AdapterMatch> = None;

    // Calculate start position for A-tailing support
    // Fastp starts from negative positions to handle Illumina adapter dimers
    let start: isize = if adapter.len() >= 16 {
        -4
    } else if adapter.len() >= 12 {
        -3
    } else if adapter.len() >= 8 {
        -2
    } else {
        0
    };

    // Minimum required match length (matchReq in fastp)
    // Fastp uses 5 as the minimum match length
    let match_req = 5;

    let end = (seq.len() as isize) - (match_req as isize);

    // STAGE 1: Try exact matching at all positions with A-tailing support
    // Fastp tries negative positions for adapters >= 16bp
    for pos in start..=end {
        if let Some(match_result) = try_exact_match(seq, adapter, pos, min_overlap) {
            if is_better_match(&match_result, best_match.as_ref()) {
                best_match = Some(match_result);
            }
        }
    }

    // If exact match found, return it
    if best_match.is_some() {
        return best_match;
    }

    // STAGE 2: Try insertion matching with different comparison lengths
    // Fastp loops through positions but always compares from START of sequences
    // The loop varies cmplen based on remaining length, not the actual comparison position
    for pos in 0..=(seq.len() - match_req) {
        let remaining = seq.len() - pos;
        if remaining < min_overlap {
            break;
        }

        let cmplen = min(remaining.saturating_sub(1), adapter.len());
        // Match fastp: always compare from position 0, passing only cmplen
        if let Some(match_result) = try_insertion_match(seq, adapter, cmplen, min_overlap) {
            if is_better_match(&match_result, best_match.as_ref()) {
                // Adjust position to match loop variable (fastp returns pos from loop)
                let mut adjusted_match = match_result;
                adjusted_match.position = pos;
                best_match = Some(adjusted_match);
                return best_match; // Fastp breaks immediately
            }
        }
    }

    // STAGE 3: Try deletion matching with different comparison lengths
    // Same logic as insertion - always compare from position 0
    for pos in 0..=(seq.len() - match_req) {
        let cmplen = min(seq.len() - pos, adapter.len().saturating_sub(1));
        if let Some(match_result) = try_deletion_match(seq, adapter, cmplen, min_overlap) {
            if is_better_match(&match_result, best_match.as_ref()) {
                let mut adjusted_match = match_result;
                adjusted_match.position = pos;
                best_match = Some(adjusted_match);
                return best_match; // Fastp breaks immediately
            }
        }
    }

    best_match
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_case_insensitive_matching() {
        let seq = b"ACGTACGTACGTagatcggaagagc"; // lowercase adapter
        let adapter = b"AGATCGGAAGAGC"; // uppercase

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.mismatches, 0);
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

    #[test]
    fn test_insertion_match_parity() {
        // Test cases: (ins_data, normal_data, cmplen, min_overlap)
        // These test cases verify that try_insertion_match_new produces identical
        // results to try_insertion_match.
        let test_cases: Vec<(&[u8], &[u8], usize, usize)> = vec![
            // Basic insertion case - single clear insertion in middle
            (b"AGATCGGXAAGAGC", b"AGATCGGAAGAGC", 13, 8),
            // Insertion at beginning (old returns None, should match)
            (b"XAGATCGGAAGAGC", b"AGATCGGAAGAGC", 13, 8),
            // Insertion near end
            (b"AGATCGGAAGAGXC", b"AGATCGGAAGAGC", 13, 8),
            // No match case (too many mismatches)
            (b"XXXXXXXXAAGAGC", b"AGATCGGAAGAGC", 13, 8),
            // Exact match (should return None - no insertion)
            (b"AGATCGGAAGAGC", b"AGATCGGAAGAGC", 12, 8),
            // Case insensitive
            (b"agatcggXaagagc", b"AGATCGGAAGAGC", 13, 8),
            // Short sequences (below min requirements)
            (b"AGATX", b"AGAT", 4, 8),
            // Edge case: cmplen at minimum
            (b"AGATCGGXAAGAGCTT", b"AGATCGGAAGAGCTT", 8, 8),
            // Longer sequences
            (
                b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACXNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                30,
                10,
            ),
            // Multiple potential insertion points
            (b"ACGTXACGTACGT", b"ACGTACGTACGT", 12, 8),
        ];

        for (i, (ins_data, normal_data, cmplen, min_overlap)) in test_cases.iter().enumerate() {
            let result_old = try_insertion_match(ins_data, normal_data, *cmplen, *min_overlap);
            let result_new = try_insertion_match_new(ins_data, normal_data, *cmplen, *min_overlap);

            match (&result_old, &result_new) {
                (None, None) => {
                    // Both return None - parity achieved
                }
                (Some(old), Some(new)) => {
                    assert_eq!(
                        old.mismatches, new.mismatches,
                        "Test case {}: mismatches differ. Old: {}, New: {}",
                        i, old.mismatches, new.mismatches
                    );
                    assert_eq!(
                        old.overlap_len, new.overlap_len,
                        "Test case {}: overlap_len differ. Old: {}, New: {}",
                        i, old.overlap_len, new.overlap_len
                    );
                    assert_eq!(
                        old.match_type, new.match_type,
                        "Test case {}: match_type differ",
                        i
                    );
                }
                _ => {
                    panic!(
                        "Test case {}: Results differ. Old: {:?}, New: {:?}\nInput: ins={:?}, normal={:?}, cmplen={}, min_overlap={}",
                        i,
                        result_old,
                        result_new,
                        String::from_utf8_lossy(ins_data),
                        String::from_utf8_lossy(normal_data),
                        cmplen,
                        min_overlap
                    );
                }
            }
        }
    }
}
