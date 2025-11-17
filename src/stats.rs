//! Statistics structures and JSON reporting
//!
//! This module provides all data structures for tracking and reporting
//! FASTQ processing statistics:
//! - `ReadStats`: Basic read statistics
//! - `QualityCurves`: Per-position quality distributions
//! - `DetailedReadStats`: Combined read + quality stats
//! - Summary: Before/after filtering summary
//! - `FilteringResult`: Filtering outcome counts
//! - `FastpReport`: Final JSON report structure
//! - `PositionStats`: Position-specific quality tracking
//! - `SimpleStats`: Simple statistics accumulator

use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

/// Statistics for a set of reads
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct ReadStats {
    pub total_reads: usize,
    pub total_bases: usize,
    pub q20_bases: usize,
    pub q30_bases: usize,
    pub q40_bases: usize,
    pub q20_rate: f64,
    pub q30_rate: f64,
    pub q40_rate: f64,
    pub read1_mean_length: usize,
    pub gc_content: f64,
}

/// Quality curves data for per-position quality statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct QualityCurves {
    #[serde(rename = "A")]
    pub a: Vec<f64>,
    #[serde(rename = "T")]
    pub t: Vec<f64>,
    #[serde(rename = "C")]
    pub c: Vec<f64>,
    #[serde(rename = "G")]
    pub g: Vec<f64>,
    pub mean: Vec<f64>,
}

/// Content curves data for per-position base content percentages
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct ContentCurves {
    #[serde(rename = "A")]
    pub a: Vec<f64>,
    #[serde(rename = "T")]
    pub t: Vec<f64>,
    #[serde(rename = "C")]
    pub c: Vec<f64>,
    #[serde(rename = "G")]
    pub g: Vec<f64>,
    #[serde(rename = "N")]
    pub n: Vec<f64>,
    #[serde(rename = "GC")]
    pub gc: Vec<f64>,
}

/// Detailed read statistics including kmer counts
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct DetailedReadStats {
    pub total_reads: usize,
    pub total_bases: usize,
    pub q20_bases: usize,
    pub q30_bases: usize,
    pub quality_curves: QualityCurves,
    pub content_curves: ContentCurves,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qual_hist: Option<IndexMap<u8, usize>>,
    pub kmer_count: IndexMap<String, usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct Summary {
    pub fastp_version: String,
    pub sequencing: String,
    pub before_filtering: ReadStats,
    pub after_filtering: ReadStats,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read1_before_filtering: Option<ReadStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_before_filtering: Option<ReadStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read1_after_filtering: Option<ReadStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_after_filtering: Option<ReadStats>,
}

#[allow(clippy::struct_field_names)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct FilteringResult {
    pub passed_filter_reads: usize,
    pub low_quality_reads: usize,
    pub low_complexity_reads: usize,
    #[serde(rename = "too_many_N_reads")]
    pub too_many_n_reads: usize,
    pub too_short_reads: usize,
    pub too_long_reads: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct DuplicationStats {
    pub rate: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct AdapterCuttingStats {
    pub adapter_trimmed_reads: usize,
    pub adapter_trimmed_bases: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct InsertSizeStats {
    pub peak: usize,
    pub unknown: usize,
    pub histogram: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct FasterpReport {
    pub summary: Summary,
    pub filtering_result: FilteringResult,
    pub read1_before_filtering: DetailedReadStats,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_before_filtering: Option<DetailedReadStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read1_after_filtering: Option<DetailedReadStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_after_filtering: Option<DetailedReadStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duplication: Option<DuplicationStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub adapter_cutting: Option<AdapterCuttingStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub insert_size: Option<InsertSizeStats>,
}

/// Accumulator for per-position quality statistics
pub(crate) struct PositionStats {
    /// Per-base quality sums: [A, T, C, G][position]
    pub base_sum: [Vec<u64>; 4],
    /// Per-base quality counts: [A, T, C, G][position]
    pub base_cnt: [Vec<u64>; 4],
    /// Total quality sum per position
    pub total_sum: Vec<u64>,
    /// Total count per position
    pub total_cnt: Vec<u64>,
    /// Quality score histogram: count of bases at each quality score (0-93, phred33)
    pub qual_hist: [usize; 94],
}

impl PositionStats {
    pub(crate) fn new() -> Self {
        // Pre-allocate for typical read length (150bp) to avoid reallocations
        const TYPICAL_READ_LEN: usize = 150;

        Self {
            base_sum: [
                Vec::with_capacity(TYPICAL_READ_LEN),
                Vec::with_capacity(TYPICAL_READ_LEN),
                Vec::with_capacity(TYPICAL_READ_LEN),
                Vec::with_capacity(TYPICAL_READ_LEN),
            ],
            base_cnt: [
                Vec::with_capacity(TYPICAL_READ_LEN),
                Vec::with_capacity(TYPICAL_READ_LEN),
                Vec::with_capacity(TYPICAL_READ_LEN),
                Vec::with_capacity(TYPICAL_READ_LEN),
            ],
            total_sum: Vec::with_capacity(TYPICAL_READ_LEN),
            total_cnt: Vec::with_capacity(TYPICAL_READ_LEN),
            qual_hist: [0; 94],
        }
    }

    /// Ensure capacity for at least `len` positions
    pub(crate) fn ensure_capacity(&mut self, len: usize) {
        if self.total_sum.len() < len {
            self.total_sum.resize(len, 0);
            self.total_cnt.resize(len, 0);
            for i in 0..4 {
                self.base_sum[i].resize(len, 0);
                self.base_cnt[i].resize(len, 0);
            }
        }
    }

    /// Convert to `QualityCurves` for JSON output
    pub(crate) fn to_quality_curves(&self) -> QualityCurves {
        let mean: Vec<f64> = self
            .total_sum
            .iter()
            .zip(&self.total_cnt)
            .map(|(&sum, &cnt)| {
                if cnt > 0 {
                    sum as f64 / cnt as f64
                } else {
                    0.0
                }
            })
            .collect();

        let mut curves = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
        for (i, curve) in curves.iter_mut().enumerate() {
            *curve = self.base_sum[i]
                .iter()
                .zip(&self.base_cnt[i])
                .enumerate()
                .map(|(pos, (&sum, &cnt))| {
                    if cnt > 0 {
                        sum as f64 / cnt as f64
                    } else {
                        mean[pos]
                    }
                })
                .collect();
        }

        QualityCurves {
            a: curves[0].clone(),
            t: curves[1].clone(),
            c: curves[2].clone(),
            g: curves[3].clone(),
            mean,
        }
    }

    /// Convert to `ContentCurves` for JSON output (base content percentages per position)
    pub(crate) fn to_content_curves(&self) -> ContentCurves {
        let mut a_pct = Vec::with_capacity(self.total_cnt.len());
        let mut t_pct = Vec::with_capacity(self.total_cnt.len());
        let mut c_pct = Vec::with_capacity(self.total_cnt.len());
        let mut g_pct = Vec::with_capacity(self.total_cnt.len());
        let mut n_pct = Vec::with_capacity(self.total_cnt.len());
        let mut gc_content_pct = Vec::with_capacity(self.total_cnt.len());

        for pos in 0..self.total_cnt.len() {
            let total = self.total_cnt[pos] as f64;
            if total > 0.0 {
                let count_a = self.base_cnt[0].get(pos).copied().unwrap_or(0) as f64;
                let count_t = self.base_cnt[1].get(pos).copied().unwrap_or(0) as f64;
                let count_c = self.base_cnt[2].get(pos).copied().unwrap_or(0) as f64;
                let count_g = self.base_cnt[3].get(pos).copied().unwrap_or(0) as f64;
                let bases_sum = count_a + count_t + count_c + count_g;
                let count_n = total - bases_sum;

                a_pct.push(count_a * 100.0 / total);
                t_pct.push(count_t * 100.0 / total);
                c_pct.push(count_c * 100.0 / total);
                g_pct.push(count_g * 100.0 / total);
                n_pct.push(count_n * 100.0 / total);
                gc_content_pct.push((count_g + count_c) * 100.0 / total);
            } else {
                a_pct.push(0.0);
                t_pct.push(0.0);
                c_pct.push(0.0);
                g_pct.push(0.0);
                n_pct.push(0.0);
                gc_content_pct.push(0.0);
            }
        }

        ContentCurves {
            a: a_pct,
            t: t_pct,
            c: c_pct,
            g: g_pct,
            n: n_pct,
            gc: gc_content_pct,
        }
    }

    /// Convert quality histogram to `IndexMap` for JSON output
    pub(crate) fn to_qual_hist(&self) -> Option<IndexMap<u8, usize>> {
        let mut hist = IndexMap::new();
        for (q, &count) in self.qual_hist.iter().enumerate() {
            if count > 0 {
                hist.insert(q as u8, count);
            }
        }
        if hist.is_empty() { None } else { Some(hist) }
    }
}

/// Simple stats accumulator
#[derive(Default)]
pub(crate) struct SimpleStats {
    pub total_reads: usize,
    pub total_bases: usize,
    pub q20_bases: usize,
    pub q30_bases: usize,
    pub q40_bases: usize,
    pub gc_bases: usize,
}

impl SimpleStats {
    pub(crate) fn add(&mut self, bases: usize, q20: usize, q30: usize, q40: usize, gc: usize) {
        self.total_reads += 1;
        self.total_bases += bases;
        self.q20_bases += q20;
        self.q30_bases += q30;
        self.q40_bases += q40;
        self.gc_bases += gc;
    }

    pub(crate) fn to_read_stats(&self) -> ReadStats {
        ReadStats {
            total_reads: self.total_reads,
            total_bases: self.total_bases,
            q20_bases: self.q20_bases,
            q30_bases: self.q30_bases,
            q40_bases: self.q40_bases,
            q20_rate: if self.total_bases > 0 {
                self.q20_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
            q30_rate: if self.total_bases > 0 {
                self.q30_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
            q40_rate: if self.total_bases > 0 {
                self.q40_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
            read1_mean_length: if self.total_reads > 0 {
                self.total_bases / self.total_reads
            } else {
                0
            },
            gc_content: if self.total_bases > 0 {
                self.gc_bases as f64 / self.total_bases as f64
            } else {
                0.0
            },
        }
    }
}

/// Calculate duplication rate from kmer counts
/// Duplication is estimated by comparing total kmers to unique kmers
pub(crate) fn calculate_duplication_rate(kmer_counts: &indexmap::IndexMap<String, usize>) -> f64 {
    let total_kmers: usize = kmer_counts.values().sum();
    let unique_kmers = kmer_counts.len();

    if total_kmers == 0 {
        return 0.0;
    }

    // Duplication rate = 1 - (unique / total)
    1.0 - (unique_kmers as f64 / total_kmers as f64)
}
