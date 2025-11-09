//! A fast and simple FASTQ preprocessor
//!
//! This tool processes FASTQ files with quality control filters including:
//! - Length filtering (minimum read length)
//! - Quality filtering (mean quality score)
//! - N-base filtering (maximum ambiguous bases)

use anyhow::{Context, Result};
use clap::Parser;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Parser, Debug)]
#[command(author, version, about = "A fast FASTQ preprocessor", long_about = None)]
struct Args {
    /// Input FASTQ file
    #[arg(short = 'i', long)]
    input: String,

    /// Output FASTQ file
    #[arg(short = 'o', long)]
    output: String,

    /// Minimum length required (default: 15)
    #[arg(short = 'l', long, default_value = "15")]
    length_required: usize,

    /// Mean quality score threshold (default: 0, disabled)
    #[arg(short = 'q', long, default_value = "0")]
    qualified_quality_phred: u8,

    /// Max number of N bases allowed (default: 5)
    #[arg(short = 'n', long, default_value = "5")]
    n_base_limit: usize,

    /// JSON report file (default: fastp.json)
    #[arg(short = 'j', long, default_value = "fastp.json")]
    json: String,
}

/// Statistics for a set of reads
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ReadStats {
    total_reads: usize,
    total_bases: usize,
    q20_bases: usize,
    q30_bases: usize,
    q20_rate: f64,
    q30_rate: f64,
    read1_mean_length: usize,
    gc_content: f64,
}

/// Quality curves data for per-position quality statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
struct QualityCurves {
    #[serde(rename = "A")]
    a: Vec<f64>,
    #[serde(rename = "T")]
    t: Vec<f64>,
    #[serde(rename = "C")]
    c: Vec<f64>,
    #[serde(rename = "G")]
    g: Vec<f64>,
    mean: Vec<f64>,
}

/// Detailed read statistics including kmer counts
#[derive(Debug, Clone, Serialize, Deserialize)]
struct DetailedReadStats {
    total_reads: usize,
    total_bases: usize,
    q20_bases: usize,
    q30_bases: usize,
    quality_curves: QualityCurves,
    kmer_count: IndexMap<String, usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Summary {
    fastp_version: String,
    sequencing: String,
    before_filtering: ReadStats,
    after_filtering: ReadStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FilteringResult {
    passed_filter_reads: usize,
    low_quality_reads: usize,
    #[serde(rename = "too_many_N_reads")]
    too_many_n_reads: usize,
    too_short_reads: usize,
    too_long_reads: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FastpReport {
    summary: Summary,
    filtering_result: FilteringResult,
    read1_before_filtering: DetailedReadStats,
}

/// Represents a single FASTQ record (4 lines)
#[derive(Debug, Clone)]
struct FastqRecord {
    header: String,
    sequence: String,
    plus: String,
    quality: String,
}

impl FastqRecord {
    /// Write this record to the output writer
    fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "{}", self.header)?;
        writeln!(writer, "{}", self.sequence)?;
        writeln!(writer, "{}", self.plus)?;
        writeln!(writer, "{}", self.quality)?;
        Ok(())
    }

    /// Calculate mean quality score (Phred+33 encoding)
    fn mean_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }

        let sum: u32 = self
            .quality
            .bytes()
            .map(|b| (b - 33) as u32) // Phred+33 encoding
            .sum();

        sum as f64 / self.quality.len() as f64
    }

    /// Count the number of N bases in the sequence
    fn count_n_bases(&self) -> usize {
        self.sequence
            .bytes()
            .filter(|&b| b == b'N' || b == b'n')
            .count()
    }

    /// Count bases with quality >= Q20
    fn count_q20_bases(&self) -> usize {
        self.quality.bytes().filter(|&b| (b - 33) >= 20).count()
    }

    /// Count bases with quality >= Q30
    fn count_q30_bases(&self) -> usize {
        self.quality.bytes().filter(|&b| (b - 33) >= 30).count()
    }

    /// Count GC bases
    fn count_gc_bases(&self) -> usize {
        self.sequence
            .bytes()
            .filter(|&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
            .count()
    }
}

/// Parse FASTQ format from a buffered reader
///
/// PERFORMANCE NOTE: This loads ALL records into memory at once.
/// - For large files (100K+ reads), this creates significant memory pressure
/// - Each record allocates 4 strings (header, sequence, plus, quality)
/// - A streaming approach would be more memory-efficient
fn read_fastq_records<R: BufRead>(reader: R) -> Result<Vec<FastqRecord>> {
    let mut records = Vec::new();
    let mut lines = reader.lines();

    while let Some(header) = lines.next() {
        let header = header?;

        // Skip empty lines
        if header.is_empty() {
            continue;
        }

        // Ensure it starts with @
        if !header.starts_with('@') {
            continue;
        }

        // Read sequence
        let sequence = match lines.next() {
            Some(s) => s?,
            None => break,
        };

        // Read plus line
        let plus = match lines.next() {
            Some(p) => p?,
            None => break,
        };

        // Read quality
        let quality = match lines.next() {
            Some(q) => q?,
            None => break,
        };

        // Validate: sequence and quality must have same length
        if sequence.len() != quality.len() {
            eprintln!(
                "WARNING: sequence and quality have different length:\n{}\n{}\n{}\n{}",
                header, sequence, plus, quality
            );
            continue; // Skip this invalid record
        }

        records.push(FastqRecord {
            header,
            sequence,
            plus,
            quality,
        });
    }

    Ok(records)
}

/// Calculate statistics for a set of records
///
/// PERFORMANCE HOTSPOT: Multiple passes through the data
/// - Pass 1: count sequence lengths
/// - Pass 2: count Q20 bases (iterates through quality string)
/// - Pass 3: count Q30 bases (iterates through quality string again)
/// - Pass 4: count GC bases (iterates through sequence string)
/// Could be optimized to a single pass
fn calculate_stats(records: &[FastqRecord]) -> ReadStats {
    let total_reads = records.len();
    let total_bases: usize = records.iter().map(|r| r.sequence.len()).sum();
    let q20_bases: usize = records.iter().map(|r| r.count_q20_bases()).sum();
    let q30_bases: usize = records.iter().map(|r| r.count_q30_bases()).sum();
    let gc_bases: usize = records.iter().map(|r| r.count_gc_bases()).sum();

    ReadStats {
        total_reads,
        total_bases,
        q20_bases,
        q30_bases,
        q20_rate: if total_bases > 0 {
            q20_bases as f64 / total_bases as f64
        } else {
            0.0
        },
        q30_rate: if total_bases > 0 {
            q30_bases as f64 / total_bases as f64
        } else {
            0.0
        },
        read1_mean_length: if total_reads > 0 {
            total_bases / total_reads
        } else {
            0
        },
        gc_content: if total_bases > 0 {
            gc_bases as f64 / total_bases as f64
        } else {
            0.0
        },
    }
}

/// Count all kmers of given size in a sequence
///
/// PERFORMANCE HOTSPOT: String allocations in tight loop
/// - Called for EVERY sequence in the dataset
/// - Allocates a new String for each kmer (e.g., 146 strings per 150bp read)
/// - For 100K reads with 150bp, this is ~14.6 million allocations
/// - Could use a fixed-size buffer or pre-allocated strings
fn count_kmers_in_sequence(sequence: &str, k: usize) -> IndexMap<String, usize> {
    let mut kmer_counts = IndexMap::new();
    let seq_bytes = sequence.as_bytes();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = &seq_bytes[i..i + k];
        // Convert to uppercase and create string
        // HOTSPOT: String allocation for every kmer
        let kmer_str: String = kmer
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        *kmer_counts.entry(kmer_str).or_insert(0) += 1;
    }

    kmer_counts
}

/// Generate all possible kmers of length k
fn generate_all_kmers(k: usize) -> IndexMap<String, usize> {
    let bases = ['A', 'C', 'G', 'T'];
    let mut kmers = IndexMap::new();

    // Generate all possible kmers recursively
    fn generate_recursive(
        k: usize,
        current: String,
        bases: &[char],
        kmers: &mut IndexMap<String, usize>,
    ) {
        if current.len() == k {
            kmers.insert(current, 0);
            return;
        }

        for &base in bases {
            let mut next = current.clone();
            next.push(base);
            generate_recursive(k, next, bases, kmers);
        }
    }

    generate_recursive(k, String::new(), &bases, &mut kmers);
    kmers
}

/// Count all kmers across all records
///
/// PERFORMANCE HOTSPOT: Major bottleneck for large datasets
/// - Calls count_kmers_in_sequence() for every record
/// - For 100K reads with 150bp:
///   - 100,000 records × ~146 kmers per read = 14.6M kmers
///   - Each kmer allocates a String = 14.6M allocations
/// - Then merges into total_counts with more HashMap operations
/// - Could be optimized by:
///   - Using integer encoding for kmers instead of strings
///   - Pre-allocating a single buffer and reusing it
///   - Using array indexing instead of HashMap (4^5 = 1024 possible 5-mers)
fn count_all_kmers(records: &[FastqRecord], k: usize) -> IndexMap<String, usize> {
    // Start with all possible kmers initialized to 0
    let mut total_counts = generate_all_kmers(k);

    // Count kmers that actually appear
    for record in records {
        let seq_kmers = count_kmers_in_sequence(&record.sequence, k);
        for (kmer, count) in seq_kmers {
            *total_counts.entry(kmer).or_insert(0) += count;
        }
    }

    total_counts
}

/// Calculate quality curves for per-position quality statistics
///
/// PERFORMANCE HOTSPOT: Nested iteration over all bases
/// - Allocates 10 vectors of length max_len (e.g., 10 * 151 = 1510 u64s)
/// - Iterates through EVERY base in EVERY record
/// - For 100K reads × 150bp = 15 million base-quality pairs
/// - Many branches in the inner loop (match statement)
/// - Could be optimized with SIMD or lookup tables
fn calculate_quality_curves(records: &[FastqRecord]) -> QualityCurves {
    if records.is_empty() {
        return QualityCurves {
            a: Vec::new(),
            t: Vec::new(),
            c: Vec::new(),
            g: Vec::new(),
            mean: Vec::new(),
        };
    }

    // Find maximum read length
    let max_len = records.iter().map(|r| r.sequence.len()).max().unwrap_or(0);

    // Initialize vectors to track quality sums and counts per position
    let mut a_sums = vec![0u64; max_len];
    let mut a_counts = vec![0u64; max_len];
    let mut t_sums = vec![0u64; max_len];
    let mut t_counts = vec![0u64; max_len];
    let mut c_sums = vec![0u64; max_len];
    let mut c_counts = vec![0u64; max_len];
    let mut g_sums = vec![0u64; max_len];
    let mut g_counts = vec![0u64; max_len];
    let mut total_sums = vec![0u64; max_len];
    let mut total_counts = vec![0u64; max_len];

    // Accumulate quality scores per position per base
    // HOTSPOT: This is the innermost loop, iterating over millions of bases
    for record in records {
        let seq_bytes = record.sequence.as_bytes();
        let qual_bytes = record.quality.as_bytes();

        for (pos, (&base, &qual)) in seq_bytes.iter().zip(qual_bytes.iter()).enumerate() {
            let quality = (qual - 33) as u64;

            total_sums[pos] += quality;
            total_counts[pos] += 1;

            match base {
                b'A' | b'a' => {
                    a_sums[pos] += quality;
                    a_counts[pos] += 1;
                }
                b'T' | b't' => {
                    t_sums[pos] += quality;
                    t_counts[pos] += 1;
                }
                b'C' | b'c' => {
                    c_sums[pos] += quality;
                    c_counts[pos] += 1;
                }
                b'G' | b'g' => {
                    g_sums[pos] += quality;
                    g_counts[pos] += 1;
                }
                _ => {} // Skip N or other bases
            }
        }
    }

    // Calculate overall mean quality per position first
    let mean: Vec<f64> = total_sums
        .iter()
        .zip(total_counts.iter())
        .map(|(&sum, &count)| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                0.0
            }
        })
        .collect();

    // Calculate mean quality per position for each base
    // If no bases of a particular type at a position, use overall mean (matches fastp behavior)
    let a_mean: Vec<f64> = a_sums
        .iter()
        .zip(a_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    let t_mean: Vec<f64> = t_sums
        .iter()
        .zip(t_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    let c_mean: Vec<f64> = c_sums
        .iter()
        .zip(c_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    let g_mean: Vec<f64> = g_sums
        .iter()
        .zip(g_counts.iter())
        .enumerate()
        .map(|(pos, (&sum, &count))| {
            if count > 0 {
                sum as f64 / count as f64
            } else {
                mean[pos]
            }
        })
        .collect();

    QualityCurves {
        a: a_mean,
        t: t_mean,
        c: c_mean,
        g: g_mean,
        mean,
    }
}

/// Calculate detailed statistics with kmer counts
fn calculate_detailed_stats(records: &[FastqRecord]) -> DetailedReadStats {
    let total_reads = records.len();
    let total_bases: usize = records.iter().map(|r| r.sequence.len()).sum();
    let q20_bases: usize = records.iter().map(|r| r.count_q20_bases()).sum();
    let q30_bases: usize = records.iter().map(|r| r.count_q30_bases()).sum();
    let quality_curves = calculate_quality_curves(records);
    let kmer_count = count_all_kmers(records, 5); // Use k=5 to match fastp

    DetailedReadStats {
        total_reads,
        total_bases,
        q20_bases,
        q30_bases,
        quality_curves,
        kmer_count,
    }
}

/// PERFORMANCE ANALYSIS:
///
/// Current implementation makes MULTIPLE PASSES through the data:
/// 1. Load all records into memory (HOTSPOT #1)
/// 2. Calculate before-filtering stats - 4 passes (HOTSPOT #2)
/// 3. Calculate detailed stats - quality_curves + kmer counting (HOTSPOT #3)
///    - quality_curves: iterates every base in every record
///    - kmer counting: millions of string allocations
/// 4. Filter records - another full pass (HOTSPOT #4)
/// 5. Calculate after-filtering stats - 4 more passes (HOTSPOT #5)
///
/// OPTIMIZATION OPPORTUNITIES:
/// - Combine multiple stat calculations into single pass
/// - Use streaming instead of loading all into memory
/// - Reduce string allocations in kmer counting
/// - Avoid cloning records when filtering
/// - Use SIMD for quality score calculations
/// - Pre-allocate buffers for kmer strings
fn main() -> Result<()> {
    let args = Args::parse();

    let input_file =
        File::open(&args.input).context(format!("Failed to open input file: {}", args.input))?;
    let reader = BufReader::new(input_file);

    let output_file = File::create(&args.output)
        .context(format!("Failed to create output file: {}", args.output))?;
    let mut writer = BufWriter::new(output_file);

    // HOTSPOT #1: Parse FASTQ records - loads entire file into memory
    let all_records = read_fastq_records(reader)?;

    // HOTSPOT #2: Calculate before-filtering statistics - 4 passes through records
    let before_stats = calculate_stats(&all_records);
    // HOTSPOT #3: Calculate detailed stats - includes quality_curves and kmer counting
    // - quality_curves: iterates over every base in every record
    // - kmer counting: 14.6M string allocations for 100K reads
    let detailed_before_stats = calculate_detailed_stats(&all_records);

    // Track filtering results
    let mut low_quality_count = 0;
    let mut too_many_n_count = 0;
    let mut too_short_count = 0;
    let mut passed_records = Vec::new();

    // Get max cycle length for sequencing description
    let max_cycle = all_records
        .iter()
        .map(|r| r.sequence.len())
        .max()
        .unwrap_or(0);

    // HOTSPOT #4: Filter records - another pass through all records
    // - mean_quality() scans the entire quality string for each record
    // - count_n_bases() scans the entire sequence string
    // - Clones records that pass (expensive for large datasets)
    for record in &all_records {
        let seq_len = record.sequence.len();
        let mean_qual = record.mean_quality();
        let n_count = record.count_n_bases();

        // Check filters in order (matching fastp priority)
        if seq_len < args.length_required {
            too_short_count += 1;
        } else if n_count > args.n_base_limit {
            too_many_n_count += 1;
        } else if args.qualified_quality_phred > 0
            && mean_qual < args.qualified_quality_phred as f64
        {
            low_quality_count += 1;
        } else {
            // Passed all filters
            record.write_to(&mut writer)?;
            // HOTSPOT: Clone creates 4 new string allocations per passing record
            passed_records.push(record.clone());
        }
    }

    writer.flush()?;

    // HOTSPOT #5: Calculate after-filtering statistics - 4 more passes through passed records
    let after_stats = calculate_stats(&passed_records);

    // Create JSON report
    let report = FastpReport {
        summary: Summary {
            fastp_version: env!("CARGO_PKG_VERSION").to_string(),
            sequencing: format!("single end ({} cycles)", max_cycle),
            before_filtering: before_stats,
            after_filtering: after_stats,
        },
        filtering_result: FilteringResult {
            passed_filter_reads: passed_records.len(),
            low_quality_reads: low_quality_count,
            too_many_n_reads: too_many_n_count,
            too_short_reads: too_short_count,
            too_long_reads: 0,
        },
        read1_before_filtering: detailed_before_stats,
    };

    // Write JSON report
    let json_file =
        File::create(&args.json).context(format!("Failed to create JSON file: {}", args.json))?;
    serde_json::to_writer_pretty(json_file, &report)?;

    Ok(())
}
