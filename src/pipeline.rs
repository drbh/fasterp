//! Multi-threaded processing pipeline
//!
//! This module provides the 3-stage parallel processing pipeline:
//! - Producer thread: Reads and parses FASTQ records into batches
//! - Worker threads: Process batches in parallel
//! - Merger thread: Writes output in order and reduces statistics
//!
//! This pipeline achieves 2-4x speedup on large datasets with multi-threading.

use anyhow::Result;
use crossbeam_channel::{Receiver, Sender};
use std::collections::BTreeMap;
use std::io::{IoSlice, Read, Write};
use std::ops::Range;
use std::sync::Arc;

use crate::kmer::{base_idx, count_k5_2bit};
use crate::processor::StreamAccumulator;
use crate::simd;
use crate::stats::{PositionStats, SimpleStats};
use crate::trimming::{TrimmingConfig, TrimmingResult, trim_read};

/// Calculate sequence complexity as percentage of bases different from next base
///
/// Complexity is defined as: (count of positions where base[i] != base[i+1]) / (length - 1) * 100
/// This matches fastp's low complexity filter algorithm.
///
/// Returns complexity percentage (0-100)
#[inline]
fn calculate_complexity(seq: &[u8]) -> usize {
    let len = seq.len();
    if len <= 1 {
        return 100; // Single base or empty is considered max complexity
    }

    let mut different_count = 0usize;

    // SAFETY: Loop bounds ensure i and i+1 are always valid indices
    unsafe {
        for i in 0..len - 1 {
            if *seq.get_unchecked(i) != *seq.get_unchecked(i + 1) {
                different_count += 1;
            }
        }
    }

    // Calculate percentage: (different_count / (len - 1)) * 100
    // Use integer math to avoid floating point
    (different_count * 100) / (len - 1)
}

/// A batch of FASTQ records parsed from a buffer
///
/// Contains raw bytes and record positions (no String allocations)
/// Each record is [`header_start`, `seq_start`, `plus_start`, `qual_start`] as byte offsets
///
/// Uses Arc to allow zero-copy sharing of buffer between workers
#[derive(Clone)]
pub(crate) struct Batch {
    pub id: u64,
    pub buf: Arc<Vec<u8>>, // Shared buffer - no copying needed
    /// Each element is [`header_start`, `seq_start`, `plus_start`, `qual_start`]
    /// Lengths are implicit: header len = `seq_start` - `header_start`, etc.
    /// Quality ends at the next record's `header_start` (or `buf.len()` for last record)
    pub recs: Vec<[usize; 4]>,
}

/// A single FASTQ record represented as ranges into a shared buffer
///
/// Zero-copy: stores ranges instead of copying bytes
pub(crate) struct RecordPiece {
    pub buf: Arc<Vec<u8>>,    // Shared buffer reference
    pub header: Range<usize>, // Header range
    pub seq: Range<usize>,    // Sequence range
    pub plus: Range<usize>,   // Plus line range
    pub qual: Range<usize>,   // Quality range
}

/// Result from a worker thread (zero-copy version)
///
/// Uses `RecordPiece` to avoid copying bytes - stores ranges instead
pub(crate) struct WorkerResult {
    pub id: u64,
    pub pieces: Vec<RecordPiece>, // Zero-copy: ranges into shared buffers
    pub before: SimpleStats,
    pub after: SimpleStats,
    pub pos: PositionStats,
    pub pos_after: PositionStats,
    pub k5: [usize; 1024],
    pub k5_after: [usize; 1024],
    // Filter counts
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
    // Adapter trimming stats
    pub adapter_trimmed_reads: usize,
    pub adapter_trimmed_bases: usize,
}

// PAIRED-END MULTI-THREADING STRUCTURES

/// A batch of paired-end FASTQ records
///
/// Contains raw bytes and record positions for both R1 and R2
/// Uses Arc for zero-copy sharing between workers
#[derive(Clone)]
pub(crate) struct PairedBatch {
    pub id: u64,
    pub buf_r1: Arc<Vec<u8>>,     // R1 shared buffer
    pub buf_r2: Arc<Vec<u8>>,     // R2 shared buffer
    pub recs_r1: Vec<[usize; 4]>, // R1 record positions
    pub recs_r2: Vec<[usize; 4]>, // R2 record positions
}

/// A paired-end record represented as ranges into shared buffers
///
/// Zero-copy: stores ranges for both R1 and R2
pub(crate) struct PairedRecordPiece {
    pub r1: RecordPiece,
    pub r2: RecordPiece,
}

/// Result from a paired-end worker thread
pub(crate) struct PairedWorkerResult {
    pub id: u64,
    pub pieces: Vec<PairedRecordPiece>,
    // R1 stats
    pub before_r1: SimpleStats,
    pub after_r1: SimpleStats,
    pub pos_r1: PositionStats,
    pub pos_r1_after: PositionStats,
    pub k5_r1: [usize; 1024],
    pub k5_r1_after: [usize; 1024],
    // R2 stats
    pub before_r2: SimpleStats,
    pub after_r2: SimpleStats,
    pub pos_r2: PositionStats,
    pub pos_r2_after: PositionStats,
    pub k5_r2: [usize; 1024],
    pub k5_r2_after: [usize; 1024],
    // Filter counts (shared - pairs pass/fail together)
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
    pub duplicated: usize,
    // Adapter trimming stats
    pub adapter_trimmed_reads: usize,
    pub adapter_trimmed_bases: usize,
    // Insert size stats
    pub insert_size_histogram: Vec<usize>,
    pub insert_size_unknown: usize,
}

/// Producer thread: read blocks and parse into batches
///
/// Reads large blocks (`batch_bytes`) from input, parses FASTQ records,
/// and emits Batch structures with NO string allocations - just byte slices.
///
/// Handles partial records at block boundaries by carrying them over to next batch.
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn producer_thread(
    input_path: String,
    batch_bytes: usize,
    sender: Sender<Option<Batch>>,
) -> Result<()> {
    let mut reader = crate::io::open_input(&input_path)?; // Handles compression automatically

    let mut batch_id = 0u64;
    let mut carryover = Vec::new();
    let mut buffer = vec![0u8; batch_bytes]; // Reuse buffer allocation

    loop {
        // Read a chunk - reusing buffer
        let bytes_read = reader.read(&mut buffer)?;

        // Prepend carryover from previous iteration
        let carryover_len = carryover.len();
        let actual_len = if carryover_len > 0 {
            // Move buffer data to make room for carryover at the start
            if carryover_len + bytes_read > buffer.len() {
                // Need to grow buffer
                buffer.resize(carryover_len + bytes_read, 0);
            }
            // Shift read data to make room
            buffer.copy_within(0..bytes_read, carryover_len);
            // Copy carryover to start
            buffer[..carryover_len].copy_from_slice(&carryover);
            carryover.clear();
            carryover_len + bytes_read
        } else {
            bytes_read
        };

        if actual_len == 0 {
            break; // EOF and no carryover
        }

        // Determine if this is the last chunk
        // Note: For compressed streams, bytes_read < batch_bytes doesn't mean EOF!
        // Decompression streams can return fewer bytes than requested even when more data exists.
        let is_eof = bytes_read == 0;

        // SIMD-accelerated newline scan using memchr (10-50x faster than naive loop)
        let mut complete_end = 0;
        let mut line_count = 0;
        let mut line_starts = vec![0];

        for newline_pos in memchr::memchr_iter(b'\n', &buffer[..actual_len]) {
            line_count += 1;
            // After every 4 lines, we have a complete record
            if line_count % 4 == 0 {
                complete_end = newline_pos + 1;
            }
            // Track line starts for parsing
            if newline_pos + 1 < actual_len {
                line_starts.push(newline_pos + 1);
            }
        }

        // On EOF, if we have remaining lines that form complete records, use them
        let complete_len = if is_eof && line_count > 0 {
            // Check if we have complete records
            if line_count % 4 == 0 {
                // All lines end with newline and form complete records
                actual_len
            } else if line_count % 4 == 3 && complete_end < actual_len {
                // We have 3 newlines, meaning 4th line exists but no trailing newline
                // This is still a complete record
                actual_len
            } else {
                complete_end
            }
        } else {
            complete_end
        };

        // Save incomplete part for next iteration (if not EOF)
        if !is_eof && complete_len < actual_len {
            carryover = buffer[complete_len..actual_len].to_vec();
        }

        // Parse complete records - use pre-computed line_starts
        if complete_len > 0 {
            // Filter line_starts to only those within complete_len
            line_starts.retain(|&pos| pos < complete_len);

            // Group into 4-line records
            let mut recs = Vec::new();
            for i in (0..line_starts.len()).step_by(4) {
                if i + 3 < line_starts.len() {
                    recs.push([
                        line_starts[i],
                        line_starts[i + 1],
                        line_starts[i + 2],
                        line_starts[i + 3],
                    ]);
                }
            }

            if !recs.is_empty() {
                let batch = Batch {
                    id: batch_id,
                    buf: Arc::new(buffer[..complete_len].to_vec()), // Wrap in Arc for sharing
                    recs,
                };
                batch_id += 1;

                if sender.send(Some(batch)).is_err() {
                    break; // Receiver disconnected
                }
            }
        }

        if is_eof {
            break;
        }
    }

    // Send sentinel
    let _ = sender.send(None);
    Ok(())
}

/// Worker thread: process batches with thread-local accumulators
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn worker_thread(
    receiver: Receiver<Option<Batch>>,
    sender: Sender<Option<WorkerResult>>,
    min_len: usize,
    n_limit: usize,
    qualified_quality_phred: u8,
    unqualified_percent_limit: usize,
    average_qual: u8,
    low_complexity_filter: bool,
    complexity_threshold: usize,
    no_kmer: bool,
    trimming_config: TrimmingConfig,
) {
    while let Ok(Some(batch)) = receiver.recv() {
        let mut before = SimpleStats::default();
        let mut after = SimpleStats::default();
        let mut pos = PositionStats::new();
        let mut pos_after = PositionStats::new();
        let mut k5 = [0usize; 1024];
        let mut k5_after = [0usize; 1024];
        let mut too_short = 0usize;
        let mut too_many_n = 0usize;
        let mut low_quality = 0usize;
        let mut low_complexity = 0usize;
        let mut invalid = 0usize;
        let mut adapter_trimmed_reads = 0usize;
        let mut adapter_trimmed_bases = 0usize;
        let mut pieces = Vec::new(); // Zero-copy: store ranges instead of bytes

        // Process each record in the batch
        for (idx, &[h_start, s_start, p_start, q_start]) in batch.recs.iter().enumerate() {
            // Calculate end positions
            let s_end = p_start - 1; // -1 to skip newline
            let p_end = q_start - 1;
            let q_end = if idx + 1 < batch.recs.len() {
                batch.recs[idx + 1][0] - 1
            } else {
                // For the last record, exclude trailing newline if present
                let buf_len = batch.buf.len();
                if buf_len > 0 && batch.buf[buf_len - 1] == b'\n' {
                    buf_len - 1
                } else {
                    buf_len
                }
            };

            // Extract slices for processing
            let seq = &batch.buf[s_start..s_end];
            let qual = &batch.buf[q_start..q_end];

            // Validate
            if seq.len() != qual.len() {
                invalid += 1;
                continue;
            }

            // Compute stats - use SIMD when available, otherwise single-pass
            let (_qsum, q20, q30, q40, _ncnt, gc) = if simd::is_simd_available() {
                // SIMD path: compute basic stats fast, then position-specific
                let stats = simd::compute_stats(seq, qual, 0); // 0 = don't count unqualified for before stats

                pos.ensure_capacity(seq.len());
                for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
                    let quality_score = u64::from(q - 33);
                    pos.total_sum[i] += quality_score;
                    pos.total_cnt[i] += 1;

                    if let Some(bi) = base_idx(b) {
                        pos.base_sum[bi][i] += quality_score;
                        pos.base_cnt[bi][i] += 1;
                    }

                    // Track quality histogram
                    if quality_score < 94 {
                        pos.qual_hist[quality_score as usize] += 1;
                    }
                }

                (
                    stats.qsum, stats.q20, stats.q30, stats.q40, stats.ncnt, stats.gc,
                )
            } else {
                // Non-SIMD path: single pass for both basic and position stats
                let mut qsum = 0u32;
                let mut q20 = 0usize;
                let mut q30 = 0usize;
                let mut q40 = 0usize;
                let mut ncnt = 0usize;
                let mut gc = 0usize;

                pos.ensure_capacity(seq.len());
                for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
                    let quality_u32 = u32::from(q - 33);
                    qsum += quality_u32;
                    if q >= 53 {
                        q20 += 1;
                    }
                    if q >= 63 {
                        q30 += 1;
                    }
                    if q >= 73 {
                        q40 += 1;
                    }
                    if b == b'N' || b == b'n' {
                        ncnt += 1;
                    }
                    if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                        gc += 1;
                    }

                    pos.total_sum[i] += u64::from(quality_u32);
                    pos.total_cnt[i] += 1;

                    if let Some(bi) = base_idx(b) {
                        pos.base_sum[bi][i] += u64::from(quality_u32);
                        pos.base_cnt[bi][i] += 1;
                    }

                    // Track quality histogram
                    if (quality_u32 as usize) < 94 {
                        pos.qual_hist[quality_u32 as usize] += 1;
                    }
                }

                (qsum, q20, q30, q40, ncnt, gc)
            };

            // K-mer counting
            if !no_kmer {
                count_k5_2bit(seq, &mut k5);
            }

            // Update before stats
            before.add(seq.len(), q20, q30, q40, gc);

            // Apply trimming if enabled
            let trimming_result = if trimming_config.is_enabled() {
                trim_read(seq, qual, &trimming_config)
            } else {
                TrimmingResult {
                    start_pos: 0,
                    end_pos: seq.len(),
                    poly_g_trimmed: 0,
                    poly_x_trimmed: 0,
                    adapter_trimmed: false,
                    adapter_bases_trimmed: 0,
                }
            };

            // Track adapter trimming stats
            if trimming_result.adapter_trimmed {
                adapter_trimmed_reads += 1;
                adapter_trimmed_bases += trimming_result.adapter_bases_trimmed;
            }

            // Get trimmed sequences
            let trimmed_seq = &seq[trimming_result.start_pos..trimming_result.end_pos];
            let trimmed_qual = &qual[trimming_result.start_pos..trimming_result.end_pos];
            let trimmed_len = trimmed_seq.len();

            // Recompute stats for trimmed read (used for filtering) - SIMD accelerated
            let trimmed_stats =
                simd::compute_stats(trimmed_seq, trimmed_qual, qualified_quality_phred);
            let trimmed_qsum = trimmed_stats.qsum;
            let trimmed_q20 = trimmed_stats.q20;
            let trimmed_q30 = trimmed_stats.q30;
            let trimmed_ncnt = trimmed_stats.ncnt;
            let trimmed_gc = trimmed_stats.gc;
            let unqualified_count = trimmed_stats.unqualified;

            // Apply filters in fastp order: quality first, then length, then N count, then complexity

            // Check unqualified percent (fastp -q/-u logic)
            // Match fastp's exact order of operations: lowQualNum > (unqualifiedPercentLimit * rlen / 100.0)
            if qualified_quality_phred > 0 && trimmed_len > 0 {
                let threshold = (unqualified_percent_limit * trimmed_len) as f64 / 100.0;
                if (unqualified_count as f64) > threshold {
                    low_quality += 1;
                    continue;
                }
            }

            // Check average quality (fastp -e logic)
            // Use integer division to match fastp exactly
            if average_qual > 0 && trimmed_len > 0 {
                let mean_qual = trimmed_qsum / trimmed_len as u32;
                if mean_qual < u32::from(average_qual) {
                    low_quality += 1;
                    continue;
                }
            }

            // Check N base limit (fastp -n logic) - fastp checks this AFTER quality but BEFORE length
            if trimmed_ncnt > n_limit {
                too_many_n += 1;
                continue;
            }

            // Length check (fastp -l logic) - AFTER quality checks to match fastp behavior
            if trimmed_len < min_len {
                too_short += 1;
                continue;
            }

            // Check low complexity (fastp -y/-Y logic)
            if low_complexity_filter && trimmed_len > 0 {
                let complexity = calculate_complexity(trimmed_seq);
                if complexity < complexity_threshold {
                    low_complexity += 1;
                    continue;
                }
            }

            // Passed - emit RANGES for TRIMMED read (zero-copy!)
            // Calculate trimmed ranges relative to original buffer
            let trimmed_seq_start = s_start + trimming_result.start_pos;
            let trimmed_seq_end = s_start + trimming_result.end_pos;
            let trimmed_qual_start = q_start + trimming_result.start_pos;
            let trimmed_qual_end = q_start + trimming_result.end_pos;

            pieces.push(RecordPiece {
                buf: Arc::clone(&batch.buf),  // Clone Arc, not the buffer!
                header: h_start..s_start - 1, // Exclude newline
                seq: trimmed_seq_start..trimmed_seq_end,
                plus: p_start..p_end,
                qual: trimmed_qual_start..trimmed_qual_end,
            });

            // Update after stats with trimmed read
            after.add(
                trimmed_len,
                trimmed_q20,
                trimmed_q30,
                trimmed_stats.q40,
                trimmed_gc,
            );

            // Track position stats for after-filtering data (including histogram)
            pos_after.ensure_capacity(trimmed_len);
            for (i, (&b, &q)) in trimmed_seq.iter().zip(trimmed_qual).enumerate() {
                let quality_score = u64::from(q - 33);
                pos_after.total_sum[i] += quality_score;
                pos_after.total_cnt[i] += 1;

                if let Some(bi) = base_idx(b) {
                    pos_after.base_sum[bi][i] += quality_score;
                    pos_after.base_cnt[bi][i] += 1;
                }

                // Track quality histogram for after-filtering
                if quality_score < 94 {
                    pos_after.qual_hist[quality_score as usize] += 1;
                }
            }

            // K-mer counting for after-filtering data
            if !no_kmer {
                count_k5_2bit(trimmed_seq, &mut k5_after);
            }
        }

        let result = WorkerResult {
            id: batch.id,
            pieces, // Zero-copy ranges instead of copied bytes
            before,
            after,
            pos,
            pos_after,
            k5,
            k5_after,
            too_short,
            too_many_n,
            low_quality,
            low_complexity,
            invalid,
            adapter_trimmed_reads,
            adapter_trimmed_bases,
        };

        if sender.send(Some(result)).is_err() {
            break; // Receiver disconnected
        }
    }

    // Send sentinel
    let _ = sender.send(None);
}

/// Merger thread: write output in order and reduce stats (zero-copy vectored I/O)
#[allow(clippy::too_many_lines)]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn merger_thread(
    receiver: Receiver<Option<WorkerResult>>,
    output_path: String,
    num_workers: usize,
    compression_level: Option<u32>,
    split_config: crate::split::SplitConfig,
    parallel_compression: bool,
) -> Result<StreamAccumulator> {
    // Create split writer (OutputWriter already has buffering, so no need for additional BufWriter)
    let compression = compression_level.unwrap_or(6);
    let mut writer = crate::split::SplitWriter::new(
        &output_path,
        split_config,
        compression,
        parallel_compression,
    )?;

    let mut acc = StreamAccumulator::new();
    let mut next_id = 0u64;
    let mut pending: BTreeMap<u64, WorkerResult> = BTreeMap::new();
    let mut workers_done = 0;

    let newline = [b'\n']; // Reusable newline for vectored I/O

    while workers_done < num_workers {
        match receiver.recv() {
            Ok(Some(result)) => {
                pending.insert(result.id, result);

                // Write all consecutive results starting from next_id
                while let Some(result) = pending.remove(&next_id) {
                    // Write output using zero-copy vectored I/O
                    for rec in &result.pieces {
                        let b = &rec.buf;
                        // Build IoSlice array for vectored write (no copying!)
                        let mut iov = [
                            IoSlice::new(&b[rec.header.clone()]),
                            IoSlice::new(&newline),
                            IoSlice::new(&b[rec.seq.clone()]),
                            IoSlice::new(&newline),
                            IoSlice::new(&b[rec.plus.clone()]),
                            IoSlice::new(&newline),
                            IoSlice::new(&b[rec.qual.clone()]),
                            IoSlice::new(&newline),
                        ];

                        // Manual implementation of write_all_vectored
                        let mut written = 0;
                        let total: usize = iov.iter().map(|s| s.len()).sum();
                        while written < total {
                            let n = writer.write_vectored(&iov)?;
                            if n == 0 {
                                return Err(std::io::Error::new(
                                    std::io::ErrorKind::WriteZero,
                                    "failed to write vectored",
                                )
                                .into());
                            }
                            written += n;
                            // Advance the slices
                            let mut skip = n;
                            for slice in &mut iov {
                                let len = slice.len();
                                if skip >= len {
                                    skip -= len;
                                    *slice = IoSlice::new(&[]);
                                } else {
                                    let data = unsafe {
                                        std::slice::from_raw_parts(
                                            slice.as_ptr().add(skip),
                                            len - skip,
                                        )
                                    };
                                    *slice = IoSlice::new(data);
                                    break;
                                }
                            }
                        }
                    }

                    // Merge stats
                    acc.before.total_reads += result.before.total_reads;
                    acc.before.total_bases += result.before.total_bases;
                    acc.before.q20_bases += result.before.q20_bases;
                    acc.before.q30_bases += result.before.q30_bases;
                    acc.before.q40_bases += result.before.q40_bases;
                    acc.before.gc_bases += result.before.gc_bases;

                    acc.after.total_reads += result.after.total_reads;
                    acc.after.total_bases += result.after.total_bases;
                    acc.after.q20_bases += result.after.q20_bases;
                    acc.after.q30_bases += result.after.q30_bases;
                    acc.after.q40_bases += result.after.q40_bases;
                    acc.after.gc_bases += result.after.gc_bases;

                    // Merge position stats
                    acc.pos.ensure_capacity(result.pos.total_sum.len());
                    for i in 0..result.pos.total_sum.len() {
                        acc.pos.total_sum[i] += result.pos.total_sum[i];
                        acc.pos.total_cnt[i] += result.pos.total_cnt[i];
                        for b in 0..4 {
                            acc.pos.base_sum[b][i] += result.pos.base_sum[b][i];
                            acc.pos.base_cnt[b][i] += result.pos.base_cnt[b][i];
                        }
                    }

                    // Merge quality histograms
                    for (q, &count) in result.pos.qual_hist.iter().enumerate() {
                        acc.pos.qual_hist[q] += count;
                    }

                    // Merge after-filtering position stats
                    acc.pos_after
                        .ensure_capacity(result.pos_after.total_sum.len());
                    for i in 0..result.pos_after.total_sum.len() {
                        acc.pos_after.total_sum[i] += result.pos_after.total_sum[i];
                        acc.pos_after.total_cnt[i] += result.pos_after.total_cnt[i];
                        for b in 0..4 {
                            acc.pos_after.base_sum[b][i] += result.pos_after.base_sum[b][i];
                            acc.pos_after.base_cnt[b][i] += result.pos_after.base_cnt[b][i];
                        }
                    }

                    // Merge after-filtering quality histograms
                    for (q, &count) in result.pos_after.qual_hist.iter().enumerate() {
                        acc.pos_after.qual_hist[q] += count;
                    }

                    // Merge k-mer counts
                    for (i, &count) in result.k5.iter().enumerate() {
                        acc.kmer_table[i] += count;
                    }

                    // Merge after-filtering k-mer counts
                    for (i, &count) in result.k5_after.iter().enumerate() {
                        acc.kmer_table_after[i] += count;
                    }

                    // Merge filter counts
                    acc.too_short += result.too_short;
                    acc.too_many_n += result.too_many_n;
                    acc.low_quality += result.low_quality;
                    acc.low_complexity += result.low_complexity;
                    acc.invalid += result.invalid;

                    // Merge adapter trimming stats
                    acc.adapter_trimmed_reads += result.adapter_trimmed_reads;
                    acc.adapter_trimmed_bases += result.adapter_trimmed_bases;

                    // Track max cycle
                    if result.pos.total_sum.len() > acc.max_cycle {
                        acc.max_cycle = result.pos.total_sum.len();
                    }

                    next_id += 1;
                }
            }
            Ok(None) => {
                workers_done += 1;
            }
            Err(_) => break,
        }
    }

    // Finish split writer
    writer.finish()?;
    Ok(acc)
}

// PAIRED-END MULTI-THREADING PIPELINE

/// Paired-end producer thread: opens and reads both R1 and R2 files
///
/// Opens files IN THIS THREAD to isolate `GzDecoders` from each other (workaround for flate2 bug).
/// Reads blocks from both input files, parses FASTQ records, and emits
/// `PairedBatch` structures with synchronized record counts.
#[allow(clippy::too_many_lines)]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn paired_producer_thread_with_paths(
    input1: String,
    input2: String,
    batch_bytes: usize,
    sender: Sender<Option<PairedBatch>>,
    _num_workers: usize,
) -> Result<()> {
    use crossbeam_channel::bounded;
    use std::io::Read;

    // Use bounded channels with larger capacity to prevent deadlock
    // Capacity of 16 provides enough buffering to handle speed differences between R1 and R2
    let (r1_sender, r1_receiver) = bounded::<(Vec<u8>, usize)>(16); // R1 decompressed chunks
    let (r2_sender, r2_receiver) = bounded::<(Vec<u8>, usize)>(16); // R2 decompressed chunks

    // Spawn dedicated decompressor thread for R1
    let input1_clone = input1.clone();
    let decomp_chunk_size = batch_bytes; // Use same size as batch_bytes
    let _decomp_thread1 = std::thread::spawn(move || -> Result<()> {
        let mut reader = crate::io::open_input(&input1_clone)?;
        // Reuse buffer to reduce allocations
        let mut buffer = vec![0u8; decomp_chunk_size];
        loop {
            // Keep reading until buffer is full or EOF
            // This is critical for gzip streams which return small chunks
            let mut total_read = 0;
            while total_read < decomp_chunk_size {
                match reader.read(&mut buffer[total_read..])? {
                    0 => break, // EOF
                    n => total_read += n,
                }
            }

            if total_read == 0 {
                break;
            }
            // Send a copy of the data (receiver owns it)
            let data = buffer[..total_read].to_vec();
            if r1_sender.send((data, total_read)).is_err() {
                break; // Receiver dropped
            }
        }
        Ok(())
    });

    // Spawn dedicated decompressor thread for R2
    let input2_clone = input2.clone();
    let _decomp_thread2 = std::thread::spawn(move || -> Result<()> {
        let mut reader = crate::io::open_input(&input2_clone)?;
        // Reuse buffer to reduce allocations
        let mut buffer = vec![0u8; decomp_chunk_size];
        loop {
            // Keep reading until buffer is full or EOF
            // This is critical for gzip streams which return small chunks
            let mut total_read = 0;
            while total_read < decomp_chunk_size {
                match reader.read(&mut buffer[total_read..])? {
                    0 => break, // EOF
                    n => total_read += n,
                }
            }

            if total_read == 0 {
                break;
            }
            // Send a copy of the data (receiver owns it)
            let data = buffer[..total_read].to_vec();
            if r2_sender.send((data, total_read)).is_err() {
                break; // Receiver dropped
            }
        }
        Ok(())
    });

    let mut batch_id = 0u64;
    let mut carryover1 = Vec::new();
    let mut carryover2 = Vec::new();

    loop {
        // With larger bounded channels (capacity 16), deadlock is unlikely
        // Simple sequential receives work fine
        let (buffer1, bytes_read1) = r1_receiver.recv().unwrap_or_default();
        let (buffer2, bytes_read2) = r2_receiver.recv().unwrap_or_default();

        // Prepend carryover for R1
        let carryover_len1 = carryover1.len();
        let (data1, actual_len1) = if carryover_len1 > 0 {
            let mut combined = Vec::with_capacity(carryover_len1 + bytes_read1);
            combined.extend_from_slice(&carryover1);
            combined.extend_from_slice(&buffer1[..bytes_read1]);
            carryover1.clear();
            let len = combined.len();
            (combined, len)
        } else {
            (buffer1, bytes_read1)
        };

        // Prepend carryover for R2
        let carryover_len2 = carryover2.len();
        let (data2, actual_len2) = if carryover_len2 > 0 {
            let mut combined = Vec::with_capacity(carryover_len2 + bytes_read2);
            combined.extend_from_slice(&carryover2);
            combined.extend_from_slice(&buffer2[..bytes_read2]);
            carryover2.clear();
            let len = combined.len();
            (combined, len)
        } else {
            (buffer2, bytes_read2)
        };

        // Check for EOF
        if actual_len1 == 0 && actual_len2 == 0 {
            break; // Both files ended and no carryover
        }

        // Verify both files have data (or both don't)
        if (actual_len1 == 0) != (actual_len2 == 0) {
            anyhow::bail!(
                "Paired-end files are not synchronized: R1 has {actual_len1} bytes, R2 has {actual_len2} bytes"
            );
        }

        let is_eof1 = bytes_read1 == 0;
        let is_eof2 = bytes_read2 == 0;

        // Parse R1
        let (complete_len1, mut recs1) = parse_fastq_buffer(&data1, actual_len1, is_eof1);

        // Parse R2
        let (complete_len2, mut recs2) = parse_fastq_buffer(&data2, actual_len2, is_eof2);

        // Handle record count mismatch due to gzip decompression boundaries
        // Take the minimum number of complete records from both files
        let num_pairs = std::cmp::min(recs1.len(), recs2.len());

        let adjusted_complete_len1;
        let adjusted_complete_len2;

        if recs1.len() == recs2.len() {
            // Same number of records - use all complete data
            adjusted_complete_len1 = complete_len1;
            adjusted_complete_len2 = complete_len2;

            // Save incomplete data to carryover
            if complete_len1 < actual_len1 {
                carryover1.extend_from_slice(&data1[complete_len1..actual_len1]);
            }
            if complete_len2 < actual_len2 {
                carryover2.extend_from_slice(&data2[complete_len2..actual_len2]);
            }
        } else {
            // Trim to same number of records
            if recs1.len() > num_pairs {
                // R1 has more records - truncate and save excess to carryover
                recs1.truncate(num_pairs);
                // Find the byte position after the last record we're keeping
                if num_pairs > 0 {
                    let last_rec = recs1[num_pairs - 1];
                    // Find end of last record (after 4th line's newline)
                    adjusted_complete_len1 = data1[last_rec[3]..complete_len1]
                        .iter()
                        .position(|&b| b == b'\n')
                        .map_or(complete_len1, |pos| last_rec[3] + pos + 1);
                } else {
                    adjusted_complete_len1 = 0;
                }
                adjusted_complete_len2 = complete_len2;
            } else {
                // R2 has more records - truncate and save excess to carryover
                recs2.truncate(num_pairs);
                if num_pairs > 0 {
                    let last_rec = recs2[num_pairs - 1];
                    adjusted_complete_len2 = data2[last_rec[3]..complete_len2]
                        .iter()
                        .position(|&b| b == b'\n')
                        .map_or(complete_len2, |pos| last_rec[3] + pos + 1);
                } else {
                    adjusted_complete_len2 = 0;
                }
                adjusted_complete_len1 = complete_len1;
            }

            // Save excess data to carryover
            if adjusted_complete_len1 < actual_len1 {
                carryover1.extend_from_slice(&data1[adjusted_complete_len1..actual_len1]);
            }
            if adjusted_complete_len2 < actual_len2 {
                carryover2.extend_from_slice(&data2[adjusted_complete_len2..actual_len2]);
            }
        }

        if !recs1.is_empty() {
            let batch = PairedBatch {
                id: batch_id,
                buf_r1: Arc::new(data1[..adjusted_complete_len1].to_vec()),
                buf_r2: Arc::new(data2[..adjusted_complete_len2].to_vec()),
                recs_r1: recs1,
                recs_r2: recs2,
            };
            batch_id += 1;

            if sender.send(Some(batch)).is_err() {
                break; // Receiver disconnected
            }
        }
    }

    // Send sentinel
    let _ = sender.send(None);
    Ok(())
}

/// Helper function to parse FASTQ records from a buffer
/// Returns (`complete_len`, records)
fn parse_fastq_buffer(buffer: &[u8], actual_len: usize, is_eof: bool) -> (usize, Vec<[usize; 4]>) {
    if actual_len == 0 {
        return (0, Vec::new());
    }

    // SIMD-accelerated newline search using memchr (10-50x faster than naive loop)
    let mut complete_end = 0;
    let mut line_count = 0;
    let mut line_starts = Vec::with_capacity(actual_len / 100); // Pre-allocate for ~100 byte lines
    line_starts.push(0);

    // Use memchr to find all newlines with SIMD
    for newline_pos in memchr::memchr_iter(b'\n', &buffer[..actual_len]) {
        line_count += 1;
        // After every 4 lines, we have a complete record
        if line_count % 4 == 0 {
            complete_end = newline_pos + 1;
        }
        // Track line starts for parsing
        if newline_pos + 1 < actual_len {
            line_starts.push(newline_pos + 1);
        }
    }

    // On EOF, if we have remaining lines that form complete records, use them
    let complete_len = if is_eof && line_count > 0 {
        // Check if we have complete records
        if line_count % 4 == 0 {
            // All lines end with newline and form complete records
            actual_len
        } else if line_count % 4 == 3 && complete_end < actual_len {
            // We have 3 newlines, meaning 4th line exists but no trailing newline
            // This is still a complete record
            actual_len
        } else {
            complete_end
        }
    } else {
        complete_end
    };

    // Parse complete records - use pre-computed line_starts
    // Pre-allocate based on line count (4 lines per record)
    let mut recs = Vec::with_capacity(line_count / 4);
    if complete_len > 0 {
        // Filter line_starts to only those within complete_len
        line_starts.retain(|&pos| pos < complete_len);

        // Group into 4-line records
        for i in (0..line_starts.len()).step_by(4) {
            if i + 3 < line_starts.len() {
                recs.push([
                    line_starts[i],
                    line_starts[i + 1],
                    line_starts[i + 2],
                    line_starts[i + 3],
                ]);
            }
        }
    }

    (complete_len, recs)
}

/// Paired-end worker thread: process read pairs with synchronization
///
/// Processes both R1 and R2 together. Read pairs pass/fail together -
/// if either read fails ANY filter, BOTH reads are discarded.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn paired_worker_thread(
    receiver: Receiver<Option<PairedBatch>>,
    sender: Sender<Option<PairedWorkerResult>>,
    min_len: usize,
    n_limit: usize,
    qualified_quality_phred: u8,
    unqualified_percent_limit: usize,
    average_qual: u8,
    low_complexity_filter: bool,
    complexity_threshold: usize,
    no_kmer: bool,
    trimming_config_r1: TrimmingConfig,
    trimming_config_r2: TrimmingConfig,
    overlap_config: Option<crate::overlap::OverlapConfig>,
    umi_config: Option<crate::umi::UmiConfig>,
    dedup_config: Option<crate::dedup::DedupConfig>,
    dup_detector: Option<std::sync::Arc<crate::bloom::DuplicateDetector>>,
) {
    // Initialize dedup tracker if deduplication is enabled
    let mut dedup_tracker = if let Some(ref cfg) = dedup_config {
        if cfg.enabled {
            Some(crate::dedup::DedupTracker::new(cfg.accuracy))
        } else {
            None
        }
    } else {
        None
    };

    while let Ok(Some(batch)) = receiver.recv() {
        let mut before_r1 = SimpleStats::default();
        let mut after_r1 = SimpleStats::default();
        let mut pos_r1 = PositionStats::new();
        let mut pos_r1_after = PositionStats::new();
        let mut k5_r1 = [0usize; 1024];
        let mut k5_r1_after = [0usize; 1024];

        let mut before_r2 = SimpleStats::default();
        let mut after_r2 = SimpleStats::default();
        let mut pos_r2 = PositionStats::new();
        let mut pos_r2_after = PositionStats::new();
        let mut k5_r2 = [0usize; 1024];
        let mut k5_r2_after = [0usize; 1024];

        let mut too_short = 0usize;
        let mut too_many_n = 0usize;
        let mut low_quality = 0usize;
        let mut low_complexity = 0usize;
        let mut invalid = 0usize;
        let mut duplicated = 0usize;
        let mut adapter_trimmed_reads = 0usize;
        let mut adapter_trimmed_bases = 0usize;
        let mut insert_size_histogram = vec![0usize; 512];
        let mut insert_size_unknown = 0usize;
        let mut pieces = Vec::new();

        // Pre-allocate reusable buffers for output (avoids allocation per read)
        // Typical read ~150bp * 4 lines + headers ≈ 700 bytes
        let mut r1_buf_reusable = Vec::with_capacity(800);
        let mut r2_buf_reusable = Vec::with_capacity(800);

        // Pre-create modified trimming configs without adapters (avoids clone per read)
        let trimming_config_no_adapter_r1 = {
            let mut cfg = trimming_config_r1.clone();
            cfg.adapter_config.adapter_seq = None;
            cfg
        };
        let trimming_config_no_adapter_r2 = {
            let mut cfg = trimming_config_r2.clone();
            cfg.adapter_config.adapter_seq = None;
            cfg.adapter_config.adapter_seq_r2 = None;
            cfg
        };

        // Process each pair
        for (idx, (&rec1, &rec2)) in batch.recs_r1.iter().zip(batch.recs_r2.iter()).enumerate() {
            // Parse R1
            let (seq1, qual1, _seq1_end, _qual1_end) =
                parse_record(&batch.buf_r1, rec1, idx, &batch.recs_r1);

            // Parse R2
            let (seq2, qual2, _seq2_end, _qual2_end) =
                parse_record(&batch.buf_r2, rec2, idx, &batch.recs_r2);

            // Validate both records
            if seq1.len() != qual1.len() || seq2.len() != qual2.len() {
                invalid += 1;
                continue;
            }

            // Track duplicates using Bloom filter (before any processing)
            // This is for statistics only, not filtering
            if let Some(ref detector) = dup_detector {
                detector.check_pair(seq1, seq2);
            }

            // Extract UMI if enabled (happens BEFORE any other processing)
            let umi_seq: Vec<u8>;
            let (final_umi_seq1, final_umi_qual1, final_umi_seq2, final_umi_qual2);
            let umi_seq1_buf;
            let umi_qual1_buf;
            let umi_seq2_buf;
            let umi_qual2_buf;
            let umi_extracted;

            if let Some(ref cfg) = umi_config {
                match cfg.location {
                    crate::umi::UmiLocation::Read1 => {
                        if let Some(extraction) = crate::umi::extract_umi(seq1, cfg) {
                            // UMI extraction succeeded
                            umi_seq = extraction.umi_seq;
                            umi_seq1_buf = seq1[extraction.end_pos..].to_vec();
                            umi_qual1_buf = qual1[extraction.end_pos..].to_vec();
                            umi_seq2_buf = seq2.to_vec();
                            umi_qual2_buf = qual2.to_vec();

                            final_umi_seq1 = &umi_seq1_buf[..];
                            final_umi_qual1 = &umi_qual1_buf[..];
                            final_umi_seq2 = &umi_seq2_buf[..];
                            final_umi_qual2 = &umi_qual2_buf[..];
                            umi_extracted = true;
                        } else {
                            // UMI extraction failed - use original
                            umi_seq = Vec::new();
                            final_umi_seq1 = seq1;
                            final_umi_qual1 = qual1;
                            final_umi_seq2 = seq2;
                            final_umi_qual2 = qual2;
                            umi_extracted = false;
                        }
                    }
                    crate::umi::UmiLocation::Read2 => {
                        if let Some(extraction) = crate::umi::extract_umi(seq2, cfg) {
                            // UMI extraction succeeded
                            umi_seq = extraction.umi_seq;
                            umi_seq1_buf = seq1.to_vec();
                            umi_qual1_buf = qual1.to_vec();
                            umi_seq2_buf = seq2[extraction.end_pos..].to_vec();
                            umi_qual2_buf = qual2[extraction.end_pos..].to_vec();

                            final_umi_seq1 = &umi_seq1_buf[..];
                            final_umi_qual1 = &umi_qual1_buf[..];
                            final_umi_seq2 = &umi_seq2_buf[..];
                            final_umi_qual2 = &umi_qual2_buf[..];
                            umi_extracted = true;
                        } else {
                            // UMI extraction failed - use original
                            umi_seq = Vec::new();
                            final_umi_seq1 = seq1;
                            final_umi_qual1 = qual1;
                            final_umi_seq2 = seq2;
                            final_umi_qual2 = qual2;
                            umi_extracted = false;
                        }
                    }
                    _ => {
                        // Other UMI locations not yet supported
                        umi_seq = Vec::new();
                        final_umi_seq1 = seq1;
                        final_umi_qual1 = qual1;
                        final_umi_seq2 = seq2;
                        final_umi_qual2 = qual2;
                        umi_extracted = false;
                    }
                }
            } else {
                // No UMI processing
                umi_seq = Vec::new();
                final_umi_seq1 = seq1;
                final_umi_qual1 = qual1;
                final_umi_seq2 = seq2;
                final_umi_qual2 = qual2;
                umi_extracted = false;
            }

            // Compute stats for R1 (after UMI removal)
            let stats1 = simd::compute_stats(final_umi_seq1, final_umi_qual1, 0);
            before_r1.add(
                final_umi_seq1.len(),
                stats1.q20,
                stats1.q30,
                stats1.q40,
                stats1.gc,
            );
            update_position_stats(&mut pos_r1, final_umi_seq1, final_umi_qual1);
            if !no_kmer {
                count_k5_2bit(final_umi_seq1, &mut k5_r1);
            }

            // Compute stats for R2 (after UMI removal)
            let stats2 = simd::compute_stats(final_umi_seq2, final_umi_qual2, 0);
            before_r2.add(
                final_umi_seq2.len(),
                stats2.q20,
                stats2.q30,
                stats2.q40,
                stats2.gc,
            );
            update_position_stats(&mut pos_r2, final_umi_seq2, final_umi_qual2);
            if !no_kmer {
                count_k5_2bit(final_umi_seq2, &mut k5_r2);
            }

            // Always detect overlap for paired-end adapter trimming (fastp does this)
            // Overlap-based trimming works even when adapters aren't detected via k-mer analysis
            let overlap_result = {
                let overlap_detect_config = crate::overlap::OverlapConfig {
                    min_overlap_len: 30, // fastp default (options.cpp:overlapRequire = 30)
                    max_diff: 5,
                    max_diff_percent: 20,
                };
                crate::overlap::detect_overlap(
                    final_umi_seq1,
                    final_umi_seq2,
                    &overlap_detect_config,
                )
            };

            // Try overlap-based adapter trimming first
            // Always use overlap-based trimming when overlap is detected (fastp does this)
            let overlap_trimmed = if let Some(ref overlap) = overlap_result {
                crate::overlap::trim_by_overlap_analysis(
                    final_umi_seq1.len(),
                    final_umi_seq2.len(),
                    overlap,
                )
            } else {
                None
            };

            // Apply trimming to R1 and R2 (after UMI removal)
            let mut trim_result1;
            let mut trim_result2;

            if let Some((overlap_trim_len1, overlap_trim_len2)) = overlap_trimmed {
                // Overlap-based adapter trimming succeeded
                // Still need to apply other trimming (poly-G, quality, etc.)
                // but skip adapter trimming since we already did it with overlap

                // Use pre-created configs without adapters (avoids clone per read)
                trim_result1 = if trimming_config_r1.is_enabled() {
                    trim_read(
                        final_umi_seq1,
                        final_umi_qual1,
                        &trimming_config_no_adapter_r1,
                    )
                } else {
                    TrimmingResult {
                        start_pos: 0,
                        end_pos: final_umi_seq1.len(),
                        poly_g_trimmed: 0,
                        poly_x_trimmed: 0,
                        adapter_trimmed: false,
                        adapter_bases_trimmed: 0,
                    }
                };

                trim_result2 = if trimming_config_r2.is_enabled() {
                    trim_read(
                        final_umi_seq2,
                        final_umi_qual2,
                        &trimming_config_no_adapter_r2,
                    )
                } else {
                    TrimmingResult {
                        start_pos: 0,
                        end_pos: final_umi_seq2.len(),
                        poly_g_trimmed: 0,
                        poly_x_trimmed: 0,
                        adapter_trimmed: false,
                        adapter_bases_trimmed: 0,
                    }
                };

                // Apply overlap-based adapter trim lengths
                // Track the end position before overlap trim to calculate adapter bases correctly
                let end_before_overlap1 = trim_result1.end_pos;
                let end_before_overlap2 = trim_result2.end_pos;

                trim_result1.end_pos = std::cmp::min(trim_result1.end_pos, overlap_trim_len1);
                trim_result2.end_pos = std::cmp::min(trim_result2.end_pos, overlap_trim_len2);

                // Mark as adapter trimmed (only count bases trimmed beyond other trims)
                if trim_result1.end_pos < end_before_overlap1 {
                    trim_result1.adapter_trimmed = true;
                    trim_result1.adapter_bases_trimmed = end_before_overlap1 - trim_result1.end_pos;
                }
                if trim_result2.end_pos < end_before_overlap2 {
                    trim_result2.adapter_trimmed = true;
                    trim_result2.adapter_bases_trimmed = end_before_overlap2 - trim_result2.end_pos;
                }
            } else {
                // No overlap-based trimming - use normal adapter trimming
                trim_result1 = if trimming_config_r1.is_enabled() {
                    trim_read(final_umi_seq1, final_umi_qual1, &trimming_config_r1)
                } else {
                    TrimmingResult {
                        start_pos: 0,
                        end_pos: final_umi_seq1.len(),
                        poly_g_trimmed: 0,
                        poly_x_trimmed: 0,
                        adapter_trimmed: false,
                        adapter_bases_trimmed: 0,
                    }
                };

                trim_result2 = if trimming_config_r2.is_enabled() {
                    trim_read(final_umi_seq2, final_umi_qual2, &trimming_config_r2)
                } else {
                    TrimmingResult {
                        start_pos: 0,
                        end_pos: final_umi_seq2.len(),
                        poly_g_trimmed: 0,
                        poly_x_trimmed: 0,
                        adapter_trimmed: false,
                        adapter_bases_trimmed: 0,
                    }
                };
            }

            // Track adapter trimming stats (count each read individually)
            if trim_result1.adapter_trimmed {
                adapter_trimmed_reads += 1;
                adapter_trimmed_bases += trim_result1.adapter_bases_trimmed;
            }
            if trim_result2.adapter_trimmed {
                adapter_trimmed_reads += 1;
                adapter_trimmed_bases += trim_result2.adapter_bases_trimmed;
            }

            // Get trimmed sequences (from UMI-processed sequences)
            let trimmed_seq1 = &final_umi_seq1[trim_result1.start_pos..trim_result1.end_pos];
            let trimmed_qual1 = &final_umi_qual1[trim_result1.start_pos..trim_result1.end_pos];
            let trimmed_seq2 = &final_umi_seq2[trim_result2.start_pos..trim_result2.end_pos];
            let trimmed_qual2 = &final_umi_qual2[trim_result2.start_pos..trim_result2.end_pos];

            // Apply base correction using overlap analysis (if enabled)
            // Reuse the overlap_result from adapter trimming detection to avoid redundant work
            let (final_seq1, final_qual1, final_seq2, final_qual2);
            let mut seq1_corrected;
            let mut qual1_corrected;
            let mut seq2_corrected;
            let mut qual2_corrected;

            // Only create mutable copies if correction is enabled AND there's an overlap to correct
            if overlap_config.is_some() && overlap_result.is_some() {
                // Create mutable copies for correction
                seq1_corrected = trimmed_seq1.to_vec();
                qual1_corrected = trimmed_qual1.to_vec();
                seq2_corrected = trimmed_seq2.to_vec();
                qual2_corrected = trimmed_qual2.to_vec();

                // Apply correction using cached overlap result
                let _correction_stats = crate::overlap::correct_by_overlap(
                    &mut seq1_corrected,
                    &mut qual1_corrected,
                    &mut seq2_corrected,
                    &mut qual2_corrected,
                    overlap_result.as_ref().unwrap(),
                );
                // TODO: Track correction statistics

                final_seq1 = &seq1_corrected[..];
                final_qual1 = &qual1_corrected[..];
                final_seq2 = &seq2_corrected[..];
                final_qual2 = &qual2_corrected[..];
            } else {
                // No correction needed - use trimmed sequences directly
                final_seq1 = trimmed_seq1;
                final_qual1 = trimmed_qual1;
                final_seq2 = trimmed_seq2;
                final_qual2 = trimmed_qual2;
            }

            // Track insert size using cached overlap result (avoids third detect_overlap call)
            if let Some(ref overlap) = overlap_result {
                if overlap.overlapped {
                    // Adjust for any trimming that occurred
                    let insert_size = final_seq1.len() + final_seq2.len()
                        - overlap.overlap_len.min(final_seq1.len() + final_seq2.len());
                    if insert_size < insert_size_histogram.len() {
                        insert_size_histogram[insert_size] += 1;
                    }
                } else {
                    insert_size_unknown += 1;
                }
            } else {
                insert_size_unknown += 1;
            }

            // Check length filter (EITHER read too short = BOTH fail)
            if final_seq1.len() < min_len || final_seq2.len() < min_len {
                too_short += 1;
                continue;
            }

            // Recompute stats for final reads (after trimming and correction)
            let trimmed_stats1 =
                simd::compute_stats(final_seq1, final_qual1, qualified_quality_phred);
            let trimmed_stats2 =
                simd::compute_stats(final_seq2, final_qual2, qualified_quality_phred);

            // Check unqualified percent FIRST to match fastp order (EITHER read fails = BOTH fail)
            let mut fail_quality = false;
            if qualified_quality_phred > 0 {
                if !final_seq1.is_empty()
                    && 100 * trimmed_stats1.unqualified
                        > unqualified_percent_limit * final_seq1.len()
                {
                    fail_quality = true;
                }
                if !final_seq2.is_empty()
                    && 100 * trimmed_stats2.unqualified
                        > unqualified_percent_limit * final_seq2.len()
                {
                    fail_quality = true;
                }
            }

            if fail_quality {
                low_quality += 1;
                continue;
            }

            // Check N limit SECOND (EITHER read too many N = BOTH fail)
            if trimmed_stats1.ncnt > n_limit || trimmed_stats2.ncnt > n_limit {
                too_many_n += 1;
                continue;
            }

            // Check average quality (EITHER read fails = BOTH fail)
            // Use integer comparison: mean >= threshold  =>  qsum >= threshold * len
            if average_qual > 0 {
                let threshold = u64::from(average_qual);
                if !final_seq1.is_empty()
                    && u64::from(trimmed_stats1.qsum) < threshold * final_seq1.len() as u64
                {
                    fail_quality = true;
                }
                if !final_seq2.is_empty()
                    && u64::from(trimmed_stats2.qsum) < threshold * final_seq2.len() as u64
                {
                    fail_quality = true;
                }
            }

            if fail_quality {
                low_quality += 1;
                continue;
            }

            // Check low complexity (EITHER read fails = BOTH fail)
            if low_complexity_filter {
                let mut fail_complexity = false;
                if !final_seq1.is_empty() {
                    let complexity1 = calculate_complexity(final_seq1);
                    if complexity1 < complexity_threshold {
                        fail_complexity = true;
                    }
                }
                if !final_seq2.is_empty() {
                    let complexity2 = calculate_complexity(final_seq2);
                    if complexity2 < complexity_threshold {
                        fail_complexity = true;
                    }
                }

                if fail_complexity {
                    low_complexity += 1;
                    continue;
                }
            }

            // Check for duplicates if deduplication is enabled
            if let Some(ref mut tracker) = dedup_tracker {
                if tracker.is_duplicate_pe(final_seq1, final_seq2) {
                    // This is a duplicate - skip it
                    duplicated += 1;
                    continue;
                }
            }

            // BOTH reads passed - emit ranges for BOTH
            let [h1_start, s1_start, p1_start, q1_start] = rec1;
            let [h2_start, s2_start, p2_start, q2_start] = rec2;

            // When correction or UMI is enabled, we need to create new buffers with modified data
            // Otherwise, use zero-copy with ranges into original buffer
            let (r1_piece, r2_piece) = if overlap_config.is_some() || umi_config.is_some() {
                // Correction or UMI was enabled - reuse pre-allocated buffers
                r1_buf_reusable.clear();
                r2_buf_reusable.clear();
                let r1_buf = &mut r1_buf_reusable;
                let r2_buf = &mut r2_buf_reusable;

                // Build R1 buffer: header + '\n' + seq + '\n' + plus + '\n' + qual + '\n'
                let r1_header_start = 0;
                // Add original header
                r1_buf.extend_from_slice(&batch.buf_r1[h1_start..s1_start - 1]);
                // Add UMI to header if extracted
                if umi_extracted {
                    if let Some(ref cfg) = umi_config {
                        r1_buf.push(b':');
                        r1_buf.extend_from_slice(cfg.prefix.as_bytes());
                        r1_buf.push(b'_');
                        r1_buf.extend_from_slice(&umi_seq);
                    }
                }
                let r1_header_end = r1_buf.len();
                r1_buf.push(b'\n');
                let r1_seq_start = r1_buf.len();
                r1_buf.extend_from_slice(final_seq1);
                let r1_seq_end = r1_buf.len();
                r1_buf.push(b'\n');
                let r1_plus_start = r1_buf.len();
                r1_buf.extend_from_slice(&batch.buf_r1[p1_start..q1_start - 1]);
                let r1_plus_end = r1_buf.len();
                r1_buf.push(b'\n');
                let r1_qual_start = r1_buf.len();
                r1_buf.extend_from_slice(final_qual1);
                let r1_qual_end = r1_buf.len();
                r1_buf.push(b'\n');

                // Build R2 buffer
                let r2_header_start = 0;
                // Add original header
                r2_buf.extend_from_slice(&batch.buf_r2[h2_start..s2_start - 1]);
                // Add UMI to header if extracted
                if umi_extracted {
                    if let Some(ref cfg) = umi_config {
                        r2_buf.push(b':');
                        r2_buf.extend_from_slice(cfg.prefix.as_bytes());
                        r2_buf.push(b'_');
                        r2_buf.extend_from_slice(&umi_seq);
                    }
                }
                let r2_header_end = r2_buf.len();
                r2_buf.push(b'\n');
                let r2_seq_start = r2_buf.len();
                r2_buf.extend_from_slice(final_seq2);
                let r2_seq_end = r2_buf.len();
                r2_buf.push(b'\n');
                let r2_plus_start = r2_buf.len();
                r2_buf.extend_from_slice(&batch.buf_r2[p2_start..q2_start - 1]);
                let r2_plus_end = r2_buf.len();
                r2_buf.push(b'\n');
                let r2_qual_start = r2_buf.len();
                r2_buf.extend_from_slice(final_qual2);
                let r2_qual_end = r2_buf.len();
                r2_buf.push(b'\n');

                let r1_arc = Arc::new(r1_buf.clone());
                let r2_arc = Arc::new(r2_buf.clone());

                (
                    RecordPiece {
                        buf: r1_arc,
                        header: r1_header_start..r1_header_end,
                        seq: r1_seq_start..r1_seq_end,
                        plus: r1_plus_start..r1_plus_end,
                        qual: r1_qual_start..r1_qual_end,
                    },
                    RecordPiece {
                        buf: r2_arc,
                        header: r2_header_start..r2_header_end,
                        seq: r2_seq_start..r2_seq_end,
                        plus: r2_plus_start..r2_plus_end,
                        qual: r2_qual_start..r2_qual_end,
                    },
                )
            } else {
                // No correction - use zero-copy with ranges into original buffer
                let trimmed_seq1_start = s1_start + trim_result1.start_pos;
                let trimmed_seq1_end = s1_start + trim_result1.end_pos;
                let trimmed_qual1_start = q1_start + trim_result1.start_pos;
                let trimmed_qual1_end = q1_start + trim_result1.end_pos;

                let trimmed_seq2_start = s2_start + trim_result2.start_pos;
                let trimmed_seq2_end = s2_start + trim_result2.end_pos;
                let trimmed_qual2_start = q2_start + trim_result2.start_pos;
                let trimmed_qual2_end = q2_start + trim_result2.end_pos;

                (
                    RecordPiece {
                        buf: batch.buf_r1.clone(),
                        header: h1_start..s1_start - 1,
                        seq: trimmed_seq1_start..trimmed_seq1_end,
                        plus: p1_start..q1_start - 1,
                        qual: trimmed_qual1_start..trimmed_qual1_end,
                    },
                    RecordPiece {
                        buf: batch.buf_r2.clone(),
                        header: h2_start..s2_start - 1,
                        seq: trimmed_seq2_start..trimmed_seq2_end,
                        plus: p2_start..q2_start - 1,
                        qual: trimmed_qual2_start..trimmed_qual2_end,
                    },
                )
            };

            pieces.push(PairedRecordPiece {
                r1: r1_piece,
                r2: r2_piece,
            });

            // Update "after" stats for BOTH reads
            after_r1.add(
                final_seq1.len(),
                trimmed_stats1.q20,
                trimmed_stats1.q30,
                trimmed_stats1.q40,
                trimmed_stats1.gc,
            );
            after_r2.add(
                final_seq2.len(),
                trimmed_stats2.q20,
                trimmed_stats2.q30,
                trimmed_stats2.q40,
                trimmed_stats2.gc,
            );

            // Track position stats for reads that passed filters
            crate::processor::track_position_stats(final_seq1, final_qual1, &mut pos_r1_after);
            crate::processor::track_position_stats(final_seq2, final_qual2, &mut pos_r2_after);

            // Track kmers for reads that passed filters
            count_k5_2bit(final_seq1, &mut k5_r1_after);
            count_k5_2bit(final_seq2, &mut k5_r2_after);
        }

        let result = PairedWorkerResult {
            id: batch.id,
            pieces,
            before_r1,
            after_r1,
            pos_r1,
            pos_r1_after,
            k5_r1,
            k5_r1_after,
            before_r2,
            after_r2,
            pos_r2,
            pos_r2_after,
            k5_r2,
            k5_r2_after,
            too_short,
            too_many_n,
            low_quality,
            low_complexity,
            invalid,
            duplicated,
            adapter_trimmed_reads,
            adapter_trimmed_bases,
            insert_size_histogram,
            insert_size_unknown,
        };

        if sender.send(Some(result)).is_err() {
            break;
        }
    }

    // Send sentinel
    let _ = sender.send(None);
}

/// Helper to parse a single record from buffer
#[inline]
fn parse_record<'a>(
    buf: &'a [u8],
    rec: [usize; 4],
    idx: usize,
    recs: &[[usize; 4]],
) -> (&'a [u8], &'a [u8], usize, usize) {
    let [_h_start, s_start, p_start, q_start] = rec;
    let s_end = p_start - 1;
    let q_end = if idx + 1 < recs.len() {
        recs[idx + 1][0] - 1
    } else {
        // For the last record, exclude trailing newline if present
        let buf_len = buf.len();
        if buf_len > 0 && buf[buf_len - 1] == b'\n' {
            buf_len - 1
        } else {
            buf_len
        }
    };

    let seq = &buf[s_start..s_end];
    let qual = &buf[q_start..q_end];

    (seq, qual, s_end, q_end)
}

/// Bases for validation: indices match (b >> 1) & 3
static BASES: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Helper to update position stats
#[inline]
fn update_position_stats(pos: &mut PositionStats, seq: &[u8], qual: &[u8]) {
    pos.ensure_capacity(seq.len());

    let len = seq.len();

    // SAFETY: ensure_capacity guarantees all arrays have at least `len` elements
    unsafe {
        for i in 0..len {
            let q = *qual.get_unchecked(i);
            let quality_score = u64::from(q.wrapping_sub(33));
            let b = *seq.get_unchecked(i);

            *pos.total_sum.get_unchecked_mut(i) += quality_score;
            *pos.total_cnt.get_unchecked_mut(i) += 1;

            // Bit trick: (b >> 1) & 3 maps A->0, C->1, T->2, G->3
            // For N (0x4E) this maps to 3 (G bucket) - acceptable since N is rare
            let bi = ((b >> 1) & 3) as usize;
            *pos.base_sum.get_unchecked_mut(bi).get_unchecked_mut(i) += quality_score;
            *pos.base_cnt.get_unchecked_mut(bi).get_unchecked_mut(i) += 1;

            // Track quality histogram (clamp to valid range)
            let hist_idx = (quality_score as usize).min(93);
            *pos.qual_hist.get_unchecked_mut(hist_idx) += 1;
        }
    }
}

/// Paired-end merger thread: write both outputs in order
#[allow(clippy::too_many_lines)]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn paired_merger_thread(
    receiver: Receiver<Option<PairedWorkerResult>>,
    output1_path: String,
    output2_path: String,
    num_workers: usize,
    compression_level: Option<u32>,
    split_config: crate::split::SplitConfig,
    parallel_compression: bool,
) -> Result<crate::processor::PairedEndAccumulator> {
    use std::io::Write;

    // Create split writers
    let compression = compression_level.unwrap_or(6);
    let mut writer1 = crate::split::SplitWriter::new(
        &output1_path,
        split_config.clone(),
        compression,
        parallel_compression,
    )?;
    let mut writer2 = crate::split::SplitWriter::new(
        &output2_path,
        split_config,
        compression,
        parallel_compression,
    )?;

    let mut acc = crate::processor::PairedEndAccumulator::new();
    let mut next_id = 0u64;
    let mut pending: BTreeMap<u64, PairedWorkerResult> = BTreeMap::new();
    let mut workers_done = 0;

    // Reusable buffers for pre-concatenating records (1 memcpy instead of 8)
    let mut r1_buf = Vec::with_capacity(1024);
    let mut r2_buf = Vec::with_capacity(1024);

    while let Ok(msg) = receiver.recv() {
        if let Some(result) = msg {
            pending.insert(result.id, result);

            // Process all consecutive batches
            while let Some(result) = pending.remove(&next_id) {
                // Write all passed records
                for piece in &result.pieces {
                    // Pre-concatenate R1 into contiguous buffer (1 memcpy instead of 8)
                    let b1 = &piece.r1.buf;
                    r1_buf.clear();
                    r1_buf.extend_from_slice(&b1[piece.r1.header.clone()]);
                    r1_buf.push(b'\n');
                    r1_buf.extend_from_slice(&b1[piece.r1.seq.clone()]);
                    r1_buf.push(b'\n');
                    r1_buf.extend_from_slice(&b1[piece.r1.plus.clone()]);
                    r1_buf.push(b'\n');
                    r1_buf.extend_from_slice(&b1[piece.r1.qual.clone()]);
                    r1_buf.push(b'\n');
                    writer1.write_all(&r1_buf)?;

                    // Pre-concatenate R2 into contiguous buffer
                    let b2 = &piece.r2.buf;
                    r2_buf.clear();
                    r2_buf.extend_from_slice(&b2[piece.r2.header.clone()]);
                    r2_buf.push(b'\n');
                    r2_buf.extend_from_slice(&b2[piece.r2.seq.clone()]);
                    r2_buf.push(b'\n');
                    r2_buf.extend_from_slice(&b2[piece.r2.plus.clone()]);
                    r2_buf.push(b'\n');
                    r2_buf.extend_from_slice(&b2[piece.r2.qual.clone()]);
                    r2_buf.push(b'\n');
                    writer2.write_all(&r2_buf)?;
                }

                // Merge stats - R1
                acc.before_r1.total_reads += result.before_r1.total_reads;
                acc.before_r1.total_bases += result.before_r1.total_bases;
                acc.before_r1.q20_bases += result.before_r1.q20_bases;
                acc.before_r1.q30_bases += result.before_r1.q30_bases;
                acc.before_r1.q40_bases += result.before_r1.q40_bases;
                acc.before_r1.gc_bases += result.before_r1.gc_bases;

                acc.after_r1.total_reads += result.after_r1.total_reads;
                acc.after_r1.total_bases += result.after_r1.total_bases;
                acc.after_r1.q20_bases += result.after_r1.q20_bases;
                acc.after_r1.q30_bases += result.after_r1.q30_bases;
                acc.after_r1.q40_bases += result.after_r1.q40_bases;
                acc.after_r1.gc_bases += result.after_r1.gc_bases;

                // Merge stats - R2
                acc.before_r2.total_reads += result.before_r2.total_reads;
                acc.before_r2.total_bases += result.before_r2.total_bases;
                acc.before_r2.q20_bases += result.before_r2.q20_bases;
                acc.before_r2.q30_bases += result.before_r2.q30_bases;
                acc.before_r2.q40_bases += result.before_r2.q40_bases;
                acc.before_r2.gc_bases += result.before_r2.gc_bases;

                acc.after_r2.total_reads += result.after_r2.total_reads;
                acc.after_r2.total_bases += result.after_r2.total_bases;
                acc.after_r2.q20_bases += result.after_r2.q20_bases;
                acc.after_r2.q30_bases += result.after_r2.q30_bases;
                acc.after_r2.q40_bases += result.after_r2.q40_bases;
                acc.after_r2.gc_bases += result.after_r2.gc_bases;

                // Merge position stats - R1
                acc.pos_r1.ensure_capacity(result.pos_r1.total_sum.len());
                for i in 0..result.pos_r1.total_sum.len() {
                    acc.pos_r1.total_sum[i] += result.pos_r1.total_sum[i];
                    acc.pos_r1.total_cnt[i] += result.pos_r1.total_cnt[i];
                    for b in 0..4 {
                        acc.pos_r1.base_sum[b][i] += result.pos_r1.base_sum[b][i];
                        acc.pos_r1.base_cnt[b][i] += result.pos_r1.base_cnt[b][i];
                    }
                }
                // Merge quality histogram - R1 before filtering
                for (i, &count) in result.pos_r1.qual_hist.iter().enumerate() {
                    acc.pos_r1.qual_hist[i] += count;
                }

                // Merge position stats - R2
                acc.pos_r2.ensure_capacity(result.pos_r2.total_sum.len());
                for i in 0..result.pos_r2.total_sum.len() {
                    acc.pos_r2.total_sum[i] += result.pos_r2.total_sum[i];
                    acc.pos_r2.total_cnt[i] += result.pos_r2.total_cnt[i];
                    for b in 0..4 {
                        acc.pos_r2.base_sum[b][i] += result.pos_r2.base_sum[b][i];
                        acc.pos_r2.base_cnt[b][i] += result.pos_r2.base_cnt[b][i];
                    }
                }
                // Merge quality histogram - R2 before filtering
                for (i, &count) in result.pos_r2.qual_hist.iter().enumerate() {
                    acc.pos_r2.qual_hist[i] += count;
                }

                // Merge position stats after filtering - R1
                acc.pos_r1_after
                    .ensure_capacity(result.pos_r1_after.total_sum.len());
                for i in 0..result.pos_r1_after.total_sum.len() {
                    acc.pos_r1_after.total_sum[i] += result.pos_r1_after.total_sum[i];
                    acc.pos_r1_after.total_cnt[i] += result.pos_r1_after.total_cnt[i];
                    for b in 0..4 {
                        acc.pos_r1_after.base_sum[b][i] += result.pos_r1_after.base_sum[b][i];
                        acc.pos_r1_after.base_cnt[b][i] += result.pos_r1_after.base_cnt[b][i];
                    }
                }

                // Merge position stats after filtering - R2
                acc.pos_r2_after
                    .ensure_capacity(result.pos_r2_after.total_sum.len());
                for i in 0..result.pos_r2_after.total_sum.len() {
                    acc.pos_r2_after.total_sum[i] += result.pos_r2_after.total_sum[i];
                    acc.pos_r2_after.total_cnt[i] += result.pos_r2_after.total_cnt[i];
                    for b in 0..4 {
                        acc.pos_r2_after.base_sum[b][i] += result.pos_r2_after.base_sum[b][i];
                        acc.pos_r2_after.base_cnt[b][i] += result.pos_r2_after.base_cnt[b][i];
                    }
                }

                // Merge quality histograms for after filtering
                for (i, &count) in result.pos_r1_after.qual_hist.iter().enumerate() {
                    acc.pos_r1_after.qual_hist[i] += count;
                }
                for (i, &count) in result.pos_r2_after.qual_hist.iter().enumerate() {
                    acc.pos_r2_after.qual_hist[i] += count;
                }

                for (i, &count) in result.k5_r1.iter().enumerate() {
                    acc.kmer_table_r1[i] += count;
                }
                for (i, &count) in result.k5_r2.iter().enumerate() {
                    acc.kmer_table_r2[i] += count;
                }

                // Merge kmer tables for after filtering
                for (i, &count) in result.k5_r1_after.iter().enumerate() {
                    acc.kmer_table_r1_after[i] += count;
                }
                for (i, &count) in result.k5_r2_after.iter().enumerate() {
                    acc.kmer_table_r2_after[i] += count;
                }

                acc.too_short += result.too_short;
                acc.too_many_n += result.too_many_n;
                acc.low_quality += result.low_quality;
                acc.low_complexity += result.low_complexity;
                acc.invalid += result.invalid;
                acc.duplicated += result.duplicated;

                // Merge adapter trimming stats
                acc.adapter_trimmed_reads += result.adapter_trimmed_reads;
                acc.adapter_trimmed_bases += result.adapter_trimmed_bases;

                // Merge insert size histograms
                for (i, &count) in result.insert_size_histogram.iter().enumerate() {
                    if i < acc.insert_size_histogram.len() {
                        acc.insert_size_histogram[i] += count;
                    }
                }
                acc.insert_size_unknown += result.insert_size_unknown;

                if result.pos_r1.total_sum.len() > acc.max_cycle_r1 {
                    acc.max_cycle_r1 = result.pos_r1.total_sum.len();
                }
                if result.pos_r2.total_sum.len() > acc.max_cycle_r2 {
                    acc.max_cycle_r2 = result.pos_r2.total_sum.len();
                }

                next_id += 1;
            }
        } else {
            workers_done += 1;
            if workers_done == num_workers {
                break;
            }
        }
    }

    writer1.finish()?;
    writer2.finish()?;

    Ok(acc)
}
