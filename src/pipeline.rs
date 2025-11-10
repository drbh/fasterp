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
use std::io::{BufWriter, IoSlice, Read, Write};
use std::ops::Range;
use std::sync::Arc;

use crate::io::open_output;
use crate::kmer::*;
use crate::processor::StreamAccumulator;
use crate::simd;
use crate::stats::*;
use crate::trimming::*;

/// Calculate sequence complexity as percentage of bases different from next base
///
/// Complexity is defined as: (count of positions where base[i] != base[i+1]) / (length - 1) * 100
/// This matches fastp's low complexity filter algorithm.
///
/// Returns complexity percentage (0-100)
#[inline]
fn calculate_complexity(seq: &[u8]) -> usize {
    if seq.len() <= 1 {
        return 100; // Single base or empty is considered max complexity
    }

    let mut different_count = 0;
    for i in 0..seq.len() - 1 {
        if seq[i] != seq[i + 1] {
            different_count += 1;
        }
    }

    // Calculate percentage: (different_count / (len - 1)) * 100
    // Use integer math to avoid floating point
    (different_count * 100) / (seq.len() - 1)
}

/// A batch of FASTQ records parsed from a buffer
///
/// Contains raw bytes and record positions (no String allocations)
/// Each record is [header_start, seq_start, plus_start, qual_start] as byte offsets
///
/// Uses Arc to allow zero-copy sharing of buffer between workers
#[derive(Clone)]
pub(crate) struct Batch {
    pub id: u64,
    pub buf: Arc<Vec<u8>>, // Shared buffer - no copying needed
    /// Each element is [header_start, seq_start, plus_start, qual_start]
    /// Lengths are implicit: header len = seq_start - header_start, etc.
    /// Quality ends at the next record's header_start (or buf.len() for last record)
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
/// Uses RecordPiece to avoid copying bytes - stores ranges instead
pub(crate) struct WorkerResult {
    pub id: u64,
    pub pieces: Vec<RecordPiece>, // Zero-copy: ranges into shared buffers
    pub before: SimpleStats,
    pub after: SimpleStats,
    pub pos: PositionStats,
    pub k5: [usize; 1024],
    // Filter counts
    pub too_short: usize,
    pub too_many_n: usize,
    pub low_quality: usize,
    pub low_complexity: usize,
    pub invalid: usize,
}

/// Producer thread: read blocks and parse into batches
///
/// Reads large blocks (batch_bytes) from input, parses FASTQ records,
/// and emits Batch structures with NO string allocations - just byte slices.
///
/// Handles partial records at block boundaries by carrying them over to next batch.
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

        // Single-pass scan: find complete records AND line starts
        let mut complete_end = 0;
        let mut line_count = 0;
        let mut line_starts = vec![0];

        for i in 0..actual_len {
            if buffer[i] == b'\n' {
                line_count += 1;
                // After every 4 lines, we have a complete record
                if line_count % 4 == 0 {
                    complete_end = i + 1;
                }
                // Track line starts for parsing
                if i + 1 < actual_len {
                    line_starts.push(i + 1);
                }
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
        let mut k5 = [0usize; 1024];
        let mut too_short = 0usize;
        let mut too_many_n = 0usize;
        let mut low_quality = 0usize;
        let mut low_complexity = 0usize;
        let mut invalid = 0usize;
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
            let (_qsum, q20, q30, _ncnt, gc) = if simd::is_simd_available() {
                // SIMD path: compute basic stats fast, then position-specific
                let stats = simd::compute_stats(seq, qual, 0); // 0 = don't count unqualified for before stats

                pos.ensure_capacity(seq.len());
                for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
                    pos.total_sum[i] += (q - 33) as u64;
                    pos.total_cnt[i] += 1;

                    if let Some(bi) = base_idx(b) {
                        pos.base_sum[bi][i] += (q - 33) as u64;
                        pos.base_cnt[bi][i] += 1;
                    }
                }

                (stats.qsum, stats.q20, stats.q30, stats.ncnt, stats.gc)
            } else {
                // Non-SIMD path: single pass for both basic and position stats
                let mut qsum = 0u32;
                let mut q20 = 0usize;
                let mut q30 = 0usize;
                let mut ncnt = 0usize;
                let mut gc = 0usize;

                pos.ensure_capacity(seq.len());
                for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
                    let qval = (q - 33) as u32;
                    qsum += qval;
                    if q >= 53 {
                        q20 += 1;
                    }
                    if q >= 63 {
                        q30 += 1;
                    }
                    if b == b'N' || b == b'n' {
                        ncnt += 1;
                    }
                    if b == b'G' || b == b'g' || b == b'C' || b == b'c' {
                        gc += 1;
                    }

                    pos.total_sum[i] += qval as u64;
                    pos.total_cnt[i] += 1;

                    if let Some(bi) = base_idx(b) {
                        pos.base_sum[bi][i] += qval as u64;
                        pos.base_cnt[bi][i] += 1;
                    }
                }

                (qsum, q20, q30, ncnt, gc)
            };

            // K-mer counting
            if !no_kmer {
                count_k5_2bit(seq, &mut k5);
            }

            // Update before stats
            before.add(seq.len(), q20, q30, gc);

            // Apply trimming if enabled
            let trimming_result = if trimming_config.is_enabled() {
                trim_read(seq, qual, &trimming_config)
            } else {
                TrimmingResult {
                    start_pos: 0,
                    end_pos: seq.len(),
                    poly_g_trimmed: 0,
                    poly_x_trimmed: 0,
                }
            };

            // Get trimmed sequences
            let trimmed_seq = &seq[trimming_result.start_pos..trimming_result.end_pos];
            let trimmed_qual = &qual[trimming_result.start_pos..trimming_result.end_pos];

            // Early length check - skip expensive stats computation for reads that are too short
            let trimmed_len = trimmed_seq.len();
            if trimmed_len < min_len {
                too_short += 1;
                continue;
            }

            // Recompute stats for trimmed read (used for filtering) - SIMD accelerated
            let trimmed_stats =
                simd::compute_stats(trimmed_seq, trimmed_qual, qualified_quality_phred);
            let trimmed_qsum = trimmed_stats.qsum;
            let trimmed_q20 = trimmed_stats.q20;
            let trimmed_q30 = trimmed_stats.q30;
            let trimmed_ncnt = trimmed_stats.ncnt;
            let trimmed_gc = trimmed_stats.gc;
            let unqualified_count = trimmed_stats.unqualified;

            // Apply remaining filters on TRIMMED read

            if trimmed_ncnt > n_limit {
                too_many_n += 1;
                continue;
            }

            // Check unqualified percent (fastp -q/-u logic)
            // unqualified_count already computed by SIMD above
            if qualified_quality_phred > 0 && trimmed_len > 0 {
                // Avoid division to prevent rounding issues: check if 100*unqualified > limit*len
                if 100 * unqualified_count > unqualified_percent_limit * trimmed_len {
                    low_quality += 1;
                    continue;
                }
            }

            // Check average quality (fastp -e logic)
            if average_qual > 0 && trimmed_len > 0 {
                let mean_qual = trimmed_qsum as f64 / trimmed_len as f64;
                if mean_qual < average_qual as f64 {
                    low_quality += 1;
                    continue;
                }
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
            after.add(trimmed_len, trimmed_q20, trimmed_q30, trimmed_gc);
        }

        let result = WorkerResult {
            id: batch.id,
            pieces, // Zero-copy ranges instead of copied bytes
            before,
            after,
            pos,
            k5,
            too_short,
            too_many_n,
            low_quality,
            low_complexity,
            invalid,
        };

        if sender.send(Some(result)).is_err() {
            break; // Receiver disconnected
        }
    }

    // Send sentinel
    let _ = sender.send(None);
}

/// Merger thread: write output in order and reduce stats (zero-copy vectored I/O)
pub(crate) fn merger_thread(
    receiver: Receiver<Option<WorkerResult>>,
    output_path: String,
    num_workers: usize,
    compression_level: Option<u32>,
) -> Result<StreamAccumulator> {
    let inner_writer = open_output(&output_path, compression_level)?;
    // Large BufWriter (16 MiB) to batch writes and reduce syscalls
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, inner_writer);

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
                    acc.before.gc_bases += result.before.gc_bases;

                    acc.after.total_reads += result.after.total_reads;
                    acc.after.total_bases += result.after.total_bases;
                    acc.after.q20_bases += result.after.q20_bases;
                    acc.after.q30_bases += result.after.q30_bases;
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

                    // Merge k-mer counts
                    for (i, &count) in result.k5.iter().enumerate() {
                        acc.kmer_table[i] += count;
                    }

                    // Merge filter counts
                    acc.too_short += result.too_short;
                    acc.too_many_n += result.too_many_n;
                    acc.low_quality += result.low_quality;
                    acc.low_complexity += result.low_complexity;
                    acc.invalid += result.invalid;

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

    // Flush BufWriter and finish inner writer
    writer.flush()?;
    writer.into_inner().map_err(|e| e.into_error())?.finish()?;
    Ok(acc)
}
