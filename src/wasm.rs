//! WebAssembly bindings for fasterp
//!
//! This module provides WASM-compatible functions for FASTQ processing.
//! All I/O is done via byte arrays instead of files.

use wasm_bindgen::prelude::*;
use std::io::Cursor;

use crate::adapter::AdapterConfig;
use crate::processor::{process_fastq_stream, process_paired_fastq_stream};
use crate::stats;
use crate::trimming::TrimmingConfig;

/// Result of processing a FASTQ file
#[wasm_bindgen]
pub struct WasmProcessResult {
    /// Processed FASTQ data
    output: Vec<u8>,
    /// JSON report
    json_report: String,
    /// Number of reads that passed filtering
    passed_reads: usize,
    /// Number of reads that failed filtering
    failed_reads: usize,
}

#[wasm_bindgen]
impl WasmProcessResult {
    /// Get the processed FASTQ output
    #[wasm_bindgen(getter)]
    pub fn output(&self) -> Vec<u8> {
        self.output.clone()
    }

    /// Get the JSON report
    #[wasm_bindgen(getter)]
    pub fn json_report(&self) -> String {
        self.json_report.clone()
    }

    /// Get the number of passed reads
    #[wasm_bindgen(getter)]
    pub fn passed_reads(&self) -> usize {
        self.passed_reads
    }

    /// Get the number of failed reads
    #[wasm_bindgen(getter)]
    pub fn failed_reads(&self) -> usize {
        self.failed_reads
    }
}

/// Process configuration for WASM
#[wasm_bindgen]
pub struct WasmConfig {
    min_length: usize,
    max_length: usize,
    average_qual: u8,
    n_base_limit: usize,
    qualified_quality_phred: u8,
    unqualified_percent_limit: usize,
    low_complexity_filter: bool,
    complexity_threshold: usize,
    trim_front: usize,
    trim_tail: usize,
    cut_front: bool,
    cut_tail: bool,
    cut_mean_quality: u8,
    cut_window_size: usize,
    trim_poly_g: bool,
    trim_poly_x: bool,
    poly_g_min_len: usize,
}

#[wasm_bindgen]
impl WasmConfig {
    /// Create a new configuration with default values
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            min_length: 15,
            max_length: 0,
            average_qual: 0,
            n_base_limit: 5,
            qualified_quality_phred: 15,
            unqualified_percent_limit: 40,
            low_complexity_filter: false,
            complexity_threshold: 30,
            trim_front: 0,
            trim_tail: 0,
            cut_front: false,
            cut_tail: true,
            cut_mean_quality: 20,
            cut_window_size: 4,
            trim_poly_g: true,
            trim_poly_x: false,
            poly_g_min_len: 10,
        }
    }

    /// Set minimum read length after processing
    #[wasm_bindgen(setter)]
    pub fn set_min_length(&mut self, value: usize) {
        self.min_length = value;
    }

    /// Set maximum read length (0 = no limit)
    #[wasm_bindgen(setter)]
    pub fn set_max_length(&mut self, value: usize) {
        self.max_length = value;
    }

    /// Set minimum average quality score
    #[wasm_bindgen(setter)]
    pub fn set_average_qual(&mut self, value: u8) {
        self.average_qual = value;
    }

    /// Set maximum N bases allowed
    #[wasm_bindgen(setter)]
    pub fn set_n_base_limit(&mut self, value: usize) {
        self.n_base_limit = value;
    }

    /// Set quality threshold for qualified bases
    #[wasm_bindgen(setter)]
    pub fn set_qualified_quality_phred(&mut self, value: u8) {
        self.qualified_quality_phred = value;
    }

    /// Set maximum percentage of unqualified bases
    #[wasm_bindgen(setter)]
    pub fn set_unqualified_percent_limit(&mut self, value: usize) {
        self.unqualified_percent_limit = value;
    }

    /// Enable/disable low complexity filter
    #[wasm_bindgen(setter)]
    pub fn set_low_complexity_filter(&mut self, value: bool) {
        self.low_complexity_filter = value;
    }

    /// Set complexity threshold percentage
    #[wasm_bindgen(setter)]
    pub fn set_complexity_threshold(&mut self, value: usize) {
        self.complexity_threshold = value;
    }

    /// Enable sliding window quality trimming from front
    #[wasm_bindgen(setter)]
    pub fn set_cut_front(&mut self, value: bool) {
        self.cut_front = value;
    }

    /// Enable sliding window quality trimming from tail
    #[wasm_bindgen(setter)]
    pub fn set_cut_tail(&mut self, value: bool) {
        self.cut_tail = value;
    }

    /// Set mean quality threshold for sliding window
    #[wasm_bindgen(setter)]
    pub fn set_cut_mean_quality(&mut self, value: u8) {
        self.cut_mean_quality = value;
    }

    /// Set sliding window size
    #[wasm_bindgen(setter)]
    pub fn set_cut_window_size(&mut self, value: usize) {
        self.cut_window_size = value;
    }

    /// Enable poly-G tail trimming
    #[wasm_bindgen(setter)]
    pub fn set_trim_poly_g(&mut self, value: bool) {
        self.trim_poly_g = value;
    }

    /// Enable poly-X tail trimming
    #[wasm_bindgen(setter)]
    pub fn set_trim_poly_x(&mut self, value: bool) {
        self.trim_poly_x = value;
    }
}

impl Default for WasmConfig {
    fn default() -> Self {
        Self::new()
    }
}

/// Process single-end FASTQ data
///
/// # Arguments
/// * `input` - Raw FASTQ data (can be gzip compressed)
/// * `config` - Processing configuration
///
/// # Returns
/// `WasmProcessResult` containing processed data and statistics
#[wasm_bindgen]
pub fn process_single_end(input: &[u8], config: &WasmConfig) -> Result<WasmProcessResult, JsValue> {
    // Detect if input is gzip compressed
    let reader: Box<dyn std::io::BufRead> = if input.len() >= 2 && input[0] == 0x1f && input[1] == 0x8b {
        // Gzip compressed
        use flate2::read::GzDecoder;
        let decoder = GzDecoder::new(Cursor::new(input));
        Box::new(std::io::BufReader::new(decoder))
    } else {
        // Plain text
        Box::new(std::io::BufReader::new(Cursor::new(input)))
    };

    let mut output = Vec::new();

    // Create adapter config (disabled for simplicity in WASM)
    let adapter_config = AdapterConfig::new();

    // Create trimming config
    let trimming_config = TrimmingConfig {
        enable_trim_front: config.cut_front && config.cut_mean_quality > 0,
        enable_trim_tail: config.cut_tail && config.cut_mean_quality > 0,
        cut_mean_quality: config.cut_mean_quality,
        cut_window_size: config.cut_window_size,
        trim_front_bases: config.trim_front,
        trim_tail_bases: config.trim_tail,
        max_len: config.max_length,
        enable_poly_g: config.trim_poly_g,
        enable_poly_x: config.trim_poly_x,
        poly_min_len: config.poly_g_min_len,
        adapter_config,
    };

    // Process
    let acc = process_fastq_stream(
        reader,
        &mut output,
        config.min_length,
        config.n_base_limit,
        config.qualified_quality_phred,
        config.unqualified_percent_limit,
        config.average_qual,
        config.low_complexity_filter,
        config.complexity_threshold,
        &trimming_config,
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;

    // Build report
    let before_stats = acc.before.to_read_stats();
    let after_stats = acc.after.to_read_stats();
    let kmer_map = acc.kmer_table_to_map();
    let duplication_rate = stats::calculate_duplication_rate(&kmer_map);

    // Create a simple JSON report
    let report = serde_json::json!({
        "summary": {
            "before_filtering": {
                "total_reads": before_stats.total_reads,
                "total_bases": before_stats.total_bases,
                "q20_bases": before_stats.q20_bases,
                "q30_bases": before_stats.q30_bases,
                "gc_content": before_stats.gc_content,
            },
            "after_filtering": {
                "total_reads": after_stats.total_reads,
                "total_bases": after_stats.total_bases,
                "q20_bases": after_stats.q20_bases,
                "q30_bases": after_stats.q30_bases,
                "gc_content": after_stats.gc_content,
            }
        },
        "filtering_result": {
            "passed_filter_reads": acc.after.total_reads,
            "low_quality_reads": acc.low_quality,
            "low_complexity_reads": acc.low_complexity,
            "too_many_n_reads": acc.too_many_n,
            "too_short_reads": acc.too_short,
        },
        "duplication": {
            "rate": duplication_rate
        }
    });

    let json_report = serde_json::to_string_pretty(&report)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

    Ok(WasmProcessResult {
        output,
        json_report,
        passed_reads: acc.after.total_reads,
        failed_reads: acc.before.total_reads - acc.after.total_reads,
    })
}

/// Process paired-end FASTQ data
///
/// # Arguments
/// * `input1` - Raw FASTQ data for read 1
/// * `input2` - Raw FASTQ data for read 2
/// * `config` - Processing configuration
///
/// # Returns
/// Tuple of (output1, output2, json_report)
#[wasm_bindgen]
pub fn process_paired_end(
    input1: &[u8],
    input2: &[u8],
    config: &WasmConfig,
) -> Result<WasmPairedResult, JsValue> {
    // Detect if inputs are gzip compressed
    let reader1: Box<dyn std::io::BufRead> = if input1.len() >= 2 && input1[0] == 0x1f && input1[1] == 0x8b {
        use flate2::read::GzDecoder;
        let decoder = GzDecoder::new(Cursor::new(input1));
        Box::new(std::io::BufReader::new(decoder))
    } else {
        Box::new(std::io::BufReader::new(Cursor::new(input1)))
    };

    let reader2: Box<dyn std::io::BufRead> = if input2.len() >= 2 && input2[0] == 0x1f && input2[1] == 0x8b {
        use flate2::read::GzDecoder;
        let decoder = GzDecoder::new(Cursor::new(input2));
        Box::new(std::io::BufReader::new(decoder))
    } else {
        Box::new(std::io::BufReader::new(Cursor::new(input2)))
    };

    let mut output1 = Vec::new();
    let mut output2 = Vec::new();

    // Create adapter config
    let adapter_config = AdapterConfig::new();

    // Create trimming configs
    let trimming_config = TrimmingConfig {
        enable_trim_front: config.cut_front && config.cut_mean_quality > 0,
        enable_trim_tail: config.cut_tail && config.cut_mean_quality > 0,
        cut_mean_quality: config.cut_mean_quality,
        cut_window_size: config.cut_window_size,
        trim_front_bases: config.trim_front,
        trim_tail_bases: config.trim_tail,
        max_len: config.max_length,
        enable_poly_g: config.trim_poly_g,
        enable_poly_x: config.trim_poly_x,
        poly_min_len: config.poly_g_min_len,
        adapter_config: adapter_config.clone(),
    };

    let trimming_config_r2 = trimming_config.clone();

    // Process
    let pe_acc = process_paired_fastq_stream(
        reader1,
        reader2,
        &mut output1,
        &mut output2,
        None::<&mut Vec<u8>>,
        false,
        false,
        config.min_length,
        config.n_base_limit,
        config.qualified_quality_phred,
        config.unqualified_percent_limit,
        config.average_qual,
        config.low_complexity_filter,
        config.complexity_threshold,
        &trimming_config,
        &trimming_config_r2,
        None,
        None,
        None,
        None,
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;

    // Build report
    let before_stats_r1 = pe_acc.before_r1.to_read_stats();
    let after_stats_r1 = pe_acc.after_r1.to_read_stats();
    let before_stats_r2 = pe_acc.before_r2.to_read_stats();
    let after_stats_r2 = pe_acc.after_r2.to_read_stats();

    let report = serde_json::json!({
        "summary": {
            "read1_before_filtering": {
                "total_reads": before_stats_r1.total_reads,
                "total_bases": before_stats_r1.total_bases,
                "gc_content": before_stats_r1.gc_content,
            },
            "read1_after_filtering": {
                "total_reads": after_stats_r1.total_reads,
                "total_bases": after_stats_r1.total_bases,
                "gc_content": after_stats_r1.gc_content,
            },
            "read2_before_filtering": {
                "total_reads": before_stats_r2.total_reads,
                "total_bases": before_stats_r2.total_bases,
                "gc_content": before_stats_r2.gc_content,
            },
            "read2_after_filtering": {
                "total_reads": after_stats_r2.total_reads,
                "total_bases": after_stats_r2.total_bases,
                "gc_content": after_stats_r2.gc_content,
            }
        },
        "filtering_result": {
            "passed_filter_reads": pe_acc.after_r1.total_reads,
            "low_quality_reads": pe_acc.low_quality,
            "low_complexity_reads": pe_acc.low_complexity,
            "too_many_n_reads": pe_acc.too_many_n,
            "too_short_reads": pe_acc.too_short,
        }
    });

    let json_report = serde_json::to_string_pretty(&report)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

    let total_before = before_stats_r1.total_reads + before_stats_r2.total_reads;
    let total_after = after_stats_r1.total_reads + after_stats_r2.total_reads;

    Ok(WasmPairedResult {
        output1,
        output2,
        json_report,
        passed_reads: total_after,
        failed_reads: total_before - total_after,
    })
}

/// Result of processing paired-end FASTQ files
#[wasm_bindgen]
pub struct WasmPairedResult {
    output1: Vec<u8>,
    output2: Vec<u8>,
    json_report: String,
    passed_reads: usize,
    failed_reads: usize,
}

#[wasm_bindgen]
impl WasmPairedResult {
    /// Get processed output for read 1
    #[wasm_bindgen(getter)]
    pub fn output1(&self) -> Vec<u8> {
        self.output1.clone()
    }

    /// Get processed output for read 2
    #[wasm_bindgen(getter)]
    pub fn output2(&self) -> Vec<u8> {
        self.output2.clone()
    }

    /// Get the JSON report
    #[wasm_bindgen(getter)]
    pub fn json_report(&self) -> String {
        self.json_report.clone()
    }

    /// Get the number of passed reads
    #[wasm_bindgen(getter)]
    pub fn passed_reads(&self) -> usize {
        self.passed_reads
    }

    /// Get the number of failed reads
    #[wasm_bindgen(getter)]
    pub fn failed_reads(&self) -> usize {
        self.failed_reads
    }
}

/// Get the version of fasterp
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}
