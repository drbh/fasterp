//! WebAssembly bindings for fasterp
//!
//! These run the same fastp 1.4.0 compatible engine as the CLI (`compat::run`)
//! over an in-memory file system: inputs are passed in as named byte arrays and
//! every file the run writes is handed back.

use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;

use crate::Args;
use clap::Parser;

/// Run the engine with command-line style `args` (without the program name)
/// over `inputs`, returning every file written: outputs, `fasterp.json` and
/// `fasterp.html` unless `-j`/`-h` name other paths.
pub fn run_files(
    args: &[String],
    inputs: Vec<(String, Vec<u8>)>,
) -> anyhow::Result<BTreeMap<String, Vec<u8>>> {
    let args =
        Args::try_parse_from(std::iter::once("fasterp".to_string()).chain(args.iter().cloned()))?;
    let (result, files) = crate::io::with_files(inputs, || crate::compat::run(args));
    result.map(|()| files)
}

/// Run fasterp exactly like the CLI. `args` is an array of strings such as
/// `["-i", "in.fq", "-o", "out.fq"]` and `files` maps input names to
/// `Uint8Array`s. Returns a `Map` of every file written (FASTQ outputs and the
/// JSON and HTML reports).
#[wasm_bindgen]
pub fn run(args: js_sys::Array, files: js_sys::Map) -> Result<js_sys::Map, JsValue> {
    let args: Vec<String> = args.iter().filter_map(|a| a.as_string()).collect();
    let mut inputs = Vec::new();
    files.for_each(&mut |value, key| {
        if let Some(name) = key.as_string() {
            inputs.push((name, js_sys::Uint8Array::new(&value).to_vec()));
        }
    });
    let outputs = run_files(&args, inputs).map_err(|e| JsValue::from_str(&format!("{e:#}")))?;
    let map = js_sys::Map::new();
    for (name, data) in outputs {
        map.set(
            &JsValue::from_str(&name),
            &js_sys::Uint8Array::from(data.as_slice()),
        );
    }
    Ok(map)
}

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

impl WasmConfig {
    /// The CLI flags for this configuration.
    fn args(&self) -> Vec<String> {
        let mut args = vec![
            format!("-l{}", self.min_length),
            format!("-b{}", self.max_length),
            format!("-e{}", self.average_qual),
            format!("-n{}", self.n_base_limit),
            format!("-q{}", self.qualified_quality_phred),
            format!("-u{}", self.unqualified_percent_limit),
            format!("-Y{}", self.complexity_threshold),
            format!("-f{}", self.trim_front),
            format!("-t{}", self.trim_tail),
            format!("-M{}", self.cut_mean_quality),
            format!("-W{}", self.cut_window_size),
            format!("--poly_g_min_len={}", self.poly_g_min_len),
        ];
        for (on, flag) in [
            (self.low_complexity_filter, "-y"),
            (self.cut_front, "-5"),
            (self.cut_tail, "-3"),
            (self.trim_poly_x, "-x"),
            (self.trim_poly_g, "-g"),
            (!self.trim_poly_g, "-G"),
        ] {
            if on {
                args.push(flag.into());
            }
        }
        args
    }
}

fn js_error(e: anyhow::Error) -> JsValue {
    JsValue::from_str(&format!("{e:#}"))
}

/// Run `args` plus the config flags and return the outputs, the JSON report
/// and the passed and failed read counts.
fn run_config(
    config: &WasmConfig,
    inputs: Vec<(String, Vec<u8>)>,
    mut args: Vec<String>,
    outputs: &[&str],
) -> anyhow::Result<(Vec<Vec<u8>>, String, usize, usize)> {
    args.extend([
        "-j".into(),
        "fasterp.json".into(),
        "-h".into(),
        "fasterp.html".into(),
    ]);
    args.extend(config.args());
    let mut files = run_files(&args, inputs)?;
    let json_report = String::from_utf8(files.remove("fasterp.json").unwrap_or_default())?;
    let report: serde_json::Value = serde_json::from_str(&json_report)?;
    let count = |v: &serde_json::Value| v.as_u64().unwrap_or(0) as usize;
    let total = count(&report["summary"]["before_filtering"]["total_reads"]);
    let passed = count(&report["filtering_result"]["passed_filter_reads"]);
    let outputs = outputs
        .iter()
        .map(|name| files.remove(*name).unwrap_or_default())
        .collect();
    Ok((outputs, json_report, passed, total - passed))
}

/// Process single-end FASTQ data (plain or gzip) with `config`.
#[wasm_bindgen]
pub fn process_single_end(input: &[u8], config: &WasmConfig) -> Result<WasmProcessResult, JsValue> {
    let args = ["-i", "in.fq", "-o", "out.fq"].map(String::from).to_vec();
    let (mut outputs, json_report, passed_reads, failed_reads) = run_config(
        config,
        vec![("in.fq".into(), input.to_vec())],
        args,
        &["out.fq"],
    )
    .map_err(js_error)?;
    Ok(WasmProcessResult {
        output: outputs.remove(0),
        json_report,
        passed_reads,
        failed_reads,
    })
}

/// Process paired-end FASTQ data (plain or gzip) with `config`.
#[wasm_bindgen]
pub fn process_paired_end(
    input1: &[u8],
    input2: &[u8],
    config: &WasmConfig,
) -> Result<WasmPairedResult, JsValue> {
    let args = [
        "-i", "in1.fq", "-I", "in2.fq", "-o", "out1.fq", "-O", "out2.fq",
    ]
    .map(String::from)
    .to_vec();
    let inputs = vec![
        ("in1.fq".into(), input1.to_vec()),
        ("in2.fq".into(), input2.to_vec()),
    ];
    let (mut outputs, json_report, passed_reads, failed_reads) =
        run_config(config, inputs, args, &["out1.fq", "out2.fq"]).map_err(js_error)?;
    Ok(WasmPairedResult {
        output2: outputs.pop().unwrap(),
        output1: outputs.pop().unwrap(),
        json_report,
        passed_reads,
        failed_reads,
    })
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;

    const READ1: &[u8] = b"@read1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    const READ2: &[u8] = b"@read2\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n";

    fn gzip_members(parts: &[&[u8]]) -> Vec<u8> {
        let mut compressed = Vec::new();
        for part in parts {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(part).unwrap();
            compressed.extend(encoder.finish().unwrap());
        }
        compressed
    }

    #[test]
    fn single_end_plain_and_single_member_gzip() {
        let expected = [READ1, READ2].concat();
        for input in [expected.clone(), gzip_members(&[&expected])] {
            let result = process_single_end(&input, &WasmConfig::new()).unwrap();
            assert_eq!(result.passed_reads(), 2);
            assert_eq!(result.failed_reads(), 0);
            assert_eq!(result.output(), expected);
        }
    }

    #[test]
    fn single_end_concatenated_gzip_members() {
        let input = gzip_members(&[READ1, b"", READ2]);
        let result = process_single_end(&input, &WasmConfig::new()).unwrap();
        assert_eq!(result.passed_reads(), 2);
        assert_eq!(result.failed_reads(), 0);
        assert_eq!(result.output(), [READ1, READ2].concat());
    }

    #[test]
    fn gzip_member_boundaries_can_split_fastq_records() {
        let expected = [READ1, READ2].concat();
        for split in 1..expected.len() {
            let input = gzip_members(&[&expected[..split], b"", &expected[split..]]);
            let result = process_single_end(&input, &WasmConfig::new()).unwrap();
            assert_eq!(result.passed_reads(), 2, "split={split}");
            assert_eq!(result.output(), expected, "split={split}");
        }
    }

    #[test]
    fn single_end_empty_first_gzip_member() {
        let input = gzip_members(&[b"", READ1, READ2]);
        let result = process_single_end(&input, &WasmConfig::new()).unwrap();
        assert_eq!(result.passed_reads(), 2);
        assert_eq!(result.failed_reads(), 0);
        assert_eq!(result.output(), [READ1, READ2].concat());
    }

    #[test]
    fn paired_end_plain_single_and_concatenated_gzip_members() {
        let expected1 = [READ1, READ2].concat();
        let expected2 = [READ2, READ1].concat();
        let inputs1 = [
            gzip_members(&[READ1, b"", READ2]),
            expected1.clone(),
            gzip_members(&[&expected1]),
        ];
        let inputs2 = [
            gzip_members(&[READ2, b"", READ1]),
            expected2.clone(),
            gzip_members(&[&expected2]),
        ];
        for input1 in &inputs1 {
            for input2 in &inputs2 {
                let result = process_paired_end(input1, input2, &WasmConfig::new()).unwrap();
                assert_eq!(result.passed_reads(), 4);
                assert_eq!(result.failed_reads(), 0);
                assert_eq!(result.output1(), expected1);
                assert_eq!(result.output2(), expected2);
            }
        }
    }
}
