//! Input/output handling with compression support
//!
//! This module provides I/O abstraction with automatic compression detection:
//! - `CompressionFormat`: Detect compression from file extension or magic bytes
//! - `open_input()`: Open input with auto-decompression (gzip support)
//! - `OutputWriter`: Wrapper for `BufWriter` with optional compression
//! - `open_output()`: Open output with optional compression
//!
//! Note: This module requires the `native` feature (file I/O, S3, HTTP).

use crate::gzip::Gzip;
use anyhow::{bail, Context, Result};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use gzp::ZWriter;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::sync::Arc;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

use aws_sdk_s3::primitives::ByteStream;
use aws_sdk_s3::Client as S3Client;
use tokio::runtime::Runtime;

/// Check if a path is a URL
fn is_url(path: &str) -> bool {
    path.starts_with("http://") || path.starts_with("https://")
}

/// Check if a path is an S3 URI
pub(crate) fn is_s3_path(path: &str) -> bool {
    path.starts_with("s3://")
}

/// Parse an S3 URI into bucket and key
/// Format: s3://bucket/key/path
pub(crate) fn parse_s3_path(path: &str) -> Result<(String, String)> {
    let path = path.strip_prefix("s3://").context("Invalid S3 path")?;
    let mut parts = path.splitn(2, '/');
    let bucket = parts
        .next()
        .context("Missing bucket in S3 path")?
        .to_string();
    let key = parts.next().context("Missing key in S3 path")?.to_string();

    if bucket.is_empty() {
        bail!("Empty bucket name in S3 path");
    }
    if key.is_empty() {
        bail!("Empty key in S3 path");
    }

    Ok((bucket, key))
}

/// S3 streaming reader that wraps GetObject response
pub struct S3Reader {
    runtime: Arc<Runtime>,
    byte_stream: Option<ByteStream>,
    buffer: Vec<u8>,
    buffer_pos: usize,
}

impl S3Reader {
    /// Create a new S3Reader from an S3 path
    pub fn new(path: &str) -> Result<Self> {
        let (bucket, key) = parse_s3_path(path)?;

        // Create tokio runtime for blocking on async operations
        let runtime = Runtime::new().context("Failed to create tokio runtime")?;

        // Load AWS config and create S3 client
        let config = runtime.block_on(aws_config::load_defaults(
            aws_config::BehaviorVersion::latest(),
        ));

        // Check for custom endpoint (e.g., for Cloudflare R2)
        let client = if let Ok(endpoint) = std::env::var("AWS_ENDPOINT_URL_S3") {
            let s3_config = aws_sdk_s3::config::Builder::from(&config)
                .endpoint_url(&endpoint)
                .force_path_style(true)
                .build();
            S3Client::from_conf(s3_config)
        } else if let Ok(endpoint) = std::env::var("AWS_ENDPOINT_URL") {
            let s3_config = aws_sdk_s3::config::Builder::from(&config)
                .endpoint_url(&endpoint)
                .force_path_style(true)
                .build();
            S3Client::from_conf(s3_config)
        } else {
            S3Client::new(&config)
        };

        // Start GetObject request
        eprintln!("Streaming from s3://{bucket}/{key}...");
        let response = runtime
            .block_on(async { client.get_object().bucket(&bucket).key(&key).send().await })
            .context(format!("Failed to get S3 object: s3://{bucket}/{key}"))?;

        Ok(Self {
            runtime: Arc::new(runtime),
            byte_stream: Some(response.body),
            buffer: Vec::new(),
            buffer_pos: 0,
        })
    }
}

impl Read for S3Reader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // If we have buffered data, return it first
        if self.buffer_pos < self.buffer.len() {
            let available = self.buffer.len() - self.buffer_pos;
            let to_copy = available.min(buf.len());
            buf[..to_copy]
                .copy_from_slice(&self.buffer[self.buffer_pos..self.buffer_pos + to_copy]);
            self.buffer_pos += to_copy;
            return Ok(to_copy);
        }

        // Get next chunk from S3 stream
        let byte_stream = match self.byte_stream.as_mut() {
            Some(stream) => stream,
            None => return Ok(0), // Stream exhausted
        };

        let chunk = self
            .runtime
            .block_on(async { byte_stream.try_next().await })
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;

        match chunk {
            Some(data) => {
                // Store chunk in buffer
                self.buffer = data.to_vec();
                self.buffer_pos = 0;

                let to_copy = self.buffer.len().min(buf.len());
                buf[..to_copy].copy_from_slice(&self.buffer[..to_copy]);
                self.buffer_pos = to_copy;
                Ok(to_copy)
            }
            None => {
                // Stream exhausted
                self.byte_stream = None;
                Ok(0)
            }
        }
    }
}

// S3Reader is Send because it owns all its data
unsafe impl Send for S3Reader {}

/// S3 streaming writer using multipart upload
pub struct S3Writer {
    runtime: Arc<Runtime>,
    client: S3Client,
    bucket: String,
    key: String,
    upload_id: String,
    buffer: Vec<u8>,
    part_number: i32,
    completed_parts: Vec<aws_sdk_s3::types::CompletedPart>,
    min_part_size: usize, // Minimum part size for multipart upload (5MB)
}

impl S3Writer {
    /// Create a new S3Writer for streaming uploads
    pub fn new(path: &str) -> Result<Self> {
        let (bucket, key) = parse_s3_path(path)?;

        // Create tokio runtime
        let runtime = Runtime::new().context("Failed to create tokio runtime")?;

        // Load AWS config and create S3 client
        let config = runtime.block_on(aws_config::load_defaults(
            aws_config::BehaviorVersion::latest(),
        ));

        // Check for custom endpoint (e.g., for Cloudflare R2)
        let client = if let Ok(endpoint) = std::env::var("AWS_ENDPOINT_URL_S3") {
            let s3_config = aws_sdk_s3::config::Builder::from(&config)
                .endpoint_url(&endpoint)
                .force_path_style(true)
                .build();
            S3Client::from_conf(s3_config)
        } else if let Ok(endpoint) = std::env::var("AWS_ENDPOINT_URL") {
            let s3_config = aws_sdk_s3::config::Builder::from(&config)
                .endpoint_url(&endpoint)
                .force_path_style(true)
                .build();
            S3Client::from_conf(s3_config)
        } else {
            S3Client::new(&config)
        };

        // Start multipart upload
        eprintln!("Starting multipart upload to s3://{bucket}/{key}...");
        let create_response = runtime
            .block_on(async {
                client
                    .create_multipart_upload()
                    .bucket(&bucket)
                    .key(&key)
                    .send()
                    .await
            })
            .context(format!(
                "Failed to create multipart upload: s3://{bucket}/{key}"
            ))?;

        let upload_id = create_response
            .upload_id()
            .context("No upload ID returned")?
            .to_string();

        Ok(Self {
            runtime: Arc::new(runtime),
            client,
            bucket,
            key,
            upload_id,
            buffer: Vec::with_capacity(8 * 1024 * 1024), // 8MB initial capacity
            part_number: 1,
            completed_parts: Vec::new(),
            min_part_size: 5 * 1024 * 1024, // 5MB minimum for S3
        })
    }

    /// Upload buffered data as a part
    /// If `force` is true, upload whatever is in the buffer (for final part)
    /// Otherwise, only upload if buffer has exactly min_part_size bytes
    fn upload_part(&mut self, force: bool) -> std::io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        // For non-final parts, only upload when we have exactly min_part_size
        // This ensures all non-trailing parts have the same size (R2 requirement)
        if !force && self.buffer.len() < self.min_part_size {
            return Ok(());
        }

        let data = if force {
            // Final part - take everything
            std::mem::take(&mut self.buffer)
        } else {
            // Take exactly min_part_size bytes, keep the rest
            let remaining = self.buffer.split_off(self.min_part_size);
            std::mem::replace(&mut self.buffer, remaining)
        };

        let part_number = self.part_number;

        let response = self
            .runtime
            .block_on(async {
                self.client
                    .upload_part()
                    .bucket(&self.bucket)
                    .key(&self.key)
                    .upload_id(&self.upload_id)
                    .part_number(part_number)
                    .body(ByteStream::from(data))
                    .send()
                    .await
            })
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;

        let e_tag = response.e_tag().unwrap_or_default().to_string();

        let completed_part = aws_sdk_s3::types::CompletedPart::builder()
            .part_number(part_number)
            .e_tag(e_tag)
            .build();

        self.completed_parts.push(completed_part);
        self.part_number += 1;

        Ok(())
    }

    /// Complete the multipart upload
    pub fn finish(mut self) -> Result<()> {
        // Upload any remaining data (force=true for final part)
        if !self.buffer.is_empty() {
            self.upload_part(true)
                .context("Failed to upload final part")?;
        }

        // Complete multipart upload
        let completed_upload = aws_sdk_s3::types::CompletedMultipartUpload::builder()
            .set_parts(Some(self.completed_parts.clone()))
            .build();

        self.runtime
            .block_on(async {
                self.client
                    .complete_multipart_upload()
                    .bucket(&self.bucket)
                    .key(&self.key)
                    .upload_id(&self.upload_id)
                    .multipart_upload(completed_upload)
                    .send()
                    .await
            })
            .context(format!(
                "Failed to complete multipart upload: s3://{}/{}",
                self.bucket, self.key
            ))?;

        eprintln!("Upload completed: s3://{}/{}", self.bucket, self.key);

        Ok(())
    }
}

impl Write for S3Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.buffer.extend_from_slice(buf);

        // Upload parts while buffer has enough data
        // Each part will be exactly min_part_size bytes
        while self.buffer.len() >= self.min_part_size {
            self.upload_part(false)?;
        }

        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        // Don't upload on flush - wait until buffer is full or finish is called
        Ok(())
    }
}

// S3Writer is Send because it owns all its data
unsafe impl Send for S3Writer {}

/// Detect compression format from file extension or magic bytes
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CompressionFormat {
    None,
    Gzip,
    Zstd,
}

impl CompressionFormat {
    /// Detect from file path
    #[allow(clippy::case_sensitive_file_extension_comparisons)]
    pub(crate) fn from_path(path: &str) -> Self {
        if path == "-" {
            return CompressionFormat::None; // stdin/stdout defaults to uncompressed
        }

        // path_lower is already lowercased, so comparison is case-insensitive
        let path_lower = path.to_lowercase();
        if path_lower.ends_with(".gz")
            || path_lower.ends_with(".gzip")
            || path_lower.ends_with(".bgz")
            || path_lower.ends_with(".bgzf")
        {
            CompressionFormat::Gzip
        } else if path_lower.ends_with(".zst") || path_lower.ends_with(".zstd") {
            CompressionFormat::Zstd
        } else {
            CompressionFormat::None
        }
    }

    /// Detect from magic bytes (for future auto-detection)
    #[allow(dead_code)]
    pub(crate) fn from_magic_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
            Some(CompressionFormat::Gzip)
        } else if bytes.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
            Some(CompressionFormat::Zstd)
        } else {
            None
        }
    }
}

/// Bounded decompression for regular local files. Paired inputs can decode concurrently.
struct PrefetchInput {
    receiver: Option<crossbeam_channel::Receiver<std::io::Result<Vec<u8>>>>,
    recycle: crossbeam_channel::Sender<Vec<u8>>,
    worker: Option<std::thread::JoinHandle<()>>,
    buffer: Vec<u8>,
    position: usize,
    finished: bool,
}

impl PrefetchInput {
    fn new(mut input: impl Read + Send + 'static) -> Self {
        let (sender, receiver) = crossbeam_channel::bounded(2);
        let (recycle, reused) = crossbeam_channel::bounded::<Vec<u8>>(2);
        let worker = std::thread::spawn(move || loop {
            let mut buffer = reused.try_recv().unwrap_or_default();
            buffer.resize(1024 * 1024, 0);
            let mut filled = 0;
            while filled < buffer.len() {
                match input.read(&mut buffer[filled..]) {
                    Ok(0) => break,
                    Ok(n) => filled += n,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    Err(e) => {
                        let _ = sender.send(Err(e));
                        return;
                    }
                }
            }
            buffer.truncate(filled);
            if sender.send(Ok(buffer)).is_err() || filled == 0 {
                break;
            }
        });
        Self {
            receiver: Some(receiver),
            recycle,
            worker: Some(worker),
            buffer: Vec::new(),
            position: 0,
            finished: false,
        }
    }
}

impl BufRead for PrefetchInput {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.position == self.buffer.len() && !self.finished {
            let next = self.receiver.as_ref().unwrap().recv().map_err(|_| {
                std::io::Error::other("Decompression worker stopped before end of stream")
            })??;
            let old = std::mem::replace(&mut self.buffer, next);
            let _ = self.recycle.try_send(old);
            self.position = 0;
            self.finished = self.buffer.is_empty();
        }
        Ok(&self.buffer[self.position..])
    }

    fn consume(&mut self, amount: usize) {
        self.position = (self.position + amount).min(self.buffer.len());
    }
}

impl Read for PrefetchInput {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }
        let available = self.fill_buf()?;
        let count = output.len().min(available.len());
        output[..count].copy_from_slice(&available[..count]);
        self.consume(count);
        Ok(count)
    }
}

impl Drop for PrefetchInput {
    fn drop(&mut self) {
        // Disconnect first: a reader stopped by --reads_to_process or an output
        // failure must not wait on a decompressor blocked in send(). Only regular
        // local files use this wrapper, so the underlying read can finish.
        drop(self.receiver.take());
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

/// Open input file, URL, S3, or stdin with automatic decompression.
pub(crate) fn open_input(path: &str) -> Result<Box<dyn BufRead + Send>> {
    let raw: Box<dyn Read + Send> = if path == "-" || path == "/dev/stdin" {
        Box::new(std::io::stdin())
    } else if is_s3_path(path) {
        Box::new(S3Reader::new(path)?)
    } else if is_url(path) {
        Box::new(
            ureq::get(path)
                .call()
                .with_context(|| format!("Failed to open URL: {path}"))?
                .into_reader(),
        )
    } else {
        Box::new(File::open(path).with_context(|| format!("Failed to open input file: {path}"))?)
    };
    let mut buffered = buffer_input(raw)?;
    let format = CompressionFormat::from_magic_bytes(buffered.fill_buf()?)
        .unwrap_or_else(|| CompressionFormat::from_path(path));
    match format {
        CompressionFormat::Gzip => {
            // Reuse the compressed-input buffer instead of adding the read
            // decoder's smaller internal buffer and an extra copy.
            let decoder = MultiGzDecoder::new(buffered);
            if std::fs::metadata(path).is_ok_and(|m| m.is_file()) {
                Ok(Box::new(PrefetchInput::new(decoder)))
            } else {
                Ok(Box::new(BufReader::with_capacity(1024 * 1024, decoder)))
            }
        }
        CompressionFormat::Zstd => Ok(Box::new(BufReader::with_capacity(
            1024 * 1024,
            ZstdDecoder::new(buffered)?,
        ))),
        CompressionFormat::None => Ok(Box::new(buffered)),
    }
}

fn buffer_input(raw: Box<dyn Read + Send>) -> std::io::Result<BufReader<Box<dyn Read + Send>>> {
    let mut buffered = BufReader::with_capacity(1024 * 1024, raw);
    let available = buffered.fill_buf()?.len();
    if available > 0 && available < 4 {
        // A pipe may return only part of the compression magic. Preserve that
        // prefix and obtain enough bytes to detect gzip/zstd before parsing it.
        // Ordinary file reads keep the existing buffer without this wrapper.
        let mut prefix = Vec::with_capacity(4);
        Read::by_ref(&mut buffered)
            .take(4)
            .read_to_end(&mut prefix)?;
        buffered = BufReader::with_capacity(
            1024 * 1024,
            Box::new(std::io::Cursor::new(prefix).chain(buffered)),
        );
    }
    Ok(buffered)
}

#[cfg(test)]
mod input_tests {
    use super::*;

    struct Fragmented(std::io::Cursor<Vec<u8>>);
    impl Read for Fragmented {
        fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
            let count = bytes.len().min(1);
            self.0.read(&mut bytes[..count])
        }
    }

    #[test]
    fn compression_magic_survives_one_byte_reads() {
        let data = b"@read\nACGT\n+\nIIII\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        let gzip = encoder.finish().unwrap();
        let zstd = zstd::stream::encode_all(&data[..], 1).unwrap();
        for (encoded, expected) in [
            (gzip, CompressionFormat::Gzip),
            (zstd, CompressionFormat::Zstd),
        ] {
            let mut input =
                buffer_input(Box::new(Fragmented(std::io::Cursor::new(encoded)))).unwrap();
            assert_eq!(
                CompressionFormat::from_magic_bytes(input.fill_buf().unwrap()),
                Some(expected)
            );
            let mut decoded = Vec::new();
            if expected == CompressionFormat::Gzip {
                MultiGzDecoder::new(input)
                    .read_to_end(&mut decoded)
                    .unwrap();
            } else {
                ZstdDecoder::new(input)
                    .unwrap()
                    .read_to_end(&mut decoded)
                    .unwrap();
            }
            assert_eq!(decoded, data);
        }
    }

    #[test]
    fn short_uncompressed_inputs_are_not_lost() {
        for len in 0..=4 {
            let data = b"@abc"[..len].to_vec();
            let mut input =
                buffer_input(Box::new(Fragmented(std::io::Cursor::new(data.clone())))).unwrap();
            let mut actual = Vec::new();
            input.read_to_end(&mut actual).unwrap();
            assert_eq!(actual, data);
        }
    }
}

/// Wrapper for output writer that ensures proper cleanup
pub enum OutputWriter {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
    Zstd(BufWriter<ZstdEncoder<'static, File>>),
    ParGzip(ParCompress<Gzip>),
    Stdout(BufWriter<std::io::Stdout>),
    S3Plain(S3Writer),
    S3Gzip(BufWriter<GzEncoder<S3Writer>>),
    S3Zstd(BufWriter<ZstdEncoder<'static, S3Writer>>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Gzip(w) => w.write(buf),
            OutputWriter::Zstd(w) => w.write(buf),
            OutputWriter::ParGzip(w) => w.write(buf),
            OutputWriter::Stdout(w) => w.write(buf),
            OutputWriter::S3Plain(w) => w.write(buf),
            OutputWriter::S3Gzip(w) => w.write(buf),
            OutputWriter::S3Zstd(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Gzip(w) => w.flush(),
            OutputWriter::Zstd(w) => w.flush(),
            OutputWriter::ParGzip(w) => w.flush(),
            OutputWriter::Stdout(w) => w.flush(),
            OutputWriter::S3Plain(w) => w.flush(),
            OutputWriter::S3Gzip(w) => w.flush(),
            OutputWriter::S3Zstd(w) => w.flush(),
        }
    }
}

impl OutputWriter {
    /// Finish writing and ensure gzip encoder is properly finalized
    pub fn finish(self) -> Result<()> {
        match self {
            OutputWriter::Plain(mut w) => {
                w.flush()?;
                Ok(())
            }
            OutputWriter::Gzip(mut w) => {
                w.flush()?;
                let encoder = w.into_inner().context("Failed to finish gzip writer")?;
                encoder.finish().context("Failed to finish gzip encoding")?;
                Ok(())
            }
            OutputWriter::Zstd(mut w) => {
                w.flush()?;
                let encoder = w.into_inner().map_err(|e| anyhow::anyhow!("{}", e))?;
                encoder.finish()?;
                Ok(())
            }
            OutputWriter::ParGzip(mut w) => {
                // MUST call finish() before ParCompress goes out of scope
                w.flush().context("Failed to flush parallel gzip writer")?;
                w.finish()
                    .context("Failed to finish parallel gzip compression")?;
                Ok(())
            }
            OutputWriter::Stdout(mut w) => {
                w.flush()?;
                Ok(())
            }
            OutputWriter::S3Plain(w) => {
                w.finish().context("Failed to finish S3 upload")?;
                Ok(())
            }
            OutputWriter::S3Gzip(mut w) => {
                w.flush()?;
                let encoder = w.into_inner().map_err(|e| {
                    anyhow::anyhow!("Failed to finish gzip writer for S3: {}", e.error())
                })?;
                let s3_writer = encoder
                    .finish()
                    .context("Failed to finish gzip encoding for S3")?;
                s3_writer.finish().context("Failed to finish S3 upload")?;
                Ok(())
            }
            OutputWriter::S3Zstd(mut w) => {
                w.flush()?;
                let encoder = w.into_inner().map_err(|e| {
                    anyhow::anyhow!("Failed to finish zstd writer for S3: {}", e.error())
                })?;
                let s3_writer = encoder.finish()?;
                s3_writer.finish().context("Failed to finish S3 upload")?;
                Ok(())
            }
        }
    }
}

/// Open output file, S3, or stdout with optional compression
pub(crate) fn open_output(
    path: &str,
    compression_level: Option<u32>,
    parallel: bool,
) -> Result<OutputWriter> {
    open_output_with_threads(
        path,
        compression_level,
        parallel,
        crate::gzip::threads_per_output(1),
    )
}

pub(crate) fn open_output_with_threads(
    path: &str,
    compression_level: Option<u32>,
    parallel: bool,
    compression_threads: usize,
) -> Result<OutputWriter> {
    if path == "-" {
        // Write to stdout
        let stdout = std::io::stdout();
        // Bound pipeline latency and detect a closed consumer promptly. Packed
        // worker output still writes large contiguous slices through this buffer.
        let writer = BufWriter::with_capacity(256 * 1024, stdout);
        Ok(OutputWriter::Stdout(writer))
    } else if is_s3_path(path) {
        // Stream to S3
        let s3_writer = S3Writer::new(path)?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                let level = compression_level.unwrap_or(6);
                let compression = Compression::new(level);
                let encoder = GzEncoder::new(s3_writer, compression);
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, encoder);
                Ok(OutputWriter::S3Gzip(writer))
            }
            CompressionFormat::Zstd => {
                // Zstd supports levels 1-22, default is 3
                let level = compression_level.unwrap_or(3).clamp(1, 22) as i32;
                let encoder = ZstdEncoder::new(s3_writer, level)?;
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, encoder);
                Ok(OutputWriter::S3Zstd(writer))
            }
            CompressionFormat::None => Ok(OutputWriter::S3Plain(s3_writer)),
        }
    } else {
        // Create file and detect compression
        let file = File::create(path).context(format!("Failed to create output file: {path}"))?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                let level = compression_level.unwrap_or(6); // Default to level 6

                if parallel {
                    let compressor = ParCompressBuilder::<Gzip>::new()
                        .compression_level(gzp::Compression::new(level))
                        .num_threads(compression_threads)
                        .context("Failed to set number of threads for parallel compression")?
                        .from_writer(file);
                    Ok(OutputWriter::ParGzip(compressor))
                } else {
                    // Use standard single-threaded
                    let compression = Compression::new(level);
                    let encoder = GzEncoder::new(file, compression);
                    let writer = BufWriter::with_capacity(16 * 1024 * 1024, encoder);
                    Ok(OutputWriter::Gzip(writer))
                }
            }
            CompressionFormat::Zstd => {
                // Zstd supports levels 1-22, default is 3
                let level = compression_level.unwrap_or(3).clamp(1, 22) as i32;
                let encoder = ZstdEncoder::new(file, level)?;
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, encoder);
                Ok(OutputWriter::Zstd(writer))
            }
            CompressionFormat::None => {
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, file);
                Ok(OutputWriter::Plain(writer))
            }
        }
    }
}
