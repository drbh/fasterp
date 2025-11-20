//! Input/output handling with compression support
//!
//! This module provides I/O abstraction with automatic compression detection:
//! - `CompressionFormat`: Detect compression from file extension or magic bytes
//! - `open_input()`: Open input with auto-decompression (gzip support)
//! - `OutputWriter`: Wrapper for `BufWriter` with optional compression
//! - `open_output()`: Open output with optional compression
//!
//! Note: This module requires the `native` feature (file I/O, S3, HTTP).

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use gzp::ZWriter;
use gzp::deflate::Gzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::Arc;

use aws_sdk_s3::Client as S3Client;
use aws_sdk_s3::primitives::ByteStream;
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

/// Download a URL to a temporary file and return the path
/// Uses cached file if it already exists
fn download_url(url: &str) -> Result<PathBuf> {
    // Get filename from URL or generate one
    let filename = url
        .split('/')
        .next_back()
        .unwrap_or("downloaded_file")
        .split('?')
        .next()
        .unwrap_or("downloaded_file");

    // Create temp directory if it doesn't exist
    let temp_dir = std::env::temp_dir().join("fasterp_downloads");
    std::fs::create_dir_all(&temp_dir).context("Failed to create temp directory for downloads")?;

    let temp_path = temp_dir.join(filename);

    // Check if file already exists (cached)
    if temp_path.exists() {
        let metadata =
            std::fs::metadata(&temp_path).context("Failed to get cached file metadata")?;
        eprintln!(
            "Using cached {} ({} bytes)",
            temp_path.display(),
            metadata.len()
        );
        return Ok(temp_path);
    }

    eprintln!("Downloading {url}...");

    let response = ureq::get(url)
        .call()
        .context(format!("Failed to download URL: {url}"))?;

    // Download to file
    let mut file =
        File::create(&temp_path).context(format!("Failed to create temp file: {temp_path:?}"))?;

    let mut reader = response.into_reader();
    let mut buffer = vec![0u8; 8 * 1024 * 1024]; // 8MB buffer
    let mut total_bytes = 0usize;

    loop {
        let bytes_read = reader
            .read(&mut buffer)
            .context("Failed to read from URL")?;
        if bytes_read == 0 {
            break;
        }
        file.write_all(&buffer[..bytes_read])
            .context("Failed to write to temp file")?;
        total_bytes += bytes_read;
    }

    eprintln!("Downloaded {total_bytes} bytes to {temp_path:?}");

    Ok(temp_path)
}

/// Detect compression format from file extension or magic bytes
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CompressionFormat {
    None,
    Gzip,
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
        if path_lower.ends_with(".gz") || path_lower.ends_with(".gzip") {
            CompressionFormat::Gzip
        } else {
            CompressionFormat::None
        }
    }

    /// Detect from magic bytes (for future auto-detection)
    #[allow(dead_code)]
    pub(crate) fn from_magic_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
            Some(CompressionFormat::Gzip)
        } else {
            None
        }
    }
}

/// Open input file, URL, S3, or stdin with automatic decompression
pub(crate) fn open_input(path: &str) -> Result<Box<dyn BufRead + Send>> {
    if path == "-" {
        // Read from stdin
        let stdin = std::io::stdin();
        let reader = BufReader::with_capacity(16 * 1024 * 1024, stdin);
        Ok(Box::new(reader))
    } else if is_s3_path(path) {
        // Stream from S3
        let s3_reader = S3Reader::new(path)?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                let buffered = BufReader::with_capacity(64 * 1024 * 1024, s3_reader);
                let decoder = GzDecoder::new(buffered);
                let reader = BufReader::with_capacity(32 * 1024 * 1024, decoder);
                Ok(Box::new(reader))
            }
            CompressionFormat::None => {
                let reader = BufReader::with_capacity(32 * 1024 * 1024, s3_reader);
                Ok(Box::new(reader))
            }
        }
    } else if is_url(path) {
        // Download URL to temp file, then open
        let temp_path = download_url(path)?;
        let temp_path_str = temp_path.to_string_lossy();

        let file = File::open(&temp_path)
            .context(format!("Failed to open downloaded file: {temp_path:?}"))?;
        let format = CompressionFormat::from_path(&temp_path_str);

        match format {
            CompressionFormat::Gzip => {
                let buffered_file = BufReader::with_capacity(64 * 1024 * 1024, file);
                let decoder = GzDecoder::new(buffered_file);
                let reader = BufReader::with_capacity(32 * 1024 * 1024, decoder);
                Ok(Box::new(reader))
            }
            CompressionFormat::None => {
                let reader = BufReader::with_capacity(16 * 1024 * 1024, file);
                Ok(Box::new(reader))
            }
        }
    } else {
        // Open file and detect compression
        let file = File::open(path).context(format!("Failed to open input file: {path}"))?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                // Use larger input buffer for compressed file (64MB) to reduce syscalls
                let buffered_file = BufReader::with_capacity(64 * 1024 * 1024, file);
                let decoder = GzDecoder::new(buffered_file);
                // Use 32MB output buffer for decompressed data
                let reader = BufReader::with_capacity(32 * 1024 * 1024, decoder);
                Ok(Box::new(reader))
            }
            CompressionFormat::None => {
                let reader = BufReader::with_capacity(16 * 1024 * 1024, file);
                Ok(Box::new(reader))
            }
        }
    }
}

/// Wrapper for output writer that ensures proper cleanup
pub enum OutputWriter {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
    ParGzip(ParCompress<Gzip>),
    Stdout(BufWriter<std::io::Stdout>),
    S3Plain(S3Writer),
    S3Gzip(BufWriter<GzEncoder<S3Writer>>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Gzip(w) => w.write(buf),
            OutputWriter::ParGzip(w) => w.write(buf),
            OutputWriter::Stdout(w) => w.write(buf),
            OutputWriter::S3Plain(w) => w.write(buf),
            OutputWriter::S3Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Gzip(w) => w.flush(),
            OutputWriter::ParGzip(w) => w.flush(),
            OutputWriter::Stdout(w) => w.flush(),
            OutputWriter::S3Plain(w) => w.flush(),
            OutputWriter::S3Gzip(w) => w.flush(),
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
        }
    }
}

/// Open output file, S3, or stdout with optional compression
pub(crate) fn open_output(
    path: &str,
    compression_level: Option<u32>,
    parallel: bool,
) -> Result<OutputWriter> {
    if path == "-" {
        // Write to stdout
        let stdout = std::io::stdout();
        let writer = BufWriter::with_capacity(16 * 1024 * 1024, stdout);
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
                    // Use parallel compression (4-8 threads)
                    let num_threads = (num_cpus::get() / 2).clamp(4, 8);
                    let compressor = ParCompressBuilder::<Gzip>::new()
                        .compression_level(gzp::Compression::new(level))
                        .num_threads(num_threads)
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
            CompressionFormat::None => {
                let writer = BufWriter::with_capacity(16 * 1024 * 1024, file);
                Ok(OutputWriter::Plain(writer))
            }
        }
    }
}
