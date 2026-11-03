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

/// The total size of the gzip member starting at `data`, from a size subfield
/// in its extra field (fasterp's `FP`, or BGZF's `BC`), or `None` when it has
/// none. Needs the member's first 12 bytes plus its extra field.
fn indexed_member_size(data: &[u8]) -> Option<usize> {
    if data.len() < 12 || data[..3] != [0x1f, 0x8b, 8] || data[3] & 4 == 0 {
        return None;
    }
    let xlen = usize::from(u16::from_le_bytes([data[10], data[11]]));
    let mut extra = data.get(12..12 + xlen)?;
    while extra.len() >= 4 {
        let len = usize::from(u16::from_le_bytes([extra[2], extra[3]]));
        let field = extra.get(4..4 + len)?;
        match (&extra[..2], len) {
            (id, 4) if id == crate::gzip::MEMBER_SIZE_ID => {
                return Some(u32::from_le_bytes(field.try_into().unwrap()) as usize)
            }
            (b"BC", 2) => return Some(usize::from(u16::from_le_bytes([field[0], field[1]])) + 1),
            _ => extra = &extra[4 + len..],
        }
    }
    None
}

enum MemberJob {
    /// Whole gzip members, back to back, and their total decompressed size.
    Members(Vec<u8>, usize),
    /// Bytes already decompressed by the sequential fallback.
    Plain(Vec<u8>),
    /// A read or decode error, delivered in order.
    Failed(std::io::Error),
}

/// Decompress a gzip file whose members record their sizes on `threads`
/// threads with libdeflate, sending output in order on `sender`. From the
/// first member without a size it falls back to the sequential decoder.
fn decompress_members(
    mut file: File,
    threads: usize,
    sender: crossbeam_channel::Sender<std::io::Result<Vec<u8>>>,
) {
    use std::collections::BTreeMap;
    const JOB_BYTES: usize = 1024 * 1024;
    let (jobs, queue) = crossbeam_channel::bounded::<(usize, MemberJob)>(threads * 2);
    let (done, results) =
        crossbeam_channel::bounded::<(usize, std::io::Result<Vec<u8>>)>(threads * 2);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let (queue, done) = (queue.clone(), done.clone());
            scope.spawn(move || {
                let mut decompressor = libdeflater::Decompressor::new();
                for (seq, job) in queue {
                    let result = match job {
                        MemberJob::Plain(bytes) => Ok(bytes),
                        MemberJob::Failed(error) => Err(error),
                        MemberJob::Members(bytes, size) => {
                            let mut out = vec![0; size];
                            let mut at = 0;
                            let mut written = 0;
                            let mut result = Ok(());
                            while at < bytes.len() {
                                let total = indexed_member_size(&bytes[at..]).unwrap();
                                let member = &bytes[at..at + total];
                                let isize =
                                    u32::from_le_bytes(member[total - 4..].try_into().unwrap())
                                        as usize;
                                match decompressor
                                    .gzip_decompress(member, &mut out[written..written + isize])
                                {
                                    Ok(n) if n == isize => written += n,
                                    other => {
                                        result = Err(std::io::Error::new(
                                            std::io::ErrorKind::InvalidData,
                                            format!("Corrupt gzip member: {other:?}"),
                                        ));
                                        break;
                                    }
                                }
                                at += total;
                            }
                            result.map(|()| out)
                        }
                    };
                    if done.send((seq, result)).is_err() {
                        break;
                    }
                }
            });
        }
        drop((queue, done));
        // Split the file into jobs of whole members.
        scope.spawn(move || {
            let mut buffer: Vec<u8> = Vec::new();
            let mut start = 0;
            let mut eof = false;
            let mut seq = 0;
            let mut fill =
                |buffer: &mut Vec<u8>, start: &mut usize, need: usize| -> std::io::Result<bool> {
                    while buffer.len() - *start < need {
                        if *start > 0 {
                            buffer.drain(..*start);
                            *start = 0;
                        }
                        let old = buffer.len();
                        buffer.resize(old + (8 << 20).max(need), 0);
                        let n = loop {
                            match file.read(&mut buffer[old..]) {
                                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                                other => break other?,
                            }
                        };
                        buffer.truncate(old + n);
                        if n == 0 {
                            return Ok(false);
                        }
                    }
                    Ok(true)
                };
            loop {
                let mut job = Vec::new();
                let mut size = 0;
                let mut fallback = false;
                while size < JOB_BYTES {
                    // Header, extra field and trailer of the next member.
                    match fill(&mut buffer, &mut start, 12) {
                        Ok(true) => {}
                        Ok(false) => {
                            eof = buffer.len() == start;
                            fallback = !eof;
                            break;
                        }
                        Err(e) => {
                            let _ = jobs.send((seq, MemberJob::Failed(e)));
                            return;
                        }
                    }
                    let xlen =
                        usize::from(u16::from_le_bytes([buffer[start + 10], buffer[start + 11]]));
                    if !matches!(fill(&mut buffer, &mut start, 12 + xlen), Ok(true)) {
                        fallback = true;
                        break;
                    }
                    let Some(total) =
                        indexed_member_size(&buffer[start..]).filter(|&t| t >= 18 + 10)
                    else {
                        fallback = true;
                        break;
                    };
                    if !matches!(fill(&mut buffer, &mut start, total), Ok(true)) {
                        fallback = true;
                        break;
                    }
                    let isize = u32::from_le_bytes(
                        buffer[start + total - 4..start + total].try_into().unwrap(),
                    );
                    job.extend_from_slice(&buffer[start..start + total]);
                    size += isize as usize;
                    start += total;
                }
                if !job.is_empty() {
                    if jobs.send((seq, MemberJob::Members(job, size))).is_err() {
                        return;
                    }
                    seq += 1;
                }
                if fallback {
                    // Decode the rest sequentially, exactly as MultiGzDecoder would.
                    let rest = std::io::Cursor::new(buffer.split_off(start)).chain(file);
                    let mut decoder = MultiGzDecoder::new(BufReader::with_capacity(1 << 20, rest));
                    loop {
                        let mut chunk = vec![0; JOB_BYTES];
                        let mut filled = 0;
                        while filled < chunk.len() {
                            match decoder.read(&mut chunk[filled..]) {
                                Ok(0) => break,
                                Ok(n) => filled += n,
                                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                                Err(e) => {
                                    let _ = jobs.send((seq, MemberJob::Failed(e)));
                                    return;
                                }
                            }
                        }
                        chunk.truncate(filled);
                        if filled == 0 || jobs.send((seq, MemberJob::Plain(chunk))).is_err() {
                            return;
                        }
                        seq += 1;
                    }
                }
                if eof {
                    return;
                }
            }
        });
        // Deliver in order.
        let mut pending = BTreeMap::new();
        let mut next = 0;
        for (seq, result) in results {
            pending.insert(seq, result);
            while let Some(result) = pending.remove(&next) {
                let failed = result.is_err();
                if sender.send(result).is_err() || failed {
                    return;
                }
                next += 1;
            }
        }
        let _ = sender.send(Ok(Vec::new()));
    });
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
static SPLIT_DECODING: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Whether large single-member gzip inputs opened from now on are decoded in
/// parallel regions (`pgzip`) rather than by one sequential decoder.
pub(crate) fn set_split_decoding(enabled: bool) {
    SPLIT_DECODING.store(enabled, std::sync::atomic::Ordering::Relaxed);
}

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
            let local = std::fs::metadata(path).is_ok_and(|m| m.is_file());
            if local && indexed_member_size(buffered.fill_buf()?).is_some() {
                // Members record their sizes: decompress them in parallel.
                drop(buffered);
                let file = File::open(path)?;
                let threads = num_cpus::get().div_ceil(4).clamp(1, 4);
                let (sender, receiver) = crossbeam_channel::bounded(2);
                let (recycle, _) = crossbeam_channel::bounded::<Vec<u8>>(1);
                let worker = std::thread::spawn(move || decompress_members(file, threads, sender));
                return Ok(Box::new(PrefetchInput {
                    receiver: Some(receiver),
                    recycle,
                    worker: Some(worker),
                    buffer: Vec::new(),
                    position: 0,
                    finished: false,
                }));
            }
            if local
                && SPLIT_DECODING.load(std::sync::atomic::Ordering::Relaxed)
                && std::fs::metadata(path).is_ok_and(|m| crate::pgzip::worth_it(m.len()))
            {
                // One large member: decode regions of it in parallel.
                drop(buffered);
                let threads = num_cpus::get().div_ceil(2).clamp(2, 8);
                let (sender, receiver) = crossbeam_channel::bounded(4);
                let (recycle, _) = crossbeam_channel::bounded::<Vec<u8>>(1);
                let path = path.to_string();
                let worker =
                    std::thread::spawn(move || crate::pgzip::decompress(&path, threads, sender));
                return Ok(Box::new(PrefetchInput {
                    receiver: Some(receiver),
                    recycle,
                    worker: Some(worker),
                    buffer: Vec::new(),
                    position: 0,
                    finished: false,
                }));
            }
            let decoder = MultiGzDecoder::new(buffered);
            if local {
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

    fn read_all(path: &std::path::Path) -> std::io::Result<Vec<u8>> {
        let mut out = Vec::new();
        open_input(path.to_str().unwrap())
            .map_err(std::io::Error::other)?
            .read_to_end(&mut out)?;
        Ok(out)
    }

    /// fasterp's own gzip output, written with the parallel compressor.
    fn fasterp_gzip(data: &[u8]) -> Vec<u8> {
        let file = tempfile::NamedTempFile::new().unwrap();
        let mut writer = ParCompressBuilder::<Gzip>::new()
            .num_threads(3)
            .unwrap()
            .from_writer(file.reopen().unwrap());
        writer.write_all(data).unwrap();
        writer.finish().unwrap();
        std::fs::read(file.path()).unwrap()
    }

    fn plain_gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::new(6));
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    /// A BGZF member: `BC` subfield holding the member size minus one.
    fn bgzf_member(data: &[u8]) -> Vec<u8> {
        let plain = plain_gzip(data);
        let total = plain.len() + 8;
        let mut out = plain[..3].to_vec();
        out.push(plain[3] | 4);
        out.extend_from_slice(&plain[4..10]);
        out.extend_from_slice(&6u16.to_le_bytes());
        out.extend_from_slice(b"BC");
        out.extend_from_slice(&2u16.to_le_bytes());
        out.extend_from_slice(&((total - 1) as u16).to_le_bytes());
        out.extend_from_slice(&plain[10..]);
        out
    }

    #[test]
    fn sized_members_decompress_in_parallel_and_in_order() {
        let data: Vec<u8> = (0..5_000_000u32)
            .map(|i| b"ACGTN\n@+"[(i.wrapping_mul(2654435761) >> 29) as usize])
            .collect();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("in.fq.gz");
        let ours = fasterp_gzip(&data);
        assert!(
            indexed_member_size(&ours).is_some(),
            "output records member sizes"
        );
        for (name, bytes, expected) in [
            ("fasterp", ours.clone(), data.clone()),
            // Sized members, then a member without a size: the rest is read
            // sequentially.
            (
                "mixed",
                [ours.clone(), plain_gzip(b"tail\n"), fasterp_gzip(b"more\n")].concat(),
                [&data[..], b"tail\n", b"more\n"].concat(),
            ),
            (
                "bgzf",
                data.chunks(60_000).flat_map(bgzf_member).collect(),
                data.clone(),
            ),
            ("empty", fasterp_gzip(b""), Vec::new()),
        ] {
            std::fs::write(&path, &bytes).unwrap();
            assert_eq!(read_all(&path).unwrap(), expected, "{name}");
        }
        // A corrupt member and a truncated file are errors, not short output.
        let mut corrupt = ours.clone();
        let middle = corrupt.len() / 2;
        corrupt[middle] ^= 0xff;
        std::fs::write(&path, &corrupt).unwrap();
        assert!(read_all(&path).is_err(), "corrupt member");
        std::fs::write(&path, &ours[..ours.len() - 100]).unwrap();
        assert!(read_all(&path).is_err(), "truncated file");
    }

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
pub(crate) use crate::gzip::threads_per_output;

pub(crate) fn create_file(path: &str) -> Result<Box<dyn Write>> {
    Ok(Box::new(File::create(path)?))
}

pub(crate) fn read_to_string(path: &str) -> Result<String> {
    Ok(std::fs::read_to_string(path)?)
}

pub(crate) fn write_file(path: &str, data: &[u8]) -> Result<()> {
    Ok(std::fs::write(path, data)?)
}

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
