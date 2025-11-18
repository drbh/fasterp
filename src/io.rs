//! Input/output handling with compression support
//!
//! This module provides I/O abstraction with automatic compression detection:
//! - `CompressionFormat`: Detect compression from file extension or magic bytes
//! - `open_input()`: Open input with auto-decompression (gzip support)
//! - `OutputWriter`: Wrapper for `BufWriter` with optional compression
//! - `open_output()`: Open output with optional compression

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use gzp::ZWriter;
use gzp::deflate::Gzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

/// Check if a path is a URL
fn is_url(path: &str) -> bool {
    path.starts_with("http://") || path.starts_with("https://")
}

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

/// Open input file, URL, or stdin with automatic decompression
pub(crate) fn open_input(path: &str) -> Result<Box<dyn BufRead + Send>> {
    if path == "-" {
        // Read from stdin
        let stdin = std::io::stdin();
        let reader = BufReader::with_capacity(16 * 1024 * 1024, stdin);
        Ok(Box::new(reader))
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
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Gzip(w) => w.write(buf),
            OutputWriter::ParGzip(w) => w.write(buf),
            OutputWriter::Stdout(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Gzip(w) => w.flush(),
            OutputWriter::ParGzip(w) => w.flush(),
            OutputWriter::Stdout(w) => w.flush(),
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
        }
    }
}

/// Open output file or stdout with optional compression
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
