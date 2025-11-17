//! Input/output handling with compression support
//!
//! This module provides I/O abstraction with automatic compression detection:
//! - CompressionFormat: Detect compression from file extension or magic bytes
//! - open_input(): Open input with auto-decompression (gzip support)
//! - OutputWriter: Wrapper for BufWriter with optional compression
//! - open_output(): Open output with optional compression

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use gzp::ZWriter;
use gzp::deflate::Gzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Detect compression format from file extension or magic bytes
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CompressionFormat {
    None,
    Gzip,
}

impl CompressionFormat {
    /// Detect from file path
    pub(crate) fn from_path(path: &str) -> Self {
        if path == "-" {
            return CompressionFormat::None; // stdin/stdout defaults to uncompressed
        }

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

/// Open input file or stdin with automatic decompression
pub(crate) fn open_input(path: &str) -> Result<Box<dyn BufRead + Send>> {
    if path == "-" {
        // Read from stdin
        let stdin = std::io::stdin();
        let reader = BufReader::with_capacity(16 * 1024 * 1024, stdin);
        Ok(Box::new(reader))
    } else {
        // Open file and detect compression
        let file = File::open(path).context(format!("Failed to open input file: {path}"))?;
        let format = CompressionFormat::from_path(path);

        match format {
            CompressionFormat::Gzip => {
                let decoder = GzDecoder::new(file);
                // Use same 16MB buffer as uncompressed - decompression is fast enough
                let reader = BufReader::with_capacity(16 * 1024 * 1024, decoder);
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
                    let num_threads = (num_cpus::get() / 2).max(4).min(8);
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
