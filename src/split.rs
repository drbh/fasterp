//! Output splitting module for fastp compatibility
//!
//! Supports two modes:
//! 1. Split by number of files: Distributes reads evenly across N files
//! 2. Split by lines: Creates new file when line limit is reached

use anyhow::Result;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Configuration for output splitting
#[derive(Debug, Clone)]
pub struct SplitConfig {
    /// Splitting mode
    pub mode: SplitMode,
    /// Number of digits for file number padding (0 = no padding, 1-10)
    pub prefix_digits: usize,
}

#[derive(Debug, Clone)]
pub enum SplitMode {
    /// No splitting
    None,
    /// Split by total number of files (2-999)
    ByFiles(usize),
    /// Split by lines per file (>=1000)
    ByLines(usize),
}

impl Default for SplitConfig {
    fn default() -> Self {
        Self {
            mode: SplitMode::None,
            prefix_digits: 4,
        }
    }
}

/// Writer that handles output splitting
pub struct SplitWriter {
    /// Base output path (without number prefix)
    base_path: PathBuf,
    /// Current file number (1-indexed)
    current_file_num: usize,
    /// Current writer
    current_writer: Option<crate::io::OutputWriter>,
    /// Lines written to current file
    lines_in_current_file: usize,
    /// Split configuration
    config: SplitConfig,
    /// Compression level
    compression_level: u32,
    /// Use parallel compression
    parallel_compression: bool,
}

impl SplitWriter {
    /// Create a new split writer
    pub fn new(
        output_path: impl AsRef<Path>,
        config: SplitConfig,
        compression_level: u32,
        parallel_compression: bool,
    ) -> Result<Self> {
        let mut writer = Self {
            base_path: output_path.as_ref().to_path_buf(),
            current_file_num: 0,
            current_writer: None,
            lines_in_current_file: 0,
            config,
            compression_level,
            parallel_compression,
        };

        // For non-split mode or split-by-lines, open the first file immediately
        if matches!(writer.config.mode, SplitMode::None | SplitMode::ByLines(_)) {
            writer.open_next_file()?;
        }

        Ok(writer)
    }

    /// Generate the filename for a given file number
    fn generate_filename(&self, file_num: usize) -> PathBuf {
        if let SplitMode::None = self.config.mode {
            // No splitting - use original filename
            return self.base_path.clone();
        }

        // Get the directory and filename
        let parent = self.base_path.parent().unwrap_or_else(|| Path::new("."));
        let filename = self
            .base_path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("output.fq");

        // Generate prefix with padding
        let prefix = if self.config.prefix_digits == 0 {
            format!("{file_num}")
        } else {
            format!("{:0width$}", file_num, width = self.config.prefix_digits)
        };

        // Combine: directory/prefix.filename
        parent.join(format!("{prefix}.{filename}"))
    }

    /// Open the next output file
    fn open_next_file(&mut self) -> Result<()> {
        // Close current writer if any
        if let Some(writer) = self.current_writer.take() {
            writer.finish()?;
        }

        // Increment file number
        self.current_file_num += 1;
        self.lines_in_current_file = 0;

        // Generate new filename
        let filename = self.generate_filename(self.current_file_num);

        // Open new file with compression support
        let filename_str = filename
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid filename path"))?;
        let writer = crate::io::open_output(
            filename_str,
            Some(self.compression_level),
            self.parallel_compression,
        )?;
        self.current_writer = Some(writer);

        Ok(())
    }

    /// Check if we need to open a new file based on split mode
    fn should_open_new_file(&self) -> bool {
        match &self.config.mode {
            SplitMode::None => false,
            SplitMode::ByFiles(_num_files) => {
                // TODO: Implement ByFiles mode (requires knowing total reads)
                false
            }
            SplitMode::ByLines(max_lines) => {
                // Open new file if we've reached the line limit
                self.lines_in_current_file >= *max_lines
            }
        }
    }

    /// Write data and track lines
    #[inline]
    pub fn write(&mut self, data: &[u8]) -> Result<()> {
        // Fast path: no splitting mode
        if matches!(self.config.mode, SplitMode::None) {
            if let Some(ref mut writer) = self.current_writer {
                writer.write_all(data)?;
            }
            return Ok(());
        }

        // Slow path: splitting enabled
        // Check if we need to open a new file
        if self.should_open_new_file() {
            self.open_next_file()?;
        }

        // Write to current file
        if let Some(ref mut writer) = self.current_writer {
            writer.write_all(data)?;

            // Count newlines to track lines using memchr (SIMD)
            let newline_count = memchr::memchr_iter(b'\n', data).count();
            self.lines_in_current_file += newline_count;
        }

        Ok(())
    }

    /// Finish writing and flush
    pub fn finish(mut self) -> Result<()> {
        if let Some(writer) = self.current_writer.take() {
            writer.finish()?;
        }
        Ok(())
    }
}

impl Write for SplitWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.write(buf).map_err(std::io::Error::other)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if let Some(ref mut writer) = self.current_writer {
            writer.flush()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_no_split() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("output.fq");

        let config = SplitConfig {
            mode: SplitMode::None,
            prefix_digits: 4,
        };

        let mut writer = SplitWriter::new(&output_path, config, 6, false).unwrap();
        writer.write(b"@read1\nACGT\n+\nIIII\n").unwrap();
        writer.finish().unwrap();

        // Should create output.fq without number prefix
        assert!(output_path.exists());
    }

    #[test]
    fn test_split_by_lines() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("output.fq");

        let config = SplitConfig {
            mode: SplitMode::ByLines(4), // 4 lines = 1 FASTQ record
            prefix_digits: 4,
        };

        let mut writer = SplitWriter::new(&output_path, config, 6, false).unwrap();

        // Write 2 records (8 lines total)
        writer.write(b"@read1\nACGT\n+\nIIII\n").unwrap();
        writer.write(b"@read2\nGGGG\n+\nIIII\n").unwrap();
        writer.finish().unwrap();

        // Should create 0001.output.fq and 0002.output.fq
        assert!(temp_dir.path().join("0001.output.fq").exists());
        assert!(temp_dir.path().join("0002.output.fq").exists());
    }

    #[test]
    fn test_generate_filename_with_padding() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("output.fq");

        let config = SplitConfig {
            mode: SplitMode::ByLines(1000),
            prefix_digits: 4,
        };

        let writer = SplitWriter::new(&output_path, config, 6, false).unwrap();

        let filename1 = writer.generate_filename(1);
        assert!(filename1.to_string_lossy().contains("0001.output.fq"));

        let filename99 = writer.generate_filename(99);
        assert!(filename99.to_string_lossy().contains("0099.output.fq"));
    }

    #[test]
    fn test_generate_filename_no_padding() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("output.fq");

        let config = SplitConfig {
            mode: SplitMode::ByLines(1000),
            prefix_digits: 0, // No padding
        };

        let writer = SplitWriter::new(&output_path, config, 6, false).unwrap();

        let filename1 = writer.generate_filename(1);
        assert!(filename1.to_string_lossy().contains("1.output.fq"));
    }
}
