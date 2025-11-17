//! UMI (Unique Molecular Identifier) processing
//!
//! This module handles extraction and processing of UMI sequences from reads.
//! UMIs are short sequences (typically 6-12bp) added during library preparation
//! to uniquely tag individual molecules for PCR duplicate removal and quantification.

/// UMI processing configuration
#[derive(Debug, Clone)]
pub struct UmiConfig {
    /// Enable UMI processing
    pub enabled: bool,

    /// Where to extract UMI from
    pub location: UmiLocation,

    /// Length of UMI sequence
    pub length: usize,

    /// Prefix to add before UMI in read name (default: "UMI")
    pub prefix: String,

    /// Skip N bases before extracting UMI
    pub skip: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UmiLocation {
    /// Extract from 5' end of read1
    Read1,
    /// Extract from 5' end of read2
    Read2,
    /// Extract from index1 read (separate file)
    Index1,
    /// Extract from index2 read (separate file)
    Index2,
    /// Extract from `per_read` locations
    PerRead,
    /// Extract from `per_index` locations
    PerIndex,
}

impl Default for UmiConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            location: UmiLocation::Read1,
            length: 0,
            prefix: String::from("UMI"),
            skip: 0,
        }
    }
}

/// Result of UMI extraction
#[derive(Debug, Clone)]
pub struct UmiExtraction {
    /// The extracted UMI sequence
    pub umi_seq: Vec<u8>,
    /// Start position in original read (after skip)
    pub start_pos: usize,
    /// End position in original read
    pub end_pos: usize,
}

/// Extract UMI from read sequence
///
/// Returns None if:
/// - UMI is disabled
/// - Read is too short for UMI extraction
pub fn extract_umi(seq: &[u8], config: &UmiConfig) -> Option<UmiExtraction> {
    if !config.enabled {
        return None;
    }

    let total_umi_start = config.skip;
    let total_umi_end = config.skip + config.length;

    // Check if read is long enough
    if seq.len() < total_umi_end {
        return None;
    }

    Some(UmiExtraction {
        umi_seq: seq[total_umi_start..total_umi_end].to_vec(),
        start_pos: total_umi_start,
        end_pos: total_umi_end,
    })
}

/// Add UMI to read header
///
/// Formats header as: @`original_header:PREFIX_UMISEQUENCE`
/// Example: @read1 -> @`read1:UMI_ACGTACGT`
pub fn add_umi_to_header(header: &[u8], umi_seq: &[u8], prefix: &str) -> Vec<u8> {
    let mut new_header = header.to_vec();
    new_header.push(b':');
    new_header.extend_from_slice(prefix.as_bytes());
    new_header.push(b'_');
    new_header.extend_from_slice(umi_seq);
    new_header
}

/// Statistics for UMI processing
#[derive(Debug, Default, Clone)]
pub struct UmiStats {
    /// Total reads processed
    pub total_reads: usize,
    /// Reads with UMI extracted successfully
    pub umi_extracted: usize,
    /// Reads where UMI extraction failed (too short)
    pub umi_failed: usize,
    /// Reads removed as duplicates (if deduplication enabled)
    pub duplicates_removed: usize,
}

impl UmiStats {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn record_extraction(&mut self, success: bool) {
        self.total_reads += 1;
        if success {
            self.umi_extracted += 1;
        } else {
            self.umi_failed += 1;
        }
    }

    pub fn record_duplicate(&mut self) {
        self.duplicates_removed += 1;
    }

    pub fn merge(&mut self, other: &UmiStats) {
        self.total_reads += other.total_reads;
        self.umi_extracted += other.umi_extracted;
        self.umi_failed += other.umi_failed;
        self.duplicates_removed += other.duplicates_removed;
    }

    pub fn extraction_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            (self.umi_extracted as f64 / self.total_reads as f64) * 100.0
        }
    }

    pub fn duplication_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            (self.duplicates_removed as f64 / self.total_reads as f64) * 100.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_umi_basic() {
        let config = UmiConfig {
            enabled: true,
            location: UmiLocation::Read1,
            length: 8,
            prefix: String::from("UMI"),
            skip: 0,
        };

        let seq = b"ACGTACGTNNNNNNNNNNNN";
        let extraction = extract_umi(seq, &config).unwrap();

        assert_eq!(extraction.umi_seq, b"ACGTACGT");
        assert_eq!(extraction.start_pos, 0);
        assert_eq!(extraction.end_pos, 8);
    }

    #[test]
    fn test_extract_umi_with_skip() {
        let config = UmiConfig {
            enabled: true,
            location: UmiLocation::Read1,
            length: 6,
            prefix: String::from("UMI"),
            skip: 2,
        };

        let seq = b"NNACGTACNNNNNN";
        let extraction = extract_umi(seq, &config).unwrap();

        assert_eq!(extraction.umi_seq, b"ACGTAC");
        assert_eq!(extraction.start_pos, 2);
        assert_eq!(extraction.end_pos, 8);
    }

    #[test]
    fn test_add_umi_to_header() {
        let header = b"@read1";
        let umi = b"ACGTACGT";
        let result = add_umi_to_header(header, umi, "UMI");

        assert_eq!(result, b"@read1:UMI_ACGTACGT");
    }

    #[test]
    fn test_add_umi_to_header_custom_prefix() {
        let header = b"@seq_12345";
        let umi = b"TTGGCC";
        let result = add_umi_to_header(header, umi, "BC");

        assert_eq!(result, b"@seq_12345:BC_TTGGCC");
    }

    #[test]
    fn test_extract_umi_too_short() {
        let config = UmiConfig {
            enabled: true,
            location: UmiLocation::Read1,
            length: 10,
            prefix: String::from("UMI"),
            skip: 0,
        };

        let seq = b"ACGT"; // Only 4 bases
        let extraction = extract_umi(seq, &config);

        assert!(extraction.is_none());
    }

    #[test]
    fn test_extract_umi_disabled() {
        let config = UmiConfig {
            enabled: false,
            location: UmiLocation::Read1,
            length: 8,
            prefix: String::from("UMI"),
            skip: 0,
        };

        let seq = b"ACGTACGTNNNNNNNNNNNN";
        let extraction = extract_umi(seq, &config);

        assert!(extraction.is_none());
    }

    #[test]
    fn test_extract_umi_with_skip_too_short() {
        let config = UmiConfig {
            enabled: true,
            location: UmiLocation::Read1,
            length: 6,
            prefix: String::from("UMI"),
            skip: 10,
        };

        let seq = b"ACGTACGT"; // Only 8 bases, need skip(10) + length(6) = 16
        let extraction = extract_umi(seq, &config);

        assert!(extraction.is_none());
    }

    #[test]
    fn test_umi_stats_basic() {
        let mut stats = UmiStats::new();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.umi_extracted, 0);
        assert_eq!(stats.umi_failed, 0);

        stats.record_extraction(true);
        assert_eq!(stats.total_reads, 1);
        assert_eq!(stats.umi_extracted, 1);
        assert_eq!(stats.umi_failed, 0);

        stats.record_extraction(false);
        assert_eq!(stats.total_reads, 2);
        assert_eq!(stats.umi_extracted, 1);
        assert_eq!(stats.umi_failed, 1);
    }

    #[test]
    fn test_umi_stats_merge() {
        let mut stats1 = UmiStats {
            total_reads: 100,
            umi_extracted: 95,
            umi_failed: 5,
            duplicates_removed: 0,
        };

        let stats2 = UmiStats {
            total_reads: 50,
            umi_extracted: 48,
            umi_failed: 2,
            duplicates_removed: 0,
        };

        stats1.merge(&stats2);

        assert_eq!(stats1.total_reads, 150);
        assert_eq!(stats1.umi_extracted, 143);
        assert_eq!(stats1.umi_failed, 7);
        assert_eq!(stats1.duplicates_removed, 0);
    }

    #[test]
    fn test_umi_stats_extraction_rate() {
        let stats = UmiStats {
            total_reads: 100,
            umi_extracted: 95,
            umi_failed: 5,
            duplicates_removed: 0,
        };

        assert_eq!(stats.extraction_rate(), 95.0);
    }

    #[test]
    fn test_default_config() {
        let config = UmiConfig::default();
        assert!(!config.enabled);
        assert_eq!(config.location, UmiLocation::Read1);
        assert_eq!(config.length, 0);
        assert_eq!(config.prefix, "UMI");
        assert_eq!(config.skip, 0);
    }
}
