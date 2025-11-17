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
    fn test_default_config() {
        let config = UmiConfig::default();
        assert!(!config.enabled);
        assert_eq!(config.location, UmiLocation::Read1);
        assert_eq!(config.length, 0);
        assert_eq!(config.prefix, "UMI");
        assert_eq!(config.skip, 0);
    }
}
