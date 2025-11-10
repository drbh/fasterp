//! Adapter trimming functionality
//!
//! This module provides adapter detection and trimming for FASTQ preprocessing:
//! - Manual adapter specification
//! - Built-in common adapter sequences
//! - Overlap-based adapter detection for paired-end reads
//! - Mismatch-tolerant matching

use std::cmp::min;

/// Configuration for adapter trimming
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// Adapter sequence for read1/single-end
    pub adapter_seq: Option<Vec<u8>>,
    /// Adapter sequence for read2 (paired-end only)
    pub adapter_seq_r2: Option<Vec<u8>>,
    /// Enable auto-detection for paired-end
    pub detect_adapter_for_pe: bool,
    /// Minimum overlap length for detection
    pub min_overlap: usize,
    /// Maximum mismatches allowed in overlap
    pub max_mismatches: usize,
}

impl AdapterConfig {
    pub fn new() -> Self {
        Self {
            adapter_seq: None,
            adapter_seq_r2: None,
            detect_adapter_for_pe: false,
            min_overlap: 5, // fastp's default for adapter trimming
            max_mismatches: 2,
        }
    }

    pub fn is_enabled(&self) -> bool {
        self.adapter_seq.is_some() || self.adapter_seq_r2.is_some() || self.detect_adapter_for_pe
    }
}

/// Built-in Illumina adapter sequences
pub mod adapters {
    use std::collections::HashMap;
    use std::sync::OnceLock;

    /// Illumina TruSeq Universal Adapter
    pub const TRUSEQ_UNIVERSAL: &[u8] = b"AGATCGGAAGAGC";

    /// Illumina TruSeq Read 1 Adapter
    pub const TRUSEQ_READ1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

    /// Illumina TruSeq Read 2 Adapter
    pub const TRUSEQ_READ2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

    /// Illumina Small RNA 3' Adapter
    pub const SMALL_RNA_3P: &[u8] = b"TGGAATTCTCGGGTGCCAAGG";

    /// Nextera Transposase Sequence
    pub const NEXTERA: &[u8] = b"CTGTCTCTTATACACATCT";

    /// Get the comprehensive known adapters database
    /// Returns a HashMap mapping adapter sequences to their descriptions
    pub fn get_known_adapters() -> &'static HashMap<&'static str, &'static str> {
        static KNOWN_ADAPTERS: OnceLock<HashMap<&'static str, &'static str>> = OnceLock::new();
        KNOWN_ADAPTERS.get_or_init(|| {
            let mut map = HashMap::new();
            map.insert("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", ">Illumina TruSeq Adapter Read 1");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", ">Illumina TruSeq Adapter Read 2");
            map.insert("GATCGTCGGACTGTAGAACTCTGAACGTGTAGA", ">Illumina Small RNA Adapter Read 2");
            map.insert("AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA", ">Illumina DpnII expression PCR Primer 2 | >Illumina NlaIII expression PCR Primer 2 | >Illumina Small RNA PCR Primer 2 | >Illumina DpnII Gex PCR Primer 2 | >Illumina NlaIII Gex PCR Primer 2");
            map.insert("AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA", ">Illumina RNA PCR Primer");
            map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1 | >Illumina Single End PCR Primer 1 | >Illumina Paried End PCR Primer 1 | >Illumina Multiplexing PCR Primer 1.01 | >TruSeq Universal Adapter | >TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1");
            map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">pcr_dimer");
            map.insert("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT", ">PCR_Primers");
            map.insert("ACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">Illumina Single End Sequencing Primer | >Illumina Paired End Adapter 1 | >Illumina Paried End Sequencing Primer 1 | >Illumina Multiplexing Adapter 2 | >Illumina Multiplexing Read1 Sequencing Primer");
            map.insert("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", ">PE2_rc | >TruSeq3_IndexedAdapter | >PE2_rc | >TruSeq3_IndexedAdapter");
            map.insert("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG", ">Reverse_adapter");
            map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", ">TruSeq2_PE_r");
            map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">PCR_Primer2_rc");
            map.insert("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA", ">PhiX_read1_adapter");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA", ">PE1_rc | >TruSeq3_UniversalAdapter | >PE1_rc | >TruSeq3_UniversalAdapter");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">PCR_Primer1_rc");
            map.insert("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA", ">PhiX_read2_adapter");
            map.insert("AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", ">TruSeq2_SE");
            map.insert("CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 35");
            map.insert("CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC", ">Illumina PCR Primer Index 10");
            map.insert("CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 10");
            map.insert("CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC", ">Illumina PCR Primer Index 2");
            map.insert("CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 2");
            map.insert("CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 38");
            map.insert("CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 27");
            map.insert("CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 25");
            map.insert("CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 31");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 44");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 37");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC", ">Illumina PCR Primer Index 6");
            map.insert("CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 6");
            map.insert("CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC", ">Illumina PCR Primer Index 5");
            map.insert("CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 5");
            map.insert("CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 23");
            map.insert("CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 30");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 21");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 42");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 33");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT", ">PrefixPE/2 | >PCR_Primer2 | >Illumina Paired End PCR Primer 2 | >PrefixPE/2 | >PCR_Primer2");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 22");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC", ">Illumina PCR Primer Index 1");
            map.insert("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 1");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 17");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC", ">Illumina PCR Primer Index 9");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 9");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 47");
            map.insert("CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 28");
            map.insert("CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 45");
            map.insert("CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC", ">Illumina PCR Primer Index 7");
            map.insert("CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 7");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 34");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC", ">Illumina PCR Primer Index 3");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 3");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 18");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 24");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 26");
            map.insert("CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 43");
            map.insert("CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 14");
            map.insert("CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 16");
            map.insert("CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 20");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC", ">Illumina PCR Primer Index 11");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 11");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 39");
            map.insert("CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 41");
            map.insert("CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC", ">Illumina PCR Primer Index 12");
            map.insert("CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 12");
            map.insert("CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 29");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC", ">Illumina PCR Primer Index 8");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 8");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 46");
            map.insert("CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 40");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 15");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 32");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 48");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC", ">Illumina PCR Primer Index 4");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 4");
            map.insert("CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 36");
            map.insert("CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 13");
            map.insert("CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", ">RNA PCR Primer, Index 19");
            map.insert("CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT", ">Illumina Single End Adapter 2 | >Illumina Single End PCR Primer 2");
            map.insert("CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT", ">ABI Solid3 Adapter B");
            map.insert("CCGACAGGTTCAGAGTTCTACAGTCCGACATG", ">Illumina NlaIII expression Sequencing Primer | >Illumina NlaIII Gex Sequencing Primer");
            map.insert("CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | >I7_Primer_Nextera_XT_Index_Kit_v2_N711 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | >I7_Primer_Nextera_XT_Index_Kit_v2_N711");
            map.insert("CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N716");
            map.insert("CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N724");
            map.insert("CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | >I7_Primer_Nextera_XT_Index_Kit_v2_N703 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | >I7_Primer_Nextera_XT_Index_Kit_v2_N703");
            map.insert("CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N715");
            map.insert("CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N722");
            map.insert("CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N708");
            map.insert("CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N726");
            map.insert("CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | >I7_Primer_Nextera_XT_Index_Kit_v2_N710 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | >I7_Primer_Nextera_XT_Index_Kit_v2_N710");
            map.insert("CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N727");
            map.insert("CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N720");
            map.insert("CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | >I7_Primer_Nextera_XT_Index_Kit_v2_N702 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | >I7_Primer_Nextera_XT_Index_Kit_v2_N702");
            map.insert("CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | >I7_Primer_Nextera_XT_Index_Kit_v2_N707 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | >I7_Primer_Nextera_XT_Index_Kit_v2_N707");
            map.insert("CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N719");
            map.insert("CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N709");
            map.insert("CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N714");
            map.insert("CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | >I7_Primer_Nextera_XT_Index_Kit_v2_N705 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | >I7_Primer_Nextera_XT_Index_Kit_v2_N705");
            map.insert("CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N718");
            map.insert("CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | >I7_Primer_Nextera_XT_Index_Kit_v2_N712 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | >I7_Primer_Nextera_XT_Index_Kit_v2_N712");
            map.insert("CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | >I7_Primer_Nextera_XT_Index_Kit_v2_N701 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | >I7_Primer_Nextera_XT_Index_Kit_v2_N701");
            map.insert("CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N721");
            map.insert("CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N723");
            map.insert("CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | >I7_Primer_Nextera_XT_Index_Kit_v2_N706 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | >I7_Primer_Nextera_XT_Index_Kit_v2_N706");
            map.insert("CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | >I7_Primer_Nextera_XT_Index_Kit_v2_N704 | >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | >I7_Primer_Nextera_XT_Index_Kit_v2_N704");
            map.insert("CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N729");
            map.insert("CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG", ">I7_Primer_Nextera_XT_Index_Kit_v2_N728");
            map.insert("CGACAGGTTCAGAGTTCTACAGTCCGACGATC", ">Illumina DpnII expression Sequencing Primer | >Illumina Small RNA Sequencing Primer | >Illumina DpnII Gex Sequencing Primer");
            map.insert("CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT", ">Illumina Paired End Sequencing Primer 2");
            map.insert("CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT", ">Clontech Universal Primer Mix Long");
            map.insert("CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG", ">I7_Adapter_Nextera_No_Barcode");
            map.insert("CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Adapter_Nextera");
            map.insert("CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG", ">ABI Solid3 Adapter A");
            map.insert("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", ">I7_Nextera_Transposase_1 | >Trans2_rc | >I7_Nextera_Transposase_1 | >Trans2_rc");
            map.insert("CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC", ">I7_Nextera_Transposase_2");
            map.insert("CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC", ">I5_Nextera_Transposase_2");
            map.insert("CTGTCTCTTATACACATCTGACGCTGCCGACGA", ">I5_Nextera_Transposase_1 | >Trans1_rc | >I5_Nextera_Transposase_1 | >Trans1_rc");
            map.insert("GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S516");
            map.insert("GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | >I5_Primer_Nextera_XT_Index_Kit_v2_S503 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | >I5_Primer_Nextera_XT_Index_Kit_v2_S503");
            map.insert("GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S515");
            map.insert("GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | >I5_Primer_Nextera_XT_Index_Kit_v2_S508 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | >I5_Primer_Nextera_XT_Index_Kit_v2_S508");
            map.insert("GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | >I5_Primer_Nextera_XT_Index_Kit_v2_S502 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | >I5_Primer_Nextera_XT_Index_Kit_v2_S502");
            map.insert("GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S520");
            map.insert("GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S510");
            map.insert("GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S511");
            map.insert("GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S513");
            map.insert("GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | >I5_Primer_Nextera_XT_Index_Kit_v2_S505 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | >I5_Primer_Nextera_XT_Index_Kit_v2_S505");
            map.insert("GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S518");
            map.insert("GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]501");
            map.insert("GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S521");
            map.insert("GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | >I5_Primer_Nextera_XT_Index_Kit_v2_S507 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | >I5_Primer_Nextera_XT_Index_Kit_v2_S507");
            map.insert("GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | >I5_Primer_Nextera_XT_Index_Kit_v2_S506 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | >I5_Primer_Nextera_XT_Index_Kit_v2_S506");
            map.insert("GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_Index_Kit_v2_S522");
            map.insert("GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]504");
            map.insert("GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT", ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | >I5_Primer_Nextera_XT_Index_Kit_v2_S517 | >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | >I5_Primer_Nextera_XT_Index_Kit_v2_S517");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", ">Nextera_LMP_Read1_External_Adapter | >Illumina Multiplexing Index Sequencing Primer");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_5 | >TruSeq Adapter, Index 5");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_25");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 25");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_8 | >TruSeq Adapter, Index 8");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_13");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 13");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_14");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 14");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_1_6 | >TruSeq Adapter, Index 1");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_15");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 15");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_27");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 27");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_7 | >TruSeq Adapter, Index 7");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 23");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_16");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 16");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_2 | >TruSeq Adapter, Index 2");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_22");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 22");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_12 | >TruSeq Adapter, Index 12");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_23");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_9 | >TruSeq Adapter, Index 9");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_6 | >TruSeq Adapter, Index 6");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_11 | >TruSeq Adapter, Index 11");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_18_7");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 18");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_19");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 19");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_20");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 20");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_21");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG", ">TruSeq Adapter, Index 21");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_10 | >TruSeq Adapter, Index 10");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_4 | >TruSeq Adapter, Index 4");
            map.insert("GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG", ">TruSeq_Adapter_Index_3 | >TruSeq Adapter, Index 3");
            map.insert("GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", ">Illumina Paired End Adapter 2");
            map.insert("GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", ">Nextera_LMP_Read2_External_Adapter");
            map.insert("GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", ">Illumina Single End Adapter 1");
            map.insert("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", ">Trans2");
            map.insert("GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", ">PrefixPE/2 | >PE2 | >Illumina Multiplexing PCR Primer 2.01 | >Illumina Multiplexing Read2 Sequencing Primer | >PrefixPE/2 | >PE2");
            map.insert("TACACTCTTTCCCTACACGACGCTCTTCCGATCT", ">PrefixPE/1 | >PE1 | >PrefixPE/1 | >PE1");
            map.insert("TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT", ">RNA_PCR_Primer_(RP1)_part_#_15013198");
            map.insert("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", ">Trans1");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_5_(RPI5)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_25_(RPI25)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_8_(RPI8)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_13_(RPI13)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_14_(RPI14)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_1_(RPI1)_2,9");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_26_(RPI26)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_15_(RPI15)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_27_(RPI27)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_28_(RPI28)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_29_(RPI29)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_30_(RPI30)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_31_(RPI31)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_32_(RPI32)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_7_(RPI7)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_33_(RPI33)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_34_(RPI34)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_35_(RPI35)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_36_(RPI36)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_16_(RPI16)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_2_(RPI2)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_37_(RPI37)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_22_(RPI22)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_38_(RPI38)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_39_(RPI39)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_40_(RPI40)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_12_(RPI12)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_41_(RPI41)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_23_(RPI23)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_9_(RPI9)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_6_(RPI6)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_11_(RPI11)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_24_(RPI24)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_17_(RPI17)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_18_(RPI18)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_19_(RPI19)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_20_(RPI20)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_21_(RPI21)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_42_(RPI42)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_43_(RPI43)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_10_(RPI10)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_44_(RPI44)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_45_(RPI45)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_46_(RPI46)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_47_(RPI47)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_48_(RPI48)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_4_(RPI4)");
            map.insert("TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG", ">RNA_PCR_Primer_Index_3_(RPI3)");
            map.insert("TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC", ">FlowCell1");
            map.insert("TTTTTTTTTTCAAGCAGAAGACGGCATACGA", ">FlowCell2");
            map.insert("AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA", ">MGI/BGI adapter (forward)");
            map.insert("AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG", ">MGI/BGI adapter (reverse)");
            map.insert("AACTGTAGGCACCATCAAT", ">QIASeq miRNA adapter");
            map
        })
    }
}

/// Type of adapter match found
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatchType {
    /// Exact match with allowed mismatches
    Exact,
    /// Match with single insertion in read sequence
    Insertion,
    /// Match with single deletion from adapter sequence
    Deletion,
}

/// Result of adapter detection
#[derive(Debug, Clone)]
pub struct AdapterMatch {
    /// Position where adapter starts
    pub position: usize,
    /// Number of matching bases
    pub matched_bases: usize,
    /// Number of mismatches
    pub mismatches: usize,
    /// Type of match found
    pub match_type: MatchType,
}

/// Try exact matching with allowed mismatches (Stage 1)
/// Supports negative start_pos for A-tailing (Illumina adapter dimer handling)
fn try_exact_match(
    seq: &[u8],
    adapter: &[u8],
    start_pos: isize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    let (seq_offset, adapter_offset, compare_len) = if start_pos < 0 {
        // Negative position: skip first abs(start_pos) bytes of adapter
        let skip = (-start_pos) as usize;
        if skip >= adapter.len() {
            return None;
        }
        let remaining_adapter = adapter.len() - skip;
        let compare_len = min(seq.len(), remaining_adapter);
        (0, skip, compare_len)
    } else {
        // Positive position: normal matching
        let start = start_pos as usize;
        if start >= seq.len() {
            return None;
        }
        let remaining = seq.len() - start;
        let compare_len = min(remaining, adapter.len());
        (start, 0, compare_len)
    };

    if compare_len < min_overlap {
        return None;
    }

    let allowed_mismatches = compare_len / 8;
    let mut matches = 0;
    let mut mismatches = 0;

    for i in 0..compare_len {
        let seq_base = seq[seq_offset + i].to_ascii_uppercase();
        let adapter_base = adapter[adapter_offset + i].to_ascii_uppercase();

        if seq_base == adapter_base {
            matches += 1;
        } else {
            mismatches += 1;
            if mismatches > allowed_mismatches {
                return None;
            }
        }
    }

    if mismatches <= allowed_mismatches && matches + mismatches >= min_overlap {
        Some(AdapterMatch {
            position: start_pos.max(0) as usize,
            matched_bases: matches,
            mismatches,
            match_type: MatchType::Exact,
        })
    } else {
        None
    }
}

/// Try matching with single insertion (Stage 2)
/// Implements fastp's Matcher::matchWithOneInsertion algorithm exactly
///
/// ins_data: sequence with suspected insertion (longer, e.g., read)
/// normal_data: reference sequence (baseline, e.g., adapter)
/// cmplen: comparison length (calculated by caller based on insertion/deletion case)
fn try_insertion_match(
    ins_data: &[u8],      // Sequence with insertion (read from start_pos)
    normal_data: &[u8],   // Reference sequence (adapter or read, depending on case)
    ins_start_pos: usize, // Starting position in ins_data
    cmplen: usize,        // Length to compare
    min_overlap: usize,
) -> Option<AdapterMatch> {
    let remaining = ins_data.len() - ins_start_pos;

    // Need at least cmplen + 1 bases to have an insertion
    if remaining < cmplen + 1 || cmplen < min_overlap {
        return None;
    }

    // Stricter threshold for indel matching: cmplen/8 - 1
    // In fastp, when cmplen < 8, this becomes negative, preventing matches
    let diff_limit = if cmplen >= 8 {
        (cmplen / 8).saturating_sub(1)
    } else {
        // For cmplen < 8, fastp gets -1, which fails all comparisons
        // We simulate this by returning early
        return None;
    };

    // Arrays of size cmplen (matching fastp exactly)
    // Initialize both to high values to prevent false matches from uncomputed positions
    let mut acc_mismatch_from_left = vec![diff_limit + 1; cmplen];
    let mut acc_mismatch_from_right = vec![diff_limit + 1; cmplen];

    // Initialize first and last elements
    acc_mismatch_from_left[0] =
        if ins_data[ins_start_pos].to_ascii_uppercase() == normal_data[0].to_ascii_uppercase() {
            0
        } else {
            1
        };

    if ins_start_pos + cmplen < ins_data.len() {
        acc_mismatch_from_right[cmplen - 1] = if ins_data[ins_start_pos + cmplen]
            .to_ascii_uppercase()
            == normal_data[cmplen - 1].to_ascii_uppercase()
        {
            0
        } else {
            1
        };
    } else {
        return None;
    }

    // Build left array with early termination
    for i in 1..cmplen {
        if ins_start_pos + i >= ins_data.len() {
            return None;
        }

        if ins_data[ins_start_pos + i].to_ascii_uppercase() != normal_data[i].to_ascii_uppercase() {
            acc_mismatch_from_left[i] = acc_mismatch_from_left[i - 1] + 1;
        } else {
            acc_mismatch_from_left[i] = acc_mismatch_from_left[i - 1];
        }

        // Early termination: if left + rightmost already exceeds limit, stop
        if acc_mismatch_from_left[i] + acc_mismatch_from_right[cmplen - 1] > diff_limit {
            break;
        }
    }

    // Build right array with early termination
    for i in (0..cmplen - 1).rev() {
        if ins_start_pos + i + 1 >= ins_data.len() {
            continue;
        }

        if ins_data[ins_start_pos + i + 1].to_ascii_uppercase()
            != normal_data[i].to_ascii_uppercase()
        {
            acc_mismatch_from_right[i] = acc_mismatch_from_right[i + 1] + 1;
        } else {
            acc_mismatch_from_right[i] = acc_mismatch_from_right[i + 1];
        }

        // Early termination: if right + leftmost exceeds limit
        if acc_mismatch_from_right[i] + acc_mismatch_from_left[0] > diff_limit {
            // Set all remaining positions to high value
            for p in 0..i {
                acc_mismatch_from_right[p] = diff_limit + 1;
            }
            break;
        }
    }

    // Check each potential skip position
    for i in 1..cmplen {
        // Early termination check
        if acc_mismatch_from_left[i - 1] + acc_mismatch_from_right[cmplen - 1] > diff_limit {
            return None;
        }

        let diff = acc_mismatch_from_left[i - 1] + acc_mismatch_from_right[i];
        if diff <= diff_limit {
            let matches = cmplen - diff;
            return Some(AdapterMatch {
                position: ins_start_pos,
                matched_bases: matches,
                mismatches: diff,
                match_type: MatchType::Insertion,
            });
        }
    }

    None
}

/// Try matching with single deletion from adapter (Stage 3)
/// Fastp handles deletion by swapping arguments to matchWithOneInsertion:
/// A deletion in adapter = insertion in adapter sequence relative to read
fn try_deletion_match(
    seq: &[u8],
    adapter: &[u8],
    start_pos: usize,
    min_overlap: usize,
) -> Option<AdapterMatch> {
    let remaining = seq.len() - start_pos;

    // Fastp uses: cmplen = min(rlen - pos, alen - 1)
    let cmplen = min(remaining, adapter.len().saturating_sub(1));

    if cmplen < min_overlap {
        return None;
    }

    // Swap arguments: treat adapter as having "insertion" relative to seq
    // Pass adapter as ins_data, seq (from start_pos) as normal_data
    if let Some(mut result) =
        try_insertion_match(adapter, &seq[start_pos..], 0, cmplen, min_overlap)
    {
        // Adjust position back to original sequence coordinates
        result.position = start_pos;
        result.match_type = MatchType::Deletion;
        Some(result)
    } else {
        None
    }
}

/// Check if new match is better than current best
fn is_better_match(new_match: &AdapterMatch, current_best: &Option<AdapterMatch>) -> bool {
    match current_best {
        None => true,
        Some(best) => {
            // Prefer earlier position
            if new_match.position < best.position {
                return true;
            }
            if new_match.position > best.position {
                return false;
            }

            // Same position: prefer exact > deletion > insertion
            use MatchType::*;
            match (&new_match.match_type, &best.match_type) {
                (Exact, Insertion) | (Exact, Deletion) => true,
                (Deletion, Insertion) => true,
                (Insertion, Exact) | (Deletion, Exact) | (Insertion, Deletion) => false,
                _ => new_match.mismatches < best.mismatches,
            }
        }
    }
}

/// Find adapter sequence in read allowing mismatches
///
/// Uses a three-stage matching approach (exact, insertion, deletion).
/// Returns the best match (earliest position with fewest mismatches).
pub fn find_adapter(
    seq: &[u8],
    adapter: &[u8],
    min_overlap: usize,
    max_mismatches: usize,
) -> Option<AdapterMatch> {
    if adapter.is_empty() || seq.len() < min_overlap {
        return None;
    }

    let mut best_match: Option<AdapterMatch> = None;

    // Calculate start position for A-tailing support
    // Fastp starts from negative positions to handle Illumina adapter dimers
    let start: isize = if adapter.len() >= 16 {
        -4
    } else if adapter.len() >= 12 {
        -3
    } else if adapter.len() >= 8 {
        -2
    } else {
        0
    };

    // Minimum required match length (matchReq in fastp)
    let match_req = 4;

    // Try all possible positions using three-stage matching
    // Fastp uses: Stage 1 (exact), Stage 2 (insertion), Stage 3 (deletion)
    let end = (seq.len() as isize) - (match_req as isize);
    for pos in start..end {
        // Stage 1: Try exact matching with mismatches
        if let Some(match_result) = try_exact_match(seq, adapter, pos, min_overlap) {
            if is_better_match(&match_result, &best_match) {
                best_match = Some(match_result);
            }
            continue; // Exact match found, skip indel stages
        }

        // Indel matching only for non-negative positions
        if pos < 0 {
            continue;
        }

        let start_pos = pos as usize;
        let remaining = seq.len() - start_pos;

        if remaining < min_overlap {
            break;
        }

        // Stage 2 & 3: Indel matching ONLY at position 0
        // Fastp only runs indel matching at the start of the sequence
        if start_pos == 0 {
            let cmplen_insertion = min(remaining.saturating_sub(1), adapter.len());
            if let Some(match_result) =
                try_insertion_match(seq, adapter, 0, cmplen_insertion, min_overlap)
            {
                if is_better_match(&match_result, &best_match) {
                    best_match = Some(match_result);
                }
                continue; // Insertion match found, skip deletion stage
            }

            if let Some(match_result) = try_deletion_match(seq, adapter, 0, min_overlap) {
                if is_better_match(&match_result, &best_match) {
                    best_match = Some(match_result);
                }
            }
        }
    }

    best_match
}

/// Trim adapter from sequence and quality
///
/// Returns (trimmed_seq, trimmed_qual) as slices
pub fn trim_adapter<'a>(
    seq: &'a [u8],
    qual: &'a [u8],
    adapter_match: &AdapterMatch,
) -> (&'a [u8], &'a [u8]) {
    // Trim everything from the adapter position onwards
    let trim_pos = adapter_match.position;
    (&seq[..trim_pos], &qual[..trim_pos])
}

/// Detect adapter using paired-end overlap information
///
/// When read1 and read2 overlap, we can detect adapter contamination
/// by checking for reverse-complement matching beyond the insert size.
pub fn detect_adapter_from_pe_overlap(
    _seq1: &[u8],
    _seq2: &[u8],
    _min_overlap: usize,
) -> (Option<AdapterMatch>, Option<AdapterMatch>) {
    // TODO: Implement paired-end overlap-based adapter detection
    // This requires:
    // 1. Find overlap between read1 and reverse-complement of read2
    // 2. Determine insert size
    // 3. Identify adapter sequences beyond the insert

    // For now, return None (not implemented)
    (None, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_adapter_match() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGC";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 12); // Adapter starts at position 12
        assert_eq!(m.mismatches, 0);
    }

    #[test]
    fn test_adapter_with_mismatch() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGX"; // X instead of C
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 12);
        assert_eq!(m.mismatches, 1);
    }

    #[test]
    fn test_partial_adapter_at_end() {
        // Only partial adapter at the end of read
        let seq = b"ACGTACGTACGTACGTAGATCGGAA"; // Only first 9 bases of adapter
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        // Should not match with min_overlap=10
        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());

        // Should match with min_overlap=8
        let result = find_adapter(seq, adapter, 8, 2);
        assert!(result.is_some());
    }

    #[test]
    fn test_no_adapter() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());
    }

    #[test]
    fn test_trim_adapter() {
        let seq = b"ACGTACGTACGTAGATCGGAAGAGC";
        let qual = b"############IIIIIIIIIIIII";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let adapter_match = find_adapter(seq, adapter, 10, 2).unwrap();
        let (trimmed_seq, trimmed_qual) = trim_adapter(seq, qual, &adapter_match);

        assert_eq!(trimmed_seq, b"ACGTACGTACGT");
        assert_eq!(trimmed_qual, b"############");
    }

    #[test]
    fn test_adapter_too_many_mismatches() {
        // 3 mismatches (XXX), but max_mismatches=2
        let seq = b"ACGTACGTACGTAGATCGGAAXXXC";
        let adapter = adapters::TRUSEQ_UNIVERSAL;

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_none());
    }

    #[test]
    fn test_case_insensitive_matching() {
        let seq = b"ACGTACGTACGTagatcggaagagc"; // lowercase adapter
        let adapter = adapters::TRUSEQ_UNIVERSAL; // uppercase

        let result = find_adapter(seq, adapter, 10, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.mismatches, 0);
    }

    #[test]
    fn test_insertion_in_read() {
        // Sequence has extra base (X) that needs to be skipped
        let seq = b"AGATCGGXAAGAGCTTTTTTTT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 5, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.position, 0);
        assert_eq!(m.match_type, MatchType::Insertion);
    }

    #[test]
    fn test_deletion_from_adapter() {
        // seq.2033 case: AGATCGGAGCTC vs AGATCGGAAGAGC (missing 'A' at position 8)
        // This has too many mismatches (diff_limit = cmplen/8 - 1 = 0 for cmplen=12)
        // Fastp does NOT match this at position 0 - verified by testing
        let seq = b"AGATCGGAGCTCACGGATCAGGTGAAT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 5, 2);

        // Fastp finds a match at position 18, not position 0
        // The beginning doesn't match due to strict deletion matching criteria
        if let Some(m) = result {
            assert_ne!(m.position, 0, "Should not match at position 0 with deletion logic");
        }
    }

    #[test]
    fn test_exact_match_preferred_over_indel() {
        // When both exact and indel match exist, prefer exact
        let seq = b"AGATCGGAAGAGCTTTTTTTT";
        let adapter = b"AGATCGGAAGAGC";

        let result = find_adapter(seq, adapter, 5, 2);
        assert!(result.is_some());

        let m = result.unwrap();
        assert_eq!(m.match_type, MatchType::Exact);
    }

    // #[test]
    // fn test_earlier_position_preferred() {
    //     // Earlier position should be preferred over later position
    //     let seq = b"AAAAAAGATCGGAAGAGCXXXXAGATCGGAAGAGC";
    //     let adapter = b"AGATCGGAAGAGC";

    //     let result = find_adapter(seq, adapter, 5, 2);
    //     assert!(result.is_some());

    //     let m = result.unwrap();
    //     assert_eq!(m.position, 6); // First occurrence
    // }
}
