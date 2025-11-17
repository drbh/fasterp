//! Unit tests for trimming functionality
//!
//! Tests cover:
//! - Sliding window quality trimming (front and tail)
//! - `PolyG` tail detection
//! - `PolyX` tail detection
//! - Full trimming integration
//! - Statistics accumulation

use super::*;
use crate::trimming::{
    TrimmingResult, TrimmingStats, detect_poly_g_tail, detect_poly_x_tail,
    trim_front_sliding_window, trim_read, trim_tail_sliding_window,
};

// Sliding Window Trimming Tests

#[test]
fn test_trim_tail_no_trimming_needed() {
    // High quality throughout
    let qual = b"IIIIIIIIII"; // Quality 40
    let result = trim_tail_sliding_window(qual, 4, 20);
    assert_eq!(result, 10, "Should not trim high quality read");
}

#[test]
fn test_trim_tail_low_quality_end() {
    // Quality drops at position 6
    let qual = b"IIIIII####"; // Quality 40 then 2
    let result = trim_tail_sliding_window(qual, 4, 20);
    // Window at pos 10 (len=10): #### mean=2 < 20
    // Window at pos 9: I### mean=11.5 < 20
    // Window at pos 8: II## mean=21 >= 20
    assert_eq!(result, 8, "Should trim low quality tail");
}

#[test]
fn test_trim_tail_entire_read_low_quality() {
    let qual = b"##########"; // Quality 2 throughout
    let result = trim_tail_sliding_window(qual, 4, 20);
    assert_eq!(result, 0, "Should trim entire low quality read");
}

#[test]
fn test_trim_tail_shorter_than_window() {
    let qual = b"III"; // Shorter than window size
    let result = trim_tail_sliding_window(qual, 4, 20);
    assert_eq!(result, 3, "Should not trim reads shorter than window");
}

#[test]
fn test_trim_tail_gradual_decline() {
    // Actual quality values: I=40, H=39, E=36, <=27, 6=21, /=14, +=10
    let qual = b"IHE<6/+";
    let result = trim_tail_sliding_window(qual, 4, 25);
    // Window at position 7: [27,21,14,10] mean=18 < 25
    // Window at position 6: [36,27,21,14] mean=24.5 < 25
    // Window at position 5: [39,36,27,21] mean=30.75 >= 25
    assert_eq!(result, 5, "Should trim at gradual quality decline");
}

#[test]
fn test_trim_front_no_trimming_needed() {
    let qual = b"IIIIIIIIII"; // Quality 40
    let result = trim_front_sliding_window(qual, 4, 20);
    assert_eq!(result, 0, "Should not trim high quality read");
}

#[test]
fn test_trim_front_low_quality_start() {
    let qual = b"####IIIIII"; // Quality 2 then 40
    let result = trim_front_sliding_window(qual, 4, 20);
    // Window at pos 0: #### mean=2 < 20
    // Window at pos 1: ###I mean=11.5 < 20
    // Window at pos 2: ##II mean=21 >= 20
    assert_eq!(result, 2, "Should trim low quality front");
}

#[test]
fn test_trim_front_entire_read_low_quality() {
    let qual = b"##########";
    let result = trim_front_sliding_window(qual, 4, 20);
    assert_eq!(result, 10, "Should trim entire low quality read");
}

#[test]
fn test_trim_front_gradual_improvement() {
    // Quality improves: 10, 15, 20, 25, 30, 35, 40
    let qual = b"+/6<EHI";
    let result = trim_front_sliding_window(qual, 4, 25);
    // Window at position 0: [10,15,20,25] mean=17.5 < 25
    // Window at position 1: [15,20,25,30] mean=22.5 < 25
    // Window at position 2: [20,25,30,35] mean=27.5 >= 25
    assert_eq!(result, 2, "Should trim until quality improves");
}

// PolyG Detection Tests

#[test]
fn test_poly_g_detection_basic() {
    let seq = b"ACGTACGTGGGGGGGGGGGG"; // 12 G's at end
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 12, "Should detect 12 G's at tail");
}

#[test]
fn test_poly_g_detection_below_threshold() {
    let seq = b"ACGTACGTGGG"; // Only 3 G's
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 0, "Should not detect polyG below threshold");
}

#[test]
fn test_poly_g_detection_exact_threshold() {
    let seq = b"ACGTGGGGGGGGGG"; // Exactly 10 G's
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 10, "Should detect polyG at exact threshold");
}

#[test]
fn test_poly_g_detection_lowercase() {
    let seq = b"ACGTgggggggggggg"; // Lowercase g's
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 12, "Should detect lowercase g's");
}

#[test]
fn test_poly_g_detection_interrupted() {
    let seq = b"ACGTGGGGAGGGGG"; // G's interrupted by A (only 5 G's at tail, below threshold)
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 0, "Should not detect polyG below threshold");
}

#[test]
fn test_poly_g_detection_entire_sequence() {
    let seq = b"GGGGGGGGGGGGGGGG"; // All G's
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 16, "Should detect entire sequence as polyG");
}

#[test]
fn test_poly_g_detection_mixed_case() {
    let seq = b"ACGTGGGGggggGGGG"; // Mixed case
    let result = detect_poly_g_tail(seq, 10);
    assert_eq!(result, 12, "Should detect mixed case G's");
}

// PolyX Detection Tests

#[test]
fn test_poly_x_detection_poly_a() {
    let seq = b"ACGTAAAAAAAAAAAA"; // 12 A's at end
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 12, "Should detect polyA");
}

#[test]
fn test_poly_x_detection_poly_t() {
    let seq = b"ACGTTTTTTTTTTT"; // 11 T's at end
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 11, "Should detect polyT");
}

#[test]
fn test_poly_x_detection_poly_c() {
    let seq = b"ACGTCCCCCCCCCCCC"; // 12 C's at end
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 12, "Should detect polyC");
}

#[test]
fn test_poly_x_detection_poly_n() {
    let seq = b"ACGTNNNNNNNNNNNN"; // 12 N's at end
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 12, "Should detect polyN");
}

#[test]
fn test_poly_x_detection_below_threshold() {
    let seq = b"ACGTAAAAA"; // Only 5 A's
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 0, "Should not detect polyX below threshold");
}

#[test]
fn test_poly_x_detection_empty_sequence() {
    let seq = b"";
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 0, "Should handle empty sequence");
}

#[test]
fn test_poly_x_detection_interrupted() {
    let seq = b"ACGTAAAACAAAA"; // A's interrupted (only 4 A's at tail, below threshold)
    let result = detect_poly_x_tail(seq, 10);
    assert_eq!(result, 0, "Should not detect polyX below threshold");
}

// Full Trimming Integration Tests

#[test]
fn test_trim_read_no_trimming() {
    let seq = b"ACGTACGTACGT";
    let qual = b"IIIIIIIIIIII"; // High quality
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    assert_eq!(result.start_pos, 0);
    assert_eq!(result.end_pos, 12);
    assert_eq!(result.trimmed_length(), 12);
}

#[test]
fn test_trim_read_fixed_trimming() {
    let seq = b"ACGTACGTACGT";
    let qual = b"IIIIIIIIIIII";
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 2,
        trim_tail_bases: 3,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    assert_eq!(result.start_pos, 2);
    assert_eq!(result.end_pos, 9); // 12 - 3
    assert_eq!(result.trimmed_length(), 7);
}

#[test]
fn test_trim_read_quality_trimming() {
    let seq = b"ACGTACGTACGT";
    let qual = b"####IIII####"; // Low quality at ends
    let config = TrimmingConfig {
        enable_trim_front: true,
        enable_trim_tail: true,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    // Front: trims to position 2 (##II has mean 21)
    // Tail: from qual[2..12] = "##IIII####", all windows with IIII have mean >= 20
    // So no tail trimming occurs
    assert_eq!(result.start_pos, 2);
    assert_eq!(result.end_pos, 10); // Length after front trim, before tail would remove more
    assert_eq!(result.trimmed_length(), 8);
}

#[test]
fn test_trim_read_poly_g_trimming() {
    let seq = b"ACGTACGTGGGGGGGGGGGG";
    let qual = b"IIIIIIIIIIIIIIIIIIII"; // High quality
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: true,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    assert_eq!(result.start_pos, 0);
    assert_eq!(result.end_pos, 8); // Trimmed 12 G's
    assert_eq!(result.poly_g_trimmed, 12);
    assert_eq!(result.trimmed_length(), 8);
}

#[test]
fn test_trim_read_poly_x_trimming() {
    let seq = b"ACGTACGTAAAAAAAAAAAA";
    let qual = b"IIIIIIIIIIIIIIIIIIII";
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: true,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    assert_eq!(result.start_pos, 0);
    assert_eq!(result.end_pos, 8); // Trimmed 12 A's
    assert_eq!(result.poly_x_trimmed, 12);
    assert_eq!(result.trimmed_length(), 8);
}

#[test]
fn test_trim_read_combined_trimming() {
    // Test all trimming types together
    let seq = b"NNACGTACGTACGTGGGGGGGGGGGGNN";
    let qual = b"##IIIIIIIIIIII##############";
    let config = TrimmingConfig {
        enable_trim_front: true,
        enable_trim_tail: true,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 2,
        trim_tail_bases: 2,
        max_len: 0,
        enable_poly_g: true,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    // Fixed: 2->26 (trim 2 from each end)
    // Quality: 4->14 (trim low quality)
    // PolyG: trim 12 G's, but they're in low quality region already removed
    assert!(result.trimmed_length() > 0);
    assert!(result.start_pos >= 2);
}

#[test]
fn test_trim_read_trim_entire_read() {
    let seq = b"GGGGGGGGGGGG";
    let qual = b"############"; // All low quality
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: true,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: true,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    let result = trim_read(seq, qual, &config);
    assert_eq!(result.trimmed_length(), 0, "Should trim entire read");
}

// TrimmingStats Tests

#[test]
fn test_trimming_stats_add() {
    let mut stats = TrimmingStats::default();
    let result = TrimmingResult {
        start_pos: 5,
        end_pos: 90,
        poly_g_trimmed: 10,
        poly_x_trimmed: 0,
        adapter_trimmed: false,
        adapter_bases_trimmed: 0,
    };

    stats.add(&result);

    assert_eq!(stats.reads_trimmed, 1);
    assert_eq!(stats.bases_trimmed_5, 5);
    assert_eq!(stats.bases_trimmed_3, 10);
    assert_eq!(stats.poly_g_trimmed_reads, 1);
    assert_eq!(stats.poly_g_trimmed_bases, 10);
}

#[test]
fn test_trimming_stats_merge() {
    let mut stats1 = TrimmingStats {
        reads_trimmed: 100,
        bases_trimmed_5: 500,
        bases_trimmed_3: 1000,
        poly_g_trimmed_reads: 50,
        poly_g_trimmed_bases: 600,
        poly_x_trimmed_reads: 20,
        poly_x_trimmed_bases: 250,
    };

    let stats2 = TrimmingStats {
        reads_trimmed: 200,
        bases_trimmed_5: 1000,
        bases_trimmed_3: 2000,
        poly_g_trimmed_reads: 100,
        poly_g_trimmed_bases: 1200,
        poly_x_trimmed_reads: 40,
        poly_x_trimmed_bases: 500,
    };

    stats1.merge(&stats2);

    assert_eq!(stats1.reads_trimmed, 300);
    assert_eq!(stats1.bases_trimmed_5, 1500);
    assert_eq!(stats1.bases_trimmed_3, 3000);
    assert_eq!(stats1.poly_g_trimmed_reads, 150);
    assert_eq!(stats1.poly_g_trimmed_bases, 1800);
}

// TrimmingConfig Tests

#[test]
fn test_trimming_config_disabled() {
    let mut adapter_config = crate::adapter::AdapterConfig::new();
    adapter_config.detect_adapter_for_pe = false; // Explicitly disable auto-detection for this test

    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 0,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config,
    };

    assert!(!config.is_enabled(), "Config should be disabled");
}

#[test]
fn test_trimming_config_enabled_quality() {
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: true,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    assert!(config.is_enabled(), "Config should be enabled");
}

#[test]
fn test_trimming_config_enabled_fixed() {
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 0,
        cut_window_size: 4,
        trim_front_bases: 5,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: false,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    assert!(config.is_enabled(), "Config should be enabled");
}

#[test]
fn test_trimming_config_enabled_poly() {
    let config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: false,
        cut_mean_quality: 0,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: true,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: crate::adapter::AdapterConfig::new(),
    };

    assert!(config.is_enabled(), "Config should be enabled");
}
