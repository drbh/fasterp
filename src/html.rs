//! HTML report generation
//!
//! Generates a clean, minimal HTML report with statistics tables.
//! No external dependencies - just pure Rust string formatting.

use crate::Args;
use crate::stats::{DetailedReadStats, FasterpReport, InsertSizeStats};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::Write;

/// Generate HTML report from processing statistics
pub fn generate_html_report(report: &FasterpReport, args: &Args, output_path: &str) -> Result<()> {
    let mut html = String::with_capacity(50_000);

    write_header(&mut html);
    write_body(&mut html, report, args);
    write_footer(&mut html);

    // Write to file
    let mut file =
        File::create(output_path).context(format!("Failed to create HTML file: {output_path}"))?;
    file.write_all(html.as_bytes())?;

    Ok(())
}

fn write_header(html: &mut String) {
    html.push_str("<!DOCTYPE html>\n<html lang='en'>\n<head>\n");
    html.push_str("<meta charset='utf-8'>\n");
    html.push_str("<meta name='viewport' content='width=device-width, initial-scale=1.0'>\n");
    html.push_str("<title>fasterp Quality Control Report</title>\n");

    // Add Plotly.js for charts
    html.push_str("<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>\n");

    // Add JetBrains Mono font for code display
    html.push_str("<link rel='preconnect' href='https://fonts.googleapis.com'>\n");
    html.push_str("<link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>\n");
    html.push_str("<link href='https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500&display=swap' rel='stylesheet'>\n");

    write_css(html);
    html.push_str("</head>\n<body>\n");
    html.push_str("<div id='container'>\n");
}

fn write_css(html: &mut String) {
    html.push_str("<style>\n");
    html.push_str("* { margin: 0; padding: 0; box-sizing: border-box; }\n");
    html.push_str(
        "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; ",
    );
    html.push_str("background: #1a1a1a; padding: 20px; line-height: 1.6; color: #e0e0e0; }\n");
    html.push_str("#container { max-width: 1000px; margin: 0 auto; background: #242424; ");
    html.push_str("border: 1px solid #3a3a3a; padding: 40px; }\n");
    html.push_str("h1 { color: #f0f0f0; font-size: 28px; font-weight: 600; ");
    html.push_str(
        "margin-bottom: 30px; border-bottom: 2px solid #4a4a4a; padding-bottom: 10px; }\n",
    );
    html.push_str("h2 { color: #c0c0c0; font-size: 18px; font-weight: 600; ");
    html.push_str("margin-top: 35px; margin-bottom: 15px; text-transform: uppercase; ");
    html.push_str("letter-spacing: 0.5px; }\n");
    html.push_str("table { width: 100%; border-collapse: collapse; margin: 15px 0 30px 0; }\n");
    html.push_str("th { background: #2a2a2a; color: #e0e0e0; padding: 10px; ");
    html.push_str("text-align: left; font-weight: 600; font-size: 13px; ");
    html.push_str("border-bottom: 2px solid #3a3a3a; }\n");
    html.push_str("td { padding: 10px; border-bottom: 1px solid #333; font-size: 14px; }\n");
    html.push_str("tr:last-child td { border-bottom: none; }\n");
    html.push_str(".col1 { font-weight: 500; color: #999; width: 250px; }\n");
    html.push_str(".value { color: #d0d0d0; }\n");
    html.push_str(".badge { display: inline-block; padding: 2px 8px; ");
    html.push_str("border-radius: 3px; font-size: 12px; ");
    html.push_str("background: #3a3a3a; color: #b0b0b0; margin-left: 6px; }\n");
    html.push_str("footer { margin-top: 50px; padding-top: 20px; ");
    html.push_str("border-top: 1px solid #3a3a3a; color: #888; font-size: 12px; }\n");
    html.push_str(".command { background: #2a2a2a; padding: 10px; ");
    html.push_str("border-left: 3px solid #4a4a4a; margin: 10px 0; ");
    html.push_str("font-family: 'JetBrains Mono', 'Courier New', monospace; font-size: 12px; ");
    html.push_str("color: #b0b0b0; word-break: break-all; }\n");
    html.push_str(".chart { width: 100%; height: 400px; margin: 20px 0; ");
    html.push_str("border: 1px solid #3a3a3a; }\n");
    html.push_str(".kmer_table { border-collapse: collapse; margin: 20px 0; ");
    html.push_str("font-family: 'JetBrains Mono', 'Courier New', monospace; font-size: 3px; }\n");
    html.push_str(
        ".kmer_table td { border: 1px solid #3a3a3a; padding: 1px; text-align: center; }\n",
    );
    html.push_str(".sub_section_tips { color: #888; font-size: 12px; ");
    html.push_str("font-style: italic; margin: 10px 0; }\n");
    html.push_str("</style>\n");
}

fn write_body(html: &mut String, report: &FasterpReport, args: &Args) {
    // Header
    html.push_str("<h1>fasterp Quality Control Report</h1>\n");

    // Summary section
    write_summary(html, report);

    // Before/After comparison
    write_before_after_comparison(html, report);

    // Filtering results
    write_filtering_results(html, report);

    // Determine if this is paired-end
    let is_paired_end = report.read2_before_filtering.is_some();

    // Insert size estimation (only for paired-end)
    if is_paired_end {
        if let Some(ref insert_size) = report.insert_size {
            write_insert_size_chart(html, insert_size);
        }
    }

    // Charts section - vertical layout
    html.push_str("<h2>Filtering Statistics</h2>\n");

    // Read1 charts
    if is_paired_end {
        html.push_str(
            "<h3 style='color: #a0a0a0; font-size: 16px; margin-top: 25px;'>Read 1</h3>\n",
        );
    }

    // Base quality charts - Read1
    write_quality_chart(
        html,
        &report.read1_before_filtering,
        "Before filtering: Base Mean Quality",
        "chart_r1_before",
    );
    if let Some(ref after_stats) = report.read1_after_filtering {
        write_quality_chart(
            html,
            after_stats,
            "After filtering: Base Mean Quality",
            "chart_r1_after",
        );
    }

    // Quality histograms - Read1
    write_quality_histogram(
        html,
        &report.read1_before_filtering,
        "Before filtering: Quality Score Histogram",
        "qual_hist_r1_before",
    );
    if let Some(ref after_stats) = report.read1_after_filtering {
        write_quality_histogram(
            html,
            after_stats,
            "After filtering: Quality Score Histogram",
            "qual_hist_r1_after",
        );
    }

    // Base contents charts - Read1
    write_base_contents_chart(
        html,
        &report.read1_before_filtering,
        "Before filtering: Base Contents",
        "contents_r1_before",
    );
    if let Some(ref after_stats) = report.read1_after_filtering {
        write_base_contents_chart(
            html,
            after_stats,
            "After filtering: Base Contents",
            "contents_r1_after",
        );
    }

    // KMER tables - Read1
    write_kmer_table(
        html,
        &report.read1_before_filtering,
        "Before filtering: KMER Counting",
    );
    if let Some(ref after_stats) = report.read1_after_filtering {
        write_kmer_table(html, after_stats, "After filtering: KMER Counting");
    }

    // Read2 charts if paired-end
    if let Some(ref r2_before) = report.read2_before_filtering {
        html.push_str(
            "<h3 style='color: #a0a0a0; font-size: 16px; margin-top: 35px;'>Read 2</h3>\n",
        );

        // Base quality charts - Read2
        write_quality_chart(
            html,
            r2_before,
            "Before filtering: Base Mean Quality",
            "chart_r2_before",
        );
        if let Some(ref after_stats) = report.read2_after_filtering {
            write_quality_chart(
                html,
                after_stats,
                "After filtering: Base Mean Quality",
                "chart_r2_after",
            );
        }

        // Quality histograms - Read2
        write_quality_histogram(
            html,
            r2_before,
            "Before filtering: Quality Score Histogram",
            "qual_hist_r2_before",
        );
        if let Some(ref after_stats) = report.read2_after_filtering {
            write_quality_histogram(
                html,
                after_stats,
                "After filtering: Quality Score Histogram",
                "qual_hist_r2_after",
            );
        }

        // Base contents charts - Read2
        write_base_contents_chart(
            html,
            r2_before,
            "Before filtering: Base Contents",
            "contents_r2_before",
        );
        if let Some(ref after_stats) = report.read2_after_filtering {
            write_base_contents_chart(
                html,
                after_stats,
                "After filtering: Base Contents",
                "contents_r2_after",
            );
        }

        // KMER tables - Read2
        write_kmer_table(html, r2_before, "Before filtering: KMER Counting");
        if let Some(ref after_stats) = report.read2_after_filtering {
            write_kmer_table(html, after_stats, "After filtering: KMER Counting");
        }
    }

    // Command info
    write_command_info(html, args, report);
}

fn write_summary(html: &mut String, report: &FasterpReport) {
    html.push_str("<h2>Summary</h2>\n");
    html.push_str("<table>\n");

    add_row(
        html,
        "fasterp version",
        &format!("{} (fasterp)", report.summary.fastp_version),
    );
    add_row(html, "sequencing", &report.summary.sequencing);

    let before = &report.summary.before_filtering;
    let after = &report.summary.after_filtering;

    add_row(
        html,
        "mean length before filtering",
        &format!("{}bp", before.read1_mean_length),
    );
    add_row(
        html,
        "mean length after filtering",
        &format!("{}bp", after.read1_mean_length),
    );

    if let Some(dup) = &report.duplication {
        add_row(html, "duplication rate", &format!("{}%", dup.rate * 100.0));
    }

    if let Some(insert_size) = &report.insert_size {
        add_row(html, "Insert size peak", &format!("{}", insert_size.peak));
    }

    html.push_str("</table>\n");
}

fn write_before_after_comparison(html: &mut String, report: &FasterpReport) {
    let before = &report.summary.before_filtering;
    let after = &report.summary.after_filtering;

    // Before Filtering section
    html.push_str("<h2>Before Filtering</h2>\n");
    html.push_str("<table>\n");
    add_row(html, "total reads", &format_number(before.total_reads));
    add_row(html, "total bases", &format_number(before.total_bases));
    add_row(
        html,
        "Q20 bases",
        &format!(
            "{} <span class='badge'>{}%</span>",
            format_number(before.q20_bases),
            before.q20_rate * 100.0
        ),
    );
    add_row(
        html,
        "Q30 bases",
        &format!(
            "{} <span class='badge'>{}%</span>",
            format_number(before.q30_bases),
            before.q30_rate * 100.0
        ),
    );
    add_row(html, "Q40 bases", "0 <span class='badge'>0%</span>");
    add_row(
        html,
        "GC content",
        &format!("{}%", before.gc_content * 100.0),
    );
    html.push_str("</table>\n");

    // After Filtering section
    html.push_str("<h2>After Filtering</h2>\n");
    html.push_str("<table>\n");
    add_row(html, "total reads", &format_number(after.total_reads));
    add_row(html, "total bases", &format_number(after.total_bases));
    add_row(
        html,
        "Q20 bases",
        &format!(
            "{} <span class='badge'>{}%</span>",
            format_number(after.q20_bases),
            after.q20_rate * 100.0
        ),
    );
    add_row(
        html,
        "Q30 bases",
        &format!(
            "{} <span class='badge'>{}%</span>",
            format_number(after.q30_bases),
            after.q30_rate * 100.0
        ),
    );
    add_row(html, "Q40 bases", "0 <span class='badge'>0%</span>");
    add_row(
        html,
        "GC content",
        &format!("{}%", after.gc_content * 100.0),
    );
    html.push_str("</table>\n");
}

fn write_filtering_results(html: &mut String, report: &FasterpReport) {
    html.push_str("<h2>Filtering Result</h2>\n");
    html.push_str("<table>\n");

    let fr = &report.filtering_result;
    let total = fr.passed_filter_reads
        + fr.low_quality_reads
        + fr.low_complexity_reads
        + fr.too_many_n_reads
        + fr.too_short_reads
        + fr.too_long_reads;

    let calc_pct = |count: usize| -> String {
        if total > 0 {
            format!("{}%", count as f64 * 100.0 / total as f64)
        } else {
            "0%".to_string()
        }
    };

    add_row(
        html,
        "reads passed filters",
        &format!(
            "{} <span class='badge'>{}</span>",
            format_number(fr.passed_filter_reads),
            calc_pct(fr.passed_filter_reads)
        ),
    );

    add_row(
        html,
        "reads with low quality",
        &format!(
            "{} <span class='badge'>{}</span>",
            format_number(fr.low_quality_reads),
            calc_pct(fr.low_quality_reads)
        ),
    );

    add_row(
        html,
        "reads with too many N",
        &format!(
            "{} <span class='badge'>{}</span>",
            format_number(fr.too_many_n_reads),
            calc_pct(fr.too_many_n_reads)
        ),
    );

    add_row(
        html,
        "reads too short",
        &format!(
            "{} <span class='badge'>{}</span>",
            format_number(fr.too_short_reads),
            calc_pct(fr.too_short_reads)
        ),
    );

    html.push_str("</table>\n");
}

fn write_base_contents_chart(
    html: &mut String,
    stats: &DetailedReadStats,
    title: &str,
    div_id: &str,
) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    html.push_str(&format!("<div id='{div_id}' class='chart'></div>\n"));
    html.push_str("<script>\n");

    // Create traces for A, T, C, G, N, GC using real data
    let bases = [
        ("A", "#8dd3c7", &stats.content_curves.a),
        ("T", "#bebada", &stats.content_curves.t),
        ("C", "#fb8072", &stats.content_curves.c),
        ("G", "#80b1d3", &stats.content_curves.g),
        ("N", "#ff0000", &stats.content_curves.n),
        ("GC", "#333", &stats.content_curves.gc),
    ];

    for (base_name, color, data) in &bases {
        html.push_str(&format!("var {}_trace = {{\n", base_name.to_lowercase()));
        html.push_str("  x: [");
        for i in 0..data.len() {
            if i > 0 {
                html.push(',');
            }
            html.push_str(&i.to_string());
        }
        html.push_str("],\n  y: [");

        for (i, val) in data.iter().enumerate() {
            if i > 0 {
                html.push(',');
            }
            html.push_str(&format!("{val:.2}"));
        }

        html.push_str("],\n");
        html.push_str(&format!("  name: '{base_name}',\n"));
        html.push_str("  type: 'scatter',\n  mode: 'lines',\n");
        html.push_str(&format!("  line: {{color: '{color}', width: 1}},\n"));
        html.push_str("  opacity: 0.8\n");
        html.push_str("};\n");
    }

    html.push_str("var data = [a_trace, t_trace, c_trace, g_trace, n_trace, gc_trace];\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Position in read (bp)', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n",
    );
    html.push_str(
        "  yaxis: {title: 'Base content (%)', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n",
    );
    html.push_str("  margin: {l: 50, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: true,\n");
    html.push_str("  legend: {x: 1, xanchor: 'right', y: 1, font: {color: '#c0c0c0'}},\n");
    html.push_str("  paper_bgcolor: '#242424',\n");
    html.push_str("  plot_bgcolor: '#2a2a2a'\n");
    html.push_str("};\n");
    html.push_str(&format!("Plotly.newPlot('{div_id}', data, layout);\n"));
    html.push_str("</script>\n");
}

fn write_kmer_table(html: &mut String, stats: &DetailedReadStats, title: &str) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    html.push_str("<div class='sub_section_tips'>Darker background means larger counts. Hover to see count.</div>\n");
    html.push_str("<table class='kmer_table' style='border-collapse: collapse; font-size: 3px; font-family: monospace;'>\n");

    // Get top 256 kmers (5-mers as in fastp)
    let mut kmer_vec: Vec<(&String, &usize)> = stats.kmer_count.iter().collect();
    kmer_vec.sort_by(|a, b| b.1.cmp(a.1));

    // Calculate mean count for normalization
    let total_count: usize = kmer_vec.iter().map(|(_, c)| **c).sum();
    let mean_count = if kmer_vec.is_empty() {
        1.0
    } else {
        total_count as f64 / kmer_vec.len() as f64
    };

    // Create 16x16 grid (256 5-mers organized by first 3 and last 2 bases)
    html.push_str("<tr><td style='font-size:3px;'></td>");
    // Header row with 2-mer suffixes
    let bases = ['A', 'T', 'C', 'G'];
    for b1 in &bases {
        for b2 in &bases {
            html.push_str(&format!("<td style='color:#666; font-weight:bold; text-align:center; padding:1px; font-size:3px;'>{b1}{b2}</td>"));
        }
    }
    html.push_str("</tr>\n");

    // Content rows with 3-mer prefixes
    for b1 in &bases {
        for b2 in &bases {
            for b3 in &bases {
                html.push_str("<tr>");
                html.push_str(&format!("<td style='color:#666; font-weight:bold; padding:1px; font-size:3px;'>{b1}{b2}{b3}</td>"));

                // For each 2-mer suffix
                for b4 in &bases {
                    for b5 in &bases {
                        let kmer = format!("{b1}{b2}{b3}{b4}{b5}");
                        let count = stats.kmer_count.get(&kmer).copied().unwrap_or(0);

                        // Calculate color intensity
                        let prop = count as f64 / mean_count;
                        let frac = if prop > 2.0 {
                            (prop - 2.0) / 20.0 + 0.5
                        } else if prop < 0.5 {
                            prop
                        } else {
                            0.5
                        };
                        let frac = frac.max(0.01).min(1.0);
                        let gray = ((1.0 - frac) * 255.0) as u8;

                        html.push_str(&format!(
                            "<td style='background:#{gray:02x}{gray:02x}{gray:02x}; text-align:center; padding:1px; cursor:help; font-size:3px;' title='{kmer}:{count} ({prop:.2}x mean)'>"
                        ));
                        html.push_str(&kmer);
                        html.push_str("</td>");
                    }
                }
                html.push_str("</tr>\n");
            }
        }
    }

    html.push_str("</table>\n");
}

fn write_quality_histogram(
    html: &mut String,
    stats: &DetailedReadStats,
    title: &str,
    div_id: &str,
) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");

    if let Some(ref qual_hist) = stats.qual_hist {
        html.push_str(&format!("<div id='{div_id}' class='chart'></div>\n"));
        html.push_str("<script>\n");

        html.push_str("var hist_trace = {\n");
        html.push_str("  x: [");
        let mut first = true;
        for (q, _) in qual_hist {
            if !first {
                html.push(',');
            }
            html.push_str(&q.to_string());
            first = false;
        }
        html.push_str("],\n  y: [");

        first = true;
        for (_, count) in qual_hist {
            if !first {
                html.push(',');
            }
            html.push_str(&count.to_string());
            first = false;
        }
        html.push_str("],\n");
        html.push_str("  type: 'bar',\n");
        html.push_str("  marker: {color: '#8dd3c7'}\n");
        html.push_str("};\n");

        html.push_str("var data = [hist_trace];\n");
        html.push_str("var layout = {\n");
        html.push_str(
            "  xaxis: {title: 'Base quality score', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n",
        );
        html.push_str("  yaxis: {title: 'Base count', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n");
        html.push_str("  margin: {l: 60, r: 30, t: 30, b: 50},\n");
        html.push_str("  showlegend: false,\n");
        html.push_str("  paper_bgcolor: '#242424',\n");
        html.push_str("  plot_bgcolor: '#2a2a2a'\n");
        html.push_str("};\n");
        html.push_str(&format!("Plotly.newPlot('{div_id}', data, layout);\n"));
        html.push_str("</script>\n");
    } else {
        html.push_str("<p>No quality histogram data available</p>\n");
    }
}

fn write_quality_chart(html: &mut String, stats: &DetailedReadStats, title: &str, div_id: &str) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    html.push_str(&format!("<div id='{div_id}' class='chart'></div>\n"));
    html.push_str("<script>\n");

    // Mean quality trace
    html.push_str("var mean_trace = {\n");
    html.push_str("  x: [");
    for i in 0..stats.quality_curves.mean.len() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&i.to_string());
    }
    html.push_str("],\n  y: [");
    for (i, val) in stats.quality_curves.mean.iter().enumerate() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&format!("{val:.2}"));
    }
    html.push_str("],\n  name: 'Mean',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#333', width: 2}\n");
    html.push_str("};\n");

    // A trace
    html.push_str("var a_trace = {\n");
    html.push_str("  x: [");
    for i in 0..stats.quality_curves.a.len() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&i.to_string());
    }
    html.push_str("],\n  y: [");
    for (i, val) in stats.quality_curves.a.iter().enumerate() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&format!("{val:.2}"));
    }
    html.push_str("],\n  name: 'A',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#8dd3c7', width: 1},\n  opacity: 0.6\n");
    html.push_str("};\n");

    // T trace
    html.push_str("var t_trace = {\n");
    html.push_str("  x: [");
    for i in 0..stats.quality_curves.t.len() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&i.to_string());
    }
    html.push_str("],\n  y: [");
    for (i, val) in stats.quality_curves.t.iter().enumerate() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&format!("{val:.2}"));
    }
    html.push_str("],\n  name: 'T',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#bebada', width: 1},\n  opacity: 0.6\n");
    html.push_str("};\n");

    // C trace
    html.push_str("var c_trace = {\n");
    html.push_str("  x: [");
    for i in 0..stats.quality_curves.c.len() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&i.to_string());
    }
    html.push_str("],\n  y: [");
    for (i, val) in stats.quality_curves.c.iter().enumerate() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&format!("{val:.2}"));
    }
    html.push_str("],\n  name: 'C',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#fb8072', width: 1},\n  opacity: 0.6\n");
    html.push_str("};\n");

    // G trace
    html.push_str("var g_trace = {\n");
    html.push_str("  x: [");
    for i in 0..stats.quality_curves.g.len() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&i.to_string());
    }
    html.push_str("],\n  y: [");
    for (i, val) in stats.quality_curves.g.iter().enumerate() {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&format!("{val:.2}"));
    }
    html.push_str("],\n  name: 'G',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#80b1d3', width: 1},\n  opacity: 0.6\n");
    html.push_str("};\n");

    // Layout and plot
    html.push_str("var data = [mean_trace, a_trace, t_trace, c_trace, g_trace];\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Position in read (bp)', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n",
    );
    html.push_str("  yaxis: {title: 'Quality score', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n");
    html.push_str("  margin: {l: 50, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: true,\n");
    html.push_str("  legend: {x: 1, xanchor: 'right', y: 1, font: {color: '#c0c0c0'}},\n");
    html.push_str("  paper_bgcolor: '#242424',\n");
    html.push_str("  plot_bgcolor: '#2a2a2a'\n");
    html.push_str("};\n");
    html.push_str(&format!("Plotly.newPlot('{div_id}', data, layout);\n"));
    html.push_str("</script>\n");
}

fn write_insert_size_chart(html: &mut String, insert_size: &InsertSizeStats) {
    html.push_str("<h2>Insert Size Estimation</h2>\n");
    html.push_str("<div id='insert_size_chart' class='chart'></div>\n");
    html.push_str("<script>\n");

    // Find the last non-zero value in the histogram
    let mut max_size = 0;
    for (size, &count) in insert_size.histogram.iter().enumerate() {
        if count > 0 {
            max_size = size;
        }
    }

    // Include all values from 0 to max_size (don't filter out zeros)
    html.push_str("var insert_sizes = [");
    for size in 0..=max_size {
        if size > 0 {
            html.push(',');
        }
        html.push_str(&size.to_string());
    }
    html.push_str("];\n");

    html.push_str("var counts = [");
    for size in 0..=max_size {
        if size > 0 {
            html.push(',');
        }
        let count = insert_size.histogram.get(size).unwrap_or(&0);
        html.push_str(&count.to_string());
    }
    html.push_str("];\n");

    html.push_str("var trace = {\n");
    html.push_str("  x: insert_sizes,\n");
    html.push_str("  y: counts,\n");
    html.push_str("  type: 'bar',\n");
    html.push_str("  marker: {color: '#fb8072'}\n");
    html.push_str("};\n");

    html.push_str("var data = [trace];\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Insert size (bp)', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n",
    );
    html.push_str("  yaxis: {title: 'Read pairs', color: '#c0c0c0', gridcolor: '#3a3a3a'},\n");
    html.push_str("  margin: {l: 60, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: false,\n");
    html.push_str("  paper_bgcolor: '#242424',\n");
    html.push_str("  plot_bgcolor: '#2a2a2a'\n");
    html.push_str("};\n");
    html.push_str("Plotly.newPlot('insert_size_chart', data, layout);\n");
    html.push_str("</script>\n");
}

fn write_command_info(html: &mut String, args: &Args, report: &FasterpReport) {
    html.push_str("<footer>\n");
    html.push_str("<p><strong class='version'>Generated by fasterp v");
    html.push_str(&report.summary.fastp_version);
    html.push_str("</strong></p>\n");

    html.push_str("<div class='command'>");
    html.push_str("fasterp -i ");
    html.push_str(&args.input);
    html.push_str(" -o ");
    html.push_str(&args.output);
    html.push_str(" -j ");
    html.push_str(&args.json);
    if let Some(ref in2) = args.input2 {
        html.push_str(" -I ");
        html.push_str(in2);
    }
    if let Some(ref out2) = args.output2 {
        html.push_str(" -O ");
        html.push_str(out2);
    }
    html.push_str("</div>\n");

    html.push_str("</footer>\n");
}

fn write_footer(html: &mut String) {
    html.push_str("</div>\n");
    html.push_str("</body>\n</html>\n");
}

// Helper functions

fn add_row(html: &mut String, label: &str, value: &str) {
    html.push_str("<tr><td class='col1'>");
    html.push_str(label);
    html.push_str("</td><td class='value'>");
    html.push_str(value);
    html.push_str("</td></tr>\n");
}

fn add_row_3col(html: &mut String, label: &str, before: &str, after: &str) {
    html.push_str("<tr><td class='col1'>");
    html.push_str(label);
    html.push_str("</td><td class='value'>");
    html.push_str(before);
    html.push_str("</td><td class='value'>");
    html.push_str(after);
    html.push_str("</td></tr>\n");
}

fn format_number(n: usize) -> String {
    if n >= 1_000_000_000 {
        format!("{:.4}B", n as f64 / 1e9)
    } else if n >= 1_000_000 {
        format!("{:.4}M", n as f64 / 1e6)
    } else if n >= 1_000 {
        format!("{:.4}K", n as f64 / 1e3)
    } else {
        n.to_string()
    }
}
