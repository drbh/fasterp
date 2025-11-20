//! HTML report generation
//!
//! Generates a clean, minimal HTML report with statistics tables.
//! No external dependencies - just pure Rust string formatting.

use crate::stats::{DetailedReadStats, FasterpReport, InsertSizeStats};
use anyhow::{Context, Result};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::Write;
use std::time::Duration;

/// Configuration for HTML report generation
pub struct HtmlReportConfig {
    pub input: String,
    pub output: String,
    pub json: String,
    pub input2: Option<String>,
    pub output2: Option<String>,
}

/// Generate HTML report from processing statistics
pub fn generate_html_report(
    report: &FasterpReport,
    args: &HtmlReportConfig,
    output_path: &str,
    elapsed: Duration,
) -> Result<()> {
    let mut html = String::with_capacity(50_000);

    write_header(&mut html);
    write_body(&mut html, report, args, elapsed);
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
    // Set theme before body renders to prevent flicker
    html.push_str("<script>\n");
    html.push_str("if (localStorage.getItem('theme') === 'light') {\n");
    html.push_str("  document.documentElement.setAttribute('data-theme', 'light');\n");
    html.push_str("}\n");
    html.push_str("function getChartColors() {\n");
    html.push_str("  var style = getComputedStyle(document.documentElement);\n");
    html.push_str("  return {\n");
    html.push_str("    paper: style.getPropertyValue('--chart-paper').trim() || '#242424',\n");
    html.push_str("    plot: style.getPropertyValue('--chart-bg').trim() || '#2a2a2a',\n");
    html.push_str("    grid: style.getPropertyValue('--chart-grid').trim() || '#3a3a3a',\n");
    html.push_str("    text: style.getPropertyValue('--text-secondary').trim() || '#c0c0c0'\n");
    html.push_str("  };\n");
    html.push_str("}\n");
    html.push_str("</script>\n");
    html.push_str("</head>\n<body>\n");
    html.push_str(
        "<button class='theme-toggle' onclick='toggleTheme()'>Toggle Light/Dark</button>\n",
    );
    html.push_str("<div id='container'>\n");
}

fn write_css(html: &mut String) {
    html.push_str("<style>\n");
    // CSS Variables for theming - clean technical style
    html.push_str(":root {\n");
    html.push_str("  --bg-primary: #1a1a1a;\n");
    html.push_str("  --bg-secondary: #242424;\n");
    html.push_str("  --bg-tertiary: #2a2a2a;\n");
    html.push_str("  --border-color: #3a3a3a;\n");
    html.push_str("  --border-accent: #4a9eff;\n");
    html.push_str("  --text-primary: #e8e8e8;\n");
    html.push_str("  --text-secondary: #c0c0c0;\n");
    html.push_str("  --text-tertiary: #909090;\n");
    html.push_str("  --text-muted: #666666;\n");
    html.push_str("  --text-value: #ffffff;\n");
    html.push_str("  --badge-bg: #2a3a4a;\n");
    html.push_str("  --badge-text: #8ac4ff;\n");
    html.push_str("  --chart-bg: #2a2a2a;\n");
    html.push_str("  --chart-paper: #242424;\n");
    html.push_str("  --chart-grid: #3a3a3a;\n");
    html.push_str("}\n");
    html.push_str(":root[data-theme='light'], [data-theme='light'] {\n");
    html.push_str("  --bg-primary: #f5f5f5;\n");
    html.push_str("  --bg-secondary: #ffffff;\n");
    html.push_str("  --bg-tertiary: #fafafa;\n");
    html.push_str("  --border-color: #e0e0e0;\n");
    html.push_str("  --border-accent: #2563eb;\n");
    html.push_str("  --text-primary: #1a1a1a;\n");
    html.push_str("  --text-secondary: #4a4a4a;\n");
    html.push_str("  --text-tertiary: #6a6a6a;\n");
    html.push_str("  --text-muted: #909090;\n");
    html.push_str("  --text-value: #000000;\n");
    html.push_str("  --badge-bg: #e8f0fe;\n");
    html.push_str("  --badge-text: #1e40af;\n");
    html.push_str("  --chart-bg: #fafafa;\n");
    html.push_str("  --chart-paper: #ffffff;\n");
    html.push_str("  --chart-grid: #e0e0e0;\n");
    html.push_str("}\n");
    html.push_str("* { margin: 0; padding: 0; box-sizing: border-box; }\n");
    html.push_str("body { font-family: 'JetBrains Mono', 'Courier New', monospace; ");
    html.push_str("background: var(--bg-primary); padding: 10px; line-height: 1.4; color: var(--text-primary); font-size: 12px; }\n");
    html.push_str("#container { width: 100%; background: var(--bg-secondary); ");
    html.push_str("border: 2px solid var(--border-color); padding: 20px; }\n");
    html.push_str(".theme-toggle { position: fixed; top: 10px; right: 10px; ");
    html.push_str("padding: 4px 8px; border: 1px solid var(--border-color); ");
    html.push_str("background: var(--bg-tertiary); color: var(--text-primary); ");
    html.push_str("cursor: pointer; font-size: 10px; font-family: 'JetBrains Mono', monospace; z-index: 1000; }\n");
    html.push_str(
        ".theme-toggle:hover { background: var(--border-accent); color: var(--bg-primary); }\n",
    );
    html.push_str("h1 { color: var(--text-primary); font-size: 18px; font-weight: 500; ");
    html.push_str("margin-bottom: 20px; border-bottom: 1px solid var(--border-accent); padding-bottom: 10px; }\n");
    html.push_str("h2 { color: var(--text-primary); font-size: 13px; font-weight: 500; ");
    html.push_str("margin-top: 25px; margin-bottom: 10px; ");
    html.push_str("border-left: 3px solid var(--border-accent); padding-left: 10px; }\n");
    html.push_str("h3 { color: var(--text-secondary); font-size: 11px; font-weight: 400; margin-top: 15px; }\n");
    html.push_str("table { width: 100%; border-collapse: collapse; margin: 10px 0 20px 0; font-size: 11px; }\n");
    html.push_str(
        "th { background: var(--bg-tertiary); color: var(--text-primary); padding: 8px 10px; ",
    );
    html.push_str("text-align: left; font-weight: 500; font-size: 11px; ");
    html.push_str("border: 1px solid var(--border-color); }\n");
    html.push_str(
        "td { padding: 6px 10px; border: 1px solid var(--border-color); font-size: 11px; }\n",
    );
    html.push_str(".col1 { font-weight: 400; color: var(--text-tertiary); width: 180px; }\n");
    html.push_str(".value { color: var(--text-value); }\n");
    html.push_str(".badge { display: inline-block; padding: 1px 4px; ");
    html.push_str("font-size: 9px; ");
    html.push_str("background: var(--badge-bg); color: var(--badge-text); margin-left: 4px; border: 1px solid var(--border-color); }\n");
    html.push_str("footer { margin-top: 30px; padding-top: 15px; ");
    html.push_str(
        "border-top: 1px solid var(--border-color); color: var(--text-muted); font-size: 9px; }\n",
    );
    html.push_str(".command { background: var(--bg-tertiary); padding: 8px; ");
    html.push_str("border: 1px solid var(--border-color); margin: 8px 0; ");
    html.push_str("font-size: 10px; color: var(--text-value); word-break: break-all; }\n");
    html.push_str(".chart { width: 100%; height: 280px; margin: 8px 0; ");
    html.push_str("border: 1px solid var(--border-color); }\n");
    html.push_str(".kmer_table { border-collapse: collapse; margin: 10px 0; font-size: 3px; }\n");
    html.push_str(
        ".kmer_table td { border: 1px solid var(--border-color); padding: 1px; text-align: center; }\n",
    );
    html.push_str(".sub_section_tips { color: var(--text-muted); font-size: 9px; ");
    html.push_str("margin: 6px 0; text-transform: lowercase; }\n");
    html.push_str(".side-by-side { display: flex; gap: 15px; }\n");
    html.push_str(".side-by-side > div { flex: 1; min-width: 0; overflow: hidden; }\n");
    html.push_str(".side-by-side table { table-layout: fixed; }\n");
    html.push_str(".side-by-side h2 { margin-top: 0; }\n");
    html.push_str("</style>\n");
}

#[allow(clippy::too_many_lines)]
fn write_body(html: &mut String, report: &FasterpReport, args: &HtmlReportConfig, elapsed: Duration) {
    // Header
    html.push_str("<h1>fasterp Quality Control Report</h1>\n");

    // Summary section
    write_summary(html, report, elapsed);

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

    // Base quality charts - Read1 (side by side)
    html.push_str("<div class='side-by-side'>\n");
    html.push_str("<div>\n");
    write_quality_chart(
        html,
        &report.read1_before_filtering,
        "Before filtering: Base Mean Quality",
        "chart_r1_before",
    );
    html.push_str("</div>\n");
    if let Some(ref after_stats) = report.read1_after_filtering {
        html.push_str("<div>\n");
        write_quality_chart(
            html,
            after_stats,
            "After filtering: Base Mean Quality",
            "chart_r1_after",
        );
        html.push_str("</div>\n");
        html.push_str("<div>\n");
        write_quality_diff_chart(
            html,
            &report.read1_before_filtering,
            after_stats,
            "Difference: Base Mean Quality",
            "chart_r1_diff",
        );
        html.push_str("</div>\n");
    }
    html.push_str("</div>\n");

    // Quality histograms - Read1 (side by side)
    html.push_str("<div class='side-by-side'>\n");
    html.push_str("<div>\n");
    write_quality_histogram(
        html,
        &report.read1_before_filtering,
        "Before filtering: Quality Score Histogram",
        "qual_hist_r1_before",
    );
    html.push_str("</div>\n");
    if let Some(ref after_stats) = report.read1_after_filtering {
        html.push_str("<div>\n");
        write_quality_histogram(
            html,
            after_stats,
            "After filtering: Quality Score Histogram",
            "qual_hist_r1_after",
        );
        html.push_str("</div>\n");
        html.push_str("<div>\n");
        write_quality_histogram_diff(
            html,
            &report.read1_before_filtering,
            after_stats,
            "Difference: Quality Score Histogram",
            "qual_hist_r1_diff",
        );
        html.push_str("</div>\n");
    }
    html.push_str("</div>\n");

    // Base contents charts - Read1 (side by side)
    html.push_str("<div class='side-by-side'>\n");
    html.push_str("<div>\n");
    write_base_contents_chart(
        html,
        &report.read1_before_filtering,
        "Before filtering: Base Contents",
        "contents_r1_before",
    );
    html.push_str("</div>\n");
    if let Some(ref after_stats) = report.read1_after_filtering {
        html.push_str("<div>\n");
        write_base_contents_chart(
            html,
            after_stats,
            "After filtering: Base Contents",
            "contents_r1_after",
        );
        html.push_str("</div>\n");
        html.push_str("<div>\n");
        write_base_contents_diff_chart(
            html,
            &report.read1_before_filtering,
            after_stats,
            "Difference: Base Contents",
            "contents_r1_diff",
        );
        html.push_str("</div>\n");
    }
    html.push_str("</div>\n");

    // KMER tables - Read1 (side by side)
    html.push_str("<div class='side-by-side'>\n");
    html.push_str("<div>\n");
    write_kmer_table(
        html,
        &report.read1_before_filtering,
        "Before filtering: KMER Counting",
    );
    html.push_str("</div>\n");
    if let Some(ref after_stats) = report.read1_after_filtering {
        html.push_str("<div>\n");
        write_kmer_table(html, after_stats, "After filtering: KMER Counting");
        html.push_str("</div>\n");
        html.push_str("<div>\n");
        write_kmer_diff_table(
            html,
            &report.read1_before_filtering,
            after_stats,
            "Difference: KMER Counting",
        );
        html.push_str("</div>\n");
    }
    html.push_str("</div>\n");

    // Read2 charts if paired-end
    if let Some(ref r2_before) = report.read2_before_filtering {
        html.push_str(
            "<h3 style='color: #a0a0a0; font-size: 16px; margin-top: 35px;'>Read 2</h3>\n",
        );

        // Base quality charts - Read2 (side by side)
        html.push_str("<div class='side-by-side'>\n");
        html.push_str("<div>\n");
        write_quality_chart(
            html,
            r2_before,
            "Before filtering: Base Mean Quality",
            "chart_r2_before",
        );
        html.push_str("</div>\n");
        if let Some(ref after_stats) = report.read2_after_filtering {
            html.push_str("<div>\n");
            write_quality_chart(
                html,
                after_stats,
                "After filtering: Base Mean Quality",
                "chart_r2_after",
            );
            html.push_str("</div>\n");
            html.push_str("<div>\n");
            write_quality_diff_chart(
                html,
                r2_before,
                after_stats,
                "Difference: Base Mean Quality",
                "chart_r2_diff",
            );
            html.push_str("</div>\n");
        }
        html.push_str("</div>\n");

        // Quality histograms - Read2 (side by side)
        html.push_str("<div class='side-by-side'>\n");
        html.push_str("<div>\n");
        write_quality_histogram(
            html,
            r2_before,
            "Before filtering: Quality Score Histogram",
            "qual_hist_r2_before",
        );
        html.push_str("</div>\n");
        if let Some(ref after_stats) = report.read2_after_filtering {
            html.push_str("<div>\n");
            write_quality_histogram(
                html,
                after_stats,
                "After filtering: Quality Score Histogram",
                "qual_hist_r2_after",
            );
            html.push_str("</div>\n");
            html.push_str("<div>\n");
            write_quality_histogram_diff(
                html,
                r2_before,
                after_stats,
                "Difference: Quality Score Histogram",
                "qual_hist_r2_diff",
            );
            html.push_str("</div>\n");
        }
        html.push_str("</div>\n");

        // Base contents charts - Read2 (side by side)
        html.push_str("<div class='side-by-side'>\n");
        html.push_str("<div>\n");
        write_base_contents_chart(
            html,
            r2_before,
            "Before filtering: Base Contents",
            "contents_r2_before",
        );
        html.push_str("</div>\n");
        if let Some(ref after_stats) = report.read2_after_filtering {
            html.push_str("<div>\n");
            write_base_contents_chart(
                html,
                after_stats,
                "After filtering: Base Contents",
                "contents_r2_after",
            );
            html.push_str("</div>\n");
            html.push_str("<div>\n");
            write_base_contents_diff_chart(
                html,
                r2_before,
                after_stats,
                "Difference: Base Contents",
                "contents_r2_diff",
            );
            html.push_str("</div>\n");
        }
        html.push_str("</div>\n");

        // KMER tables - Read2 (side by side)
        html.push_str("<div class='side-by-side'>\n");
        html.push_str("<div>\n");
        write_kmer_table(html, r2_before, "Before filtering: KMER Counting");
        html.push_str("</div>\n");
        if let Some(ref after_stats) = report.read2_after_filtering {
            html.push_str("<div>\n");
            write_kmer_table(html, after_stats, "After filtering: KMER Counting");
            html.push_str("</div>\n");
            html.push_str("<div>\n");
            write_kmer_diff_table(html, r2_before, after_stats, "Difference: KMER Counting");
            html.push_str("</div>\n");
        }
        html.push_str("</div>\n");
    }

    // Command info
    write_command_info(html, args, report);
}

fn write_summary(html: &mut String, report: &FasterpReport, elapsed: Duration) {
    html.push_str("<h2>Summary</h2>\n");
    html.push_str("<table>\n");

    add_row(
        html,
        "fasterp version",
        &format!("{} (fasterp)", report.summary.fastp_version),
    );
    add_row(html, "sequencing", &report.summary.sequencing);

    // Format runtime
    let runtime_str = if elapsed.as_secs() >= 1 {
        format!("{:.3} seconds", elapsed.as_secs_f64())
    } else {
        format!("{} ms", elapsed.as_millis())
    };
    add_row(html, "runtime", &runtime_str);

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
    let is_paired = report.read2_before_filtering.is_some();

    // Helper to format difference with sign
    fn format_diff(diff: i64) -> String {
        if diff > 0 {
            format!("+{}", format_number(diff as usize))
        } else if diff < 0 {
            format!("-{}", format_number((-diff) as usize))
        } else {
            "0".to_string()
        }
    }

    fn format_pct_diff(diff: f64) -> String {
        if diff > 0.0 {
            format!("+{:.2}%", diff)
        } else if diff < 0.0 {
            format!("{:.2}%", diff)
        } else {
            "0%".to_string()
        }
    }

    if is_paired {
        // Paired-end: show Read1 and Read2 in separate tables
        let r1_before = &report.read1_before_filtering;
        let r2_before = report.read2_before_filtering.as_ref().unwrap();
        let r1_after = report.read1_after_filtering.as_ref();
        let r2_after = report.read2_after_filtering.as_ref();
        let reads_per_file = before.total_reads / 2;
        let after_reads_per_file = after.total_reads / 2;

        // Read1 table
        html.push_str("<h2>Read 1: Before / After</h2>\n");
        html.push_str("<table>\n");
        html.push_str("<tr><th></th><th>Before</th><th>After</th><th>Diff</th></tr>\n");

        // Total reads
        let reads_diff = after_reads_per_file as i64 - reads_per_file as i64;
        html.push_str("<tr><td class='col1'>total reads</td>");
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(reads_per_file)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(after_reads_per_file)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_diff(reads_diff)
        ));

        // Total bases
        html.push_str("<tr><td class='col1'>total bases</td>");
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(r1_before.total_bases)
        ));
        if let Some(r1a) = r1_after {
            let bases_diff = r1a.total_bases as i64 - r1_before.total_bases as i64;
            html.push_str(&format!(
                "<td class='value'>{}</td>",
                format_number(r1a.total_bases)
            ));
            html.push_str(&format!(
                "<td class='value'>{}</td></tr>\n",
                format_diff(bases_diff)
            ));
        } else {
            html.push_str("<td class='value'>-</td><td class='value'>-</td></tr>\n");
        }

        // Q20 bases
        let r1_before_q20_pct = if r1_before.total_bases > 0 {
            r1_before.q20_bases as f64 * 100.0 / r1_before.total_bases as f64
        } else {
            0.0
        };
        html.push_str("<tr><td class='col1'>Q20 bases</td>");
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(r1_before.q20_bases),
            r1_before_q20_pct
        ));
        if let Some(r1a) = r1_after {
            let r1_after_q20_pct = if r1a.total_bases > 0 {
                r1a.q20_bases as f64 * 100.0 / r1a.total_bases as f64
            } else {
                0.0
            };
            let pct_diff = r1_after_q20_pct - r1_before_q20_pct;
            html.push_str(&format!(
                "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
                format_number(r1a.q20_bases),
                r1_after_q20_pct
            ));
            html.push_str(&format!(
                "<td class='value'>{}</td></tr>\n",
                format_pct_diff(pct_diff)
            ));
        } else {
            html.push_str("<td class='value'>-</td><td class='value'>-</td></tr>\n");
        }

        // Q30 bases
        let r1_before_q30_pct = if r1_before.total_bases > 0 {
            r1_before.q30_bases as f64 * 100.0 / r1_before.total_bases as f64
        } else {
            0.0
        };
        html.push_str("<tr><td class='col1'>Q30 bases</td>");
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(r1_before.q30_bases),
            r1_before_q30_pct
        ));
        if let Some(r1a) = r1_after {
            let r1_after_q30_pct = if r1a.total_bases > 0 {
                r1a.q30_bases as f64 * 100.0 / r1a.total_bases as f64
            } else {
                0.0
            };
            let pct_diff = r1_after_q30_pct - r1_before_q30_pct;
            html.push_str(&format!(
                "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
                format_number(r1a.q30_bases),
                r1_after_q30_pct
            ));
            html.push_str(&format!(
                "<td class='value'>{}</td></tr>\n",
                format_pct_diff(pct_diff)
            ));
        } else {
            html.push_str("<td class='value'>-</td><td class='value'>-</td></tr>\n");
        }

        html.push_str("</table>\n");

        // Read2 table
        html.push_str("<h2>Read 2: Before / After</h2>\n");
        html.push_str("<table>\n");
        html.push_str("<tr><th></th><th>Before</th><th>After</th><th>Diff</th></tr>\n");

        // Total reads
        html.push_str("<tr><td class='col1'>total reads</td>");
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(reads_per_file)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(after_reads_per_file)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_diff(reads_diff)
        ));

        // Total bases
        html.push_str("<tr><td class='col1'>total bases</td>");
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(r2_before.total_bases)
        ));
        if let Some(r2a) = r2_after {
            let bases_diff = r2a.total_bases as i64 - r2_before.total_bases as i64;
            html.push_str(&format!(
                "<td class='value'>{}</td>",
                format_number(r2a.total_bases)
            ));
            html.push_str(&format!(
                "<td class='value'>{}</td></tr>\n",
                format_diff(bases_diff)
            ));
        } else {
            html.push_str("<td class='value'>-</td><td class='value'>-</td></tr>\n");
        }

        // Q20 bases
        let r2_before_q20_pct = if r2_before.total_bases > 0 {
            r2_before.q20_bases as f64 * 100.0 / r2_before.total_bases as f64
        } else {
            0.0
        };
        html.push_str("<tr><td class='col1'>Q20 bases</td>");
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(r2_before.q20_bases),
            r2_before_q20_pct
        ));
        if let Some(r2a) = r2_after {
            let r2_after_q20_pct = if r2a.total_bases > 0 {
                r2a.q20_bases as f64 * 100.0 / r2a.total_bases as f64
            } else {
                0.0
            };
            let pct_diff = r2_after_q20_pct - r2_before_q20_pct;
            html.push_str(&format!(
                "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
                format_number(r2a.q20_bases),
                r2_after_q20_pct
            ));
            html.push_str(&format!(
                "<td class='value'>{}</td></tr>\n",
                format_pct_diff(pct_diff)
            ));
        } else {
            html.push_str("<td class='value'>-</td><td class='value'>-</td></tr>\n");
        }

        // Q30 bases
        let r2_before_q30_pct = if r2_before.total_bases > 0 {
            r2_before.q30_bases as f64 * 100.0 / r2_before.total_bases as f64
        } else {
            0.0
        };
        html.push_str("<tr><td class='col1'>Q30 bases</td>");
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(r2_before.q30_bases),
            r2_before_q30_pct
        ));
        if let Some(r2a) = r2_after {
            let r2_after_q30_pct = if r2a.total_bases > 0 {
                r2a.q30_bases as f64 * 100.0 / r2a.total_bases as f64
            } else {
                0.0
            };
            let pct_diff = r2_after_q30_pct - r2_before_q30_pct;
            html.push_str(&format!(
                "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
                format_number(r2a.q30_bases),
                r2_after_q30_pct
            ));
            html.push_str(&format!(
                "<td class='value'>{}</td></tr>\n",
                format_pct_diff(pct_diff)
            ));
        } else {
            html.push_str("<td class='value'>-</td><td class='value'>-</td></tr>\n");
        }

        html.push_str("</table>\n");
    } else {
        // Single-end: show simple before/after with diff
        html.push_str("<h2>Before / After Filtering</h2>\n");
        html.push_str("<table>\n");
        html.push_str("<tr><th></th><th>Before</th><th>After</th><th>Diff</th></tr>\n");

        // Total reads
        let reads_diff = after.total_reads as i64 - before.total_reads as i64;
        html.push_str("<tr><td class='col1'>total reads</td>");
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(before.total_reads)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(after.total_reads)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_diff(reads_diff)
        ));

        // Total bases
        let bases_diff = after.total_bases as i64 - before.total_bases as i64;
        html.push_str("<tr><td class='col1'>total bases</td>");
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(before.total_bases)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td>",
            format_number(after.total_bases)
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_diff(bases_diff)
        ));

        // Q20 bases
        let q20_pct_diff = (after.q20_rate - before.q20_rate) * 100.0;
        html.push_str("<tr><td class='col1'>Q20 bases</td>");
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(before.q20_bases),
            before.q20_rate * 100.0
        ));
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(after.q20_bases),
            after.q20_rate * 100.0
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_pct_diff(q20_pct_diff)
        ));

        // Q30 bases
        let q30_pct_diff = (after.q30_rate - before.q30_rate) * 100.0;
        html.push_str("<tr><td class='col1'>Q30 bases</td>");
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(before.q30_bases),
            before.q30_rate * 100.0
        ));
        html.push_str(&format!(
            "<td class='value'>{} <span class='badge'>{:.2}%</span></td>",
            format_number(after.q30_bases),
            after.q30_rate * 100.0
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_pct_diff(q30_pct_diff)
        ));

        // GC content
        let gc_diff = (after.gc_content - before.gc_content) * 100.0;
        html.push_str("<tr><td class='col1'>GC content</td>");
        html.push_str(&format!(
            "<td class='value'>{:.2}%</td>",
            before.gc_content * 100.0
        ));
        html.push_str(&format!(
            "<td class='value'>{:.2}%</td>",
            after.gc_content * 100.0
        ));
        html.push_str(&format!(
            "<td class='value'>{}</td></tr>\n",
            format_pct_diff(gc_diff)
        ));

        html.push_str("</table>\n");
    }
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
    let _ = writeln!(html, "<div id='{div_id}' class='chart'></div>");
    html.push_str("<script>\n");

    // Create traces for A, T, C, G, N, GC using real data
    let bases = [
        ("A", "#8dd3c7", &stats.content_curves.a),
        ("T", "#bebada", &stats.content_curves.t),
        ("C", "#5dade2", &stats.content_curves.c),
        ("G", "#80b1d3", &stats.content_curves.g),
        ("N", "#ff0000", &stats.content_curves.n),
        ("GC", "#333", &stats.content_curves.gc),
    ];

    for (base_name, color, data) in &bases {
        let _ = writeln!(html, "var {}_trace = {{", base_name.to_lowercase());
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
            let _ = write!(html, "{val:.2}");
        }

        html.push_str("],\n");
        let _ = writeln!(html, "  name: '{base_name}',");
        html.push_str("  type: 'scatter',\n  mode: 'lines',\n");
        let _ = writeln!(html, "  line: {{color: '{color}', width: 1}},");
        html.push_str("  opacity: 0.8\n");
        html.push_str("};\n");
    }

    html.push_str("var data = [a_trace, t_trace, c_trace, g_trace, n_trace, gc_trace];\n");
    html.push_str("var colors = getChartColors();\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Position in read (bp)', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str(
        "  yaxis: {title: 'Base content (%)', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str("  margin: {l: 50, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: true,\n");
    html.push_str("  legend: {x: 1, xanchor: 'right', y: 1, font: {color: colors.text}},\n");
    html.push_str("  paper_bgcolor: colors.paper,\n");
    html.push_str("  plot_bgcolor: colors.plot\n");
    html.push_str("};\n");
    let _ = writeln!(
        html,
        "Plotly.newPlot('{div_id}', data, layout, {{responsive: true}});"
    );
    html.push_str("</script>\n");
}

fn write_kmer_table(html: &mut String, stats: &DetailedReadStats, title: &str) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    html.push_str("<div class='sub_section_tips'>Blue intensity shows count relative to mean. Hover to see count.</div>\n");
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
            let _ = write!(
                html,
                "<td style='color:var(--text-muted); font-weight:bold; text-align:center; padding:1px; font-size:3px;'>{b1}{b2}</td>"
            );
        }
    }
    html.push_str("</tr>\n");

    // Content rows with 3-mer prefixes
    for b1 in &bases {
        for b2 in &bases {
            for b3 in &bases {
                html.push_str("<tr>");
                let _ = write!(
                    html,
                    "<td style='color:var(--text-muted); font-weight:bold; padding:1px; font-size:3px;'>{b1}{b2}{b3}</td>"
                );

                // For each 2-mer suffix
                for b4 in &bases {
                    for b5 in &bases {
                        let kmer = format!("{b1}{b2}{b3}{b4}{b5}");
                        let count = stats.kmer_count.get(&kmer).copied().unwrap_or(0);

                        // Calculate color intensity - use subtle blue with transparency
                        let prop = count as f64 / mean_count;
                        let frac = if prop > 2.0 {
                            (prop - 2.0) / 20.0 + 0.5
                        } else if prop < 0.5 {
                            prop
                        } else {
                            0.5
                        };
                        let frac = frac.clamp(0.01, 1.0);

                        // Use rgba with transparency - blue tint that works in both themes
                        let alpha = frac * 0.4;

                        let _ = write!(
                            html,
                            "<td style='background:rgba(74,158,255,{alpha:.2}); color:var(--text-primary); text-align:center; padding:1px; cursor:help; font-size:3px;' title='{kmer}:{count} ({prop:.2}x mean)'>"
                        );
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

fn write_kmer_diff_table(
    html: &mut String,
    before: &DetailedReadStats,
    after: &DetailedReadStats,
    title: &str,
) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    html.push_str("<div class='sub_section_tips'>Orange = decreased, Teal = increased. Hover to see change.</div>\n");
    html.push_str("<table class='kmer_table' style='border-collapse: collapse; font-size: 3px; font-family: monospace;'>\n");

    // Calculate mean counts for reference
    let before_total: usize = before.kmer_count.values().sum();
    let before_mean = if before.kmer_count.is_empty() {
        1.0
    } else {
        before_total as f64 / before.kmer_count.len() as f64
    };

    // Create 16x16 grid (256 5-mers organized by first 3 and last 2 bases)
    html.push_str("<tr><td style='font-size:3px;'></td>");
    // Header row with 2-mer suffixes
    let bases = ['A', 'T', 'C', 'G'];
    for b1 in &bases {
        for b2 in &bases {
            let _ = write!(
                html,
                "<td style='color:var(--text-muted); font-weight:bold; text-align:center; padding:1px; font-size:3px;'>{b1}{b2}</td>"
            );
        }
    }
    html.push_str("</tr>\n");

    // Content rows with 3-mer prefixes
    for b1 in &bases {
        for b2 in &bases {
            for b3 in &bases {
                html.push_str("<tr>");
                let _ = write!(
                    html,
                    "<td style='color:var(--text-muted); font-weight:bold; padding:1px; font-size:3px;'>{b1}{b2}{b3}</td>"
                );

                // For each 2-mer suffix
                for b4 in &bases {
                    for b5 in &bases {
                        let kmer = format!("{b1}{b2}{b3}{b4}{b5}");
                        let before_count =
                            before.kmer_count.get(&kmer).copied().unwrap_or(0) as i64;
                        let after_count = after.kmer_count.get(&kmer).copied().unwrap_or(0) as i64;
                        let diff = after_count - before_count;

                        // Calculate color based on difference using transparency
                        // Teal for increase, orange for decrease
                        let (r, g, b, alpha) = if diff > 0 {
                            // Teal with transparency for increase
                            let intensity = ((diff as f64 / before_mean).min(2.0) * 0.5).min(1.0);
                            (0, 200, 150, intensity * 0.5)
                        } else if diff < 0 {
                            // Orange with transparency for decrease
                            let intensity = ((-diff as f64 / before_mean).min(2.0) * 0.5).min(1.0);
                            (255, 150, 50, intensity * 0.5)
                        } else {
                            // Neutral - no tint
                            (128, 128, 128, 0.1)
                        };

                        let _ = write!(
                            html,
                            "<td style='background:rgba({r},{g},{b},{alpha:.2}); color:var(--text-primary); text-align:center; padding:1px; cursor:help; font-size:3px;' title='{kmer}: {diff:+} (before:{before_count}, after:{after_count})'>"
                        );
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
        let _ = writeln!(html, "<div id='{div_id}' class='chart'></div>");
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
        html.push_str("var colors = getChartColors();\n");
        html.push_str("var layout = {\n");
        html.push_str(
            "  xaxis: {title: 'Base quality score', color: colors.text, gridcolor: colors.grid},\n",
        );
        html.push_str(
            "  yaxis: {title: 'Base count', color: colors.text, gridcolor: colors.grid},\n",
        );
        html.push_str("  margin: {l: 60, r: 30, t: 30, b: 50},\n");
        html.push_str("  showlegend: false,\n");
        html.push_str("  paper_bgcolor: colors.paper,\n");
        html.push_str("  plot_bgcolor: colors.plot\n");
        html.push_str("};\n");
        let _ = writeln!(
            html,
            "Plotly.newPlot('{div_id}', data, layout, {{responsive: true}});"
        );
        html.push_str("</script>\n");
    } else {
        html.push_str("<p>No quality histogram data available</p>\n");
    }
}

#[allow(clippy::too_many_lines)]
fn write_quality_chart(html: &mut String, stats: &DetailedReadStats, title: &str, div_id: &str) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    let _ = writeln!(html, "<div id='{div_id}' class='chart'></div>");
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
        let _ = write!(html, "{val:.2}");
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
        let _ = write!(html, "{val:.2}");
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
        let _ = write!(html, "{val:.2}");
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
        let _ = write!(html, "{val:.2}");
    }
    html.push_str("],\n  name: 'C',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#5dade2', width: 1},\n  opacity: 0.6\n");
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
        let _ = write!(html, "{val:.2}");
    }
    html.push_str("],\n  name: 'G',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  line: {color: '#80b1d3', width: 1},\n  opacity: 0.6\n");
    html.push_str("};\n");

    // Layout and plot
    html.push_str("var data = [mean_trace, a_trace, t_trace, c_trace, g_trace];\n");
    html.push_str("var colors = getChartColors();\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Position in read (bp)', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str(
        "  yaxis: {title: 'Quality score', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str("  margin: {l: 50, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: true,\n");
    html.push_str("  legend: {x: 1, xanchor: 'right', y: 1, font: {color: colors.text}},\n");
    html.push_str("  paper_bgcolor: colors.paper,\n");
    html.push_str("  plot_bgcolor: colors.plot\n");
    html.push_str("};\n");
    let _ = writeln!(
        html,
        "Plotly.newPlot('{div_id}', data, layout, {{responsive: true}});"
    );
    html.push_str("</script>\n");
}

fn write_quality_diff_chart(
    html: &mut String,
    before: &DetailedReadStats,
    after: &DetailedReadStats,
    title: &str,
    div_id: &str,
) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    let _ = writeln!(html, "<div id='{div_id}' class='chart'></div>");
    html.push_str("<script>\n");

    // Calculate difference (after - before)
    let len = before
        .quality_curves
        .mean
        .len()
        .min(after.quality_curves.mean.len());

    html.push_str("var diff_trace = {\n");
    html.push_str("  x: [");
    for i in 0..len {
        if i > 0 {
            html.push(',');
        }
        html.push_str(&i.to_string());
    }
    html.push_str("],\n  y: [");
    for i in 0..len {
        if i > 0 {
            html.push(',');
        }
        let diff = after.quality_curves.mean[i] - before.quality_curves.mean[i];
        let _ = write!(html, "{diff:.2}");
    }
    html.push_str("],\n  name: 'Difference',\n  type: 'scatter',\n  mode: 'lines',\n");
    html.push_str("  fill: 'tozeroy',\n");
    html.push_str("  line: {color: '#5dade2', width: 1}\n");
    html.push_str("};\n");

    html.push_str("var data = [diff_trace];\n");
    html.push_str("var colors = getChartColors();\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Position in read (bp)', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str(
        "  yaxis: {title: 'Quality change', color: colors.text, gridcolor: colors.grid, zeroline: true, zerolinecolor: '#666'},\n",
    );
    html.push_str("  margin: {l: 50, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: false,\n");
    html.push_str("  paper_bgcolor: colors.paper,\n");
    html.push_str("  plot_bgcolor: colors.plot\n");
    html.push_str("};\n");
    let _ = writeln!(
        html,
        "Plotly.newPlot('{div_id}', data, layout, {{responsive: true}});"
    );
    html.push_str("</script>\n");
}

fn write_quality_histogram_diff(
    html: &mut String,
    before: &DetailedReadStats,
    after: &DetailedReadStats,
    title: &str,
    div_id: &str,
) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");

    if let (Some(before_hist), Some(after_hist)) = (&before.qual_hist, &after.qual_hist) {
        let _ = writeln!(html, "<div id='{div_id}' class='chart'></div>");
        html.push_str("<script>\n");

        // Build a map of quality scores to counts
        let mut diff_map = std::collections::BTreeMap::new();
        for (q, count) in before_hist {
            *diff_map.entry(*q).or_insert(0i64) -= *count as i64;
        }
        for (q, count) in after_hist {
            *diff_map.entry(*q).or_insert(0i64) += *count as i64;
        }

        html.push_str("var diff_trace = {\n");
        html.push_str("  x: [");
        let mut first = true;
        for (q, _) in &diff_map {
            if !first {
                html.push(',');
            }
            html.push_str(&q.to_string());
            first = false;
        }
        html.push_str("],\n  y: [");

        first = true;
        for (_, diff) in &diff_map {
            if !first {
                html.push(',');
            }
            html.push_str(&diff.to_string());
            first = false;
        }
        html.push_str("],\n");
        html.push_str("  type: 'bar',\n");
        html.push_str("  marker: {color: '#5dade2'}\n");
        html.push_str("};\n");

        html.push_str("var data = [diff_trace];\n");
        html.push_str("var colors = getChartColors();\n");
        html.push_str("var layout = {\n");
        html.push_str(
            "  xaxis: {title: 'Base quality score', color: colors.text, gridcolor: colors.grid},\n",
        );
        html.push_str(
            "  yaxis: {title: 'Count change', color: colors.text, gridcolor: colors.grid},\n",
        );
        html.push_str("  margin: {l: 60, r: 30, t: 30, b: 50},\n");
        html.push_str("  showlegend: false,\n");
        html.push_str("  paper_bgcolor: colors.paper,\n");
        html.push_str("  plot_bgcolor: colors.plot\n");
        html.push_str("};\n");
        let _ = writeln!(
            html,
            "Plotly.newPlot('{div_id}', data, layout, {{responsive: true}});"
        );
        html.push_str("</script>\n");
    } else {
        html.push_str("<p>No quality histogram data available</p>\n");
    }
}

fn write_base_contents_diff_chart(
    html: &mut String,
    before: &DetailedReadStats,
    after: &DetailedReadStats,
    title: &str,
    div_id: &str,
) {
    html.push_str("<h2>");
    html.push_str(title);
    html.push_str("</h2>\n");
    let _ = writeln!(html, "<div id='{div_id}' class='chart'></div>");
    html.push_str("<script>\n");

    let len = before
        .content_curves
        .a
        .len()
        .min(after.content_curves.a.len());

    // Create diff traces for A, T, C, G, GC
    let bases = [
        (
            "A",
            "#8dd3c7",
            &before.content_curves.a,
            &after.content_curves.a,
        ),
        (
            "T",
            "#bebada",
            &before.content_curves.t,
            &after.content_curves.t,
        ),
        (
            "C",
            "#5dade2",
            &before.content_curves.c,
            &after.content_curves.c,
        ),
        (
            "G",
            "#80b1d3",
            &before.content_curves.g,
            &after.content_curves.g,
        ),
        (
            "GC",
            "#333",
            &before.content_curves.gc,
            &after.content_curves.gc,
        ),
    ];

    for (base_name, color, before_data, after_data) in &bases {
        let _ = writeln!(html, "var {}_trace = {{", base_name.to_lowercase());
        html.push_str("  x: [");
        for i in 0..len {
            if i > 0 {
                html.push(',');
            }
            html.push_str(&i.to_string());
        }
        html.push_str("],\n  y: [");

        for i in 0..len {
            if i > 0 {
                html.push(',');
            }
            let diff = after_data[i] - before_data[i];
            let _ = write!(html, "{diff:.2}");
        }

        html.push_str("],\n");
        let _ = writeln!(html, "  name: '{base_name}',");
        html.push_str("  type: 'scatter',\n  mode: 'lines',\n");
        let _ = writeln!(html, "  line: {{color: '{color}', width: 1}},");
        html.push_str("  opacity: 0.8\n");
        html.push_str("};\n");
    }

    html.push_str("var data = [a_trace, t_trace, c_trace, g_trace, gc_trace];\n");
    html.push_str("var colors = getChartColors();\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Position in read (bp)', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str(
        "  yaxis: {title: 'Content change (%)', color: colors.text, gridcolor: colors.grid, zeroline: true, zerolinecolor: '#666'},\n",
    );
    html.push_str("  margin: {l: 50, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: true,\n");
    html.push_str("  legend: {x: 1, xanchor: 'right', y: 1, font: {color: colors.text}},\n");
    html.push_str("  paper_bgcolor: colors.paper,\n");
    html.push_str("  plot_bgcolor: colors.plot\n");
    html.push_str("};\n");
    let _ = writeln!(
        html,
        "Plotly.newPlot('{div_id}', data, layout, {{responsive: true}});"
    );
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
    html.push_str("  marker: {color: '#5dade2'}\n");
    html.push_str("};\n");

    html.push_str("var data = [trace];\n");
    html.push_str("var colors = getChartColors();\n");
    html.push_str("var layout = {\n");
    html.push_str(
        "  xaxis: {title: 'Insert size (bp)', color: colors.text, gridcolor: colors.grid},\n",
    );
    html.push_str("  yaxis: {title: 'Read pairs', color: colors.text, gridcolor: colors.grid},\n");
    html.push_str("  margin: {l: 60, r: 30, t: 30, b: 50},\n");
    html.push_str("  showlegend: false,\n");
    html.push_str("  paper_bgcolor: colors.paper,\n");
    html.push_str("  plot_bgcolor: colors.plot\n");
    html.push_str("};\n");
    html.push_str("Plotly.newPlot('insert_size_chart', data, layout, {responsive: true});\n");
    html.push_str("</script>\n");
}

fn write_command_info(html: &mut String, args: &HtmlReportConfig, report: &FasterpReport) {
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
    // Theme toggle and Plotly resize scripts
    html.push_str("<script>\n");
    html.push_str("function updateChartColors() {\n");
    html.push_str("  var colors = getChartColors();\n");
    html.push_str("  var charts = document.querySelectorAll('.chart');\n");
    html.push_str("  charts.forEach(function(chart) {\n");
    html.push_str("    if (chart.data) {\n");
    html.push_str("      Plotly.relayout(chart, {\n");
    html.push_str("        'paper_bgcolor': colors.paper,\n");
    html.push_str("        'plot_bgcolor': colors.plot,\n");
    html.push_str("        'xaxis.gridcolor': colors.grid,\n");
    html.push_str("        'xaxis.color': colors.text,\n");
    html.push_str("        'yaxis.gridcolor': colors.grid,\n");
    html.push_str("        'yaxis.color': colors.text,\n");
    html.push_str("        'legend.font.color': colors.text\n");
    html.push_str("      });\n");
    html.push_str("    }\n");
    html.push_str("  });\n");
    html.push_str("}\n");
    html.push_str("function toggleTheme() {\n");
    html.push_str("  var root = document.documentElement;\n");
    html.push_str("  if (root.getAttribute('data-theme') === 'light') {\n");
    html.push_str("    root.removeAttribute('data-theme');\n");
    html.push_str("    localStorage.setItem('theme', 'dark');\n");
    html.push_str("  } else {\n");
    html.push_str("    root.setAttribute('data-theme', 'light');\n");
    html.push_str("    localStorage.setItem('theme', 'light');\n");
    html.push_str("  }\n");
    html.push_str("  setTimeout(updateChartColors, 50);\n");
    html.push_str("}\n");
    html.push_str("window.addEventListener('load', function() {\n");
    html.push_str("  window.dispatchEvent(new Event('resize'));\n");
    html.push_str("  if (localStorage.getItem('theme') === 'light') {\n");
    html.push_str("    setTimeout(updateChartColors, 100);\n");
    html.push_str("  }\n");
    html.push_str("});\n");
    html.push_str("</script>\n");
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
