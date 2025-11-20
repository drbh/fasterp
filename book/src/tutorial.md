# Tutorial: Compare fasterp and fastp

This tutorial walks through downloading real sequencing data, running both **fasterp** and **fastp**, and comparing their outputs and reports.

## Prerequisites

- **fasterp** installed (`cargo install --git https://github.com/drbh/fasterp.git`)
- **fastp** installed ([installation guide](https://github.com/OpenGene/fastp#get-fastp))
- **curl** for downloading files
- **jq** (optional) for inspecting JSON reports

## Step 1: Download Test Data

Download paired-end FASTQ files from [SRX2987343](https://www.ncbi.nlm.nih.gov/sra/SRX2987343) - a chromatin accessibility study on mouse hair follicle stem cells, examining how DNA packaging changes during differentiation into hair-related cell types.

```bash
# Create a working directory
mkdir -p fasterp_tutorial && cd fasterp_tutorial

# Download R1 (forward reads) - 365 MB
curl -LO https://genedata.dholtz.com/SRX2987343/SRR5808766_1.fastq.gz

# Download R2 (reverse reads) - 364 MB
curl -LO https://genedata.dholtz.com/SRX2987343/SRR5808766_2.fastq.gz

# Verify downloads
ls -lh *.gz
```

Expected output:
```
-rw-r--r--  1 user  staff   365M  SRR5808766_1.fastq.gz
-rw-r--r--  1 user  staff   364M  SRR5808766_2.fastq.gz
```

## Step 2: Run fastp

Process the paired-end data with fastp:

```bash
fastp \
  -i SRR5808766_1.fastq.gz -I SRR5808766_2.fastq.gz \
  -o fastp_out_R1.fq.gz -O fastp_out_R2.fq.gz \
  -j fastp_report.json \
  -h fastp_report.html
```

Expected output:
```
Read1 before filtering:
total reads: 9799076
total bases: 499752876
Q20 bases: 494840792(99.0171%)
Q30 bases: 490032403(98.0549%)
Q40 bases: 333003632(66.6337%)

Read2 before filtering:
total reads: 9799076
total bases: 499752876
Q20 bases: 491032411(98.255%)
Q30 bases: 485393847(97.1268%)
Q40 bases: 335965361(67.2263%)

Read1 after filtering:
total reads: 9761787
total bases: 485985544
Q20 bases: 481444894(99.0657%)
Q30 bases: 476939326(98.1386%)
Q40 bases: 324904111(66.8547%)

Read2 after filtering:
total reads: 9761787
total bases: 485985544
Q20 bases: 479080210(98.5791%)
Q30 bases: 473663787(97.4646%)
Q40 bases: 327920585(67.4754%)

Filtering result:
reads passed filter: 19523574
reads failed due to low quality: 74384
reads failed due to too many N: 194
reads failed due to too short: 0
reads with adapter trimmed: 3113628
bases trimmed due to adapters: 23734734

Duplication rate: 28.3885%

Insert size peak (evaluated by paired-end reads): 51

JSON report: fastp_report.json
HTML report: fastp_report.html

fastp v1.0.1, time used: 11 seconds
```

## Step 3: Run fasterp

Process the same data with fasterp:

```bash
fasterp \
  -i SRR5808766_1.fastq.gz -I SRR5808766_2.fastq.gz \
  -o fasterp_out_R1.fq.gz -O fasterp_out_R2.fq.gz \
  -j fasterp_report.json \
  --html fasterp_report.html
```

Expected output:
```
Detecting adapter sequence for read1...
No adapter detected for read1

Detecting adapter sequence for read2...
No adapter detected for read2

Read1 before filtering:
total reads: 9799076
total bases: 499752876
Q20 bases: 494840792(99.0171%)
Q30 bases: 490032403(98.0549%)
Q40 bases: 333003632(66.6337%)

Read2 before filtering:
total reads: 9799076
total bases: 499752876
Q20 bases: 491032411(98.2550%)
Q30 bases: 485393847(97.1268%)
Q40 bases: 335965361(67.2263%)

Read1 after filtering:
total reads: 9761787
total bases: 485985544
Q20 bases: 481444894(99.0657%)
Q30 bases: 476939326(98.1386%)
Q40 bases: 324904111(66.8547%)

Read2 after filtering:
total reads: 9761787
total bases: 485985544
Q20 bases: 479080210(98.5791%)
Q30 bases: 473663787(97.4646%)
Q40 bases: 327920585(67.4754%)

Filtering result:
reads passed filter: 19523574
reads failed due to low quality: 74384
reads failed due to too many N: 194
reads failed due to too short: 0
reads with adapter trimmed: 3113628
bases trimmed due to adapters: 23734734

Duplication rate: 28.3885%

Insert size peak (evaluated by paired-end reads): 51

JSON report: fasterp_report.json
HTML report: fasterp_report.html

fasterp v0.1.0, time used: 6 seconds
```

Note: fasterp processed ~10 million read pairs in 6 seconds vs fastp's 11 seconds.

## Step 4: Compare Outputs

### Verify identical FASTQ output

When comparing gzipped files, the compressed hashes may differ due to compression metadata. Compare the decompressed content:

```bash
# Compare decompressed R1 content
gunzip -c fastp_out_R1.fq.gz | shasum -a 256
gunzip -c fasterp_out_R1.fq.gz | shasum -a 256

# Compare decompressed R2 content
gunzip -c fastp_out_R2.fq.gz | shasum -a 256
gunzip -c fasterp_out_R2.fq.gz | shasum -a 256
```

Expected output:
```
716d6bc9b5aa075e5f7ff527b1638de1ea56b67439b7a7646bfe25fb14d132e7  -
716d6bc9b5aa075e5f7ff527b1638de1ea56b67439b7a7646bfe25fb14d132e7  -

c28e4e14c58ded4f4a7b776f6fb9f92cabc83056253f0ff4e268e5db1653b656  -
c28e4e14c58ded4f4a7b776f6fb9f92cabc83056253f0ff4e268e5db1653b656  -
```

The hashes match - both tools produced identical output.

### Compare key statistics

```bash
# Total reads before/after filtering
echo "=== fastp ==="
jq '{
  before: .summary.before_filtering.total_reads,
  after: .summary.after_filtering.total_reads,
  passed_rate: .summary.after_filtering.total_reads / .summary.before_filtering.total_reads * 100
}' fastp_report.json

echo "=== fasterp ==="
jq '{
  before: .summary.before_filtering.total_reads,
  after: .summary.after_filtering.total_reads,
  passed_rate: .summary.after_filtering.total_reads / .summary.before_filtering.total_reads * 100
}' fasterp_report.json
```

Expected output:
```
=== fastp ===
{
  "before": 19598152,
  "after": 19523574,
  "passed_rate": 99.61946412090282
}

=== fasterp ===
{
  "before": 19598152,
  "after": 19523574,
  "passed_rate": 99.61946412090282
}
```

### Compare k-mer counts

```bash
# Top k-mers should be identical
echo "=== fastp k-mers ==="
jq '.read1_before_filtering.kmer_count | to_entries | sort_by(-.value) | .[0:5]' fastp_report.json

echo "=== fasterp k-mers ==="
jq '.read1_before_filtering.kmer_count | to_entries | sort_by(-.value) | .[0:5]' fasterp_report.json
```

Expected output:
```
=== fastp k-mers ===
[
  { "key": "CTGTC", "value": 1651445 },
  { "key": "TGTCT", "value": 1651009 },
  { "key": "TCTCT", "value": 1576849 },
  { "key": "GTCTC", "value": 1358458 },
  { "key": "CTCTT", "value": 1314738 }
]

=== fasterp k-mers ===
[
  { "key": "CTGTC", "value": 1651445 },
  { "key": "TGTCT", "value": 1651009 },
  { "key": "TCTCT", "value": 1576849 },
  { "key": "GTCTC", "value": 1358458 },
  { "key": "CTCTT", "value": 1314738 }
]
```

## Step 5: View Reports

Open the HTML reports in your browser to compare visualizations:

```bash
# macOS
open fastp_report.html fasterp_report.html

# Linux
xdg-open fastp_report.html && xdg-open fasterp_report.html

# Or start a local server
python3 -m http.server 8000
# Then visit http://localhost:8000
```

### Example Reports

Here are the actual reports generated from this tutorial:

| Tool | Report |
|------|--------|
| **fastp** | [fastp_report.html](./fastp_report.html) |
| **fasterp** | [fasterp_report.html](./fasterp_report.html) |

<div style="display: flex; gap: 1em; margin: 1em 0;">
  <div style="flex: 1; text-align: center;">
    <strong>fastp</strong>
    <a href="./fastp_report.html" target="_blank" style="display: block; cursor: pointer;">
      <iframe src="./fastp_report.html" style="width: 100%; height: 300px; border: 1px solid #ccc; border-radius: 4px; pointer-events: none;"></iframe>
    </a>
  </div>
  <div style="flex: 1; text-align: center;">
    <strong>fasterp</strong>
    <a href="./fasterp_report.html" target="_blank" style="display: block; cursor: pointer;">
      <iframe src="./fasterp_report.html" style="width: 100%; height: 300px; border: 1px solid #ccc; border-radius: 4px; pointer-events: none;"></iframe>
    </a>
  </div>
</div>

Compare these side-by-side to see identical quality metrics, base content graphs, and k-mer distributions.

## Cleanup

```bash
# Remove generated files
rm -f *.gz fastp_report.* fasterp_report.*
cd .. && rmdir fasterp_tutorial
```

## Next Steps

- Try with your own FASTQ data
- Experiment with different filtering parameters
- Check the [CLI Reference](./README.md#cli-reference) for all available options
