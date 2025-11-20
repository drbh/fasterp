# fasterp

High-performance FASTQ preprocessing in Rust. A drop-in replacement for [fastp](https://github.com/OpenGene/fastp) with blazing fast performance.

## Installation

```bash
cargo install fasterp
```

Or build from source:

```bash
git clone https://github.com/drbh/fasterp
cd fasterp
cargo build --release
```

## Quick Start

```bash
# Single-end
fasterp -i input.fq -o output.fq

# Paired-end
fasterp -i R1.fq -I R2.fq -o out1.fq -O out2.fq

# With quality filtering
fasterp -i input.fq -o output.fq -q 20 -l 50

# With adapter trimming
fasterp -i input.fq -o output.fq -a AGATCGGAAGAGC

# Multi-threaded with memory limit
fasterp -i large.fq -o output.fq -w 8 --max-memory 4096
```

## Quick Reference

Copy-paste commands for common tasks:

| Task                        | Command                                                                   |
| --------------------------- | ------------------------------------------------------------------------- |
| **Basic cleanup**           | `fasterp -i in.fq -o out.fq`                                              |
| **Strict quality filter**   | `fasterp -i in.fq -o out.fq -q 20 -u 30 -l 50`                            |
| **Aggressive trimming**     | `fasterp -i in.fq -o out.fq --cut-front --cut-tail --cut-mean-quality 20` |
| **Remove adapters only**    | `fasterp -i in.fq -o out.fq -a AGATCGGAAGAGC`                             |
| **Paired-end basic**        | `fasterp -i R1.fq -I R2.fq -o o1.fq -O o2.fq`                             |
| **Paired-end + merge**      | `fasterp -i R1.fq -I R2.fq -o o1.fq -O o2.fq -m --merged-out merged.fq`   |
| **Paired-end + correction** | `fasterp -i R1.fq -I R2.fq -o o1.fq -O o2.fq -c`                          |
| **Fast (16 threads)**       | `fasterp -i in.fq -o out.fq -w 16`                                        |
| **Memory-limited (4GB)**    | `fasterp -i in.fq -o out.fq --max-memory 4096`                            |
| **Compressed in/out**       | `fasterp -i in.fq.gz -o out.fq.gz`                                        |
| **Split into 10 files**     | `fasterp -i in.fq -o out.fq -s 10`                                        |
| **Skip all filtering**      | `fasterp -i in.fq -o out.fq -A -G -Q -L`                                  |
| **With UMI extraction**     | `fasterp -i in.fq -o out.fq --umi --umi-len 8`                            |
| **Deduplicate**             | `fasterp -i in.fq -o out.fq -D`                                           |
| **Custom report name**      | `fasterp -i in.fq -o out.fq -j report.json --html report.html`            |

**Pro tip**: Combine flags freely. Example for production pipeline:
```bash
fasterp -i R1.fq.gz -I R2.fq.gz -o o1.fq.gz -O o2.fq.gz \
  -q 20 -l 50 --cut-tail --cut-mean-quality 20 \
  -c -w 16 --max-memory 8192
```

## CLI Reference

### Input/Output

FASTQ files contain sequence data as text. Each "read" is a short string of letters (A, C, G, T) representing a fragment of genetic material, plus a confidence score for each letter. fasterp reads these files, processes them, and writes cleaned versions.

| Flag | Long               | Default | Description                                                               |
| ---- | ------------------ | ------- | ------------------------------------------------------------------------- |
| `-i` | `--in1`            |         | Input FASTQ file or URL (use '-' for stdin, supports .gz and http(s)://)  |
| `-o` | `--out1`           |         | Output FASTQ file (use '-' for stdout, .gz extension enables compression) |
| `-I` | `--in2`            |         | Read2 input file or URL (paired-end mode, supports .gz and http(s)://)    |
| `-O` | `--out2`           |         | Read2 output file (paired-end mode)                                       |
|      | `--unpaired1`      |         | Output file for unpaired read1                                            |
|      | `--unpaired2`      |         | Output file for unpaired read2                                            |
|      | `--interleaved-in` |         | Input is interleaved paired-end                                           |

<svg viewBox="0 0 700 90" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:700px;font-family:system-ui,sans-serif;font-size:14px">
  <defs>
    <marker id="a1" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/>
    </marker>
  </defs>
  <rect x="20" y="30" width="120" height="40" rx="6" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="80" y="56" text-anchor="middle" font-weight="600" fill="#1e40af">input.fq.gz</text>
  <line x1="140" y1="50" x2="195" y2="50" stroke="#333" stroke-width="2" marker-end="url(#a1)"/>
  <rect x="200" y="30" width="120" height="40" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="260" y="56" text-anchor="middle" font-weight="600" fill="#92400e">Process</text>
  <line x1="320" y1="50" x2="375" y2="50" stroke="#333" stroke-width="2" marker-end="url(#a1)"/>
  <rect x="380" y="30" width="120" height="40" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="440" y="56" text-anchor="middle" font-weight="600" fill="#065f46">output.fq.gz</text>
  <line x1="500" y1="50" x2="555" y2="50" stroke="#333" stroke-width="2" marker-end="url(#a1)"/>
  <rect x="560" y="30" width="120" height="40" rx="6" fill="#ede9fe" stroke="#7c3aed" stroke-width="2"/>
  <text x="620" y="56" text-anchor="middle" font-weight="600" fill="#5b21b6">report.json</text>
</svg>

### Length Filtering

Reads that are too short often come from errors or failed sequencing. Think of it like filtering out incomplete sentences - if a read is shorter than your minimum threshold, it's probably not useful data and gets discarded.

| Flag | Long                      | Default | Description                                              |
| ---- | ------------------------- | ------- | -------------------------------------------------------- |
| `-l` | `--length-required`       | 15      | Minimum length required                                  |
| `-L` | `--disable-length-filter` |         | Disable length filtering (fastp compatibility)           |
| `-b` | `--max-len1`              | 0       | Maximum length for read1 - trim if longer (0 = disabled) |
| `-B` | `--max-len2`              | 0       | Maximum length for read2 - trim if longer (0 = disabled) |

<svg viewBox="0 0 460 170" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:460px;font-family:system-ui,sans-serif;font-size:13px">
  <defs>
    <marker id="a2" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/>
    </marker>
  </defs>
  <!-- Read input -->
  <rect x="20" y="50" width="100" height="40" rx="6" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="70" y="76" text-anchor="middle" font-weight="600" fill="#1e40af">Read</text>
  <line x1="120" y1="70" x2="165" y2="70" stroke="#333" stroke-width="2" marker-end="url(#a2)"/>
  <!-- Min check -->
  <rect x="170" y="50" width="120" height="40" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="230" y="76" text-anchor="middle" font-weight="600" fill="#92400e">len ≥ 15?</text>
  <!-- Yes path -->
  <line x1="290" y1="70" x2="335" y2="70" stroke="#059669" stroke-width="2" marker-end="url(#a2)"/>
  <text x="312" y="62" fill="#059669" font-weight="600">yes</text>
  <!-- Keep -->
  <rect x="340" y="50" width="100" height="40" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="390" y="76" text-anchor="middle" font-weight="600" fill="#065f46">Keep</text>
  <!-- No path -->
  <line x1="230" y1="90" x2="230" y2="115" stroke="#dc2626" stroke-width="2" marker-end="url(#a2)"/>
  <text x="245" y="108" fill="#dc2626" font-weight="600">no</text>
  <!-- Discard -->
  <rect x="170" y="120" width="120" height="35" rx="6" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="230" y="143" text-anchor="middle" font-weight="600" fill="#991b1b">Discard</text>
</svg>

### Quality Filtering

Each letter in a read has a confidence score (0-40+) indicating how certain the machine was about that letter. Higher is better. Quality filtering removes reads with too many low-confidence letters, too many unknown letters (N), or repetitive patterns that suggest errors. It's like spell-checking - removing text that's too garbled to trust.

| Flag | Long                          | Default | Description                                                       |
| ---- | ----------------------------- | ------- | ----------------------------------------------------------------- |
| `-q` | `--qualified-quality-phred`   | 15      | Quality value that a base is qualified (phred >= this)            |
| `-u` | `--unqualified-percent-limit` | 40      | Percent of bases allowed to be unqualified (0-100)                |
| `-e` | `--average-qual`              | 0       | Average quality threshold - discard if mean < this (0 = disabled) |
| `-n` | `--n-base-limit`              | 5       | Max number of N bases allowed                                     |
| `-y` | `--low-complexity-filter`     |         | Enable low complexity filter                                      |
| `-Y` | `--complexity-threshold`      | 30      | Complexity threshold (0-100), 30 = 30% complexity required        |

<svg viewBox="0 0 500 220" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:500px;font-family:system-ui,sans-serif;font-size:13px">
  <defs>
    <marker id="a3" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/>
    </marker>
  </defs>
  <!-- Read -->
  <rect x="180" y="10" width="100" height="35" rx="6" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="230" y="33" text-anchor="middle" font-weight="600" fill="#1e40af">Read</text>
  <line x1="230" y1="45" x2="230" y2="60" stroke="#333" stroke-width="2" marker-end="url(#a3)"/>
  <!-- N check -->
  <rect x="160" y="65" width="140" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="230" y="88" text-anchor="middle" font-weight="600" fill="#92400e">N bases ≤ 5?</text>
  <line x1="230" y1="100" x2="230" y2="115" stroke="#059669" stroke-width="2" marker-end="url(#a3)"/>
  <!-- Quality check -->
  <rect x="140" y="120" width="180" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="230" y="143" text-anchor="middle" font-weight="600" fill="#92400e">Mean quality ≥ threshold?</text>
  <line x1="230" y1="155" x2="230" y2="170" stroke="#059669" stroke-width="2" marker-end="url(#a3)"/>
  <!-- Pass -->
  <rect x="170" y="175" width="120" height="35" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="230" y="198" text-anchor="middle" font-weight="600" fill="#065f46">Pass Filter</text>
  <!-- Discard arrows -->
  <line x1="300" y1="82" x2="380" y2="82" stroke="#dc2626" stroke-width="2" marker-end="url(#a3)"/>
  <line x1="320" y1="137" x2="380" y2="137" stroke="#dc2626" stroke-width="2" marker-end="url(#a3)"/>
  <rect x="385" y="95" width="90" height="35" rx="6" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="430" y="118" text-anchor="middle" font-weight="600" fill="#991b1b">Discard</text>
</svg>

### Trimming - Fixed Position

Sometimes the beginning or end of reads contain known bad data (like primer sequences or low-quality regions). Fixed trimming removes a set number of characters from the start or end of every read - like cutting off the first and last few words of every sentence because you know they're always wrong.

| Flag | Long            | Default | Description                       |
| ---- | --------------- | ------- | --------------------------------- |
|      | `--trim-front`  | 0       | Trim N bases from 5' (front) end  |
|      | `--trim-tail`   | 0       | Trim N bases from 3' (tail) end   |
| `-f` | `--trim-front1` | 0       | Trim N bases from front for read1 |
| `-t` | `--trim-tail1`  | 0       | Trim N bases from tail for read1  |
| `-F` | `--trim-front2` | 0       | Trim N bases from front for read2 |
| `-T` | `--trim-tail2`  | 0       | Trim N bases from tail for read2  |

<svg viewBox="0 0 500 130" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:500px;font-family:system-ui,sans-serif;font-size:14px">
  <!-- Before -->
  <text x="20" y="30" font-weight="700" fill="#374151">Before:</text>
  <text x="100" y="30" font-family="monospace" fill="#6b7280">5'─</text>
  <rect x="130" y="15" width="40" height="24" rx="4" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="150" y="32" text-anchor="middle" font-family="monospace" font-weight="600" fill="#991b1b">NNN</text>
  <text x="175" y="30" font-family="monospace" fill="#059669" font-weight="600">ACGTACGTACGT</text>
  <rect x="315" y="15" width="30" height="24" rx="4" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="330" y="32" text-anchor="middle" font-family="monospace" font-weight="600" fill="#991b1b">NN</text>
  <text x="350" y="30" font-family="monospace" fill="#6b7280">─3'</text>
  <!-- Arrow -->
  <text x="230" y="65" text-anchor="middle" font-weight="600" fill="#d97706">-f 3 -t 2</text>
  <line x1="230" y1="70" x2="230" y2="85" stroke="#d97706" stroke-width="2"/>
  <path d="M 225 80 L 230 90 L 235 80" fill="#d97706"/>
  <!-- After -->
  <text x="20" y="110" font-weight="700" fill="#374151">After:</text>
  <text x="100" y="110" font-family="monospace" fill="#6b7280">5'─</text>
  <rect x="130" y="95" width="145" height="24" rx="4" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="202" y="112" text-anchor="middle" font-family="monospace" font-weight="600" fill="#065f46">ACGTACGTACGT</text>
  <text x="280" y="110" font-family="monospace" fill="#6b7280">─3'</text>
</svg>

### Trimming - Quality Based

Instead of cutting a fixed amount, this looks at actual quality scores and trims where the data gets bad. A sliding window moves along the read, checking average quality. When it drops below your threshold, everything past that point is cut off. It's like editing a document and stopping where the text becomes unreadable.

| Flag | Long                         | Default | Description                                               |
| ---- | ---------------------------- | ------- | --------------------------------------------------------- |
|      | `--cut-mean-quality`         | 0       | Quality cutoff for sliding-window trimming (0 = disabled) |
|      | `--cut-window-size`          | 4       | Sliding window size for quality trimming                  |
|      | `--cut-front`                |         | Enable quality trimming at 5' end (front)                 |
|      | `--cut-tail`                 |         | Enable quality trimming at 3' end (tail)                  |
|      | `--disable-quality-trimming` |         | Disable polyG/quality tail trimming                       |

<svg viewBox="0 0 550 130" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:550px;font-family:system-ui,sans-serif;font-size:13px">
  <!-- Quality scores -->
  <text x="20" y="25" font-weight="700" fill="#374151">Quality scores:</text>
  <text x="140" y="25" font-family="monospace" fill="#dc2626" font-weight="600">10 12 14</text>
  <text x="220" y="25" font-family="monospace" fill="#059669" font-weight="600">35 36 38 37 35</text>
  <text x="360" y="25" font-family="monospace" fill="#dc2626" font-weight="600">20 15 10 8</text>
  <!-- Highlight good region -->
  <rect x="215" y="8" width="140" height="24" rx="4" fill="#d1fae5" stroke="#059669" stroke-width="2" fill-opacity="0.5"/>
  <!-- Window -->
  <rect x="355" y="35" width="80" height="30" rx="4" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="395" y="55" text-anchor="middle" font-size="11" fill="#92400e">window=4</text>
  <!-- Cut line -->
  <line x1="355" y1="8" x2="355" y2="75" stroke="#dc2626" stroke-width="2" stroke-dasharray="4"/>
  <text x="365" y="85" fill="#dc2626" font-weight="600" font-size="12">cut here</text>
  <!-- Result -->
  <text x="20" y="115" font-weight="700" fill="#374151">Result:</text>
  <text x="100" y="115" fill="#374151">Keep region where sliding window mean ≥ 20</text>
</svg>

### PolyG/PolyX Trimming

Some sequencing machines produce false runs of the same letter (like "GGGGGGGG") at the end of reads when they run out of real signal. These aren't real data - they're artifacts. This feature detects and removes these repetitive tails, like removing a stuck key's output from a document.

| Flag | Long                    | Default | Description                                          |
| ---- | ----------------------- | ------- | ---------------------------------------------------- |
|      | `--trim-poly-g`         |         | Enable polyG tail trimming                           |
| `-G` | `--disable-trim-poly-g` |         | Disable polyG tail trimming                          |
|      | `--trim-poly-x`         |         | Enable generic polyX tail trimming (any homopolymer) |
|      | `--poly-g-min-len`      | 10      | Minimum length for polyG/polyX detection             |

<svg viewBox="0 0 500 130" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:500px;font-family:system-ui,sans-serif;font-size:14px">
  <!-- Before -->
  <text x="20" y="30" font-weight="700" fill="#374151">Before:</text>
  <text x="100" y="30" font-family="monospace" fill="#059669" font-weight="600">ACGTACGT</text>
  <rect x="195" y="15" width="150" height="24" rx="4" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="270" y="32" text-anchor="middle" font-family="monospace" font-weight="600" fill="#991b1b">GGGGGGGGGGGG</text>
  <!-- Label -->
  <text x="360" y="30" font-size="11" fill="#6b7280">(polyG tail)</text>
  <!-- Arrow -->
  <line x1="200" y1="50" x2="200" y2="75" stroke="#d97706" stroke-width="2"/>
  <path d="M 195 70 L 200 80 L 205 70" fill="#d97706"/>
  <!-- After -->
  <text x="20" y="105" font-weight="700" fill="#374151">After:</text>
  <rect x="100" y="90" width="95" height="24" rx="4" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="147" y="107" text-anchor="middle" font-family="monospace" font-weight="600" fill="#065f46">ACGTACGT</text>
  <text x="210" y="105" font-size="11" fill="#6b7280">(NovaSeq artifact removed)</text>
</svg>

### Adapter Trimming

Adapters are short synthetic sequences added during sample preparation - they're not part of the actual sample. They appear at the ends of reads and must be removed, like stripping metadata headers from a file. fasterp can auto-detect common adapters or you can specify them manually.

| Flag | Long                          | Default | Description                                        |
| ---- | ----------------------------- | ------- | -------------------------------------------------- |
| `-A` | `--disable-adapter-trimming`  |         | Disable adapter trimming                           |
| `-a` | `--adapter-sequence`          |         | Adapter for read1 (auto-detected if not specified) |
|      | `--adapter-sequence-r2`       |         | Adapter sequence for read2                         |
|      | `--disable-adapter-detection` |         | Disable adapter auto-detection                     |

<svg viewBox="0 0 550 150" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:550px;font-family:system-ui,sans-serif;font-size:14px">
  <!-- Auto-detect label -->
  <rect x="20" y="10" width="500" height="30" rx="6" fill="#f3f4f6" stroke="#9ca3af" stroke-width="1"/>
  <text x="270" y="30" text-anchor="middle" font-size="12" fill="#4b5563">Auto-detect: Sample reads → K-mer analysis → Find adapter sequence</text>
  <!-- Before -->
  <text x="20" y="70" font-weight="700" fill="#374151">Before:</text>
  <text x="100" y="70" font-family="monospace" fill="#059669" font-weight="600">ACGTACGT</text>
  <rect x="195" y="55" width="160" height="24" rx="4" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="275" y="72" text-anchor="middle" font-family="monospace" font-weight="600" fill="#991b1b">AGATCGGAAGAGC</text>
  <text x="370" y="70" font-size="11" fill="#6b7280">(adapter)</text>
  <!-- Arrow -->
  <line x1="200" y1="90" x2="200" y2="105" stroke="#d97706" stroke-width="2"/>
  <path d="M 195 100 L 200 110 L 205 100" fill="#d97706"/>
  <!-- After -->
  <text x="20" y="130" font-weight="700" fill="#374151">After:</text>
  <rect x="100" y="115" width="95" height="24" rx="4" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="147" y="132" text-anchor="middle" font-family="monospace" font-weight="600" fill="#065f46">ACGTACGT</text>
</svg>

### Paired-End Options

Paired-end sequencing reads the same fragment from both ends, giving you two reads that overlap in the middle. When they overlap, you can compare them: if one says "A" with high confidence and the other says "T" with low confidence, you trust the "A". You can also merge overlapping pairs into a single, longer read. It's like having two people transcribe the same audio and combining their best parts.

| Flag | Long                           | Default | Description                                          |
| ---- | ------------------------------ | ------- | ---------------------------------------------------- |
| `-c` | `--correction`                 |         | Enable base correction using overlap analysis        |
|      | `--overlap-len-require`        | 30      | Minimum overlap length required for correction       |
|      | `--overlap-diff-limit`         | 5       | Maximum allowed differences in overlap region        |
|      | `--overlap-diff-percent-limit` | 20      | Maximum allowed difference percentage (0-100)        |
| `-m` | `--merge`                      |         | Merge overlapping paired-end reads into single reads |
|      | `--merged-out`                 |         | Output file for merged reads                         |
|      | `--include-unmerged`           |         | Write unmerged reads to merged output file           |

<svg viewBox="0 0 550 210" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:550px;font-family:system-ui,sans-serif;font-size:13px">
  <!-- Overlap section -->
  <text x="20" y="25" font-weight="700" fill="#374151">Overlap Detection:</text>
  <!-- R1 -->
  <text x="20" y="55" font-family="monospace" fill="#374151">R1: 5'─</text>
  <text x="80" y="55" font-family="monospace" fill="#059669" font-weight="600">ACGTACGT</text>
  <rect x="165" y="40" width="80" height="24" rx="4" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="205" y="57" text-anchor="middle" font-family="monospace" fill="#1e40af" font-weight="600">NNNNNN</text>
  <text x="250" y="55" font-family="monospace" fill="#374151">─3'</text>
  <!-- R2 -->
  <text x="20" y="85" font-family="monospace" fill="#374151">R2: 3'─────</text>
  <rect x="165" y="70" width="80" height="24" rx="4" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="205" y="87" text-anchor="middle" font-family="monospace" fill="#1e40af" font-weight="600">NNNNNN</text>
  <text x="250" y="85" font-family="monospace" fill="#059669" font-weight="600">TGCATGCA</text>
  <text x="340" y="85" font-family="monospace" fill="#374151">─5'</text>
  <!-- Overlap label -->
  <text x="205" y="110" text-anchor="middle" font-size="11" fill="#2563eb" font-weight="600">overlap region</text>
  <!-- Correction -->
  <text x="20" y="145" font-weight="700" fill="#374151">Base Correction (-c):</text>
  <text x="180" y="145" font-size="12" fill="#374151">R1: <tspan fill="#059669" font-weight="600">A</tspan> (Q30) vs R2: <tspan fill="#dc2626" font-weight="600">T</tspan> (Q15) → Use <tspan fill="#059669" font-weight="600">A</tspan></text>
  <!-- Merged -->
  <text x="20" y="180" font-weight="700" fill="#374151">Merged (-m):</text>
  <rect x="120" y="165" width="280" height="24" rx="4" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="260" y="182" text-anchor="middle" font-family="monospace" font-size="12" fill="#065f46" font-weight="600">5'─ACGTACGT──────TGCATGCA─3'</text>
</svg>

### UMI Processing

A UMI (Unique Molecular Identifier) is a short random tag added to each original molecule before copying. Since copies of the same molecule have the same UMI, you can later identify and remove duplicates that came from the same original. Think of it like adding a unique serial number to each document before photocopying - you can then tell which copies came from which original.

| Flag | Long           | Default | Description                                                |
| ---- | -------------- | ------- | ---------------------------------------------------------- |
|      | `--umi`        |         | Enable UMI preprocessing                                   |
|      | `--umi-loc`    | read1   | UMI location: read1/read2/index1/index2/per_read/per_index |
|      | `--umi-len`    | 0       | UMI length (required when --umi is enabled)                |
|      | `--umi-prefix` | UMI     | Prefix in read name for UMI                                |
|      | `--umi-skip`   | 0       | Skip N bases before UMI                                    |

<svg viewBox="0 0 500 130" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:500px;font-family:system-ui,sans-serif;font-size:14px">
  <!-- Before -->
  <text x="20" y="30" font-weight="700" fill="#374151">Before:</text>
  <rect x="100" y="15" width="90" height="24" rx="4" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="145" y="32" text-anchor="middle" font-family="monospace" font-weight="600" fill="#1e40af">ACGTACGT</text>
  <text x="200" y="30" font-family="monospace" fill="#6b7280">NNNNNNNNNNN...</text>
  <text x="100" y="50" font-size="11" fill="#2563eb" font-weight="600">UMI (8bp)</text>
  <!-- Arrow -->
  <line x1="200" y1="60" x2="200" y2="75" stroke="#d97706" stroke-width="2"/>
  <path d="M 195 70 L 200 80 L 205 70" fill="#d97706"/>
  <!-- After -->
  <text x="20" y="105" font-weight="700" fill="#374151">After:</text>
  <text x="100" y="105" font-family="monospace" fill="#374151">@READ_NAME</text>
  <rect x="220" y="90" width="130" height="24" rx="4" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="285" y="107" text-anchor="middle" font-family="monospace" font-weight="600" fill="#1e40af">UMI:ACGTACGT</text>
</svg>

### Deduplication

During sample preparation, molecules get copied many times. Identical reads are likely copies of the same original, which can skew your analysis (like counting the same vote multiple times). Deduplication detects and removes these duplicates using a memory-efficient probabilistic filter. It estimates duplication rates without needing to store every read in memory.

| Flag | Long                      | Default | Description                                         |
| ---- | ------------------------- | ------- | --------------------------------------------------- |
| `-D` | `--dedup`                 |         | Enable deduplication to remove duplicate reads      |
|      | `--dup-calc-accuracy`     | 0       | Deduplication accuracy (1-6), higher = more memory  |
|      | `--dont-eval-duplication` |         | Don't evaluate duplication rate (saves time/memory) |

<svg viewBox="0 0 550 150" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:550px;font-family:system-ui,sans-serif;font-size:13px">
  <defs>
    <marker id="a4" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/>
    </marker>
  </defs>
  <!-- Read -->
  <rect x="20" y="30" width="80" height="35" rx="6" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="60" y="53" text-anchor="middle" font-weight="600" fill="#1e40af">Read</text>
  <line x1="100" y1="47" x2="140" y2="47" stroke="#333" stroke-width="2" marker-end="url(#a4)"/>
  <!-- Hash -->
  <rect x="145" y="30" width="80" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="185" y="53" text-anchor="middle" font-weight="600" fill="#92400e">Hash</text>
  <line x1="225" y1="47" x2="265" y2="47" stroke="#333" stroke-width="2" marker-end="url(#a4)"/>
  <!-- Bloom check -->
  <rect x="270" y="30" width="120" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="330" y="53" text-anchor="middle" font-weight="600" fill="#92400e">In Bloom?</text>
  <!-- Yes = discard -->
  <line x1="390" y1="47" x2="440" y2="47" stroke="#dc2626" stroke-width="2" marker-end="url(#a4)"/>
  <text x="415" y="40" fill="#dc2626" font-weight="600">yes</text>
  <rect x="445" y="30" width="80" height="35" rx="6" fill="#fee2e2" stroke="#dc2626" stroke-width="2"/>
  <text x="485" y="53" text-anchor="middle" font-weight="600" fill="#991b1b">Discard</text>
  <!-- No = add and keep -->
  <line x1="330" y1="65" x2="330" y2="95" stroke="#059669" stroke-width="2" marker-end="url(#a4)"/>
  <text x="345" y="85" fill="#059669" font-weight="600">no</text>
  <rect x="270" y="100" width="120" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="330" y="123" text-anchor="middle" font-weight="600" fill="#92400e">Add to filter</text>
  <line x1="390" y1="117" x2="440" y2="117" stroke="#333" stroke-width="2" marker-end="url(#a4)"/>
  <rect x="445" y="100" width="80" height="35" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="485" y="123" text-anchor="middle" font-weight="600" fill="#065f46">Keep</text>
</svg>

### Output Splitting

Large datasets can be split into multiple smaller files for parallel downstream processing. You can split by total file count (distribute reads evenly across N files) or by lines per file (each file gets up to N lines). This is useful when you want to process chunks in parallel on different machines.

| Flag | Long                    | Default | Description                                            |
| ---- | ----------------------- | ------- | ------------------------------------------------------ |
| `-s` | `--split`               |         | Split output by limiting total number of files (2-999) |
| `-S` | `--split-by-lines`      |         | Split output by limiting lines per file (>=1000)       |
| `-d` | `--split-prefix-digits` | 4       | Digits for file number padding (1-10)                  |

<svg viewBox="0 0 320 120" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:320px;font-family:system-ui,sans-serif;font-size:13px">
  <defs>
    <marker id="a5" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/>
    </marker>
  </defs>
  <!-- Input -->
  <rect x="20" y="40" width="100" height="40" rx="6" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="70" y="66" text-anchor="middle" font-weight="600" fill="#1e40af">input.fq</text>
  <!-- Arrows -->
  <line x1="120" y1="50" x2="180" y2="25" stroke="#333" stroke-width="2" marker-end="url(#a5)"/>
  <line x1="120" y1="60" x2="180" y2="60" stroke="#333" stroke-width="2" marker-end="url(#a5)"/>
  <line x1="120" y1="70" x2="180" y2="95" stroke="#333" stroke-width="2" marker-end="url(#a5)"/>
  <!-- Output files -->
  <rect x="185" y="10" width="120" height="28" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="245" y="29" text-anchor="middle" font-weight="600" font-size="12" fill="#065f46">out.0001.fq</text>
  <rect x="185" y="45" width="120" height="28" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="245" y="64" text-anchor="middle" font-weight="600" font-size="12" fill="#065f46">out.0002.fq</text>
  <rect x="185" y="80" width="120" height="28" rx="6" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="245" y="99" text-anchor="middle" font-weight="600" font-size="12" fill="#065f46">out.0003.fq</text>
</svg>

### Reporting

fasterp generates detailed statistics about your data before and after processing: total reads, quality distributions, how many were filtered out and why. The JSON report is machine-readable for pipelines; the HTML report has interactive charts for visual inspection.

| Flag | Long             | Default      | Description                                         |
| ---- | ---------------- | ------------ | --------------------------------------------------- |
| `-j` | `--json`         | fasterp.json | JSON report file                                    |
|      | `--html`         | fasterp.html | HTML report file                                    |
|      | `--stats-format` | compact      | Stats output format: compact, pretty, off, or jsonl |

### Performance

Control how fasterp uses your system resources. More threads = faster processing but more CPU/memory usage. The batch size affects memory consumption - larger batches are more efficient but use more RAM. Use `--max-memory` to set a limit and let fasterp auto-tune the other settings.

| Flag | Long                     | Default   | Description                              |
| ---- | ------------------------ | --------- | ---------------------------------------- |
| `-w` | `--thread`               | auto      | Number of worker threads                 |
| `-z` | `--compression`          | 6         | Compression level for gzip output (0-9)  |
|      | `--parallel-compression` | true      | Use parallel compression for gzip output |
|      | `--batch-bytes`          | 33554432  | Batch size in bytes (32 MiB)             |
|      | `--max-backlog`          | threads+1 | Maximum backlog of batches               |
|      | `--max-memory`           |           | Maximum memory usage in MB               |
|      | `--skip-kmer-counting`   |           | Skip k-mer counting for performance      |

## Architecture

fasterp processes data in a pipeline with three stages running in parallel. The Producer reads chunks from disk. Multiple Workers process those chunks simultaneously (filtering, trimming, counting). The Merger collects results in order and writes output. This design keeps all CPU cores busy while maintaining correct output order.

<svg viewBox="0 0 650 160" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:1em auto;max-width:650px;font-family:system-ui,sans-serif;font-size:12px">
  <defs>
    <marker id="a6" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/>
    </marker>
  </defs>
  <!-- Producer -->
  <rect x="20" y="35" width="120" height="80" rx="8" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>
  <text x="80" y="60" text-anchor="middle" font-weight="700" fill="#1e40af">Producer</text>
  <text x="80" y="80" text-anchor="middle" font-size="10" fill="#1e40af">Read &amp; Parse</text>
  <text x="80" y="95" text-anchor="middle" font-size="10" fill="#1e40af">32MB batches</text>
  <!-- Arrows to workers -->
  <line x1="140" y1="55" x2="190" y2="35" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <line x1="140" y1="75" x2="190" y2="75" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <line x1="140" y1="95" x2="190" y2="115" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <!-- Workers -->
  <rect x="195" y="15" width="100" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="245" y="38" text-anchor="middle" font-weight="600" fill="#92400e">Worker 1</text>
  <rect x="195" y="57" width="100" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="245" y="80" text-anchor="middle" font-weight="600" fill="#92400e">Worker 2</text>
  <rect x="195" y="100" width="100" height="35" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="2"/>
  <text x="245" y="123" text-anchor="middle" font-weight="600" fill="#92400e">Worker N</text>
  <!-- Arrows to merger -->
  <line x1="295" y1="32" x2="345" y2="55" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <line x1="295" y1="75" x2="345" y2="75" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <line x1="295" y1="117" x2="345" y2="95" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <!-- Merger -->
  <rect x="350" y="35" width="120" height="80" rx="8" fill="#d1fae5" stroke="#059669" stroke-width="2"/>
  <text x="410" y="60" text-anchor="middle" font-weight="700" fill="#065f46">Merger</text>
  <text x="410" y="80" text-anchor="middle" font-size="10" fill="#065f46">Reorder &amp; Write</text>
  <text x="410" y="95" text-anchor="middle" font-size="10" fill="#065f46">Aggregate stats</text>
  <!-- Output -->
  <line x1="470" y1="60" x2="520" y2="60" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <line x1="470" y1="90" x2="520" y2="90" stroke="#333" stroke-width="2" marker-end="url(#a6)"/>
  <rect x="525" y="45" width="100" height="28" rx="6" fill="#ede9fe" stroke="#7c3aed" stroke-width="2"/>
  <text x="575" y="64" text-anchor="middle" font-weight="600" font-size="11" fill="#5b21b6">output.fq</text>
  <rect x="525" y="78" width="100" height="28" rx="6" fill="#ede9fe" stroke="#7c3aed" stroke-width="2"/>
  <text x="575" y="97" text-anchor="middle" font-weight="600" font-size="11" fill="#5b21b6">report.json</text>
</svg>

fasterp uses a 3-stage pipeline for multi-threaded processing:

1. **Producer**: Reads input in 32MB batches, parses FASTQ records
2. **Workers**: Process batches in parallel (filtering, trimming, stats)
3. **Merger**: Writes output in order, aggregates statistics

Key optimizations:
- SIMD acceleration (AVX2/NEON) for quality statistics
- Zero-copy buffers with `Arc<Vec<u8>>`
- Lookup tables for base/quality checks
- Bloom filter for duplication detection

## Examples

### Quality Control Only

```bash
fasterp -i input.fq -o output.fq -A -G
```

### Aggressive Trimming

```bash
fasterp -i input.fq -o output.fq \
  --cut-front --cut-tail --cut-mean-quality 20 \
  -q 20 -u 30 -l 50
```

### Paired-End with Merging

```bash
fasterp -i R1.fq -I R2.fq -o out1.fq -O out2.fq \
  -m --merged-out merged.fq -c
```

### High-Throughput Pipeline

```bash
fasterp -i large.fq.gz -o output.fq.gz \
  -w 16 --max-memory 8192 -z 4
```

### With UMI

```bash
fasterp -i input.fq -o output.fq \
  --umi --umi-loc read1 --umi-len 12
```

## Comparison with fastp

fasterp produces byte-exact output for standard parameters:

```bash
fastp -i input.fq -o fastp_out.fq
fasterp -i input.fq -o fasterp_out.fq
sha256sum fastp_out.fq fasterp_out.fq  # hashes match
```
