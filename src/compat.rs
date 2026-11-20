//! Native fastp 1.4.0 processing semantics.
//!
//! Batches are processed in bounded waves. Results are emitted in input order;
//! an arbitrarily slow first batch cannot accumulate an unbounded reorder queue.
use crate::util::Instant;
use crate::{adapter, auto_adapter, html, io, overlap, simd, split, stats, Args};
use anyhow::{ensure, Context, Result};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::{json, Value};
use std::borrow::Cow;
use std::collections::{BTreeMap, VecDeque};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

mod block;
#[cfg(feature = "native")]
mod direct_writer;
#[cfg(target_arch = "aarch64")]
mod duplication_neon;
#[cfg(feature = "native")]
mod paired_gzip;

pub const FASTP_VERSION: &str = "1.4.0";

#[derive(Clone, Debug, Default)]
struct Read<'a> {
    // Canonical LF-delimited input, while all four fields still match it.
    verbatim: Option<&'a [u8]>,
    name: Cow<'a, [u8]>,
    seq: Cow<'a, [u8]>,
    plus: Cow<'a, [u8]>,
    qual: Cow<'a, [u8]>,
}

impl<'a> Read<'a> {
    fn validate(&mut self, record: usize, phred64: bool) -> Result<()> {
        let [name, seq, plus, qual] = [&self.name, &self.seq, &self.plus, &self.qual];
        ensure!(
            name.starts_with(b"@"),
            "Invalid FASTQ header at record {}",
            record
        );
        ensure!(
            plus.starts_with(b"+"),
            "Invalid FASTQ separator at record {}",
            record
        );
        ensure!(
            seq.len() == qual.len(),
            "Sequence/quality length mismatch at record {}",
            record
        );
        let offset = if phred64 { 64 } else { 33 };
        ensure!(
            valid_quality(&qual, offset),
            "Invalid quality encoding at record {}",
            record
        );
        if phred64 {
            self.verbatim = None;
            self.qual.to_mut().iter_mut().for_each(|q| *q -= 31);
        }
        Ok(())
    }

    fn name_mut(&mut self) -> &mut Vec<u8> {
        self.verbatim = None;
        self.name.to_mut()
    }

    fn into_owned(self) -> Read<'static> {
        Read {
            verbatim: None,
            name: Cow::Owned(self.name.into_owned()),
            seq: Cow::Owned(self.seq.into_owned()),
            plus: Cow::Owned(self.plus.into_owned()),
            qual: Cow::Owned(self.qual.into_owned()),
        }
    }

    fn capacity(&self) -> usize {
        [&self.name, &self.seq, &self.plus, &self.qual]
            .iter()
            .map(|field| match field {
                Cow::Owned(bytes) => bytes.capacity(),
                Cow::Borrowed(_) => 0,
            })
            .sum()
    }
    fn resize(&mut self, len: usize) {
        if len < self.seq.len() {
            self.verbatim = None;
        }
        for field in [&mut self.seq, &mut self.qual] {
            match field {
                Cow::Borrowed(bytes) => *bytes = &bytes[..len.min(bytes.len())],
                Cow::Owned(bytes) => bytes.truncate(len),
            }
        }
    }

    fn trim_front(&mut self, len: usize) {
        let len = len.min(self.seq.len());
        if len > 0 {
            self.verbatim = None;
        }
        for field in [&mut self.seq, &mut self.qual] {
            match field {
                Cow::Borrowed(bytes) => *bytes = &bytes[len..],
                Cow::Owned(bytes) => {
                    bytes.drain(..len);
                }
            }
        }
    }

    fn borrowed_bases(&self) -> Option<(&'a [u8], &'a [u8])> {
        match (&self.seq, &self.qual) {
            (Cow::Borrowed(seq), Cow::Borrowed(qual)) => Some((seq, qual)),
            _ => None,
        }
    }

    fn write(&self, writer: &mut impl Write, tag: Option<&str>) -> Result<()> {
        if let (Some(bytes), None) = (self.verbatim, tag) {
            writer.write_all(bytes)?;
            return Ok(());
        }
        writer.write_all(&self.name)?;
        if let Some(tag) = tag {
            writer.write_all(b" ")?;
            writer.write_all(tag.as_bytes())?;
        }
        for line in [&b""[..], &self.seq, &self.plus, &self.qual] {
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
        }
        Ok(())
    }
}

/// Strict four-line reader: report malformed/truncated data instead of silently
/// discarding it or letting sequence/quality length mismatches reach SIMD code.
struct Reader {
    inner: Box<dyn BufRead + Send>,
    record: usize,
    phred64: bool,
}

impl Reader {
    fn open(path: &str, phred64: bool) -> Result<Self> {
        Ok(Self {
            inner: io::open_input(path)?,
            record: 0,
            phred64,
        })
    }

    fn next(&mut self) -> Result<Option<Read<'static>>> {
        self.next_reusing(Read::default())
    }

    fn next_reusing(&mut self, mut read: Read<'static>) -> Result<Option<Read<'static>>> {
        Ok(self.read_into(&mut read)?.then_some(read))
    }

    /// fastp skips lines that cannot start a record. Skip only blank and
    /// whitespace-only ones, and only where a header is expected; blank lines
    /// inside a record still fail (fastp silently drops all reads from there).
    fn skip_blank_lines(&mut self) -> Result<()> {
        loop {
            let buffer = self.inner.fill_buf()?;
            if buffer.first().is_none_or(|&b| b == b'@') {
                return Ok(());
            }
            match buffer.iter().position(|&b| !is_space(b)) {
                Some(i) if buffer[i] == b'\n' => self.inner.consume(i + 1),
                Some(0) => return Ok(()),
                Some(_) => anyhow::bail!("Invalid FASTQ header at record {}", self.record + 1),
                None => {
                    // A whitespace run reaching the buffer end: read the rest
                    // of the line to see whether it is blank.
                    let mut line = Vec::new();
                    self.inner.read_until(b'\n', &mut line)?;
                    ensure!(
                        line.iter().all(|&b| is_space(b) || b == b'\n'),
                        "Invalid FASTQ header at record {}",
                        self.record + 1
                    );
                }
            }
        }
    }

    fn read_into(&mut self, read: &mut Read) -> Result<bool> {
        read.verbatim = None;
        let mut lines = [
            read.name.to_mut(),
            read.seq.to_mut(),
            read.plus.to_mut(),
            read.qual.to_mut(),
        ];
        // Ordinary records fit in the existing input buffer. Find all four
        // boundaries with one buffer borrow instead of four virtual reads.
        let mut buffer = self.inner.fill_buf()?;
        if buffer.first() != Some(&b'@') {
            self.skip_blank_lines()?;
            buffer = self.inner.fill_buf()?;
        }
        let mut newlines = memchr::memchr_iter(b'\n', buffer);
        let ends = std::array::from_fn::<_, 4, _>(|_| newlines.next());
        if let [Some(a), Some(b), Some(c), Some(d)] = ends {
            let mut start = 0;
            for (line, end) in lines.iter_mut().zip([a, b, c, d]) {
                let content = &buffer[start..end];
                line.clear();
                line.extend_from_slice(content.strip_suffix(b"\r").unwrap_or(content));
                start = end + 1;
            }
            self.inner.consume(d + 1);
        } else {
            // Preserve records spanning buffers and a final quality line without
            // a newline, including strict truncation diagnostics.
            for (i, line) in lines.iter_mut().enumerate() {
                line.clear();
                let n = self.inner.read_until(b'\n', line)?;
                if n == 0 && i == 0 {
                    return Ok(false);
                }
                ensure!(n > 0, "Truncated FASTQ record {}", self.record + 1);
                if line.last() == Some(&b'\n') {
                    line.pop();
                }
                if line.last() == Some(&b'\r') {
                    line.pop();
                }
            }
        }
        self.record += 1;
        read.validate(self.record, self.phred64)?;
        Ok(true)
    }
}

/// Whitespace other than LF, as C `isspace`.
fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\r' | b'\x0b' | b'\x0c')
}

// A reduction lets the compiler validate a vector of bytes at a time, without
// a branch per quality value. Wrapping subtraction also rejects bytes < offset.
fn valid_quality(qual: &[u8], offset: u8) -> bool {
    qual.iter().fold(0u8, |bad, &q| {
        bad | u8::from(q.wrapping_sub(offset) > 126 - offset)
    }) == 0
}

// A, C, G, T, and all other bytes. One observation touches only one adjacent
// pair of counters; totals including ambiguous bases are derived at reporting.
#[derive(Clone, Copy, Default)]
struct BaseStats {
    count: u64,
    quality_sum: u64,
}
type Cycle = [BaseStats; 5];

#[derive(Default)]
struct KmerWindow {
    code: usize,
    valid: usize,
}

impl KmerWindow {
    fn add(&mut self, base: u8, counts: &mut [usize; 1024]) {
        if base < 4 {
            self.code = ((self.code << 2) | usize::from(base)) & 1023;
            self.valid += 1;
            if self.valid >= 5 {
                counts[self.code] += 1;
            }
        } else {
            self.valid = 0;
            self.code = 0;
        }
    }

    fn add_four(&mut self, bases: [u8; 4], counts: &mut [usize; 1024]) {
        let [a, b, c, d] = bases;
        if self.valid >= 4 && (a | b | c | d) < 4 {
            // The previous four bases and these four contain four overlapping
            // 5-mers. Extract each independently from the packed eight bases.
            let next = (a << 6) | (b << 4) | (c << 2) | d;
            self.code = ((self.code & 255) << 8) | usize::from(next);
            for shift in [6, 4, 2, 0] {
                counts[(self.code >> shift) & 1023] += 1;
            }
        } else {
            for base in bases {
                self.add(base, counts);
            }
        }
    }
}

#[derive(Clone)]
struct Sample {
    simple: stats::SimpleStats,
    cycles: Vec<Cycle>,
    // Separate adjacent updates so long runs of one quality score do not form
    // a dependency chain on a single histogram counter.
    qual_hist: [[usize; 94]; 4],
    kmers: [usize; 1024],
    max_len: usize,
    overrep: Option<Overrep>,
}

impl Default for Sample {
    fn default() -> Self {
        Self {
            simple: stats::SimpleStats::default(),
            cycles: Vec::new(),
            qual_hist: [[0; 94]; 4],
            kmers: [0; 1024],
            max_len: 0,
            overrep: None,
        }
    }
}

impl Sample {
    fn empty_like(&self) -> Self {
        Self {
            overrep: self.overrep.as_ref().map(|ora| Overrep {
                counts: ora.counts.keys().cloned().map(|seq| (seq, 0)).collect(),
                steps: ora.steps,
                sampling: ora.sampling,
            }),
            ..Self::default()
        }
    }

    fn merge(&mut self, other: Self) {
        self.simple.total_reads += other.simple.total_reads;
        self.simple.total_bases += other.simple.total_bases;
        self.simple.q20_bases += other.simple.q20_bases;
        self.simple.q30_bases += other.simple.q30_bases;
        self.simple.q40_bases += other.simple.q40_bases;
        self.simple.gc_bases += other.simple.gc_bases;
        self.max_len = self.max_len.max(other.max_len);
        for (dst, src) in self.kmers.iter_mut().zip(other.kmers) {
            *dst += src;
        }
        if self.cycles.len() < other.cycles.len() {
            self.cycles.resize(other.cycles.len(), Cycle::default());
        }
        for (dst, src) in self.cycles.iter_mut().zip(other.cycles) {
            for (dst, src) in dst.iter_mut().zip(src) {
                dst.count += src.count;
                dst.quality_sum += src.quality_sum;
            }
        }
        for (dst, src) in self.qual_hist.iter_mut().zip(other.qual_hist) {
            for (dst, src) in dst.iter_mut().zip(src) {
                *dst += src;
            }
        }
        if let (Some(dst), Some(src)) = (&mut self.overrep, other.overrep) {
            for (seq, n) in src.counts {
                *dst.counts
                    .get_mut(&seq)
                    .expect("shared overrepresentation candidates") += n;
            }
        }
    }

    fn add_at(&mut self, read: &Read, no_kmer: bool, index: usize) {
        self.add_parts(&read.seq, &read.qual, no_kmer, index);
    }

    fn add_parts(&mut self, seq: &[u8], qual: &[u8], no_kmer: bool, index: usize) {
        if let Some(ora) = self.overrep.as_mut() {
            if index % ora.sampling == 0 {
                ora.add(seq);
            }
        }
        self.simple.total_reads += 1;
        self.simple.total_bases += seq.len();
        self.max_len = self.max_len.max(seq.len());
        if self.cycles.len() < seq.len() {
            self.cycles.resize(seq.len(), Cycle::default());
        }
        let mut window = KmerWindow::default();
        let (bases, remaining_bases) = seq.as_chunks::<4>();
        let (qualities, remaining_qualities) = qual.as_chunks::<4>();
        let (cycles, remaining_cycles) = self.cycles[..seq.len()].as_chunks_mut::<4>();
        let base_index = |b| crate::kmer::base_to_2bit(b).unwrap_or(4) as u8;

        for ((bases, qualities), cycles) in bases.iter().zip(qualities).zip(cycles) {
            let codes = bases.map(base_index);
            for i in 0..4 {
                let q = usize::from(qualities[i] - 33);
                let base = &mut cycles[i][usize::from(codes[i])];
                base.count += 1;
                base.quality_sum += q as u64;
                self.qual_hist[i][q] += 1;
            }
            if !no_kmer {
                window.add_four(codes, &mut self.kmers);
            }
        }
        for (i, ((&b, &q), cycle)) in remaining_bases
            .iter()
            .zip(remaining_qualities)
            .zip(remaining_cycles)
            .enumerate()
        {
            let code = base_index(b);
            let q = usize::from(q - 33);
            let base = &mut cycle[usize::from(code)];
            base.count += 1;
            base.quality_sum += q as u64;
            self.qual_hist[i][q] += 1;
            if !no_kmer {
                window.add(code, &mut self.kmers);
            }
        }
    }

    // A passing prefix contributes to both reports via the shared sample.
    // Only the removed suffix still needs to be counted before filtering.
    fn add_removed_tail(&mut self, seq: &[u8], qual: &[u8], kept: usize, no_kmer: bool) {
        self.simple.total_bases += seq.len() - kept;
        self.max_len = self.max_len.max(seq.len());
        if self.cycles.len() < seq.len() {
            self.cycles.resize(seq.len(), Cycle::default());
        }
        let mut window = KmerWindow::default();
        let code = |b| crate::kmer::base_to_2bit(b).unwrap_or(4) as u8;
        if !no_kmer {
            // Seed the window without counting any 5-mer from the prefix.
            for &b in &seq[kept.saturating_sub(4)..kept] {
                window.add(code(b), &mut self.kmers);
            }
        }
        for pos in kept..seq.len() {
            let b = code(seq[pos]);
            let q = usize::from(qual[pos] - 33);
            let base = &mut self.cycles[pos][usize::from(b)];
            base.count += 1;
            base.quality_sum += q as u64;
            self.qual_hist[pos & 3][q] += 1;
            if !no_kmer {
                window.add(b, &mut self.kmers);
            }
        }
    }

    fn finish_summary(&mut self) {
        // Derive totals from the histograms instead of rescanning reads.
        self.simple.q20_bases = self
            .qual_hist
            .iter()
            .map(|h| h[20..].iter().sum::<usize>())
            .sum();
        self.simple.q30_bases = self
            .qual_hist
            .iter()
            .map(|h| h[30..].iter().sum::<usize>())
            .sum();
        self.simple.q40_bases = self
            .qual_hist
            .iter()
            .map(|h| h[40..].iter().sum::<usize>())
            .sum();
        self.simple.gc_bases = self
            .cycles
            .iter()
            .map(|cycle| (cycle[1].count + cycle[2].count) as usize)
            .sum();
    }

    fn details(&self) -> Value {
        let mut p = stats::PositionStats::new();
        p.ensure_capacity(self.cycles.len());
        p.qual_hist = std::array::from_fn(|q| self.qual_hist.iter().map(|h| h[q]).sum());
        for (i, cycle) in self.cycles.iter().enumerate() {
            p.total_sum[i] = cycle.iter().map(|base| base.quality_sum).sum();
            p.total_cnt[i] = cycle.iter().map(|base| base.count).sum();
            // PositionStats uses A, T, C, G order.
            for (dst, src) in [0, 3, 1, 2].into_iter().enumerate() {
                p.base_sum[dst][i] = cycle[src].quality_sum;
                p.base_cnt[dst][i] = cycle[src].count;
            }
        }
        let mut content = p.to_content_curves();
        // fastp JSON uses fractions. The existing HTML renderer uses percent.
        for curve in [
            &mut content.a,
            &mut content.t,
            &mut content.c,
            &mut content.g,
            &mut content.n,
            &mut content.gc,
        ] {
            for value in curve {
                *value /= 100.0;
            }
        }
        let mut kmers = BTreeMap::new();
        let alphabet = b"ACGT";
        for (i, &count) in self.kmers.iter().enumerate() {
            let key: String = (0..5)
                .rev()
                .map(|s| alphabet[(i >> (s * 2)) & 3] as char)
                .collect();
            kmers.insert(key, count);
        }
        json!({"total_reads":self.simple.total_reads, "total_bases":self.simple.total_bases,
            "q20_bases":self.simple.q20_bases, "q30_bases":self.simple.q30_bases,
            "q40_bases":self.simple.q40_bases, "total_cycles":self.max_len,
            "quality_curves":p.to_quality_curves(), "content_curves":content,
            "qual_hist":p.to_qual_hist(), "kmer_count":kmers,
            "overrepresented_sequences":self.overrep.as_ref().map(Overrep::report).unwrap_or_default()})
    }
}

#[derive(Clone)]
struct Overrep {
    counts: BTreeMap<Vec<u8>, usize>,
    steps: [usize; 5],
    sampling: usize,
}

impl Overrep {
    fn candidates(reads: &[&Read], len: usize) -> BTreeMap<Vec<u8>, usize> {
        let mut counts = BTreeMap::new();
        let mut bases = 0;
        let mut sampled = 0;
        for r in reads {
            if bases >= 151 * 10000 {
                break;
            }
            bases += r.seq.len();
            sampled += 1;
        }
        // The evaluation reads outlive this scan. Borrow each distinct window
        // from them and keep only one length's frequency table in memory.
        let mut steps = [10, 20, 40, 100, 150.min(len.saturating_sub(2)).max(1)];
        steps.sort_unstable();
        for group in steps.chunk_by(|a, b| a == b) {
            let step = group[0];
            let threshold: usize = match step {
                n if n >= len.saturating_sub(1) => 3,
                100.. => 5,
                40.. => 20,
                20.. => 100,
                10.. => 500,
                _ => continue,
            };
            // fastp excludes the window ending exactly at the read end.
            let windows = || {
                reads[..sampled]
                    .iter()
                    .flat_map(|read| read.seq[..read.seq.len().saturating_sub(1)].windows(step))
            };
            let window_count: usize = reads[..sampled]
                .iter()
                .map(|read| read.seq.len().saturating_sub(step))
                .sum();
            let mut frequencies = FxHashMap::default();
            if step >= 8 && window_count >= 4096 {
                // Every occurrence of a sequence increments the same bucket.
                // Collisions can only admit extra windows to the exact count.
                let buckets = window_count.min(1 << 20).next_power_of_two() * 2;
                let shift = 64 - buckets.trailing_zeros();
                let bucket = |w: &[u8]| {
                    let prefix = u64::from_le_bytes(w[..8].try_into().unwrap());
                    let suffix = u64::from_le_bytes(w[w.len() - 8..].try_into().unwrap());
                    ((prefix ^ suffix.rotate_left(29)).wrapping_mul(0x9e3779b97f4a7c15) >> shift)
                        as usize
                };
                let mut estimates = vec![0u16; buckets];
                let cutoff = threshold.div_ceil(group.len()) as u16;
                for window in windows() {
                    let count = &mut estimates[bucket(window)];
                    *count = count.saturating_add(1);
                }
                for window in windows().filter(|w| estimates[bucket(w)] >= cutoff) {
                    *frequencies.entry(window).or_insert(0usize) += group.len();
                }
            } else {
                for window in windows() {
                    *frequencies.entry(window).or_insert(0usize) += group.len();
                }
            }
            for (seq, n) in frequencies {
                if n >= threshold {
                    counts.insert(seq.to_vec(), n);
                }
            }
        }
        // Match upstream's lexicographic removal of redundant substrings.
        let keys: Vec<_> = counts.keys().cloned().collect();
        for seq in keys {
            let Some(&count) = counts.get(&seq) else {
                continue;
            };
            if counts.iter().any(|(long, &n)| {
                *long != seq
                    && long.len() >= seq.len()
                    && long.windows(seq.len()).any(|w| w == seq)
                    && count / n < 10
            }) {
                counts.remove(&seq);
            }
        }
        counts
    }
    fn new(candidates: &BTreeMap<Vec<u8>, usize>, len: usize, sampling: usize) -> Self {
        Self {
            counts: candidates.keys().cloned().map(|s| (s, 0)).collect(),
            steps: [10, 20, 40, 100, 150.min(len.saturating_sub(2)).max(1)].map(|step| {
                if candidates.keys().any(|seq| seq.len() == step) {
                    step
                } else {
                    0
                }
            }),
            sampling,
        }
    }
    fn add(&mut self, seq: &[u8]) {
        for step in self.steps.into_iter().filter(|&step| step != 0) {
            let mut i = 0;
            while i < seq.len().saturating_sub(step) {
                if let Some(n) = self.counts.get_mut(&seq[i..i + step]) {
                    *n += 1;
                    i += step;
                }
                i += 1;
            }
        }
    }
    fn report(&self) -> BTreeMap<String, usize> {
        self.counts
            .iter()
            .filter(|(seq, n)| {
                self.sampling * *n
                    > match seq.len() {
                        10 => 500,
                        20 => 200,
                        40 => 100,
                        100 => 50,
                        _ => 20,
                    }
            })
            .map(|(seq, &n)| (String::from_utf8_lossy(seq).into_owned(), n))
            .collect()
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Failure {
    Pass,
    N,
    Short,
    Long,
    Quality,
    Complexity,
    Dimer,
}

impl Failure {
    fn tag(self) -> &'static str {
        match self {
            Self::Pass => "paired_read_is_failing",
            Self::Quality => "failed_quality_filter",
            Self::N => "failed_too_many_n_bases",
            Self::Short => "failed_too_short",
            Self::Long => "failed_too_long",
            Self::Complexity => "failed_low_complexity",
            Self::Dimer => "failed_adapter_dimer",
        }
    }
}

fn filter(r: &Read, a: &Args) -> Failure {
    if r.seq.is_empty() {
        return Failure::Short;
    }
    let m = simd::compute_stats(&r.seq, &r.qual, a.qualified_quality_phred);
    if !a.disable_quality_filtering {
        if m.unqualified as f64 > a.unqualified_percent_limit as f64 * r.seq.len() as f64 / 100.0
            || (a.average_qual > 0 && m.qsum as usize / r.seq.len() < a.average_qual as usize)
        {
            return Failure::Quality;
        }
        if m.ncnt > a.n_base_limit {
            return Failure::N;
        }
    }
    if !a.disable_length_filtering {
        if r.seq.len() < a.length_required {
            return Failure::Short;
        }
        if a.length_limit > 0 && r.seq.len() > a.length_limit {
            return Failure::Long;
        }
    }
    if a.low_complexity_filter
        && (r.seq.len() <= 1
            || r.seq.windows(2).filter(|w| w[0] != w[1]).count() as f64
                / ((r.seq.len() - 1) as f64)
                < a.complexity_threshold as f64 / 100.0)
    {
        return Failure::Complexity;
    }
    Failure::Pass
}

/// Port of fastp's trimAndCut: front/tail windows remove the whole bad window,
/// right cutting keeps high quality bases at the start of the first bad window.
fn cut(r: &mut Read, a: &Args, r2: bool) -> usize {
    let l = r.seq.len();
    let mut front = if r2 {
        a.trim_front2.unwrap_or(a.trim_front1.max(a.trim_front))
    } else {
        a.trim_front1.max(a.trim_front)
    };
    let tail = if r2 {
        a.trim_tail2.unwrap_or(a.trim_tail1.max(a.trim_tail))
    } else {
        a.trim_tail1.max(a.trim_tail)
    };
    if front.saturating_add(tail) > l {
        r.resize(0);
        return 0;
    }
    let mut end = l - tail;
    if a.cut_front {
        let w = a.cut_front_window_size.unwrap_or(a.cut_window_size);
        let q = a.cut_front_mean_quality.unwrap_or(a.cut_mean_quality) as usize + 33;
        if end.saturating_sub(front) <= w {
            r.resize(0);
            return 0;
        }
        let mut s = front;
        while s + w < end {
            if r.qual[s..s + w].iter().map(|&q| q as usize).sum::<usize>() >= w * q {
                break;
            }
            s += 1;
        }
        if s > 0 {
            s += w - 1;
        }
        while s < l && r.seq[s] == b'N' {
            s += 1;
        }
        front = s;
    }
    if a.cut_right {
        let w = a.cut_right_window_size.unwrap_or(a.cut_window_size);
        let q = a.cut_right_mean_quality.unwrap_or(a.cut_mean_quality) as usize + 33;
        if end.saturating_sub(front) <= w {
            r.resize(0);
            return 0;
        }
        for s in front..end - w {
            if r.qual[s..s + w].iter().map(|&q| q as usize).sum::<usize>() < w * q {
                end = s;
                while end < l - 1 && r.qual[end] as usize >= q {
                    end += 1;
                }
                break;
            }
        }
    } else if a.cut_tail && !a.disable_trim_tail {
        let w = a.cut_tail_window_size.unwrap_or(a.cut_window_size);
        let q = a.cut_tail_mean_quality.unwrap_or(a.cut_mean_quality) as usize + 33;
        if end.saturating_sub(front) <= w {
            r.resize(0);
            return 0;
        }
        let mut t = end - 1;
        while t >= front + w {
            if r.qual[t + 1 - w..=t]
                .iter()
                .map(|&q| q as usize)
                .sum::<usize>()
                >= w * q
            {
                break;
            }
            t -= 1;
        }
        if t < l - 1 {
            t = t + 1 - w;
        }
        while r.seq[t] == b'N' {
            if t == 0 {
                r.resize(0);
                return 0;
            }
            t -= 1;
        }
        end = t + 1;
    }
    if front >= end || ((a.cut_front || a.cut_tail || a.cut_right) && front >= l.saturating_sub(1))
    {
        r.resize(0);
        return 0;
    }
    r.resize(end);
    r.trim_front(front);
    front
}

fn poly_g(r: &mut Read, min: usize) {
    let mut mismatch = 0;
    let mut first = r.seq.len().saturating_sub(1);
    let mut compared = 0;
    for i in 0..r.seq.len() {
        if r.seq[r.seq.len() - i - 1] != b'G' {
            mismatch += 1;
        } else {
            first = r.seq.len() - i - 1;
        }
        if mismatch > 5 || (mismatch > (i + 1) / 8 && i + 1 >= min) {
            break;
        }
        compared += 1;
    }
    if compared >= min {
        r.resize(first);
    }
}

fn poly_x(r: &mut Read, min: usize) -> Option<(usize, usize)> {
    if r.seq.is_empty() {
        return None;
    }
    let bases = b"ATCG";
    let mut counts = [0usize; 4];
    let mut pos = 0;
    while pos < r.seq.len() {
        let b = r.seq[r.seq.len() - pos - 1];
        for i in 0..4 {
            if b == bases[i] || b == b'N' {
                counts[i] += 1;
            }
        }
        let cmp = pos + 1;
        if counts.iter().all(|&n| cmp - n > 5.min(cmp / 8)) && (pos >= 8 || pos + 2 >= min) {
            break;
        }
        pos += 1;
    }
    if pos + 1 < min {
        return None;
    }
    let base = (0..4).max_by_key(|&i| (counts[i], 4 - i)).unwrap();
    // Upstream scans through the first mismatch before choosing the trim point.
    let mut n = (pos + 1).min(r.seq.len());
    while n > 0 && r.seq[r.seq.len() - n] != bases[base] {
        n -= 1;
    }
    if n == 0 {
        return None;
    }
    r.resize(r.seq.len() - n);
    Some((base, n))
}

fn index(name: &[u8], last: bool) -> Vec<u8> {
    if name.len() < 5 {
        return Vec::new();
    }
    let mut end = name.len();
    for i in (0..=name.len() - 3).rev() {
        if last && (name[i] == b':' || name[i] == b'+') {
            return name[i + 1..].to_vec();
        }
        if !last && name[i] == b'+' {
            end = i;
        }
        if !last && name[i] == b':' {
            return name[i + 1..end].to_vec();
        }
    }
    Vec::new()
}

fn add_umi(r: &mut Read, umi: &[u8], a: &Args) {
    let pos = r
        .name
        .iter()
        .position(|&c| c == b' ')
        .unwrap_or(r.name.len());
    let mut tag = a.umi_delim.as_bytes().to_vec();
    if !a.umi_prefix.is_empty() {
        tag.extend_from_slice(a.umi_prefix.as_bytes());
        tag.push(b'_');
    }
    tag.extend_from_slice(umi);
    r.name_mut().splice(pos..pos, tag);
}

fn umi(r1: &mut Read, r2: Option<&mut Read>, a: &Args) {
    if !a.umi {
        return;
    }
    let mut r2 = r2;
    let mut value = Vec::new();
    match a.umi_loc.as_str() {
        "index1" => value = index(&r1.name, false),
        "index2" => {
            if let Some(r2) = r2.as_ref() {
                value = index(&r2.name, true);
            }
        }
        "read1" | "per_read" => {
            value = r1.seq[..a.umi_len.min(r1.seq.len())].to_vec();
            r1.trim_front(
                value
                    .len()
                    .saturating_add(a.umi_skip)
                    .min(r1.seq.len().saturating_sub(1)),
            );
            if a.umi_loc == "per_read" {
                if let Some(r2) = r2.as_mut() {
                    let n = a.umi_len.min(r2.seq.len());
                    value.push(b'_');
                    value.extend_from_slice(&r2.seq[..n]);
                    r2.trim_front(
                        n.saturating_add(a.umi_skip)
                            .min(r2.seq.len().saturating_sub(1)),
                    );
                }
            }
        }
        "read2" => {
            if let Some(r2) = r2.as_mut() {
                value = r2.seq[..a.umi_len.min(r2.seq.len())].to_vec();
                r2.trim_front(
                    value
                        .len()
                        .saturating_add(a.umi_skip)
                        .min(r2.seq.len().saturating_sub(1)),
                );
            }
        }
        "per_index" => {
            value = index(&r1.name, false);
            if let Some(r2) = r2.as_ref() {
                value.push(b'_');
                value.extend(index(&r2.name, true));
            }
        }
        _ => unreachable!(),
    }
    if !value.is_empty() {
        add_umi(r1, &value, a);
        if let Some(r2) = r2 {
            add_umi(r2, &value, a);
        }
    }
}

/// Lazily allocated Bloom pages use the same hashes, bit addresses, and
/// accuracy levels as fastp, without eagerly touching gigabytes of zeroes.
struct Duplication {
    buffers: Vec<BloomPages>,
    primes: Vec<u64>,
    bits: u64,
}

enum BloomPages {
    Dense(BucketTable),
    Sparse(SparseBuckets),
}

// Keep small sets of exact offsets inline. Spill to a short vector for
// moderately populated pages, then use indexed lookups as membership scans
// grow more expensive. A full bitmap bounds the storage for a dense page.
enum BloomPage {
    Small { len: u8, bits: [u16; 24] },
    Offsets(Vec<u16>),
    Sparse(FxHashSet<u16>),
    Bitmap(Box<[u8; 4096]>),
}

impl Default for BloomPage {
    fn default() -> Self {
        Self::Small {
            len: 0,
            bits: [u16::MAX; 24],
        }
    }
}

impl BloomPage {
    fn check_and_set(&mut self, bit: u16) -> bool {
        match self {
            Self::Small { len, bits } => {
                // Unused slots contain a value outside the 15-bit page.
                // Searching a fixed-size array lets the compiler use SIMD.
                if bits.contains(&bit) {
                    return true;
                }
                if usize::from(*len) < bits.len() {
                    bits[usize::from(*len)] = bit;
                    *len += 1;
                } else {
                    let mut offsets = Vec::with_capacity(32);
                    offsets.extend_from_slice(bits);
                    offsets.push(bit);
                    *self = Self::Offsets(offsets);
                }
                false
            }
            Self::Offsets(offsets) => {
                if offsets.contains(&bit) {
                    return true;
                }
                if offsets.len() < 64 {
                    offsets.push(bit);
                } else {
                    *self = Self::Sparse(
                        offsets
                            .iter()
                            .copied()
                            .chain(std::iter::once(bit))
                            .collect(),
                    );
                }
                false
            }
            Self::Sparse(offsets) => {
                if !offsets.insert(bit) {
                    return true;
                }
                // At 896 entries the next hash table allocation would exceed
                // the size of the complete bitmap (u16 keys plus control bytes).
                if offsets.len() == 896 {
                    let mut bitmap = Box::new([0; 4096]);
                    for &offset in offsets.iter() {
                        bitmap[usize::from(offset) / 8] |= 1 << (offset % 8);
                    }
                    *self = Self::Bitmap(bitmap);
                }
                false
            }
            Self::Bitmap(bitmap) => {
                let byte = &mut bitmap[usize::from(bit) / 8];
                let mask = 1 << (bit % 8);
                let duplicate = *byte & mask != 0;
                *byte |= mask;
                duplicate
            }
        }
    }
}

// Pack fixed-width offsets in a word, using zero as the empty-slot marker.
// The top bit records membership of offset zero.
#[derive(Default)]
struct Bucket<const WIDTH: u32 = 11>(u64);

impl<const WIDTH: u32> Bucket<WIDTH> {
    const MASK: u64 = (1 << WIDTH) - 1;
    const LANES: u32 = 63 / WIDTH;
    const ONES: u64 = ((1 << (WIDTH * Self::LANES)) - 1) / Self::MASK;

    /// A full bucket returns None for an unseen nonzero offset.
    fn check_and_set(&mut self, offset: u16) -> Option<bool> {
        const ZERO: u64 = 1 << 63;
        if offset == 0 {
            let seen = self.0 & ZERO != 0;
            self.0 |= ZERO;
            return Some(seen);
        }
        let value = u64::from(offset);
        let diff = self.0 ^ (value * Self::ONES);
        // Subtract one in every lane. A zero lane borrows into its high bit.
        if diff.wrapping_sub(Self::ONES) & !diff & (Self::ONES << (WIDTH - 1)) != 0 {
            return Some(true);
        }
        if self.0 & (Self::MASK << (WIDTH * (Self::LANES - 1))) != 0 {
            return None;
        }
        self.0 = (self.0 & ZERO) | (self.0 << WIDTH) | value;
        Some(false)
    }
}

struct BucketTable {
    buckets: Vec<Bucket>,
    overflow: Vec<Option<Box<BloomPage>>>,
}

impl BucketTable {
    // Default-level hashes use 32 address bits. Reversing the bytes spreads
    // clustered hashes across buckets while preserving every bit identity.
    #[inline]
    fn address(pos: u64) -> u64 {
        u64::from((pos as u32).swap_bytes())
    }

    fn new(pages: usize) -> Self {
        Self {
            buckets: std::iter::repeat_with(Bucket::default)
                .take(pages * 16)
                .collect(),
            overflow: std::iter::repeat_with(|| None).take(pages).collect(),
        }
    }

    #[inline]
    fn prefetch(&self, bit: u64) {
        let ptr = self.buckets.as_ptr().wrapping_add((bit >> 11) as usize);
        #[cfg(target_arch = "aarch64")]
        // SAFETY: the hint does not dereference or modify memory.
        unsafe {
            std::arch::asm!("prfm pldl1keep, [{ptr}]", ptr = in(reg) ptr,
                options(nostack, preserves_flags, readonly));
        }
        #[cfg(target_arch = "x86_64")]
        // SAFETY: the pointer addresses an allocated bucket and this is a cache hint.
        unsafe {
            std::arch::x86_64::_mm_prefetch::<{ std::arch::x86_64::_MM_HINT_T0 }>(ptr.cast());
        }
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        let _ = ptr;
    }

    fn check_and_set(&mut self, bit: u64) -> bool {
        if let Some(seen) = self.buckets[(bit >> 11) as usize].check_and_set((bit & 2047) as u16) {
            return seen;
        }
        self.overflow[(bit >> 15) as usize]
            .get_or_insert_with(Default::default)
            .check_and_set((bit & 32767) as u16)
    }
}

// A sparse directory across each larger address space. Allocate the 256-byte
// groups on first access. Four 13-bit offsets fit in a bucket, so most probes
// finish there. Crowded buckets share an overflow page per 32768 bits.
struct SparseBuckets {
    groups: Vec<Option<Box<[Bucket<13>; 32]>>>,
    overflow: FxHashMap<u64, BloomPage>,
}

impl SparseBuckets {
    fn new(bits: u64) -> Self {
        Self {
            groups: std::iter::repeat_with(|| None)
                .take((bits >> 18) as usize)
                .collect(),
            overflow: FxHashMap::default(),
        }
    }

    fn check_and_set(&mut self, bit: u64) -> bool {
        let group = self.groups[(bit >> 18) as usize]
            .get_or_insert_with(|| Box::new(std::array::from_fn(|_| Bucket::default())));
        if let Some(seen) = group[((bit >> 13) & 31) as usize].check_and_set((bit & 8191) as u16) {
            return seen;
        }
        self.overflow
            .entry(bit >> 15)
            .or_default()
            .check_and_set((bit & 32767) as u16)
    }
}

impl BloomPages {
    fn check_and_set(&mut self, bit: u64) -> bool {
        match self {
            Self::Dense(table) => table.check_and_set(BucketTable::address(bit)),
            Self::Sparse(table) => table.check_and_set(bit),
        }
    }
}

impl Duplication {
    fn new(level: u8) -> Self {
        let (scale, n) = match level {
            2 => (2, 2),
            3 => (2, 4),
            4 => (4, 4),
            5 => (8, 4),
            6 => (8, 8),
            _ => (1, 2),
        };
        let mut primes = Vec::new();
        let mut number = 10000u64;
        while primes.len() < n * 512 {
            number += 1;
            if number % 2 == 0 {
                continue;
            }
            let mut d = 3;
            while d * d <= number && number % d != 0 {
                d += 2;
            }
            if d * d > number {
                primes.push(number);
                number += 10000;
            }
        }
        let bits = (1u64 << 32) * scale;
        let pages = (bits / 8 / 4096) as usize;
        // Direct lookup for the default two directories. Higher accuracies
        // allocate groups of buckets on demand.
        let dense = pages * n <= 2 * (1 << 17);
        Self {
            buffers: (0..n)
                .map(|_| {
                    if dense {
                        BloomPages::Dense(BucketTable::new(pages))
                    } else {
                        BloomPages::Sparse(SparseBuckets::new(bits))
                    }
                })
                .collect(),
            primes,
            bits,
        }
    }
    fn hasher(&self) -> DuplicationHash {
        DuplicationHash {
            buffers: self.buffers.len(),
            primes: self.primes.clone(),
            #[cfg(target_arch = "aarch64")]
            pair: (self.buffers.len() == 2).then(|| duplication_neon::PairHash::new(&self.primes)),
        }
    }

    fn count_packed(&mut self, positions: &[u64]) -> usize {
        // The two default hash directories fit in a flat table. Bring later
        // buckets into cache while preserving the original check-and-set order.
        if let [BloomPages::Dense(a), BloomPages::Dense(b)] = self.buffers.as_mut_slice() {
            let (pairs, tail) = positions.as_chunks::<2>();
            debug_assert!(tail.is_empty());
            let check = |a: &mut BucketTable, b: &mut BucketTable, pair: &[u64; 2]| {
                usize::from(
                    a.check_and_set(BucketTable::address(pair[0]))
                        & b.check_and_set(BucketTable::address(pair[1])),
                )
            };
            const AHEAD: usize = 16;
            let mut count = 0;
            for (pair, later) in pairs.iter().zip(pairs.iter().skip(AHEAD)) {
                a.prefetch(BucketTable::address(later[0]));
                b.prefetch(BucketTable::address(later[1]));
                count += check(a, b, pair);
            }
            for pair in &pairs[pairs.len().saturating_sub(AHEAD)..] {
                count += check(a, b, pair);
            }
            return count;
        }
        positions
            .chunks_exact(self.buffers.len())
            .map(|pair| usize::from(self.check_positions(pair)))
            .sum()
    }

    fn check_positions(&mut self, positions: &[u64]) -> bool {
        let mut duplicate = true;
        for (buf, pos) in self.buffers.iter_mut().zip(positions) {
            let bit = pos & (self.bits - 1);
            duplicate &= buf.check_and_set(bit);
        }
        duplicate
    }
}

struct DuplicationHash {
    buffers: usize,
    primes: Vec<u64>,
    #[cfg(target_arch = "aarch64")]
    pair: Option<Box<duplication_neon::PairHash>>,
}

impl DuplicationHash {
    fn positions(&self, r1: &[u8], r2: Option<&[u8]>) -> [u64; 8] {
        #[cfg(target_arch = "aarch64")]
        if let Some(table) = &self.pair {
            if r1.len() + r2.map_or(0, <[u8]>::len) <= 512 {
                return table.positions(r1, r2, &self.primes);
            }
        }
        match self.buffers {
            2 => self.positions_with::<2>(r1, r2),
            4 => self.positions_with::<4>(r1, r2),
            8 => self.positions_with::<8>(r1, r2),
            _ => unreachable!("validated duplication accuracy"),
        }
    }

    fn positions_with<const N: usize>(&self, r1: &[u8], r2: Option<&[u8]>) -> [u64; 8] {
        let mut positions = [0u64; 8];
        let mut offset = 0;
        for r in std::iter::once(r1).chain(r2) {
            for (p, &b) in r.iter().enumerate() {
                let base = match b {
                    b'A' => 7,
                    b'T' => 222,
                    b'C' => 74,
                    b'G' => 31,
                    _ => 13,
                };
                for (i, pos) in positions.iter_mut().enumerate().take(N) {
                    let prime = self.primes[((p + offset) * N + i) & (N * 512 - 1)];
                    *pos = pos.wrapping_add(prime.wrapping_mul(base + (p + offset) as u64));
                }
            }
            offset += r.len();
        }
        positions
    }
}

struct Config {
    packed_output: bool,
    args: Args,
    adapters: [Option<adapter::AdapterMatcher>; 2],
    fasta: Vec<adapter::AdapterMatcher>,
    blacklist: [Vec<Vec<u8>>; 2],
    poly_g: bool,
}

impl Config {
    fn dedup_before_filtering(&self) -> bool {
        self.packed_output && self.args.dedup
    }

    fn profile_after_in_worker(&self) -> bool {
        self.packed_output || (!self.args.dedup && !self.args.overrepresentation_analysis)
    }
}

fn blacklisted(r: &Read, c: &Config, last: bool) -> bool {
    if c.blacklist[usize::from(last)].is_empty() {
        return false;
    }
    let target = index(&r.name, last);
    c.blacklist[usize::from(last)].iter().any(|seq| {
        seq.iter().zip(&target).filter(|(a, b)| a != b).count() <= c.args.filter_by_index_threshold
    })
}

struct Processed<'a> {
    duplicate_positions: [u64; 8],
    duplicate: bool,
    reads: [Option<Read<'a>>; 2],
    failed_reads: [Option<Read<'a>>; 2],
    merged: Option<Read<'a>>,
    overlapped: Option<Read<'a>>,
    failures: [Failure; 2],
    adapter_reads: usize,
    adapter_bases: usize,
    poly_x: Vec<(usize, usize)>,
    corrected_reads: usize,
    corrected_bases: usize,
    insert: Option<usize>,
    skipped: bool,
}

#[derive(Default)]
struct ReadPool {
    reads: Vec<(Read<'static>, Option<Read<'static>>)>,
    bytes: usize,
    limit: usize,
}

impl ReadPool {
    fn absorb(&mut self, mut other: Self) {
        if other.bytes <= self.limit.saturating_sub(self.bytes) {
            self.bytes += other.bytes;
            if self.reads.is_empty() {
                self.reads = other.reads;
            } else {
                self.reads.append(&mut other.reads);
            }
            return;
        }
        for pair in other.reads {
            self.put(pair);
        }
    }

    fn put(&mut self, pair: (Read<'static>, Option<Read<'static>>)) {
        let bytes = pair.0.capacity() + pair.1.as_ref().map_or(0, Read::capacity);
        if bytes <= self.limit.saturating_sub(self.bytes) {
            self.bytes += bytes;
            self.reads.push(pair);
        }
    }
}

// Return buffers even when output handling exits early for a filtered, merged,
// or duplicate record. The pool caps retained capacity, including long reads.
struct RecycleRead<'a> {
    processed: Processed<'static>,
    pool: &'a mut ReadPool,
}

impl std::ops::Deref for RecycleRead<'_> {
    type Target = Processed<'static>;
    fn deref(&self) -> &Self::Target {
        &self.processed
    }
}

impl Drop for RecycleRead<'_> {
    fn drop(&mut self) {
        if let Some(first) = self.processed.reads[0].take() {
            self.pool.put((first, self.processed.reads[1].take()));
        }
    }
}

impl Processed<'_> {
    /// Reads contributing to after-filter statistics, in emission order.
    fn passing_reads(&self, a: &Args) -> [Option<(usize, &Read<'_>)>; 2] {
        if self.skipped || (a.dedup && self.duplicate) {
            return [None, None];
        }
        if let Some(r) = &self.merged {
            return [(self.failures[0] == Failure::Pass).then_some((0, r)), None];
        }
        if a.merge && a.include_unmerged {
            return std::array::from_fn(|i| {
                self.reads[i]
                    .as_ref()
                    .filter(|_| self.failures[i] == Failure::Pass)
                    .map(|r| (0, r))
            });
        }
        if !a.merge && self.failures.iter().all(|f| *f == Failure::Pass) {
            return std::array::from_fn(|i| self.reads[i].as_ref().map(|r| (i, r)));
        }
        [None, None]
    }
}

struct Chunk<'a> {
    before: [Sample; 2],
    after: [Sample; 2],
    after_offset: [usize; 2],
    processed: Vec<Processed<'a>>,
    packed: Option<PackedChunk>,
}

struct PackedChunk {
    output: [Vec<u8>; 3],
    positions: Vec<u64>,
    duplicates: usize,
    recycled: ReadPool,
    count: usize,
    merged_pairs: usize,
    failures: [usize; 7],
    adapter_reads: usize,
    adapter_bases: usize,
    corrected_reads: usize,
    corrected_bases: usize,
    poly_reads: [usize; 4],
    poly_bases: [usize; 4],
    histogram: Vec<usize>,
    unknown: usize,
}

/// Output buffers of written chunks, reused by later chunks for the same
/// reason as batch buffers: freed large blocks are faulted in again. At most
/// a few batches of chunks are in flight.
static SPARE_OUTPUT: std::sync::Mutex<Vec<Vec<u8>>> = std::sync::Mutex::new(Vec::new());
const SPARE_OUTPUTS: usize = 256;

fn spare_output() -> Vec<u8> {
    let spare = SPARE_OUTPUT.lock().unwrap_or_else(|e| e.into_inner()).pop();
    spare.unwrap_or_default()
}

/// Return a written output buffer for a later chunk to reuse.
fn recycle_output(mut out: Vec<u8>) {
    if out.capacity() > 0 {
        let mut spare = SPARE_OUTPUT.lock().unwrap_or_else(|e| e.into_inner());
        if spare.len() < SPARE_OUTPUTS {
            out.clear();
            spare.push(out);
        }
    }
}

impl Drop for PackedChunk {
    fn drop(&mut self) {
        for out in &mut self.output {
            recycle_output(std::mem::take(out));
        }
    }
}

impl PackedChunk {
    fn new(records: usize, recycle_limit: usize, hashes: usize) -> Self {
        Self {
            output: std::array::from_fn(|_| spare_output()),
            positions: Vec::with_capacity(records * hashes),
            duplicates: 0,
            recycled: ReadPool {
                reads: Vec::with_capacity(if recycle_limit == 0 { 0 } else { records }),
                limit: recycle_limit,
                ..ReadPool::default()
            },
            count: 0,
            merged_pairs: 0,
            failures: [0; 7],
            adapter_reads: 0,
            adapter_bases: 0,
            corrected_reads: 0,
            corrected_bases: 0,
            poly_reads: [0; 4],
            poly_bases: [0; 4],
            histogram: vec![0; 513],
            unknown: 0,
        }
    }

    fn add(&mut self, mut p: Processed, a: &Args) {
        self.count += 1;
        self.duplicates += usize::from(p.duplicate);
        if !p.skipped {
            self.adapter_reads += p.adapter_reads;
            self.adapter_bases += p.adapter_bases;
            self.corrected_reads += p.corrected_reads;
            self.corrected_bases += p.corrected_bases;
            for &(base, count) in &p.poly_x {
                self.poly_reads[base] += 1;
                self.poly_bases[base] += count;
            }
            let paired = p.reads[1].is_some();
            if paired {
                if let Some(size) = p.insert {
                    self.histogram[size.min(512)] += 1;
                } else {
                    self.unknown += 1;
                }
            }
            if let Some(merged) = &p.merged {
                self.failures[p.failures[0] as usize] += 2;
                if p.failures[0] == Failure::Pass && !p.duplicate {
                    self.merged_pairs += 1;
                    if a.stdout || a.merged_out.is_some() {
                        merged
                            .write(&mut self.output[2], None)
                            .expect("writing to Vec cannot fail");
                    }
                }
            } else if a.merge && a.include_unmerged {
                for (i, read) in p.reads.iter().enumerate() {
                    if let Some(read) = read {
                        self.failures[p.failures[i] as usize] += 1;
                        if p.failures[i] == Failure::Pass
                            && !p.duplicate
                            && (a.stdout || a.merged_out.is_some())
                        {
                            read.write(&mut self.output[2], None)
                                .expect("writing to Vec cannot fail");
                        }
                    }
                }
            } else {
                let failure = if paired {
                    p.failures[0].max(p.failures[1])
                } else {
                    p.failures[0]
                };
                self.failures[failure as usize] += if paired { 2 } else { 1 };
                if failure == Failure::Pass && !p.duplicate {
                    for (i, read) in p.reads.iter().enumerate() {
                        if let Some(read) = read {
                            if (a.stdout && !a.merge)
                                || (i == 0 && !a.output.is_empty())
                                || (i == 1 && a.output2.is_some())
                            {
                                // Both stdout mates share one buffer so paired
                                // reads remain interleaved within every chunk.
                                let output = if a.stdout { 0 } else { i };
                                read.write(&mut self.output[output], None)
                                    .expect("writing to Vec cannot fail");
                            }
                        }
                    }
                }
            }
        }
        if self.recycled.limit > 0 {
            if let Some(first) = p.reads[0].take() {
                let pair = (first, p.reads[1].take());
                self.recycled
                    .put((pair.0.into_owned(), pair.1.map(Read::into_owned)));
            }
        }
    }
}

fn pack_main_output(a: &Args) -> bool {
    !a.overrepresentation_analysis
        && a.output != "-"
        && a.output2.as_deref() != Some("-")
        // Shared destinations (including the same S3 key) need interleaved
        // record emission rather than one packed buffer per mate.
        && a.output2.as_deref() != Some(a.output.as_str())
        && a.merged_out.as_deref().is_none_or(|path| {
            path != "-" && path != a.output && Some(path) != a.output2.as_deref()
        })
        && a.split.is_none()
        && a.split_by_lines.is_none()
        && a.failed_out.is_none()
        && a.unpaired1.is_none()
        && a.unpaired2.is_none()
        && a.overlapped_out.is_none()
}

fn mark_duplicates(chunks: &mut [Chunk], duplication: Option<&mut Duplication>) {
    if let Some(duplication) = duplication {
        for chunk in chunks {
            if let Some(packed) = &mut chunk.packed {
                packed.duplicates += duplication.count_packed(&packed.positions);
            } else {
                for p in &mut chunk.processed {
                    p.duplicate = duplication.check_positions(&p.duplicate_positions);
                }
            }
        }
    }
}

fn trim_adapter(r: &mut Read, seq: &adapter::AdapterMatcher, min: usize) -> Option<usize> {
    let before = r.seq.len();
    let found = seq.find(&r.seq, min)?;
    r.resize(found.position);
    Some(if found.is_atailing {
        found.overlap_len
    } else {
        before - r.seq.len()
    })
}

fn gap_overlap(r1: &[u8], r2: &[u8], c: &overlap::OverlapConfig) -> Option<overlap::OverlapResult> {
    let rc = overlap::reverse_complement(r2);
    let forward = (0..r1.len().saturating_sub(c.min_overlap_len)).map(|i| i as isize);
    let reverse = (0..r2.len().saturating_sub(c.min_overlap_len)).map(|i| -(i as isize));
    for offset in forward.chain(reverse) {
        let a = &r1[offset.max(0) as usize..];
        let b = &rc[(-offset).max(0) as usize..];
        let len = a.len().min(b.len());
        let limit = c.max_diff.min(len * c.max_diff_percent / 100);
        if len < 2 {
            continue;
        }
        if let Some(diff) = adapter::one_gap_diff(a, b, len - 1, limit)
            .or_else(|| adapter::one_gap_diff(b, a, len - 1, limit))
        {
            return Some(overlap::OverlapResult {
                overlapped: true,
                offset,
                overlap_len: len,
                diff,
            });
        }
    }
    None
}

fn process<'a>(
    mut r1: Read<'a>,
    mut r2: Option<Read<'a>>,
    c: &Config,
    overlap_scratch: &mut Vec<u8>,
) -> Processed<'a> {
    let a = &c.args;
    let mut p = Processed {
        duplicate_positions: [0; 8],
        duplicate: false,
        reads: [None, None],
        failed_reads: [None, None],
        merged: None,
        overlapped: None,
        failures: [Failure::Pass; 2],
        adapter_reads: 0,
        adapter_bases: 0,
        poly_x: Vec::new(),
        corrected_reads: 0,
        corrected_bases: 0,
        insert: None,
        skipped: false,
    };
    if blacklisted(&r1, c, false) || r2.as_ref().is_some_and(|r| blacklisted(r, c, true)) {
        p.skipped = true;
        p.reads = [Some(r1), r2];
        return p;
    }
    if a.fix_mgi_id {
        for r in std::iter::once(&mut r1).chain(r2.iter_mut()) {
            if r.name.ends_with(b"/1") || r.name.ends_with(b"/2") {
                let pos = r.name.len() - 2;
                r.name_mut().insert(pos, b' ');
            }
        }
    }
    umi(&mut r1, r2.as_mut(), a);
    let before_cut = if a.failed_out.is_some() {
        [Some(r1.clone()), r2.clone()]
    } else {
        [None, None]
    };
    let front1 = cut(&mut r1, a, false);
    let front2 = r2.as_mut().map_or(0, |r| cut(r, a, true));
    // trimAndCut returns null without changing the original when no read
    // remains; failed output still uses that original (after UMI handling).
    for (i, r) in std::iter::once(&r1).chain(r2.iter()).enumerate() {
        if r.seq.is_empty() {
            p.failed_reads[i] = before_cut[i].clone();
        }
    }
    if c.poly_g {
        poly_g(&mut r1, a.poly_g_min_len);
        if let Some(r2) = r2.as_mut() {
            poly_g(r2, a.poly_g_min_len);
        }
    }
    let ovconfig = overlap::OverlapConfig {
        min_overlap_len: a.overlap_len_require,
        max_diff: a.overlap_diff_limit,
        max_diff_percent: a.overlap_diff_percent_limit,
    };
    let overlap_lengths = (r1.seq.len(), r2.as_ref().map_or(0, |r| r.seq.len()));
    let ov = r2
        .as_ref()
        .and_then(|r| overlap::detect_overlap_reusing(&r1.seq, &r.seq, &ovconfig, overlap_scratch));
    let gap_ov = if a.allow_gap_overlap_trimming && ov.is_none() {
        r2.as_ref()
            .and_then(|r| gap_overlap(&r1.seq, &r.seq, &ovconfig))
    } else {
        None
    };
    if let (Some(r2), Some(ov)) = (r2.as_ref(), ov.as_ref()) {
        p.insert = Some(if ov.offset > 0 {
            r1.seq.len() + r2.seq.len() - ov.overlap_len + front1 + front2
        } else {
            ov.overlap_len + front1 + front2
        });
    }
    // Exact overlaps need no correction or mutable copies.
    if a.correction && ov.as_ref().is_some_and(|o| o.diff > 0) {
        if let (Some(r2), Some(ov)) = (r2.as_mut(), ov.as_ref()) {
            r1.verbatim = None;
            r2.verbatim = None;
            let seq1 = r1.seq.clone();
            let seq2 = r2.seq.clone();
            let correction = overlap::correct_by_overlap(
                r1.seq.to_mut(),
                r1.qual.to_mut(),
                r2.seq.to_mut(),
                r2.qual.to_mut(),
                ov,
            );
            p.corrected_reads = usize::from(seq1 != r1.seq) + usize::from(seq2 != r2.seq);
            p.corrected_bases = correction.corrected;
        }
    }
    let mut trimmed = [false; 2];
    if !a.disable_adapter_trimming {
        let mut used_overlap = false;
        if let (Some(r2), Some(ov)) = (r2.as_mut(), ov.as_ref().or(gap_ov.as_ref())) {
            if ov.offset < 0 {
                let n1 = r1.seq.len().min(ov.overlap_len + front2);
                let n2 = r2.seq.len().min(ov.overlap_len + front1);
                p.adapter_bases += r1.seq.len() - n1 + r2.seq.len() - n2;
                r1.resize(n1);
                r2.resize(n2);
                trimmed = [true, true];
                used_overlap = true;
            }
        }
        for (i, r) in std::iter::once(&mut r1).chain(r2.iter_mut()).enumerate() {
            if !used_overlap {
                if let Some(seq) = &c.adapters[i] {
                    if let Some(n) = trim_adapter(r, seq, 4) {
                        p.adapter_bases += n;
                        trimmed[i] = true;
                    }
                }
            }
            let min = if c.fasta.len() > 256 {
                6
            } else if c.fasta.len() > 16 {
                5
            } else {
                4
            };
            for seq in &c.fasta {
                if let Some(n) = trim_adapter(r, seq, min) {
                    p.adapter_bases += n;
                    trimmed[i] = true;
                }
            }
        }
        p.adapter_reads = trimmed.iter().filter(|&&b| b).count();
    }
    let dimer = !a.disable_adapter_trimming
        && trimmed.iter().any(|&b| b)
        && r1.seq.len() <= a.dimer_max_len
        && r2.as_ref().is_none_or(|r| r.seq.len() <= a.dimer_max_len);
    if a.overlapped_out.is_some() {
        if let Some(r2) = &r2 {
            let config = overlap::OverlapConfig {
                max_diff_percent: 0,
                ..ovconfig.clone()
            };
            if let Some(ov) = overlap::detect_overlap(&r1.seq, &r2.seq, &config) {
                if ov.diff == 0 {
                    let mut r = r1.clone();
                    let start = ov.offset.max(0) as usize;
                    // Match v1.4.0 read construction: its second string-constructor
                    // argument is an offset into the suffix, not a substring length.
                    r.trim_front(start + ov.overlap_len);
                    p.overlapped = Some(r);
                }
            }
        }
    }
    for (i, r) in std::iter::once(&mut r1).chain(r2.iter_mut()).enumerate() {
        if a.trim_poly_x {
            if let Some(x) = poly_x(r, a.poly_x_min_len) {
                p.poly_x.push(x);
            }
        }
        let max_len = if i == 1 {
            a.max_len2.unwrap_or(a.max_len1)
        } else {
            a.max_len1
        };
        if max_len > 0 {
            r.resize(max_len);
        }
        p.failures[i] = if dimer { Failure::Dimer } else { filter(r, a) };
    }
    if a.merge {
        if let Some(r2) = &r2 {
            // After the first search, only correction or tail truncation can
            // change sequence content. Reuse even a failed search when neither
            // happened; otherwise honor fastp's post-trimming overlap semantics.
            let merge_overlap =
                if overlap_lengths == (r1.seq.len(), r2.seq.len()) && p.corrected_bases == 0 {
                    ov
                } else {
                    overlap::detect_overlap_reusing(&r1.seq, &r2.seq, &ovconfig, overlap_scratch)
                };
            if let Some(ov) = merge_overlap {
                let (name, seq, qual) =
                    overlap::merge_reads(&r1.seq, &r1.qual, &r1.name, &r2.seq, &r2.qual, &ov);
                let suffix = format!(
                    " merged_{}_{}",
                    ov.overlap_len + ov.offset.max(0) as usize,
                    if ov.offset > 0 {
                        r2.seq.len() - ov.overlap_len
                    } else {
                        0
                    }
                );
                let mut plus = r1.plus.clone();
                if plus.as_ref() != b"+" {
                    plus.to_mut().extend_from_slice(suffix.as_bytes());
                }
                let merged = Read {
                    name: name.into(),
                    seq: seq.into(),
                    plus,
                    qual: qual.into(),
                    ..Read::default()
                };
                let failure = filter(&merged, a);
                p.failures = [failure; 2];
                p.merged = Some(merged);
            }
        }
    }
    p.reads = [Some(r1), r2];
    p
}

struct ReadBatch {
    duplicate_flags: Option<Vec<bool>>,
    raw: Option<[block::Records; 2]>,
    reads: Vec<(Read<'static>, Option<Read<'static>>)>,
    start: usize,
}

impl ReadBatch {
    fn mark_duplicates(&mut self, shared: &Shared, filter: &mut Duplication) {
        let raw = self
            .raw
            .as_ref()
            .expect("packed output uses raw FASTQ batches");
        let hasher = shared
            .hasher
            .as_ref()
            .expect("deduplication requires hashes");
        let interleaved = shared.config.args.interleaved_in;
        let stride = if interleaved { 2 } else { 1 };
        // Hash original sequences in parallel, then mark in input order so
        // the first occurrence still wins across worker and batch boundaries.
        let positions: Vec<_> = shared.pool().install(|| {
            (0..raw[0].len() / stride)
                .into_par_iter()
                .map(|i| {
                    let mate = if interleaved {
                        Some(raw[0].sequence(2 * i + 1))
                    } else if raw[1].len() > 0 {
                        Some(raw[1].sequence(i))
                    } else {
                        None
                    };
                    hasher.positions(raw[0].sequence(i * stride), mate)
                })
                .collect()
        });
        self.duplicate_flags = Some(
            positions
                .iter()
                .map(|pos| filter.check_positions(pos))
                .collect(),
        );
    }
}

struct BatchReader {
    raw: bool,
    reader1: Reader,
    reader2: Option<Reader>,
    sampled: std::vec::IntoIter<(Read<'static>, Option<Read<'static>>)>,
    sampled_raw: VecDeque<ReadBatch>,
    pool: ReadPool,
    interleaved: bool,
    max_bytes: usize,
    limit: usize,
    count: usize,
    eof: bool,
}

impl BatchReader {
    fn next_batch(&mut self) -> Result<Option<ReadBatch>> {
        if self.eof || (self.limit != 0 && self.count >= self.limit) {
            return Ok(None);
        }
        let stride = if self.interleaved { 2 } else { 1 };
        if let Some(mut batch) = self.sampled_raw.pop_front() {
            let raw = batch.raw.as_mut().unwrap();
            if self.limit != 0 {
                let remaining = self.limit - self.count;
                raw[0].truncate(remaining.saturating_mul(stride));
                raw[1].truncate(remaining);
            }
            batch.start = self.count;
            self.count += raw[0].len() / stride;
            return Ok(Some(batch));
        }
        if self.raw && self.sampled.len() == 0 {
            let limit = if self.limit == 0 {
                usize::MAX
            } else {
                self.limit - self.count
            };
            let raw = if self.interleaved {
                block::read_interleaved_batch(&mut self.reader1, self.max_bytes, limit)?
            } else {
                block::read_batch(
                    &mut self.reader1,
                    self.reader2.as_mut(),
                    self.max_bytes,
                    limit,
                )?
            };
            let Some(raw) = raw else {
                self.eof = true;
                return Ok(None);
            };
            let start = self.count;
            self.count += raw[0].len() / stride;
            return Ok(Some(ReadBatch {
                duplicate_flags: None,
                raw: Some(raw),
                reads: Vec::new(),
                start,
            }));
        }
        // Fill recycled slots in place to avoid moving Vec descriptors per pair.
        let mut batch = std::mem::take(&mut self.pool.reads);
        self.pool.bytes = 0;
        let mut records = 0;
        let mut bytes = 0;
        while bytes < self.max_bytes
            && !self.eof
            && (self.limit == 0 || self.count + records < self.limit)
        {
            if let Some(pair) = self.sampled.next() {
                if records == batch.len() {
                    batch.push(pair);
                } else {
                    batch[records] = pair;
                }
            } else {
                if self.raw {
                    // Drain the retained discovery sample, then switch to
                    // borrowed blocks instead of filling one extra owned batch.
                    break;
                }
                if records == batch.len() {
                    batch.push((Read::default(), None));
                }
                let pair = &mut batch[records];
                if !self.reader1.read_into(&mut pair.0)? {
                    if let Some(reader2) = self.reader2.as_mut() {
                        ensure!(
                            reader2.next()?.is_none(),
                            "Paired input files have different read counts"
                        );
                    }
                    self.eof = true;
                    break;
                }
                if self.interleaved {
                    ensure!(
                        self.reader1
                            .read_into(pair.1.get_or_insert_with(Read::default))?,
                        "Odd read count in interleaved FASTQ"
                    );
                } else if let Some(reader2) = self.reader2.as_mut() {
                    ensure!(
                        reader2.read_into(pair.1.get_or_insert_with(Read::default))?,
                        "Paired input files have different read counts"
                    );
                }
            }
            let (r1, r2) = &batch[records];
            bytes += r1.seq.len() * 2
                + r1.name.len()
                + r2.as_ref().map_or(0, |r| r.seq.len() * 2 + r.name.len());
            records += 1;
        }
        batch.truncate(records);
        if batch.is_empty() {
            return Ok(None);
        }
        let start = self.count;
        self.count += batch.len();
        Ok(Some(ReadBatch {
            duplicate_flags: None,
            raw: None,
            reads: batch,
            start,
        }))
    }
}

fn process_chunk<'a>(
    pairs: impl Iterator<Item = Result<(Read<'a>, Option<Read<'a>>)>>,
    records: usize,
    start: usize,
    config: &Config,
    before: &[Sample; 2],
    after: &[Sample; 2],
    duplication: Option<&DuplicationHash>,
    duplicate_flags: Option<&[bool]>,
    recycle_limit: usize,
) -> Result<Chunk<'a>> {
    let a = &config.args;
    let profile_after_in_worker = config.profile_after_in_worker();
    let packed_output = config.packed_output;
    let shared_stats = packed_output && !a.merge;
    let mut shared: [Sample; 2] = std::array::from_fn(|_| Sample::default());
    let mut chunk = Chunk {
        before: std::array::from_fn(|i| before[i].empty_like()),
        after: std::array::from_fn(|i| after[i].empty_like()),
        after_offset: [0; 2],
        processed: Vec::with_capacity(if packed_output { 0 } else { records }),
        packed: packed_output.then(|| {
            PackedChunk::new(records, recycle_limit, duplication.map_or(0, |d| d.buffers))
        }),
    };
    let mut overlap_scratch = Vec::new();
    for (index, pair) in pairs.enumerate() {
        let (r1, r2) = pair?;
        let ordinal = start + index;
        let mut original = [None; 2];
        for (i, r) in [Some(&r1), r2.as_ref()].into_iter().enumerate() {
            if let Some(r) = r {
                if shared_stats {
                    original[i] = r.borrowed_bases();
                }
                if original[i].is_none() {
                    chunk.before[i].add_at(r, a.no_kmer, ordinal);
                }
            }
        }
        let positions = duplication
            .as_ref()
            .map(|d| d.positions(&r1.seq, r2.as_ref().map(|r| r.seq.as_ref())))
            .unwrap_or([0; 8]);
        let mut p = process(r1, r2, &config, &mut overlap_scratch);
        p.duplicate = duplicate_flags.is_some_and(|flags| flags[index]);
        if shared_stats {
            let passing =
                !p.duplicate && !p.skipped && p.failures.iter().all(|f| *f == Failure::Pass);
            for (i, r) in p
                .reads
                .iter()
                .enumerate()
                .filter_map(|(i, r)| r.as_ref().map(|r| (i, r)))
            {
                if let Some((seq, qual)) = original[i] {
                    // A borrowed slice with the same start is an unchanged
                    // prefix, even when trimming has shortened its tail.
                    if passing && r.seq.as_ptr() == seq.as_ptr() && r.qual.as_ptr() == qual.as_ptr()
                    {
                        shared[i].add_at(r, a.no_kmer, ordinal);
                        if r.seq.len() < seq.len() {
                            chunk.before[i].add_removed_tail(seq, qual, r.seq.len(), a.no_kmer);
                        }
                        continue;
                    }
                    chunk.before[i].add_parts(seq, qual, a.no_kmer, ordinal);
                }
                if passing {
                    chunk.after[i].add_at(r, a.no_kmer, 0);
                }
            }
        } else if profile_after_in_worker {
            for (i, r) in p.passing_reads(a).into_iter().flatten() {
                chunk.after[i].add_at(r, a.no_kmer, 0);
            }
        }
        if let Some(packed) = &mut chunk.packed {
            if let Some(d) = duplication {
                packed.positions.extend_from_slice(&positions[..d.buffers]);
            }
            packed.add(p, a);
        } else {
            p.duplicate_positions = positions;
            chunk.processed.push(p);
        }
    }
    if shared_stats {
        for (i, sample) in shared.into_iter().enumerate() {
            chunk.before[i].merge(sample.clone());
            chunk.after[i].merge(sample);
        }
    }
    Ok(chunk)
}

fn process_batch(
    batch: ReadBatch,
    pool: &rayon::ThreadPool,
    config: &Config,
    before: &[Sample; 2],
    after: &[Sample; 2],
    duplication: Option<&DuplicationHash>,
    recycle_limit: usize,
) -> Result<Vec<Chunk<'static>>> {
    let ReadBatch {
        raw,
        reads,
        start,
        duplicate_flags: flags,
    } = batch;
    let threads = pool.current_num_threads();
    let interleaved = config.args.interleaved_in;
    let stride = if interleaved { 2 } else { 1 };
    let count = raw.as_ref().map_or(reads.len(), |r| r[0].len() / stride);
    // One chunk per worker leaves the batch waiting on its slowest chunk,
    // and with few workers that is often one on an efficiency core. Split
    // into chunks of about 4000 records, up to four per worker, so faster
    // workers take more of the batch; more and smaller chunks only add
    // per-chunk overhead.
    let chunks = (count / 4000).clamp(threads, 4 * threads);
    let chunk_size = count.div_ceil(chunks).max(1);
    let duplication = if flags.is_some() { None } else { duplication };
    let results: Vec<Result<Chunk<'static>>> = pool.install(|| {
        if let Some(raw) = raw {
            (0..count.div_ceil(chunk_size))
                .into_par_iter()
                .map(|chunk| {
                    let offset = chunk * chunk_size;
                    let end = (offset + chunk_size).min(count);
                    let pairs = (offset..end).map(|i| {
                        let ordinal = (start + i) * stride;
                        let first = raw[0].read(i * stride, ordinal, config.args.phred64)?;
                        let second = if interleaved {
                            Some(raw[0].read(i * stride + 1, ordinal + 1, config.args.phred64)?)
                        } else if raw[1].len() > 0 {
                            Some(raw[1].read(i, start + i, config.args.phred64)?)
                        } else {
                            None
                        };
                        Ok((first, second))
                    });
                    let chunk = process_chunk(
                        pairs,
                        end - offset,
                        start + offset,
                        config,
                        before,
                        after,
                        duplication,
                        flags.as_ref().map(|flags| &flags[offset..end]),
                        0,
                    )?;
                    // Borrowed input is used only for packed output. All output
                    // bytes and statistics own their data before this block drops.
                    assert!(chunk.processed.is_empty());
                    Ok(Chunk {
                        before: chunk.before,
                        after: chunk.after,
                        after_offset: chunk.after_offset,
                        processed: Vec::new(),
                        packed: chunk.packed,
                    })
                })
                .collect()
        } else {
            reads
                .into_par_iter()
                .chunks(chunk_size)
                .enumerate()
                .map(|(chunk, pairs)| {
                    let records = pairs.len();
                    process_chunk(
                        pairs.into_iter().map(Ok),
                        records,
                        start + chunk * chunk_size,
                        config,
                        before,
                        after,
                        duplication,
                        None,
                        recycle_limit / threads,
                    )
                })
                .collect()
        }
    });
    // Collect in input order so malformed input always reports the earliest
    // failing worker chunk, independent of scheduling.
    results.into_iter().collect()
}

struct Output {
    #[cfg(feature = "native")]
    paired_gzip: Option<paired_gzip::PairedGzip>,
    /// Uncompressed local files that are not split, written on their own threads.
    #[cfg(feature = "native")]
    direct: BTreeMap<String, direct_writer::DirectWriter>,
    writers: BTreeMap<String, split::SplitWriter>,
    scratch: Vec<u8>,
    compression_threads: usize,
}

impl Output {
    fn new(a: &Args) -> Self {
        let gzip_outputs = [
            Some(a.output.as_str()),
            a.output2.as_deref(),
            a.merged_out.as_deref(),
            a.unpaired1.as_deref(),
            a.unpaired2.as_deref(),
            a.failed_out.as_deref(),
            a.overlapped_out.as_deref(),
        ]
        .into_iter()
        .flatten()
        .filter(|path| {
            !crate::io::is_s3_path(path)
                && crate::io::CompressionFormat::from_path(path)
                    == crate::io::CompressionFormat::Gzip
        })
        .collect::<std::collections::BTreeSet<_>>()
        .len();
        Self {
            #[cfg(feature = "native")]
            paired_gzip: None,
            #[cfg(feature = "native")]
            direct: BTreeMap::new(),
            writers: BTreeMap::new(),
            scratch: Vec::new(),
            compression_threads: crate::io::threads_per_output(gzip_outputs),
        }
    }
    fn open_main(&mut self, a: &Args, split: &split::SplitConfig) -> Result<()> {
        #[cfg(feature = "native")]
        if let Some(second) = a.output2.as_deref() {
            let paths = [a.output.as_str(), second];
            if a.parallel_compression
                && !a.merge
                && pack_main_output(a)
                && paths.iter().all(|p| {
                    !io::is_s3_path(p)
                        && io::CompressionFormat::from_path(p) == io::CompressionFormat::Gzip
                })
            {
                self.paired_gzip = Some(paired_gzip::PairedGzip::new(
                    paths,
                    a.compression_level.unwrap_or(4),
                    2 * self.compression_threads,
                )?);
                return Ok(());
            }
        }
        for path in [Some(a.output.as_str()), a.output2.as_deref()]
            .into_iter()
            .flatten()
        {
            self.open(path, a, split.clone())?;
        }
        Ok(())
    }
    fn open(&mut self, path: &str, a: &Args, config: split::SplitConfig) -> Result<()> {
        #[cfg(feature = "native")]
        if self.direct.contains_key(path) {
            return Ok(());
        }
        #[cfg(feature = "native")]
        if !path.is_empty()
            && path != "-"
            && !self.writers.contains_key(path)
            && matches!(config.mode, split::SplitMode::None)
            && !io::is_s3_path(path)
            && io::CompressionFormat::from_path(path) == io::CompressionFormat::None
        {
            let writer = direct_writer::DirectWriter::create(path, recycle_output)?;
            self.direct.insert(path.into(), writer);
            return Ok(());
        }
        if !path.is_empty() && !self.writers.contains_key(path) {
            self.writers.insert(
                path.into(),
                split::SplitWriter::new_with_compression_threads(
                    path,
                    config,
                    a.compression_level.unwrap_or(4),
                    a.parallel_compression,
                    self.compression_threads,
                )?,
            );
        }
        Ok(())
    }
    fn write(&mut self, path: Option<&str>, read: &Read, tag: Option<&str>) -> Result<()> {
        #[cfg(feature = "native")]
        if let Some(writer) = path.and_then(|p| self.direct.get_mut(p)) {
            self.scratch.clear();
            read.write(&mut self.scratch, tag)?;
            return writer.write(&self.scratch);
        }
        if let Some(writer) = path.and_then(|p| self.writers.get_mut(p)) {
            self.scratch.clear();
            read.write(&mut self.scratch, tag)?;
            writer.write_all(&self.scratch)?;
        }
        Ok(())
    }
    fn write_packed_mates(
        &mut self,
        paths: [Option<&str>; 2],
        data: [&mut Vec<u8>; 2],
    ) -> Result<()> {
        #[cfg(feature = "native")]
        if let Some(pair) = &mut self.paired_gzip {
            pair.write(0, data[0])?;
            return pair.write(1, data[1]);
        }
        let [first, second] = data;
        self.write_packed(paths[0], first)?;
        self.write_packed(paths[1], second)
    }
    /// A direct writer takes the buffer itself; the chunk recycles an empty one.
    fn write_packed(&mut self, path: Option<&str>, bytes: &mut Vec<u8>) -> Result<()> {
        #[cfg(feature = "native")]
        if let Some(writer) = path.and_then(|p| self.direct.get_mut(p)) {
            return writer.write_owned(std::mem::take(bytes));
        }
        if let Some(writer) = path.and_then(|p| self.writers.get_mut(p)) {
            writer.write_all(bytes)?;
        }
        Ok(())
    }
    fn finish(self) -> Result<()> {
        #[cfg(feature = "native")]
        if let Some(pair) = self.paired_gzip {
            pair.finish()?;
        }
        #[cfg(feature = "native")]
        for (_, w) in self.direct {
            w.finish()?;
        }
        for (_, w) in self.writers {
            w.finish()?;
        }
        Ok(())
    }
}

/// fastp defaults to 3 workers. Output does not depend on the worker count, so
/// fasterp uses every CPU by default, up to the 16 workers fastp allows.
/// wasm32 cannot spawn threads, so its one worker is the calling thread. A
/// thread can join only one pool, so build it once and reuse it across runs.
#[cfg(target_arch = "wasm32")]
fn wasm_pool() -> Result<&'static rayon::ThreadPool> {
    static POOL: std::sync::OnceLock<rayon::ThreadPool> = std::sync::OnceLock::new();
    if let Some(pool) = POOL.get() {
        return Ok(pool);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .use_current_thread()
        .build()?;
    Ok(POOL.get_or_init(|| pool))
}

/// Stack for each worker. A worker that waits for part of a batch runs other
/// chunks meanwhile, on the same stack, and with two batches in flight it
/// can start the next batch there too. Each level holds a chunk's statistics,
/// about 240 KiB; eleven levels deep was seen, past the 2 MiB a thread gets
/// by default. This covers one level for every chunk of both batches. Only
/// the pages a worker touches are ever backed by memory.
#[cfg(not(target_arch = "wasm32"))]
const WORKER_STACK: usize = 64 << 20;

#[cfg(not(target_arch = "wasm32"))]
fn worker_pool(threads: usize) -> Result<rayon::ThreadPool> {
    Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(WORKER_STACK)
        .build()?)
}

fn worker_threads(a: &Args) -> usize {
    // wasm32-unknown-unknown cannot spawn threads.
    if cfg!(target_arch = "wasm32") {
        return 1;
    }
    a.threads.unwrap_or_else(|| num_cpus::get().clamp(1, 16))
}

fn validate(a: &Args) -> Result<()> {
    let pe = a.input2.is_some() || a.interleaved_in;
    ensure!(
        !a.input.is_empty() || a.stdin,
        "Input is required: specify --in1/-i or --stdin"
    );
    ensure!(
        !(a.interleaved_in && a.input2.is_some()),
        "--interleaved_in cannot be combined with --in2"
    );
    ensure!(
        !(a.input == "-" && a.input2.as_deref() == Some("-")),
        "Use --interleaved_in for paired input from stdin"
    );
    ensure!(
        pe || a.output2.is_none(),
        "--out2 requires paired-end input"
    );
    ensure!(
        a.stdout || !pe || (a.output.is_empty() == a.output2.is_none()),
        "Paired output requires both --out1 and --out2"
    );
    ensure!(
        !(a.merge && (a.split.is_some() || a.split_by_lines.is_some())),
        "Cannot combine merging and splitting"
    );
    ensure!(!a.merge || pe, "Merging requires paired-end input");
    ensure!(
        !a.merge || a.merged_out.is_some() || a.stdout,
        "Merging requires --merged_out or --stdout"
    );
    ensure!(
        !a.include_unmerged || a.merge,
        "--include_unmerged requires --merge"
    );
    ensure!(
        !(a.trim_poly_g && a.disable_trim_poly_g),
        "Cannot enable and disable polyG trimming together"
    );
    ensure!(
        (1..=16).contains(&worker_threads(a)),
        "Worker thread count must be between 1 and 16"
    );
    ensure!(
        a.batch_bytes > 0 && a.max_backlog != Some(0),
        "Batch size and backlog must be positive"
    );
    ensure!(
        a.unqualified_percent_limit <= 100 && a.complexity_threshold <= 100,
        "Percentage thresholds must be between 0 and 100"
    );
    ensure!(
        a.overlap_diff_percent_limit <= 100,
        "Overlap difference percentage must be between 0 and 100"
    );
    ensure!(
        a.qualified_quality_phred <= 93 && a.average_qual <= 93,
        "Quality thresholds must be between 0 and 93"
    );
    ensure!(
        a.dup_calc_accuracy <= 6,
        "Duplication accuracy must be between 1 and 6"
    );
    for w in [
        Some(a.cut_window_size),
        a.cut_front_window_size,
        a.cut_tail_window_size,
        a.cut_right_window_size,
    ]
    .into_iter()
    .flatten()
    {
        ensure!(
            (1..=1000).contains(&w),
            "Quality cutting window must be between 1 and 1000"
        );
    }
    for q in [
        Some(a.cut_mean_quality),
        a.cut_front_mean_quality,
        a.cut_tail_mean_quality,
        a.cut_right_mean_quality,
    ]
    .into_iter()
    .flatten()
    {
        ensure!(
            (1..=36).contains(&q),
            "Quality cutting threshold must be between 1 and 36"
        );
    }
    ensure!(
        a.poly_g_min_len > 0 && a.poly_x_min_len > 0,
        "PolyG/polyX minimum lengths must be positive"
    );
    ensure!(a.overlap_len_require > 0, "Overlap length must be positive");
    ensure!(
        (1..=10000).contains(&a.overrepresentation_sampling),
        "Overrepresentation sampling must be between 1 and 10000"
    );
    ensure!(
        ["compact", "pretty", "off", "jsonl"].contains(&a.stats_format.as_str()),
        "Invalid stats format"
    );
    if a.umi {
        ensure!(
            [
                "read1",
                "read2",
                "index1",
                "index2",
                "per_read",
                "per_index"
            ]
            .contains(&a.umi_loc.as_str()),
            "--umi requires --umi_loc (read1/read2/index1/index2/per_read/per_index)"
        );
        ensure!(
            pe || !["read2", "index2"].contains(&a.umi_loc.as_str()),
            "UMI location requires paired-end input"
        );
        ensure!(
            !["read1", "read2", "per_read"].contains(&a.umi_loc.as_str()) || a.umi_len > 0,
            "Read UMI extraction requires --umi_len"
        );
    }
    if let Some(n) = a.split {
        ensure!(
            (2..=999).contains(&n),
            "Split file count must be between 2 and 999"
        );
        ensure!(a.input != "-", "Cannot split stdin by file count");
    }
    if let Some(n) = a.split_by_lines {
        ensure!(
            n >= 1000 && n % 4 == 0,
            "Split line count must be >=1000 and a multiple of 4"
        );
    }
    ensure!(
        a.split_prefix_digits <= 10,
        "Split prefix digits must be between 0 and 10"
    );
    ensure!(
        !(a.stdout && (a.split.is_some() || a.split_by_lines.is_some())),
        "Cannot split stdout"
    );
    Ok(())
}

fn normalized_path(path: &str) -> Result<PathBuf> {
    let mut p = PathBuf::from(path);
    // Resolve dangling output symlinks too: their target can be another output
    // that has not been created yet. Bound traversal to reject symlink cycles.
    for _ in 0..40 {
        if p.exists() {
            return Ok(p.canonicalize()?);
        }
        if p.symlink_metadata()
            .is_ok_and(|m| m.file_type().is_symlink())
        {
            let target = std::fs::read_link(&p)?;
            p = if target.is_absolute() {
                target
            } else {
                p.parent().unwrap_or(Path::new(".")).join(target)
            };
            continue;
        }
        return output_entry_path(&p);
    }
    anyhow::bail!("Too many output symlinks: {path}")
}

fn output_entry_path(p: &Path) -> Result<PathBuf> {
    let parent = p
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    Ok(parent
        .canonicalize()?
        .join(p.file_name().context("Invalid output path")?))
}

fn same_file(a: &Path, b: &Path) -> bool {
    if a == b {
        return true;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(a), Ok(b)) = (a.metadata(), b.metadata()) {
            return a.dev() == b.dev() && a.ino() == b.ino();
        }
    }
    false
}

fn is_split_path(path: &Path, base: &Path, count: Option<usize>) -> bool {
    let Some(name) = base.file_name().and_then(|s| s.to_str()) else {
        return false;
    };
    let suffix = format!(".{name}");
    path.parent() == base.parent()
        && path
            .file_name()
            .and_then(|s| s.to_str())
            .and_then(|name| name.strip_suffix(&suffix))
            .is_some_and(|prefix| {
                !prefix.is_empty()
                    && prefix.bytes().all(|b| b.is_ascii_digit())
                    && prefix
                        .parse::<usize>()
                        .is_ok_and(|n| n > 0 && count.is_none_or(|max| n <= max))
            })
}

/// Check the entire output set before opening any writer, including reports and
/// numbered split files. In-place processing would truncate the source FASTQ.
fn validate_paths(a: &Args) -> Result<()> {
    let inputs: Vec<PathBuf> = std::iter::once(a.input.as_str())
        .chain(a.input2.as_deref())
        .filter(|s| !s.is_empty() && *s != "-" && !s.contains("://"))
        .map(normalized_path)
        .collect::<Result<_>>()?;
    let mut seen: BTreeMap<PathBuf, usize> = BTreeMap::new();
    let outputs = [
        Some(a.output.as_str()),
        a.output2.as_deref(),
        a.merged_out.as_deref(),
        a.unpaired1.as_deref(),
        a.unpaired2.as_deref(),
        a.failed_out.as_deref(),
        a.overlapped_out.as_deref(),
        (a.stats_format != "off").then_some(a.json.as_str()),
        Some(a.html.as_str()),
    ];
    for (i, path) in outputs
        .into_iter()
        .enumerate()
        .filter_map(|(i, p)| p.map(|p| (i, p)))
        .filter(|(_, p)| !p.is_empty() && *p != "-" && !p.starts_with("s3://"))
    {
        let p = normalized_path(path)?;
        ensure!(!p.is_dir(), "Output path is a directory: {path}");
        ensure!(
            !inputs.iter().any(|i| same_file(i, &p)),
            "Output would overwrite input: {path}"
        );
        if let Some((_, &previous)) = seen.iter().find(|(other, _)| same_file(other, &p)) {
            ensure!(
                previous == 3 && i == 4,
                "Output/report paths must be distinct: {path}"
            );
        }
        seen.insert(p, i);
        ensure!(
            !a.dont_overwrite || !Path::new(path).exists(),
            "Output already exists: {path}"
        );
    }
    if a.split.is_some() || a.split_by_lines.is_some() {
        let mut split_paths = Vec::<PathBuf>::new();
        let split_bases = [Some(a.output.as_str()), a.output2.as_deref()]
            .into_iter()
            .flatten()
            .filter(|s| !s.is_empty())
            .map(|base| output_entry_path(Path::new(base)))
            .collect::<Result<Vec<_>>>()?;
        for full in &split_bases {
            let parent = full.parent().context("Invalid split output directory")?;
            let matches = |p: &Path| is_split_path(p, full, a.split);
            for input in &inputs {
                ensure!(
                    !matches(input),
                    "Numbered split output would overwrite input: {}",
                    input.display()
                );
            }
            for path in seen.keys() {
                ensure!(
                    !matches(path),
                    "Numbered split output conflicts with another output/report: {}",
                    path.display()
                );
            }
            for entry in std::fs::read_dir(parent)? {
                let entry = entry?;
                if matches(&entry.path()) {
                    ensure!(
                        !a.dont_overwrite,
                        "Split output already exists: {}",
                        entry.path().display()
                    );
                    let path = normalized_path(
                        entry.path().to_str().context("Invalid split output path")?,
                    )?;
                    ensure!(
                        !inputs.iter().any(|input| same_file(input, &path)),
                        "Split output resolves to an input file: {}",
                        entry.path().display()
                    );
                    ensure!(
                        !seen
                            .keys()
                            .chain(split_paths.iter())
                            .any(|other| same_file(other, &path)),
                        "Split output resolves to another output/report: {}",
                        entry.path().display()
                    );
                    ensure!(
                        path == entry.path()
                            || !split_bases
                                .iter()
                                .any(|base| is_split_path(&path, base, a.split)),
                        "Split output resolves to another numbered output: {}",
                        entry.path().display()
                    );
                    split_paths.push(path);
                }
            }
        }
    }
    Ok(())
}

fn load_sequences(path: Option<&str>, fasta: bool) -> Result<Vec<Vec<u8>>> {
    let Some(path) = path else {
        return Ok(Vec::new());
    };
    let text = crate::io::read_to_string(path).with_context(|| format!("Cannot read {path}"))?;
    let mut seqs = Vec::<Vec<u8>>::new();
    for line in text.lines().map(str::trim).filter(|s| !s.is_empty()) {
        if fasta && line.starts_with('>') {
            seqs.push(Vec::new());
            continue;
        }
        if fasta {
            ensure!(!seqs.is_empty(), "Adapter FASTA must start with a > header");
            seqs.last_mut()
                .unwrap()
                .extend(line.as_bytes().iter().map(u8::to_ascii_uppercase));
        } else {
            seqs.push(line.as_bytes().to_vec());
        }
    }
    if fasta {
        ensure!(
            seqs.iter()
                .all(|s| s.len() >= 6 && s.iter().all(|b| b"ATCG".contains(b))),
            "FASTA adapters must contain at least six A/T/C/G bases"
        );
    }
    Ok(seqs)
}

pub fn run(a: Args) -> Result<()> {
    let start = Instant::now();
    let mut engine = Engine::new(a)?;
    while engine.step()? {}
    let done = engine.finish()?;
    write_reports(&done.args, &done.report, start)?;
    eprintln!(
        "Processed {} {}; passed {} reads. Compatibility target: fastp {FASTP_VERSION}",
        done.count,
        if done.pe { "pairs" } else { "reads" },
        done.passed
    );
    Ok(())
}

struct Finished {
    args: Args,
    report: Value,
    count: usize,
    passed: usize,
    pe: bool,
}

/// Where the main and merged reads go; `-` is stdout.
struct Targets {
    main1: Option<String>,
    main2: Option<String>,
    merged: Option<String>,
}

/// Counts and profiles accumulated in input order.
struct Stats {
    before: [Sample; 2],
    after: [Sample; 2],
    failures: [usize; 7],
    adapter_reads: usize,
    adapter_bases: usize,
    corrected_reads: usize,
    corrected_bases: usize,
    poly_reads: [usize; 4],
    poly_bases: [usize; 4],
    histogram: Vec<usize>,
    unknown: usize,
    duplicates: usize,
    count: usize,
    merged_pairs: usize,
}

impl Default for Stats {
    fn default() -> Self {
        Self {
            before: Default::default(),
            after: Default::default(),
            failures: [0; 7],
            adapter_reads: 0,
            adapter_bases: 0,
            corrected_reads: 0,
            corrected_bases: 0,
            poly_reads: [0; 4],
            poly_bases: [0; 4],
            histogram: vec![0; 513],
            unknown: 0,
            duplicates: 0,
            count: 0,
            merged_pairs: 0,
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
type PoolHandle = rayon::ThreadPool;
#[cfg(target_arch = "wasm32")]
type PoolHandle = &'static rayon::ThreadPool;

/// Read-only state shared with the pipeline threads.
struct Shared {
    config: Config,
    pool: PoolHandle,
    /// Empty before/after profiles that each chunk starts from.
    templates: [[Sample; 2]; 2],
    hasher: Option<DuplicationHash>,
    recycle_limit: usize,
}

impl Shared {
    fn pool(&self) -> &rayon::ThreadPool {
        &self.pool
    }

    fn process(&self, batch: ReadBatch) -> ProcessedBatch {
        let [before, after] = &self.templates;
        process_batch(
            batch,
            self.pool(),
            &self.config,
            before,
            after,
            self.hasher.as_ref(),
            self.recycle_limit,
        )
    }
}

type ProcessedBatch = Result<Vec<Chunk<'static>>>;

/// Where processed batches come from: this thread, or the read-ahead,
/// processing and duplicate-marking threads.
enum Source {
    Serial(Box<BatchReader>),
    Pipelined(Pipeline),
}

/// Pipeline threads and their channels. Dropping it stops the threads: the
/// processing thread finds its output closed, and the reader stops after its
/// current read.
struct Pipeline {
    processed: Option<crossbeam_channel::Receiver<ProcessedBatch>>,
    recycle: crossbeam_channel::Sender<Vec<ReadPool>>,
    threads: Vec<std::thread::JoinHandle<()>>,
}

impl Pipeline {
    fn start(
        mut reader: BatchReader,
        shared: &Arc<Shared>,
        mut duplication: Option<Duplication>,
    ) -> Self {
        // One queued batch per stage bounds both memory and reordering.
        let (recycle, recycle_rx) = crossbeam_channel::bounded::<Vec<ReadPool>>(2);
        let (input_tx, input_rx) = crossbeam_channel::bounded(1);
        let (processed_tx, processed) = crossbeam_channel::bounded::<ProcessedBatch>(1);
        let mut threads = vec![std::thread::spawn(move || loop {
            while let Ok(recycled) = recycle_rx.try_recv() {
                for recycled in recycled {
                    reader.pool.absorb(recycled);
                }
            }
            match reader.next_batch() {
                Ok(Some(batch)) => {
                    if input_tx.send(Ok(batch)).is_err() {
                        break;
                    }
                }
                Ok(None) => break,
                Err(error) => {
                    let _ = input_tx.send(Err(error));
                    break;
                }
            }
        })];
        let input_rx = if shared.config.dedup_before_filtering() {
            let mut filter = duplication.take().expect("deduplication requires a filter");
            let shared = Arc::clone(shared);
            let (ready, prepared) = crossbeam_channel::bounded(1);
            threads.push(std::thread::spawn(move || {
                for batch in input_rx {
                    let result = batch.map(|mut b| {
                        b.mark_duplicates(&shared, &mut filter);
                        b
                    });
                    let failed = result.is_err();
                    if ready.send(result).is_err() || failed {
                        break;
                    }
                }
            }));
            prepared
        } else {
            input_rx
        };
        let processed_tx = if duplication.is_some() {
            let (duplicate_tx, duplicate_rx) = crossbeam_channel::bounded::<ProcessedBatch>(1);
            threads.push(std::thread::spawn(move || {
                for result in duplicate_rx {
                    let result = result.map(|mut chunks| {
                        mark_duplicates(&mut chunks, duplication.as_mut());
                        chunks
                    });
                    let failed = result.is_err();
                    if processed_tx.send(result).is_err() || failed {
                        break;
                    }
                }
            }));
            duplicate_tx
        } else {
            processed_tx
        };
        let shared = Arc::clone(shared);
        threads.push(std::thread::spawn(move || {
            // A batch waits for its slowest chunk, often one on an efficiency
            // core. Starting the next batch meanwhile gives idle workers its
            // chunks. Results still leave in input order. A few workers are
            // never idle for long, and a second batch only costs memory.
            let in_flight = if shared.pool().current_num_threads() > 4 {
                2
            } else {
                1
            };
            type Running = crossbeam_channel::Receiver<std::thread::Result<ProcessedBatch>>;
            let mut running: VecDeque<Running> = VecDeque::with_capacity(in_flight);
            let mut input = input_rx.into_iter();
            let mut reading = true;
            loop {
                while reading && running.len() < in_flight {
                    let Some(batch) = input.next() else {
                        reading = false;
                        break;
                    };
                    let (done, result) = crossbeam_channel::bounded(1);
                    match batch {
                        Ok(batch) => {
                            let job = Arc::clone(&shared);
                            shared.pool().spawn(move || {
                                let result =
                                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                        job.process(batch)
                                    }));
                                let _ = done.send(result);
                            });
                        }
                        Err(error) => {
                            reading = false;
                            let _ = done.send(Ok(Err(error)));
                        }
                    }
                    running.push_back(result);
                }
                let Some(next) = running.pop_front() else {
                    break;
                };
                let result = match next.recv() {
                    Ok(Ok(result)) => result,
                    Ok(Err(panic)) => std::panic::resume_unwind(panic),
                    Err(_) => unreachable!("a processing job always reports"),
                };
                let failed = result.is_err();
                if processed_tx.send(result).is_err() || failed {
                    break;
                }
            }
            // Let started jobs finish before the pipeline reports it stopped.
            for next in running {
                let _ = next.recv();
            }
        }));
        Self {
            processed: Some(processed),
            recycle,
            threads,
        }
    }

    /// The next processed batch, or `None` once the threads have finished.
    fn next(&mut self) -> Option<ProcessedBatch> {
        let received = self.processed.as_ref()?.recv().ok();
        if received.is_none() {
            self.processed = None;
            // A thread that panicked closed its channel early; surface the
            // panic as `thread::scope` would instead of ending quietly.
            for thread in self.threads.drain(..) {
                if let Err(panic) = thread.join() {
                    std::panic::resume_unwind(panic);
                }
            }
        }
        received
    }
}

impl Drop for Pipeline {
    fn drop(&mut self) {
        self.processed = None;
        for thread in self.threads.drain(..) {
            let _ = thread.join();
        }
    }
}

/// The processing engine: set up by `new`, advanced one batch at a time by
/// `step`, and turned into the report by `finish`.
struct Engine {
    shared: Arc<Shared>,
    source: Source,
    /// Duplicate marking on this thread when not pipelined.
    duplication: Option<Duplication>,
    out: Output,
    targets: Targets,
    stats: Stats,
    pe: bool,
    max_bytes: usize,
    done: bool,
    /// A step failed, so the counts are incomplete.
    failed: bool,
}

impl Engine {
    fn new(mut a: Args) -> Result<Self> {
        if a.stdin || a.input == "/dev/stdin" {
            a.input = "-".into();
        }
        if a.output == "/dev/stdout" {
            a.output = "-".into();
        }
        if a.output == "-" {
            a.stdout = true;
        }
        if let (Some(first), Some(second)) = (&a.unpaired1, &a.unpaired2) {
            if !first.contains("://")
                && !second.contains("://")
                && same_file(&normalized_path(first)?, &normalized_path(second)?)
            {
                a.unpaired2 = Some(first.clone());
            }
        }
        a.umi_loc.make_ascii_lowercase();
        if a.merge {
            a.correction = true;
            if a.merged_out.is_none() && !a.stdout && !a.output.is_empty() && a.output2.is_none() {
                a.merged_out = Some(std::mem::take(&mut a.output));
            }
            if a.include_unmerged {
                a.output.clear();
                a.output2 = None;
                a.unpaired1 = None;
                a.unpaired2 = None;
            }
        }
        if a.stdout {
            a.output.clear();
            a.output2 = None;
        }
        if a.input2.is_none() && !a.interleaved_in {
            a.unpaired1 = None;
            a.unpaired2 = None;
            a.overlapped_out = None;
        }
        if !a.merge {
            a.merged_out = None;
        }
        validate(&a)?;
        // The in-memory file system has no symlinks or directories to check.
        #[cfg(feature = "native")]
        validate_paths(&a)?;
        let pe = a.input2.is_some() || a.interleaved_in;
        let total_for_split = if a.split.is_some() {
            let mut reader = Reader::open(&a.input, a.phred64)?;
            let mut n = 0;
            while reader.next()?.is_some() {
                n += 1;
            }
            if a.interleaved_in {
                ensure!(n % 2 == 0, "Odd read count in interleaved FASTQ");
                n /= 2;
            }
            if a.reads_to_process > 0 {
                n = n.min(a.reads_to_process);
            }
            Some(n)
        } else {
            None
        };
        // Decoding one large member in parallel costs about 2.5 times zlib's CPU
        // per byte and holds the regions in flight in memory, so it pays off
        // only when decoding limits the run. For a single input that is output
        // that is not gzip, or binned qualities, which make compressing the
        // output cheap. Paired inputs share the decoding threads and need
        // uncompressed output: any compression limits the run before two
        // sequential decoders do.
        let output = io::CompressionFormat::from_path(&a.output);
        let split = if a.input2.is_none() {
            output != io::CompressionFormat::Gzip || io::binned_gzip_input(&a.input)
        } else {
            output == io::CompressionFormat::None
        };
        io::set_split_decoding(if split {
            1 + usize::from(a.input2.is_some())
        } else {
            0
        });
        let mut reader1 = Reader::open(&a.input, a.phred64)?;
        let mut reader2 = a
            .input2
            .as_ref()
            .map(|s| Reader::open(s, a.phred64))
            .transpose()?;
        let mut sampled = Vec::new();
        let mut sampled_raw = VecDeque::new();
        // Compact-output jobs retain evaluation data in the same blocks used
        // by processing. The sample is read once without allocating four fields per
        // record, and each block is released as soon as processing consumes it.
        let pipelined = cfg!(feature = "native")
            && a.max_backlog != Some(1)
            && std::iter::once(a.input.as_str())
                .chain(a.input2.as_deref())
                .all(|path| std::fs::metadata(path).is_ok_and(|m| m.is_file()));
        // Borrowed parsing also works synchronously. Only reading ahead requires
        // regular files; pipes must stay on this thread for prompt cancellation.
        let raw_input = pack_main_output(&a);
        // Read once and retain the sample: this works for pipes, HTTP, and S3 too.
        let auto = |seq: &Option<String>| seq.as_deref().is_none_or(|s| s == "auto");
        let needs_adapter_sample = !a.disable_adapter_trimming
            && !a.disable_adapter_detection
            && if pe {
                a.detect_adapter_for_pe
                    && (auto(&a.adapter_sequence) || auto(&a.adapter_sequence_r2))
            } else {
                auto(&a.adapter_sequence)
            };
        // Otherwise only the first header is needed for instrument/polyG detection.
        let sample_limit = if a.input == "-" {
            0
        } else if needs_adapter_sample || a.overrepresentation_analysis {
            256 * 1024
        } else {
            1
        };
        let mut sample_bases = 0;
        let mut sample_count = 0;
        if raw_input {
            while sample_count < sample_limit && sample_bases < 151 * 256 * 1024 {
                let budget = a.batch_bytes.min(2 * 1024 * 1024);
                let remaining = sample_limit - sample_count;
                let base_limit = 151 * 256 * 1024 - sample_bases;
                let block = if a.interleaved_in {
                    block::read_sample_batch::<true>(&mut reader1, None, budget, remaining, base_limit)?
                } else {
                    block::read_sample_batch::<false>(
                        &mut reader1,
                        reader2.as_mut(),
                        budget,
                        remaining,
                        base_limit,
                    )?
                };
                let Some(blocks) = block else { break };
                let stride = if a.interleaved_in { 2 } else { 1 };
                // Validate the entire evaluation sample even when a read limit
                // later selects only its prefix. The base budget counts R1 only.
                for index in 0..blocks[0].len() {
                    blocks[0].read(index, sample_count * stride + index, a.phred64)?;
                    if blocks[1].len() > 0 {
                        blocks[1].read(index, sample_count + index, a.phred64)?;
                    }
                    if index % stride == 0 {
                        sample_bases += blocks[0].sequence(index).len();
                    }
                }
                sample_count += blocks[0].len() / stride;
                sampled_raw.push_back(ReadBatch {
                    duplicate_flags: None,
                    raw: Some(blocks),
                    reads: Vec::new(),
                    start: 0,
                });
            }
        } else {
            for _ in 0..sample_limit {
                let Some(r1) = reader1.next()? else {
                    break;
                };
                let r2 = if a.interleaved_in {
                    Some(
                        reader1
                            .next()?
                            .context("Odd read count in interleaved FASTQ")?,
                    )
                } else if let Some(reader2) = reader2.as_mut() {
                    Some(
                        reader2
                            .next()?
                            .context("Paired input files have different read counts")?,
                    )
                } else {
                    None
                };
                sample_bases += r1.seq.len();
                sampled.push((r1, r2));
                if sample_bases >= 151 * 256 * 1024 {
                    break;
                }
            }
            sample_count = sampled.len();
        }
        let manual = |s: &Option<String>| {
            s.as_ref().filter(|s| s.as_str() != "auto").map(|s| {
                s.as_bytes()
                    .iter()
                    .map(u8::to_ascii_uppercase)
                    .collect::<Vec<_>>()
            })
        };
        let mut adapters = [manual(&a.adapter_sequence), manual(&a.adapter_sequence_r2)];
        for seq in adapters.iter().flatten() {
            ensure!(
                seq.len() >= 4 && seq.iter().all(|b| b"ATCG".contains(b)),
                "Adapter sequences must contain at least four A/T/C/G bases"
            );
        }
        let sample_sequences = |mate: usize| -> Vec<&[u8]> {
            if raw_input {
                sampled_raw
                    .par_iter()
                    .flat_map_iter(|b| {
                        let blocks = b.raw.as_ref().unwrap();
                        let (records, start, stride) = if a.interleaved_in {
                            (&blocks[0], mate, 2)
                        } else {
                            (&blocks[mate], 0, 1)
                        };
                        (start..records.len())
                            .step_by(stride)
                            .map(|i| records.sequence(i))
                    })
                    .collect()
            } else {
                sampled
                    .iter()
                    .filter_map(|(r1, r2)| {
                        if mate == 0 {
                            Some(r1.seq.as_ref())
                        } else {
                            r2.as_ref().map(|r| r.seq.as_ref())
                        }
                    })
                    .collect()
            }
        };
        let threads = worker_threads(&a);
        #[cfg(not(target_arch = "wasm32"))]
        let owned_pool = worker_pool(threads)?;
        #[cfg(not(target_arch = "wasm32"))]
        let pool = &owned_pool;
        #[cfg(target_arch = "wasm32")]
        let owned_pool = wasm_pool()?;
        #[cfg(target_arch = "wasm32")]
        let pool = owned_pool;
        // Detection scans the whole sample; it runs on the worker pool.
        let detect = |mate: usize| {
            pool.install(|| auto_adapter::detect(&sample_sequences(mate), a.trim_tail1))
        };
        if !a.disable_adapter_trimming && !a.disable_adapter_detection && sample_count >= 10000 {
            if !pe && adapters[0].is_none() {
                adapters[0] = detect(0);
            } else if pe && a.detect_adapter_for_pe {
                if adapters[0].is_none() {
                    adapters[0] = detect(0);
                }
                if adapters[1].is_none() {
                    adapters[1] = detect(1);
                }
            }
        }
        if pe && adapters[1].is_none() && !a.detect_adapter_for_pe {
            adapters[1] = adapters[0].clone();
        }
        let first_name = sampled.first().map(|(r, _)| r.name.as_ref()).or_else(|| {
            sampled_raw
                .front()
                .map(|batch| batch.raw.as_ref().unwrap()[0].name(0))
        });
        let poly_g = !a.disable_trim_poly_g
            && (a.trim_poly_g
                || first_name.is_some_and(|name| {
                    [
                        "@FS", "@MN", "@SH", "@NS", "@NB", "@NDX", "@VL", "@VH", "@A", "@NA", "@LH",
                    ]
                    .iter()
                    .any(|s| name.starts_with(s.as_bytes()))
                }));
        let fasta = load_sequences(a.adapter_fasta.as_deref(), true)?;
        let blacklist = [
            load_sequences(a.filter_by_index1.as_deref(), false)?,
            load_sequences(a.filter_by_index2.as_deref(), false)?,
        ];
        let max_bytes = if let Some(mb) = a.max_memory {
            ensure!(mb >= 16, "--max_memory must be at least 16 MB");
            a.batch_bytes.min(mb.saturating_mul(1024 * 1024) / 8).max(1)
        } else {
            a.batch_bytes
        };
        let config = Config {
            packed_output: raw_input,
            args: a,
            adapters: adapters.map(|a| a.map(adapter::AdapterMatcher::new)),
            fasta: fasta
                .into_iter()
                .map(adapter::AdapterMatcher::new)
                .collect(),
            blacklist,
            poly_g,
        };
        let a = &config.args;
        let split = split::SplitConfig {
            mode: if let Some(n) = a.split {
                split::SplitMode::ByFiles(n)
            } else if let Some(n) = a.split_by_lines {
                split::SplitMode::ByLines(n)
            } else {
                split::SplitMode::None
            },
            prefix_digits: a.split_prefix_digits,
        };
        let mut out = Output::new(a);
        out.open_main(a, &split)?;
        for path in [
            a.merged_out.as_deref(),
            a.unpaired1.as_deref(),
            a.unpaired2.as_deref(),
            a.failed_out.as_deref(),
            a.overlapped_out.as_deref(),
        ]
        .into_iter()
        .flatten()
        {
            out.open(path, a, split::SplitConfig::default())?;
        }
        if a.stdout || a.output == "-" {
            out.open("-", a, split::SplitConfig::default())?;
        }
        if let Some(total) = total_for_split {
            for w in out.writers.values_mut() {
                w.set_total_reads(total)?;
            }
        }
        let to_stdout = (a.stdout || a.output == "-") && !a.merge;
        let targets = Targets {
            main1: Some(if to_stdout {
                "-".into()
            } else {
                a.output.clone()
            }),
            main2: if to_stdout {
                Some("-".into())
            } else {
                a.output2.clone()
            },
            merged: if a.stdout {
                Some("-".into())
            } else {
                a.merged_out.clone()
            },
        };
        let mut stats = Stats::default();
        let mut duplication = if a.dedup || !a.dont_eval_duplication {
            Some(Duplication::new(if a.dup_calc_accuracy > 0 {
                a.dup_calc_accuracy
            } else if a.dedup {
                3
            } else {
                1
            }))
        } else {
            None
        };
        if a.overrepresentation_analysis {
            for i in 0..if pe { 2 } else { 1 } {
                let reads: Vec<_> = sampled
                    .iter()
                    .filter_map(|(r1, r2)| if i == 0 { Some(r1) } else { r2.as_ref() })
                    .collect();
                let len = reads
                    .iter()
                    .take(1000)
                    .map(|r| r.seq.len())
                    .max()
                    .unwrap_or(0);
                let candidates = Overrep::candidates(&reads, len);
                stats.before[i].overrep = Some(Overrep::new(
                    &candidates,
                    len,
                    a.overrepresentation_sampling,
                ));
                stats.after[i].overrep = Some(Overrep::new(
                    &candidates,
                    len,
                    a.overrepresentation_sampling,
                ));
            }
        }
        // Only regular local files are read ahead. Pipes and remote streams retain
        // synchronous reads so a failed output cannot leave a reader blocked while
        // the pipeline threads are being joined.
        // Compact synchronous jobs use the same batch ceiling as the pipeline,
        // limiting serial parsing work before each parallel processing step.
        // Smaller batches leave the worker pool waiting between them.
        let max_bytes = if pipelined || config.packed_output {
            max_bytes.min(8 * 1024 * 1024)
        } else {
            max_bytes
        };
        // Capacity includes separators and modest Vec growth beyond FASTQ payload.
        let recycle_limit = max_bytes.saturating_add(max_bytes / 4);
        let batch_reader = BatchReader {
            raw: raw_input,
            reader1,
            reader2,
            sampled: sampled.into_iter(),
            sampled_raw,
            pool: ReadPool {
                limit: recycle_limit,
                ..ReadPool::default()
            },
            interleaved: a.interleaved_in,
            max_bytes,
            limit: a.reads_to_process,
            count: 0,
            eof: false,
        };
        let templates = [
            std::array::from_fn(|i| stats.before[i].empty_like()),
            std::array::from_fn(|i| stats.after[i].empty_like()),
        ];
        let hasher = duplication.as_ref().map(Duplication::hasher);
        let shared = Arc::new(Shared {
            config,
            pool: owned_pool,
            templates,
            hasher,
            recycle_limit,
        });
        let source = if pipelined {
            Source::Pipelined(Pipeline::start(batch_reader, &shared, duplication.take()))
        } else {
            Source::Serial(Box::new(batch_reader))
        };
        Ok(Self {
            shared,
            source,
            duplication,
            out,
            targets,
            stats,
            pe,
            max_bytes,
            done: false,
            failed: false,
        })
    }

    /// Process the next batch. Returns `false` once the input is exhausted;
    /// after an error it returns the error once and then `false`.
    fn step(&mut self) -> Result<bool> {
        if self.done {
            return Ok(false);
        }
        let result = self.process_next();
        self.done = !matches!(result, Ok(true));
        self.failed = result.is_err();
        result
    }

    fn process_next(&mut self) -> Result<bool> {
        let shared = Arc::clone(&self.shared);
        let a = &shared.config.args;
        let pe = self.pe;
        let max_bytes = self.max_bytes;
        let mut chunks = match &mut self.source {
            Source::Pipelined(pipeline) => match pipeline.next() {
                Some(result) => result?,
                None => return Ok(false),
            },
            Source::Serial(reader) => match reader.next_batch()? {
                Some(mut batch) => {
                    let mut filter = self.duplication.as_mut();
                    if shared.config.dedup_before_filtering() {
                        batch.mark_duplicates(
                            &shared,
                            filter.take().expect("deduplication requires a filter"),
                        );
                    }
                    let mut chunks = shared.process(batch)?;
                    mark_duplicates(&mut chunks, filter);
                    chunks
                }
                None => return Ok(false),
            },
        };
        let Self {
            out,
            stats,
            targets,
            ..
        } = self;
        let (main1, main2, merged) = (
            targets.main1.as_deref(),
            targets.main2.as_deref(),
            targets.merged.as_deref(),
        );
        let mut read_pool = ReadPool {
            reads: Vec::with_capacity(chunks.iter().map(|c| c.processed.len()).sum()),
            limit: max_bytes,
            ..ReadPool::default()
        };
        let mut read_pools = Vec::with_capacity(chunks.len() + 1);
        let profile_after_in_worker = shared.config.profile_after_in_worker();
        // Bloom decisions remain ordered. Determine global after-filter indices
        // before parallel profiling so sampling is independent of chunk boundaries.
        let mut offsets = std::array::from_fn(|i| stats.after[i].simple.total_reads);
        for chunk in &mut chunks {
            chunk.after_offset = offsets;
            for p in &mut chunk.processed {
                stats.duplicates += usize::from(p.duplicate);
                if !profile_after_in_worker {
                    for (i, _) in p.passing_reads(a).into_iter().flatten() {
                        offsets[i] += 1;
                    }
                }
            }
        }
        if !profile_after_in_worker {
            shared.pool().install(|| {
                chunks.par_iter_mut().for_each(|chunk| {
                    let mut offsets = chunk.after_offset;
                    for p in &chunk.processed {
                        for (i, r) in p.passing_reads(a).into_iter().flatten() {
                            chunk.after[i].add_at(r, a.no_kmer, offsets[i]);
                            offsets[i] += 1;
                        }
                    }
                })
            });
        }
        for chunk in chunks {
            for (dst, src) in stats.before.iter_mut().zip(chunk.before) {
                dst.merge(src);
            }
            for (dst, src) in stats.after.iter_mut().zip(chunk.after) {
                dst.merge(src);
            }
            if let Some(mut packed) = chunk.packed {
                let [first, second, third] = &mut packed.output;
                out.write_packed_mates([main1, main2], [first, second])?;
                out.write_packed(merged, third)?;
                stats.count += packed.count;
                stats.merged_pairs += packed.merged_pairs;
                stats.duplicates += packed.duplicates;
                stats.adapter_reads += packed.adapter_reads;
                stats.adapter_bases += packed.adapter_bases;
                stats.corrected_reads += packed.corrected_reads;
                stats.corrected_bases += packed.corrected_bases;
                stats.unknown += packed.unknown;
                for (dst, src) in stats.failures.iter_mut().zip(packed.failures) {
                    *dst += src;
                }
                for (dst, src) in stats.histogram.iter_mut().zip(&packed.histogram) {
                    *dst += src;
                }
                for (dst, src) in stats.poly_reads.iter_mut().zip(packed.poly_reads) {
                    *dst += src;
                }
                for (dst, src) in stats.poly_bases.iter_mut().zip(packed.poly_bases) {
                    *dst += src;
                }
                read_pools.push(std::mem::take(&mut packed.recycled));
                continue;
            }
            for p in chunk.processed {
                let p = RecycleRead {
                    processed: p,
                    pool: &mut read_pool,
                };
                stats.count += 1;
                let duplicate = p.duplicate;
                if p.skipped {
                    continue;
                }
                stats.corrected_reads += p.corrected_reads;
                stats.corrected_bases += p.corrected_bases;
                stats.adapter_reads += p.adapter_reads;
                stats.adapter_bases += p.adapter_bases;
                for &(base, n) in &p.poly_x {
                    stats.poly_reads[base] += 1;
                    stats.poly_bases[base] += n;
                }
                if pe {
                    if let Some(size) = p.insert {
                        stats.histogram[size.min(512)] += 1;
                    } else {
                        stats.unknown += 1;
                    }
                }
                if let Some(r) = &p.overlapped {
                    out.write(a.overlapped_out.as_deref(), r, None)?;
                }
                if let Some(r) = &p.merged {
                    stats.failures[p.failures[0] as usize] += 2;
                    if p.failures[0] == Failure::Pass && !(a.dedup && duplicate) {
                        out.write(merged, r, None)?;
                        stats.merged_pairs += 1;
                    }
                    continue;
                }
                if a.merge && a.include_unmerged {
                    for i in 0..2 {
                        if let Some(r) = &p.reads[i] {
                            stats.failures[p.failures[i] as usize] += 1;
                            if p.failures[i] == Failure::Pass && !(a.dedup && duplicate) {
                                out.write(merged, r, None)?;
                            }
                        }
                    }
                    continue;
                }
                let f = if pe {
                    p.failures[0].max(p.failures[1])
                } else {
                    p.failures[0]
                };
                stats.failures[f as usize] += if pe { 2 } else { 1 };
                if a.dedup && duplicate {
                    continue;
                }
                if f == Failure::Pass {
                    if let Some(r) = &p.reads[0] {
                        out.write(main1, r, None)?;
                    }
                    if let Some(r) = &p.reads[1] {
                        out.write(main2, r, None)?;
                    }
                } else {
                    // fastp writes failed paired reads only when the mate passes.
                    if pe && p.failures.iter().all(|f| *f != Failure::Pass) {
                        continue;
                    }
                    for i in 0..if pe { 2 } else { 1 } {
                        let Some(r) = &p.reads[i] else {
                            continue;
                        };
                        let unpaired = if i == 0 {
                            a.unpaired1.as_deref()
                        } else {
                            a.unpaired2.as_deref().or(a.unpaired1.as_deref())
                        };
                        if pe && p.failures[i] == Failure::Pass && unpaired.is_some() {
                            out.write(unpaired, r, None)?;
                        } else {
                            out.write(
                                a.failed_out.as_deref(),
                                p.failed_reads[i].as_ref().unwrap_or(r),
                                Some(p.failures[i].tag()),
                            )?;
                        }
                    }
                }
            }
        }
        if a.verbose {
            eprintln!(
                "Processed {} {}",
                stats.count,
                if pe { "read pairs" } else { "reads" }
            );
        }
        read_pools.push(read_pool);
        match &mut self.source {
            Source::Pipelined(pipeline) => {
                let _ = pipeline.recycle.try_send(read_pools);
            }
            Source::Serial(reader) => {
                for read_pool in read_pools {
                    reader.pool.absorb(read_pool);
                }
            }
        }
        Ok(true)
    }

    /// Process any remaining input, close the outputs and build the report.
    /// Fails if a step failed: the report would cover only part of the input.
    fn finish(mut self) -> Result<Finished> {
        while self.step()? {}
        ensure!(!self.failed, "Processing stopped after an error");
        drop(self.source);
        let Self {
            shared,
            out,
            mut stats,
            pe,
            ..
        } = self;
        let a = &shared.config.args;
        out.finish()?;
        for sample in stats.before.iter_mut().chain(&mut stats.after) {
            sample.finish_summary();
        }
        let mut report = make_report(&stats.before, &stats.after, pe, a.merge);
        report["filtering_result"] = json!({"passed_filter_reads":stats.failures[0],"low_quality_reads":stats.failures[4],"too_many_N_reads":stats.failures[1],"too_short_reads":stats.failures[2],"too_long_reads":stats.failures[3],"low_complexity_reads":stats.failures[5],"adapter_dimer_reads":stats.failures[6]});
        if a.correction {
            report["filtering_result"]["corrected_reads"] = json!(stats.corrected_reads);
            report["filtering_result"]["corrected_bases"] = json!(stats.corrected_bases);
        }
        if !a.dont_eval_duplication || a.dedup {
            report["duplication"] = json!({"rate":if stats.count>0 {stats.duplicates as f64/stats.count as f64} else {0.0}});
        }
        if !a.disable_adapter_trimming {
            report["adapter_cutting"] = json!({"adapter_trimmed_reads":stats.adapter_reads,"adapter_trimmed_bases":stats.adapter_bases});
        }
        if a.trim_poly_x {
            report["polyx_trimming"] = json!({"total_polyx_trimmed_reads":stats.poly_reads.iter().sum::<usize>(),"total_polyx_trimmed_bases":stats.poly_bases.iter().sum::<usize>(),"polyx_trimmed_reads":stats.poly_reads,"polyx_trimmed_bases":stats.poly_bases});
        }
        if pe {
            let peak = (0..512)
                .max_by_key(|&i| (stats.histogram[i], 512 - i))
                .unwrap();
            report["insert_size"] = json!({"peak":peak,"unknown":stats.unknown+stats.histogram[512],"histogram":&stats.histogram[..512]});
        }
        if a.merge {
            report["merged_and_filtered"] = stats.after[0].details();
            report["merged_pairs"] = json!(stats.merged_pairs);
        }
        let passed = stats.failures[0];
        let count = stats.count;
        drop(stats);
        let shared = Arc::into_inner(shared).context("pipeline threads still running")?;
        Ok(Finished {
            report,
            count,
            passed,
            pe,
            args: shared.config.args,
        })
    }
}

fn make_report(before: &[Sample; 2], after: &[Sample; 2], pe: bool, merge: bool) -> Value {
    fn summary(samples: &[Sample; 2], pe: bool) -> Value {
        let mut total = stats::SimpleStats::default();
        for s in samples {
            total.total_reads += s.simple.total_reads;
            total.total_bases += s.simple.total_bases;
            total.q20_bases += s.simple.q20_bases;
            total.q30_bases += s.simple.q30_bases;
            total.q40_bases += s.simple.q40_bases;
            total.gc_bases += s.simple.gc_bases;
        }
        let mut value = serde_json::to_value(total.to_read_stats()).unwrap();
        value["read1_mean_length"] = json!(samples[0].simple.to_read_stats().read1_mean_length);
        if pe {
            value["read2_mean_length"] = json!(samples[1].simple.to_read_stats().read1_mean_length);
        }
        value
    }
    let mut report = json!({"summary":{"fasterp_version":env!("CARGO_PKG_VERSION"),"fastp_compatibility_version":FASTP_VERSION,"sequencing":if pe {format!("paired end ({} cycles + {} cycles)",before[0].max_len,before[1].max_len)} else {format!("single end ({} cycles)",before[0].max_len)},"before_filtering":summary(before,pe),"after_filtering":summary(after,pe && !merge)},"filtering_result":{},"read1_before_filtering":before[0].details(),"read1_after_filtering":after[0].details()});
    if pe {
        report["read2_before_filtering"] = before[1].details();
        report["read2_after_filtering"] = after[1].details();
    }
    report
}

fn write_reports(a: &Args, report: &Value, start: Instant) -> Result<()> {
    if a.stats_format != "off" {
        let mut writer = BufWriter::new(crate::io::create_file(&a.json)?);
        if a.stats_format == "pretty" {
            serde_json::to_writer_pretty(&mut writer, report)?;
        } else {
            serde_json::to_writer(&mut writer, report)?;
        }
        if a.stats_format == "jsonl" {
            writeln!(writer)?;
        }
        writer.flush()?;
    }
    let mut typed: stats::FasterpReport = serde_json::from_value(report.clone())?;
    for details in std::iter::once(&mut typed.read1_before_filtering)
        .chain(typed.read2_before_filtering.iter_mut())
        .chain(typed.read1_after_filtering.iter_mut())
        .chain(typed.read2_after_filtering.iter_mut())
    {
        let c = &mut details.content_curves;
        for curve in [&mut c.a, &mut c.t, &mut c.c, &mut c.g, &mut c.n, &mut c.gc] {
            for value in curve {
                *value *= 100.0;
            }
        }
    }
    let config = html::HtmlReportConfig {
        input: a.input.clone(),
        output: a.output.clone(),
        json: a.json.clone(),
        input2: a.input2.clone(),
        output2: a.output2.clone(),
    };
    html::generate_html_report(&typed, &config, &a.html, start.elapsed())?;
    let mut text = crate::io::read_to_string(&a.html)?;
    if a.report_title != "fasterp report" {
        text = text.replace(
            "fasterp Quality Control Report",
            &escape_html(&a.report_title),
        );
    }
    let mut extra=format!("<section><h2>Additional processing statistics</h2><p>fastp compatibility target: {FASTP_VERSION}</p><table>");
    for field in ["adapter_dimer_reads", "corrected_reads", "corrected_bases"] {
        if let Some(value) = report["filtering_result"].get(field) {
            extra.push_str(&format!("<tr><td>{field}</td><td>{value}</td></tr>"));
        }
    }
    if let Some(value) = report.get("merged_pairs") {
        extra.push_str(&format!("<tr><td>merged_pairs</td><td>{value}</td></tr>"));
    }
    if let Some(poly) = report.get("polyx_trimming") {
        for field in ["total_polyx_trimmed_reads", "total_polyx_trimmed_bases"] {
            extra.push_str(&format!(
                "<tr><td>{field}</td><td>{}</td></tr>",
                poly[field]
            ));
        }
    }
    extra.push_str("</table></section>");
    if a.overrepresentation_analysis {
        extra.push_str(&format!(
            "<section><h2>Overrepresented sequences</h2><p>Sampling: 1 in {} reads</p>",
            a.overrepresentation_sampling
        ));
        for stage in [
            "read1_before_filtering",
            "read2_before_filtering",
            "read1_after_filtering",
            "read2_after_filtering",
        ] {
            let Some(sequences) = report[stage]["overrepresented_sequences"].as_object() else {
                continue;
            };
            extra.push_str(&format!(
                "<h3>{stage}</h3><table><tr><th>Sequence</th><th>Sampled occurrences</th></tr>"
            ));
            for (seq, n) in sequences {
                extra.push_str(&format!(
                    "<tr><td>{}</td><td>{n}</td></tr>",
                    escape_html(seq)
                ));
            }
            extra.push_str("</table>");
        }
        extra.push_str("</section>");
    }
    crate::io::write_file(
        &a.html,
        text.replace("</body>", &format!("{extra}</body>"))
            .as_bytes(),
    )?;
    Ok(())
}

fn escape_html(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn workers_have_room_for_nested_batches() {
        #[inline(never)]
        fn descend(levels: usize) -> usize {
            // One level's worth of chunk statistics.
            let frame = std::hint::black_box([1u8; 240 << 10]);
            if levels == 0 {
                return frame[0] as usize;
            }
            descend(levels - 1) + frame[frame.len() - 1] as usize
        }
        // Four times the deepest nesting seen with two batches in flight.
        let pool = worker_pool(2).unwrap();
        assert_eq!(pool.install(|| descend(44)), 45);
    }

    #[test]
    fn overrep_counts_repeated_window_lengths_and_excludes_final_window() {
        let read = Read {
            seq: Cow::Borrowed(b"AAAAAAAAAACC"),
            ..Read::default()
        };
        let reads = vec![&read; 250];
        // For twelve cycles the length-dependent step and fixed ten-base
        // step are equal. Both contribute toward the discovery threshold.
        assert_eq!(
            Overrep::candidates(&reads, 12),
            BTreeMap::from([(b"AAAAAAAAAA".to_vec(), 500), (b"AAAAAAAAAC".to_vec(), 500)])
        );
    }

    #[test]
    fn overrep_bucket_collisions_still_require_exact_counts() {
        // The 20-base keys share both eight-byte ends and therefore a bucket.
        // Only the first key reaches the threshold of 100. The trailing byte
        // gives each read exactly one eligible 20-base window.
        let reads: Vec<_> = (0..5000)
            .map(|i| {
                let mut seq = b"ACGTACGTAAAAACGTACGTN".to_vec();
                if i >= 100 {
                    let key = 1 + (i - 100) % 255;
                    for (j, base) in seq[8..12].iter_mut().enumerate() {
                        *base = b"ACGT"[(key >> (2 * j)) & 3];
                    }
                }
                Read {
                    seq: Cow::Owned(seq),
                    ..Read::default()
                }
            })
            .collect();
        let found = Overrep::candidates(&reads.iter().collect::<Vec<_>>(), 24);
        let exact: Vec<_> = found.iter().filter(|(seq, _)| seq.len() == 20).collect();
        assert_eq!(exact, vec![(&reads[0].seq[..20].to_vec(), &100)]);
    }

    #[test]
    fn quality_validation_matches_scalar_at_every_byte_and_vector_boundary() {
        for offset in [33, 64] {
            for len in [0, 1, 15, 16, 17, 31, 32, 33, 51, 64, 65, 148, 513] {
                let mut quality = vec![offset; len];
                assert!(valid_quality(&quality, offset));
                for i in 0..len {
                    for byte in 0..=255 {
                        quality[i] = byte;
                        assert_eq!(
                            valid_quality(&quality, offset),
                            byte >= offset && byte <= 126,
                            "offset={offset}, len={len}, index={i}, byte={byte}"
                        );
                    }
                    quality[i] = offset;
                }
            }
        }
    }

    #[test]
    fn vectorized_sample_matches_per_base_reference() {
        // Lengths up to 300 cover every tail of the four-base loop, and
        // one-read samples merged in exercise flush and merge.
        let mut fast = Sample::default();
        let mut merged = Sample::default();
        let mut kmers = [0; 1024];
        let mut cycles = vec![[(0u64, 0u64); 5]; 301];
        let mut state = 140u64;
        for ordinal in 0..2000 {
            let mut seq = Vec::new();
            let mut qual = Vec::new();
            for _ in 0..(ordinal * 7) % 301 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                seq.push(b"ATCGNatcgnX."[(state >> 32) as usize % 12]);
                qual.push(33 + ((state >> 24) % 94) as u8);
            }
            fast.add_parts(&seq, &qual, false, ordinal);
            let mut part = Sample::default();
            part.add_parts(&seq, &qual, false, ordinal);
            merged.merge(part);
            crate::kmer::count_k5_2bit(&seq, &mut kmers);
            for (i, (&b, &q)) in seq.iter().zip(&qual).enumerate() {
                let base = b"ACGT"
                    .iter()
                    .position(|&c| c == b.to_ascii_uppercase())
                    .unwrap_or(4);
                cycles[i][base].0 += 1;
                cycles[i][base].1 += u64::from(q - 33);
            }
        }
        for sample in [&mut fast, &mut merged] {
            sample.finish_summary();
        }
        assert_eq!(merged.details(), fast.details());
        for sample in [&fast, &merged] {
            assert_eq!(sample.kmers, kmers);
            for (i, expected) in cycles.iter().enumerate() {
                let actual = sample
                    .cycles
                    .get(i)
                    .map_or([(0, 0); 5], |cycle| cycle.map(|b| (b.count, b.quality_sum)));
                assert_eq!(actual, *expected, "cycle {i}");
            }
        }
    }

    #[test]
    fn shared_prefix_and_removed_tail_reconstruct_full_read_statistics() {
        let seq: Vec<_> = b"gggACGTcacgtNXaTACGA"
            .iter()
            .copied()
            .cycle()
            .take(105)
            .collect();
        let qual: Vec<_> = (33..=126).cycle().take(seq.len()).collect();
        for no_kmer in [false, true] {
            for kept in 0..=seq.len() {
                let mut expected = Sample::default();
                expected.add_parts(&seq, &qual, no_kmer, 0);
                let mut prefix = Sample::default();
                prefix.add_parts(&seq[..kept], &qual[..kept], no_kmer, 0);
                let mut before = Sample::default();
                before.add_removed_tail(&seq, &qual, kept, no_kmer);
                before.merge(prefix);
                assert_eq!(before.simple.total_reads, expected.simple.total_reads);
                assert_eq!(before.simple.total_bases, expected.simple.total_bases);
                assert_eq!(before.max_len, expected.max_len);
                assert_eq!(before.kmers, expected.kmers);
                assert_eq!(before.qual_hist, expected.qual_hist);
                for (actual, expected) in before.cycles.iter().zip(&expected.cycles) {
                    assert_eq!(
                        actual.map(|b| (b.count, b.quality_sum)),
                        expected.map(|b| (b.count, b.quality_sum))
                    );
                }
            }
        }
    }

    #[test]
    fn fused_statistics_match_separate_quality_and_kmer_passes() {
        let mut sample = Sample::default();
        let mut simple = stats::SimpleStats::default();
        let mut kmers = [0; 1024];
        let mut qualities = [0; 94];
        let mut cycles = vec![[(0u64, 0u64); 5]; 256];
        let mut state = 140u64;
        for ordinal in 0..1000 {
            let mut read = Read::default();
            for _ in 0..ordinal % 257 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                read.seq
                    .to_mut()
                    .push(b"ATCGNatcgnX"[(state >> 32) as usize % 11]);
                read.qual.to_mut().push(33 + ((state >> 24) % 94) as u8);
            }
            sample.add_at(&read, false, ordinal);
            let metrics = simd::compute_stats(&read.seq, &read.qual, 0);
            simple.add(
                read.seq.len(),
                metrics.q20,
                metrics.q30,
                metrics.q40,
                metrics.gc,
            );
            crate::kmer::count_k5_2bit(&read.seq, &mut kmers);
            for (i, (&b, &q)) in read.seq.iter().zip(read.qual.iter()).enumerate() {
                qualities[usize::from(q - 33)] += 1;
                let base = b"ACGT"
                    .iter()
                    .position(|&c| c == b.to_ascii_uppercase())
                    .unwrap_or(4);
                cycles[i][base].0 += 1;
                cycles[i][base].1 += u64::from(q - 33);
            }
        }
        sample.finish_summary();
        assert_eq!(
            serde_json::to_value(sample.simple.to_read_stats()).unwrap(),
            serde_json::to_value(simple.to_read_stats()).unwrap()
        );
        assert_eq!(sample.kmers, kmers);
        for (cycle, expected) in sample.cycles.iter().zip(cycles) {
            assert_eq!(cycle.map(|base| (base.count, base.quality_sum)), expected);
        }
        for (q, expected) in qualities.into_iter().enumerate() {
            assert_eq!(
                sample.qual_hist.iter().map(|h| h[q]).sum::<usize>(),
                expected
            );
        }
    }

    #[test]
    fn specialized_duplication_hash_matches_scalar_reference() {
        for level in [1, 2, 3, 6] {
            let duplication = Duplication::new(level);
            let hash = duplication.hasher();
            for length in [
                0, 1, 15, 16, 17, 31, 32, 33, 51, 102, 148, 255, 256, 296, 511, 512, 513, 1024,
                2049,
            ] {
                let alphabet = b"ATCGNnatcg\0\xff";
                let seq: Vec<_> = (0..length)
                    .map(|i| {
                        if length == 512 {
                            i as u8
                        } else {
                            alphabet[(i * 17 + length) % alphabet.len()]
                        }
                    })
                    .collect();
                for split in [
                    length,
                    length / 3,
                    length / 2,
                    length.min(15),
                    length.min(17),
                ] {
                    let (first, second) = seq.split_at(split);
                    let mut expected = [0u64; 8];
                    let n = duplication.buffers.len();
                    for (p, &b) in seq.iter().enumerate() {
                        let base = match b {
                            b'A' => 7,
                            b'T' => 222,
                            b'C' => 74,
                            b'G' => 31,
                            _ => 13,
                        };
                        for (i, value) in expected.iter_mut().take(n).enumerate() {
                            *value = value.wrapping_add(
                                duplication.primes[(p * n + i) % (n * 512)]
                                    .wrapping_mul(base + p as u64),
                            );
                        }
                    }
                    assert_eq!(
                        hash.positions(first, Some(second)),
                        expected,
                        "level={level}, length={length}, split={split}"
                    );
                    #[cfg(target_arch = "aarch64")]
                    if let Some(table) = hash.pair.as_deref().filter(|_| length <= 512) {
                        assert_eq!(
                            table.positions_neon(first, Some(second), &hash.primes),
                            expected
                        );
                    }
                    if split == length {
                        assert_eq!(hash.positions(first, None), expected);
                    }
                }
            }
        }
    }

    #[test]
    fn buffered_reader_preserves_lines_at_every_buffer_boundary() {
        let input =
            b"@first description\r\nACGTN\r\n+first\r\n@ABCD\r\n@empty\n\n+\n\n@last\nAC\n+\nHH";
        for capacity in 1..=input.len() + 1 {
            let mut reader = Reader {
                inner: Box::new(std::io::BufReader::with_capacity(
                    capacity,
                    std::io::Cursor::new(input),
                )),
                record: 0,
                phred64: true,
            };
            let first = reader.next().unwrap().unwrap();
            assert_eq!(first.name.as_ref(), b"@first description");
            assert_eq!(first.plus.as_ref(), b"+first");
            assert_eq!(first.seq.as_ref(), b"ACGTN");
            assert_eq!(first.qual.as_ref(), b"!\"#$%");
            let empty = reader.next_reusing(first).unwrap().unwrap();
            assert_eq!(empty.name.as_ref(), b"@empty");
            assert!(empty.seq.is_empty() && empty.qual.is_empty());
            let last = reader.next_reusing(empty).unwrap().unwrap();
            assert_eq!(last.name.as_ref(), b"@last");
            assert_eq!(last.seq.as_ref(), b"AC");
            assert_eq!(last.qual.as_ref(), b"))");
            assert!(reader.next_reusing(last).unwrap().is_none());
        }
    }

    fn buffered(input: &[u8], capacity: usize) -> Reader {
        Reader {
            inner: Box::new(std::io::BufReader::with_capacity(
                capacity,
                std::io::Cursor::new(input.to_vec()),
            )),
            record: 0,
            phred64: false,
        }
    }

    #[test]
    fn reader_skips_blank_lines_between_records_like_fastp() {
        let spaces = " ".repeat(40);
        let input = format!(
            "\n \t\r\n@a\nAC\n+\nII\n\n\n@b\nT\n+\nI\n{spaces}\n\x0b\x0c\r\n@c\nG\n+\nH\n\n\t "
        );
        for capacity in 1..=input.len() + 1 {
            let mut reader = buffered(input.as_bytes(), capacity);
            let mut names = Vec::new();
            while let Some(read) = reader.next().unwrap() {
                names.push(read.name.into_owned());
            }
            assert_eq!(names, [&b"@a"[..], b"@b", b"@c"], "capacity={capacity}");
            assert_eq!(reader.record, 3);
        }
        for input in [&b""[..], b"\n", b"  \n\t\r\n  "] {
            assert!(buffered(input, 1).next().unwrap().is_none());
        }
    }

    #[test]
    fn reader_rejects_blank_lines_inside_records_at_the_right_record() {
        for (input, error) in [
            (
                &b"\n@a\nA\n+\nI\n\n@b\n\nA\n+\nI\n"[..],
                "Invalid FASTQ separator at record 2",
            ),
            (
                b"@a\nA\n+\nI\n@b\nA\n\n+\nI\n",
                "Invalid FASTQ separator at record 2",
            ),
            (
                b"@a\nA\n+\nI\n@b\nA\n+\n\nI\n",
                "Sequence/quality length mismatch at record 2",
            ),
            (
                b"@a\nA\n+\nI\n\n  @b\nA\n+\nI\n",
                "Invalid FASTQ header at record 2",
            ),
            (b"@a\nA\n+\nI\n\n@b\nA\n+\n", "Truncated FASTQ record 2"),
        ] {
            for capacity in 1..=input.len() + 1 {
                let mut reader = buffered(input, capacity);
                assert!(reader.next().unwrap().is_some());
                let actual = reader.next().unwrap_err().to_string();
                assert_eq!(actual, error, "capacity={capacity}, input={input:?}");
            }
        }
    }

    #[test]
    fn packed_offsets_match_each_lane_and_spill_in_input_order() {
        for stored in [[1, 2, 3, 4, 5], [2047, 1024, 1023, 16, 15]] {
            let mut bucket = Bucket::<11>::default();
            for used in 0..=stored.len() {
                for offset in 0..2048 {
                    let mut probe = Bucket::<11>(bucket.0);
                    let expected = if stored[..used].contains(&offset) {
                        Some(true)
                    } else if used < stored.len() || offset == 0 {
                        Some(false)
                    } else {
                        None
                    };
                    assert_eq!(probe.check_and_set(offset), expected);
                }
                if let Some(&next) = stored.get(used) {
                    assert_eq!(bucket.check_and_set(next), Some(false));
                }
            }
            assert_eq!(bucket.check_and_set(0), Some(false));
            assert_eq!(bucket.check_and_set(0), Some(true));
        }

        let mut table = BucketTable::new(2);
        let mut reference = [0u8; 8192];
        for i in 0..200_000u64 {
            let bit = ((i * 40503) ^ (i / 7)) as usize & 65535;
            let mask = 1 << (bit % 8);
            let byte = &mut reference[bit / 8];
            assert_eq!(table.check_and_set(bit as u64), *byte & mask != 0);
            *byte |= mask;
        }
        for bit in 0..65536 {
            assert_eq!(
                table.check_and_set(bit as u64),
                reference[bit / 8] & (1 << (bit % 8)) != 0
            );
        }
    }

    #[test]
    fn adaptive_bloom_page_preserves_every_bit_through_promotion() {
        let mut page = BloomPage::default();
        let mut reference = [0u8; 4096];
        let mut state = 140u64;
        for index in 0..100_000 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let bit = if index < 66 {
                index * 991 % 32768
            } else {
                (state >> 32) as usize % 32768
            };
            let expected = reference[bit / 8] & (1 << (bit % 8)) != 0;
            assert_eq!(
                page.check_and_set(bit as u16),
                expected,
                "insertion {index}, bit {bit}"
            );
            reference[bit / 8] |= 1 << (bit % 8);
            if index == 23 {
                assert!(matches!(page, BloomPage::Small { len: 24, .. }));
            }
            if index == 24 {
                assert!(matches!(page, BloomPage::Offsets(_)));
            }
            if index == 64 {
                assert!(matches!(page, BloomPage::Sparse(_)));
            }
        }
        assert!(matches!(page, BloomPage::Bitmap(_)));
        for bit in 0..32768 {
            assert_eq!(
                page.check_and_set(bit as u16),
                reference[bit / 8] & (1 << (bit % 8)) != 0
            );
        }
    }

    #[test]
    fn packed_duplicate_counts_match_ordered_bit_sets_across_chunks() {
        for level in 1..=6 {
            let mut duplication = Duplication::new(level);
            let width = duplication.buffers.len();
            let mut seen = vec![FxHashSet::default(); width];
            let mut start: usize = 0;
            for n in [0, 1, 15, 16, 17, 33, 257, 4096] {
                let mut packed = Vec::new();
                let mut expected = 0;
                for i in start..start + n {
                    let mut duplicate = true;
                    for (j, set) in seen.iter_mut().enumerate() {
                        // Independent repeats, collisions between the filters,
                        // and high address bits that wrap at lower accuracies.
                        let k = (i / (j + 1)) % (53 + j);
                        let pos = ((k as u64) << 32) | ((k as u64 & 15) << 24);
                        duplicate &= !set.insert(pos & (duplication.bits - 1));
                        packed.push(pos);
                    }
                    expected += usize::from(duplicate);
                }
                assert_eq!(
                    duplication.count_packed(&packed),
                    expected,
                    "level {level}, batch {n}"
                );
                start += n;
            }
        }
    }

    #[test]
    fn sparse_buckets_preserve_offsets_across_groups_and_overflow() {
        let bits = 1u64 << 35;
        let mut buckets = SparseBuckets::new(bits);
        let mut seen = FxHashSet::default();
        // The first four nonzero offsets stay inline. The rest spill into
        // pages shared between buckets, including a dense bitmap in each group.
        for &base in &[0, 1 << 18, 1 << 32, bits - (1 << 18)] {
            for bucket in [0, 1, 3, 4, 31] {
                let start = base + (bucket << 13);
                for offset in (0..1024).chain((0..8192).rev()) {
                    let bit = start + offset;
                    assert_eq!(buckets.check_and_set(bit), !seen.insert(bit), "{bit}");
                }
            }
        }
        for bit in seen {
            assert!(buckets.check_and_set(bit), "{bit}");
        }
    }

    #[test]
    fn bloom_pages_preserve_bit_addresses_and_order() {
        for level in [1, 3, 6] {
            let mut duplication = Duplication::new(level);
            for bit in [0, 1, 7, 8, 32767, 32768, duplication.bits - 1] {
                let positions = [bit; 8];
                assert!(
                    !duplication.check_positions(&positions),
                    "level={level}, bit={bit}"
                );
                assert!(duplication.check_positions(&positions));
            }
            // Addresses wrap modulo the same upstream bit count.
            assert!(duplication.check_positions(&[duplication.bits; 8]));
        }
    }

    const RECORD: &[u8] =
        b"@r\nACGTACGTACGTACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    /// Run `f` on an engine over `fastq`, read in small batches so a run takes
    /// many steps. Natively the input is a regular file, which is read ahead.
    fn with_engine(fastq: Vec<u8>, f: impl FnOnce(Engine)) {
        use clap::Parser;
        let engine = |input: &str| {
            let args = Args::parse_from([
                "fasterp",
                "-i",
                input,
                "--disable_adapter_trimming",
                "--batch_bytes",
                "4096",
            ]);
            Engine::new(args).unwrap()
        };
        #[cfg(feature = "native")]
        {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("in.fq");
            std::fs::write(&path, fastq).unwrap();
            f(engine(path.to_str().unwrap()));
        }
        #[cfg(not(feature = "native"))]
        io::with_files(vec![("in.fq".into(), fastq)], || {
            f(engine("in.fq"));
            Ok(())
        })
        .0
        .unwrap();
    }

    /// Dropping an engine mid-run stops and joins the read-ahead threads,
    /// which then release the shared state.
    #[test]
    fn dropping_a_running_engine_joins_its_threads() {
        // About 2 MB, so one step covers well under 1% of the input.
        with_engine(RECORD.repeat(32_000), |mut engine| {
            assert_eq!(
                matches!(engine.source, Source::Pipelined(_)),
                cfg!(feature = "native")
            );
            assert!(engine.step().unwrap());
            assert!((1..1000).contains(&engine.stats.count));
            let shared = Arc::clone(&engine.shared);
            drop(engine);
            assert_eq!(Arc::strong_count(&shared), 1);
        });
    }

    /// A failed step reports its error once; the engine then stops, refuses
    /// to report partial counts, and drops cleanly.
    #[test]
    fn a_failed_step_stops_the_engine() {
        let mut fastq = RECORD.repeat(1000);
        fastq.extend_from_slice(b"r\nACGT\n+\nIIII\n");
        with_engine(fastq, |mut engine| {
            let mut steps = 0;
            let error = loop {
                match engine.step() {
                    Ok(more) => assert!(more, "input ended without its error"),
                    Err(error) => break error,
                }
                steps += 1;
            };
            assert!(steps > 1);
            assert!(
                format!("{error:#}").contains("Invalid FASTQ header at record 1001"),
                "{error:#}"
            );
            assert!(!engine.step().unwrap());
            let shared = Arc::clone(&engine.shared);
            assert!(engine.finish().is_err());
            assert_eq!(Arc::strong_count(&shared), 1);
        });
    }
}
