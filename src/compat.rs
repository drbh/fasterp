//! Native fastp 1.4.0 processing semantics.
//!
//! Batches are processed in bounded waves. Results are emitted in input order;
//! an arbitrarily slow first batch cannot accumulate an unbounded reorder queue.
use crate::{adapter, auto_adapter, html, io, kmer, overlap, simd, split, stats, Args};
use anyhow::{ensure, Context, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde_json::{json, Value};
use std::borrow::Cow;
use std::collections::{BTreeMap, VecDeque};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

mod block;

pub const FASTP_VERSION: &str = "1.4.0";

#[derive(Clone, Debug, Default)]
struct Read<'a> {
    name: Cow<'a, [u8]>,
    seq: Cow<'a, [u8]>,
    plus: Cow<'a, [u8]>,
    qual: Cow<'a, [u8]>,
    original: Option<Box<(Cow<'a, [u8]>, Cow<'a, [u8]>)>>,
    track_changes: bool,
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
            self.qual.to_mut().iter_mut().for_each(|q| *q -= 31);
        }
        Ok(())
    }

    fn into_owned(self) -> Read<'static> {
        Read {
            name: Cow::Owned(self.name.into_owned()),
            seq: Cow::Owned(self.seq.into_owned()),
            plus: Cow::Owned(self.plus.into_owned()),
            qual: Cow::Owned(self.qual.into_owned()),
            original: self.original.map(|b| {
                let (seq, qual) = *b;
                Box::new((Cow::Owned(seq.into_owned()), Cow::Owned(qual.into_owned())))
            }),
            track_changes: self.track_changes,
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
            self.save_original();
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
            self.save_original();
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

    fn save_original(&mut self) {
        if self.track_changes && self.original.is_none() {
            self.original = Some(Box::new((self.seq.clone(), self.qual.clone())));
        }
    }

    fn write(&self, writer: &mut impl Write, tag: Option<&str>) -> Result<()> {
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

    fn read_into(&mut self, read: &mut Read) -> Result<bool> {
        read.original = None;
        read.track_changes = false;
        let mut lines = [
            read.name.to_mut(),
            read.seq.to_mut(),
            read.plus.to_mut(),
            read.qual.to_mut(),
        ];
        // Ordinary records fit in the existing input buffer. Find all four
        // boundaries with one buffer borrow instead of four virtual reads.
        let buffer = self.inner.fill_buf()?;
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

// A reduction lets the compiler validate a vector of bytes at a time, without
// a branch per quality value. Wrapping subtraction also rejects bytes < offset.
fn valid_quality(qual: &[u8], offset: u8) -> bool {
    qual.iter().fold(0u8, |bad, &q| {
        bad | u8::from(q.wrapping_sub(offset) > 126 - offset)
    }) == 0
}

// Keep the counters touched for one cycle together. The reporting structure
// stores each counter in a separate vector, which is expensive in this hot loop.
#[derive(Clone, Default)]
struct Cycle {
    total_sum: u64,
    total_cnt: u64,
    base_sum: [u64; 4],
    base_cnt: [u64; 4],
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
            dst.total_sum += src.total_sum;
            dst.total_cnt += src.total_cnt;
            for base in 0..4 {
                dst.base_sum[base] += src.base_sum[base];
                dst.base_cnt[base] += src.base_cnt[base];
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
        let mut kmer_code = 0usize;
        let mut valid = 0usize;
        for (index, ((&b, &q), cycle)) in seq.iter().zip(qual).zip(&mut self.cycles).enumerate() {
            let q = usize::from(q - 33);
            cycle.total_sum += q as u64;
            cycle.total_cnt += 1;
            self.qual_hist[index & 3][q] += 1;
            if let Some(base) = kmer::base_idx(b) {
                cycle.base_sum[base] += q as u64;
                cycle.base_cnt[base] += 1;
                if !no_kmer {
                    kmer_code = ((kmer_code << 2) | [0, 3, 1, 2][base]) & 1023;
                    valid += 1;
                    if valid >= 5 {
                        self.kmers[kmer_code] += 1;
                    }
                }
            } else {
                kmer_code = 0;
                valid = 0;
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
            .map(|c| (c.base_cnt[2] + c.base_cnt[3]) as usize)
            .sum();
    }

    fn details(&self) -> Value {
        let mut p = stats::PositionStats::new();
        p.ensure_capacity(self.cycles.len());
        p.qual_hist = std::array::from_fn(|q| self.qual_hist.iter().map(|h| h[q]).sum());
        for (i, cycle) in self.cycles.iter().enumerate() {
            p.total_sum[i] = cycle.total_sum;
            p.total_cnt[i] = cycle.total_cnt;
            for base in 0..4 {
                p.base_sum[base][i] = cycle.base_sum[base];
                p.base_cnt[base][i] = cycle.base_cnt[base];
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
    fn candidates<'a>(
        reads: impl Iterator<Item = &'a Read<'a>>,
        len: usize,
    ) -> BTreeMap<Vec<u8>, usize> {
        let mut counts = BTreeMap::new();
        let mut bases = 0;
        let steps = [10, 20, 40, 100, 150.min(len.saturating_sub(2)).max(1)];
        for r in reads {
            if bases >= 151 * 10000 {
                break;
            }
            bases += r.seq.len();
            for step in steps {
                for i in 0..r.seq.len().saturating_sub(step) {
                    *counts.entry(r.seq[i..i + step].to_vec()).or_insert(0usize) += 1;
                }
            }
        }
        counts.retain(|seq, count| {
            *count
                >= if seq.len() >= len.saturating_sub(1) {
                    3
                } else if seq.len() >= 100 {
                    5
                } else if seq.len() >= 40 {
                    20
                } else if seq.len() >= 20 {
                    100
                } else if seq.len() >= 10 {
                    500
                } else {
                    usize::MAX
                }
        });
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
            steps: [10, 20, 40, 100, 150.min(len.saturating_sub(2)).max(1)],
            sampling,
        }
    }
    fn add(&mut self, seq: &[u8]) {
        for step in self.steps {
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
    r.name.to_mut().splice(pos..pos, tag);
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
    Dense(Vec<BloomPage>),
    Sparse(FxHashMap<u64, BloomPage>),
}

// Most 4 KiB pages contain only a few set bits on ordinary FASTQ inputs.
// Store those exact bit offsets inline, promoting a populated page to the
// original bitmap. This changes storage, never hashes or membership decisions.
enum BloomPage {
    Small {
        len: u8,
        bits: [u16; 32],
        summary: u64,
    },
    Bitmap(Box<[u8; 4096]>),
}

impl Default for BloomPage {
    fn default() -> Self {
        Self::Small {
            len: 0,
            bits: [0; 32],
            summary: 0,
        }
    }
}

impl BloomPage {
    fn check_and_set(&mut self, bit: u16) -> bool {
        match self {
            Self::Small { len, bits, summary } => {
                let mask = 1u64 << (bit % 64);
                // An unset summary bit proves absence. A hit still checks the
                // exact offsets, so this shortcut cannot create false positives.
                if *summary & mask != 0 && bits[..usize::from(*len)].contains(&bit) {
                    return true;
                }
                if usize::from(*len) < bits.len() {
                    bits[usize::from(*len)] = bit;
                    *len += 1;
                    *summary |= mask;
                } else {
                    let mut bitmap = Box::new([0; 4096]);
                    for &offset in bits.iter().chain(std::iter::once(&bit)) {
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

impl BloomPages {
    fn page(&mut self, index: u64) -> &mut BloomPage {
        match self {
            Self::Dense(pages) => &mut pages[index as usize],
            Self::Sparse(pages) => pages.entry(index).or_default(),
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
        // Direct lookup for ordinary accuracy levels, without committing the
        // Bloom data itself. Keep high-accuracy directories sparse and bounded.
        let dense = pages * n * std::mem::size_of::<BloomPage>() <= 32 * 1024 * 1024;
        Self {
            buffers: (0..n)
                .map(|_| {
                    if dense {
                        BloomPages::Dense(
                            std::iter::repeat_with(BloomPage::default)
                                .take(pages)
                                .collect(),
                        )
                    } else {
                        BloomPages::Sparse(FxHashMap::default())
                    }
                })
                .collect(),
            primes,
            bits,
        }
    }
    fn hasher(&self) -> DuplicationHash {
        #[cfg(target_arch = "aarch64")]
        let (pair_primes, pair_bias) = if self.buffers.len() == 2 {
            // For two hashes every prime * nucleotide code fits in u32. The
            // position-only term is independent of the read and is added once.
            let primes: Vec<u32> = self.primes.iter().map(|&p| p as u32).collect();
            assert!(primes.iter().all(|&p| p <= u32::MAX / 222));
            let mut bias = vec![[0u64; 8]; 513];
            for pos in 0..512 {
                bias[pos + 1] = bias[pos];
                for i in 0..2 {
                    bias[pos + 1][i] += self.primes[pos * 2 + i] * pos as u64;
                }
            }
            (primes, bias)
        } else {
            (Vec::new(), Vec::new())
        };
        DuplicationHash {
            buffers: self.buffers.len(),
            primes: self.primes.clone(),
            #[cfg(target_arch = "aarch64")]
            pair_primes,
            #[cfg(target_arch = "aarch64")]
            pair_bias,
        }
    }

    fn check_positions(&mut self, positions: [u64; 8]) -> bool {
        let mut duplicate = true;
        for (buf, pos) in self.buffers.iter_mut().zip(positions) {
            let bit = pos % self.bits;
            let page = buf.page(bit / 32768);
            duplicate &= page.check_and_set((bit % 32768) as u16);
        }
        duplicate
    }
}

struct DuplicationHash {
    buffers: usize,
    primes: Vec<u64>,
    #[cfg(target_arch = "aarch64")]
    pair_primes: Vec<u32>,
    #[cfg(target_arch = "aarch64")]
    pair_bias: Vec<[u64; 8]>,
}

impl DuplicationHash {
    fn positions(&self, r1: &Read, r2: Option<&Read>) -> [u64; 8] {
        #[cfg(target_arch = "aarch64")]
        if self.buffers == 2 && r1.seq.len() + r2.map_or(0, |r| r.seq.len()) <= 512 {
            // NEON is part of the aarch64 baseline. Longer reads retain the
            // wrapping scalar algorithm, including its repeated prime table.
            return unsafe { self.positions_two_neon(r1, r2) };
        }
        match self.buffers {
            2 => self.positions_with::<2>(r1, r2),
            4 => self.positions_with::<4>(r1, r2),
            8 => self.positions_with::<8>(r1, r2),
            _ => unreachable!("validated duplication accuracy"),
        }
    }

    fn positions_with<const N: usize>(&self, r1: &Read, r2: Option<&Read>) -> [u64; 8] {
        let mut positions = [0u64; 8];
        let mut offset = 0;
        for r in std::iter::once(r1).chain(r2) {
            for (p, &b) in r.seq.iter().enumerate() {
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
            offset += r.seq.len();
        }
        positions
    }

    #[cfg(target_arch = "aarch64")]
    unsafe fn positions_two_neon(&self, r1: &Read, r2: Option<&Read>) -> [u64; 8] {
        use std::arch::aarch64::*;
        // Each load is bounded by a full 16-byte sequence chunk; the prime
        // table covers both reads combined (at most 512 bases). Widening sums
        // preserve the exact u64 hash rather than wrapping intermediate u32s.
        unsafe {
            let mut sums = [vdupq_n_u64(0); 2];
            let mut positions = self.pair_bias[r1.seq.len() + r2.map_or(0, |r| r.seq.len())];
            let mut offset = 0;
            for read in std::iter::once(r1).chain(r2) {
                let mut p = 0;
                while p + 16 <= read.seq.len() {
                    let bytes = vld1q_u8(read.seq.as_ptr().add(p));
                    let mut bases = vdupq_n_u8(13);
                    for (letter, value) in [(b'A', 7), (b'T', 222), (b'C', 74), (b'G', 31)] {
                        bases = vbslq_u8(
                            vceqq_u8(bytes, vdupq_n_u8(letter)),
                            vdupq_n_u8(value),
                            bases,
                        );
                    }
                    let low = vmovl_u8(vget_low_u8(bases));
                    let high = vmovl_high_u8(bases);
                    let groups = [
                        vmovl_u16(vget_low_u16(low)),
                        vmovl_high_u16(low),
                        vmovl_u16(vget_low_u16(high)),
                        vmovl_high_u16(high),
                    ];
                    for (group, codes) in groups.into_iter().enumerate() {
                        let primes =
                            vld2q_u32(self.pair_primes.as_ptr().add((offset + p + group * 4) * 2));
                        sums[0] = vpadalq_u32(sums[0], vmulq_u32(codes, primes.0));
                        sums[1] = vpadalq_u32(sums[1], vmulq_u32(codes, primes.1));
                    }
                    p += 16;
                }
                for (p, &base) in read.seq.iter().enumerate().skip(p) {
                    let base = match base {
                        b'A' => 7,
                        b'T' => 222,
                        b'C' => 74,
                        b'G' => 31,
                        _ => 13,
                    };
                    for i in 0..2 {
                        positions[i] += u64::from(self.pair_primes[(offset + p) * 2 + i]) * base;
                    }
                }
                offset += read.seq.len();
            }
            for i in 0..2 {
                positions[i] = positions[i].wrapping_add(vaddvq_u64(sums[i]));
            }
            positions
        }
    }
}

struct Config {
    args: Args,
    adapters: [Option<Vec<u8>>; 2],
    fasta: Vec<Vec<u8>>,
    blacklist: [Vec<Vec<u8>>; 2],
    poly_g: bool,
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
    positions: Vec<[u64; 8]>,
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

impl PackedChunk {
    fn new(records: usize, recycle_limit: usize, duplication: bool) -> Self {
        Self {
            output: std::array::from_fn(|_| Vec::new()),
            positions: Vec::with_capacity(if duplication { records } else { 0 }),
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

    fn add(&mut self, mut p: Processed, a: &Args, duplication: bool) {
        self.count += 1;
        if duplication {
            self.positions.push(p.duplicate_positions);
        }
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
                if p.failures[0] == Failure::Pass {
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
                        if p.failures[i] == Failure::Pass && (a.stdout || a.merged_out.is_some()) {
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
                if failure == Failure::Pass {
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
                let mut pair = (first, p.reads[1].take());
                for read in std::iter::once(&mut pair.0).chain(pair.1.iter_mut()) {
                    read.original = None;
                    read.track_changes = false;
                }
                self.recycled
                    .put((pair.0.into_owned(), pair.1.map(Read::into_owned)));
            }
        }
    }
}

fn pack_main_output(a: &Args) -> bool {
    !a.dedup
        && !a.overrepresentation_analysis
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
                for &positions in &packed.positions {
                    packed.duplicates += usize::from(duplication.check_positions(positions));
                }
            } else {
                for p in &mut chunk.processed {
                    p.duplicate = duplication.check_positions(p.duplicate_positions);
                }
            }
        }
    }
}

fn trim_adapter(r: &mut Read, seq: &[u8], min: usize) -> Option<usize> {
    let before = r.seq.len();
    let found = adapter::find_adapter(&r.seq, seq, min, usize::MAX)?;
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
                r.name.to_mut().insert(pos, b' ');
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
            r1.save_original();
            r2.save_original();
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
    raw: Option<[block::Records; 2]>,
    reads: Vec<(Read<'static>, Option<Read<'static>>)>,
    start: usize,
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
        if let Some(mut batch) = self.sampled_raw.pop_front() {
            let raw = batch.raw.as_mut().unwrap();
            if self.limit != 0 {
                for mate in raw.iter_mut() {
                    mate.truncate(self.limit - self.count);
                }
            }
            batch.start = self.count;
            self.count += raw[0].len();
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
            self.count += raw[0].len() / if self.interleaved { 2 } else { 1 };
            return Ok(Some(ReadBatch {
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
    recycle_limit: usize,
) -> Result<Chunk<'a>> {
    let a = &config.args;
    let profile_after_in_worker = !a.dedup && !a.overrepresentation_analysis;
    let packed_output = pack_main_output(a);
    let shared_stats = packed_output && !a.merge;
    let mut shared: [Sample; 2] = std::array::from_fn(|_| Sample::default());
    let mut chunk = Chunk {
        before: std::array::from_fn(|i| before[i].empty_like()),
        after: std::array::from_fn(|i| after[i].empty_like()),
        after_offset: [0; 2],
        processed: Vec::with_capacity(if packed_output { 0 } else { records }),
        packed: packed_output
            .then(|| PackedChunk::new(records, recycle_limit, duplication.is_some())),
    };
    let mut overlap_scratch = Vec::new();
    for (index, pair) in pairs.enumerate() {
        let (mut r1, mut r2) = pair?;
        let ordinal = start + index;
        if shared_stats {
            r1.track_changes = true;
            if let Some(r2) = &mut r2 {
                r2.track_changes = true;
            }
        } else {
            chunk.before[0].add_at(&r1, a.no_kmer, ordinal);
            if let Some(r2) = &r2 {
                chunk.before[1].add_at(r2, a.no_kmer, ordinal);
            }
        }
        let positions = duplication
            .as_ref()
            .map(|d| d.positions(&r1, r2.as_ref()))
            .unwrap_or([0; 8]);
        let mut p = process(r1, r2, &config, &mut overlap_scratch);
        p.duplicate_positions = positions;
        if shared_stats {
            let passing = !p.skipped && p.failures.iter().all(|f| *f == Failure::Pass);
            for (i, r) in p
                .reads
                .iter()
                .enumerate()
                .filter_map(|(i, r)| r.as_ref().map(|r| (i, r)))
            {
                if passing && r.original.is_none() {
                    // The same bases and qualities contribute to both
                    // reports. Count them once, then merge the totals.
                    shared[i].add_at(r, a.no_kmer, ordinal);
                } else {
                    if let Some(original) = &r.original {
                        chunk.before[i].add_parts(&original.0, &original.1, a.no_kmer, ordinal);
                    } else {
                        chunk.before[i].add_at(r, a.no_kmer, ordinal);
                    }
                    if passing {
                        chunk.after[i].add_at(r, a.no_kmer, 0);
                    }
                }
            }
        } else if profile_after_in_worker {
            for (i, r) in p.passing_reads(a).into_iter().flatten() {
                chunk.after[i].add_at(r, a.no_kmer, 0);
            }
        }
        if let Some(packed) = &mut chunk.packed {
            packed.add(p, a, duplication.is_some());
        } else {
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
    let ReadBatch { raw, reads, start } = batch;
    let threads = pool.current_num_threads();
    let interleaved = config.args.interleaved_in;
    let stride = if interleaved { 2 } else { 1 };
    let count = raw.as_ref().map_or(reads.len(), |r| r[0].len() / stride);
    let chunk_size = count.div_ceil(threads).max(1);
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
            writers: BTreeMap::new(),
            scratch: Vec::new(),
            compression_threads: crate::gzip::threads_per_output(gzip_outputs),
        }
    }
    fn open(&mut self, path: &str, a: &Args, config: split::SplitConfig) -> Result<()> {
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
        if let Some(writer) = path.and_then(|p| self.writers.get_mut(p)) {
            self.scratch.clear();
            read.write(&mut self.scratch, tag)?;
            writer.write_all(&self.scratch)?;
        }
        Ok(())
    }
    fn write_packed(&mut self, path: Option<&str>, bytes: &[u8]) -> Result<()> {
        if let Some(writer) = path.and_then(|p| self.writers.get_mut(p)) {
            writer.write_all(bytes)?;
        }
        Ok(())
    }
    fn finish(self) -> Result<()> {
        for (_, w) in self.writers {
            w.finish()?;
        }
        Ok(())
    }
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
        (1..=16).contains(&a.threads.unwrap_or(3)),
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
    let text = std::fs::read_to_string(path).with_context(|| format!("Cannot read {path}"))?;
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

pub fn run(mut a: Args) -> Result<()> {
    let start = Instant::now();
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
    let pipelined = a.max_backlog != Some(1)
        && std::iter::once(a.input.as_str())
            .chain(a.input2.as_deref())
            .all(|path| std::fs::metadata(path).is_ok_and(|m| m.is_file()));
    // Borrowed parsing also works synchronously. Only reading ahead requires
    // regular files; pipes must stay on this thread for prompt cancellation.
    // Interleaved stdin has no evaluation prefix. Keep local interleaved
    // sampling on the owned path until that separate path is generalized.
    let raw_input = pack_main_output(&a) && (!a.interleaved_in || a.input == "-");
    // Read once and retain the sample: this works for pipes, HTTP, and S3 too.
    let auto = |seq: &Option<String>| seq.as_deref().is_none_or(|s| s == "auto");
    let needs_adapter_sample = !a.disable_adapter_trimming
        && !a.disable_adapter_detection
        && if pe {
            a.detect_adapter_for_pe && (auto(&a.adapter_sequence) || auto(&a.adapter_sequence_r2))
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
            let Some(blocks) = block::read_sample_batch(
                &mut reader1,
                reader2.as_mut(),
                a.batch_bytes.min(2 * 1024 * 1024),
                sample_limit - sample_count,
                151 * 256 * 1024 - sample_bases,
            )?
            else {
                break;
            };
            for index in 0..blocks[0].len() {
                // Preserve validation of the entire evaluation sample, including
                // when --reads_to_process subsequently processes only its prefix.
                for mate in &blocks {
                    if mate.len() > 0 {
                        mate.read(index, sample_count + index, a.phred64)?;
                    }
                }
                sample_bases += blocks[0].sequence(index).len();
            }
            sample_count += blocks[0].len();
            sampled_raw.push_back(ReadBatch {
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
                .iter()
                .flat_map(|b| {
                    let records = &b.raw.as_ref().unwrap()[mate];
                    (0..records.len()).map(|i| records.sequence(i))
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
    if !a.disable_adapter_trimming && !a.disable_adapter_detection && sample_count >= 10000 {
        if !pe && adapters[0].is_none() {
            adapters[0] = auto_adapter::detect(&sample_sequences(0), a.trim_tail1);
        } else if pe && a.detect_adapter_for_pe {
            if adapters[0].is_none() {
                adapters[0] = auto_adapter::detect(&sample_sequences(0), a.trim_tail1);
            }
            if adapters[1].is_none() {
                adapters[1] = auto_adapter::detect(&sample_sequences(1), a.trim_tail1);
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
    let threads = a.threads.unwrap_or(3);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    let max_bytes = if let Some(mb) = a.max_memory {
        ensure!(mb >= 16, "--max_memory must be at least 16 MB");
        a.batch_bytes.min(mb.saturating_mul(1024 * 1024) / 8).max(1)
    } else {
        a.batch_bytes
    };
    let config = Config {
        args: a,
        adapters,
        fasta,
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
    for path in [Some(a.output.as_str()), a.output2.as_deref()]
        .into_iter()
        .flatten()
    {
        out.open(path, a, split.clone())?;
    }
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
    let main1 = if (a.stdout || a.output == "-") && !a.merge {
        Some("-")
    } else {
        Some(a.output.as_str())
    };
    let main2 = if (a.stdout || a.output == "-") && !a.merge {
        Some("-")
    } else {
        a.output2.as_deref()
    };
    let merged_out = if a.stdout {
        Some("-")
    } else {
        a.merged_out.as_deref()
    };
    let mut before = [Sample::default(), Sample::default()];
    let mut after = [Sample::default(), Sample::default()];
    let mut failures = [0usize; 7];
    let mut adapter_reads = 0;
    let mut adapter_bases = 0;
    let mut corrected_reads = 0;
    let mut corrected_bases = 0;
    let mut poly_reads = [0usize; 4];
    let mut poly_bases = [0usize; 4];
    let mut histogram = vec![0usize; 513];
    let mut unknown = 0;
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
    let mut duplicates = 0;
    let mut count = 0;
    let mut merged_pairs = 0;
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
            let candidates = Overrep::candidates(reads.iter().copied(), len);
            before[i].overrep = Some(Overrep::new(
                &candidates,
                len,
                a.overrepresentation_sampling,
            ));
            after[i].overrep = Some(Overrep::new(
                &candidates,
                len,
                a.overrepresentation_sampling,
            ));
        }
    }
    // Only regular local files are read ahead. Pipes and remote streams retain
    // synchronous reads so a failed output cannot leave a reader blocked while
    // scoped threads are being joined.
    // Compact synchronous jobs use the same batch ceiling as the pipeline,
    // limiting serial parsing work before each parallel processing step.
    let max_bytes = if pipelined || pack_main_output(a) {
        max_bytes.min(2 * 1024 * 1024)
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
    let before_templates = std::array::from_fn(|i| before[i].empty_like());
    let after_templates = std::array::from_fn(|i| after[i].empty_like());
    let hasher = duplication.as_ref().map(Duplication::hasher);
    std::thread::scope(|scope| -> Result<()> {
        let mut serial_reader = Some(batch_reader);
        let (recycle_tx, recycle_rx) = crossbeam_channel::bounded::<Vec<ReadPool>>(2);
        let processed_rx = if pipelined {
            // One queued batch per stage bounds both memory and reordering.
            let (input_tx, input_rx) = crossbeam_channel::bounded(1);
            let (processed_tx, processed_rx) = crossbeam_channel::bounded(1);
            let mut reader = serial_reader.take().unwrap();
            scope.spawn(move || loop {
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
            });
            let pool = &pool;
            let config = &config;
            let before = &before_templates;
            let after = &after_templates;
            let hasher = hasher.as_ref();
            let mut worker_duplication = duplication.take();
            let processed_tx = if worker_duplication.is_some() {
                let (duplicate_tx, duplicate_rx) =
                    crossbeam_channel::bounded::<Result<Vec<Chunk<'static>>>>(1);
                scope.spawn(move || {
                    for result in duplicate_rx {
                        let result = result.map(|mut chunks| {
                            mark_duplicates(&mut chunks, worker_duplication.as_mut());
                            chunks
                        });
                        let failed = result.is_err();
                        if processed_tx.send(result).is_err() || failed {
                            break;
                        }
                    }
                });
                duplicate_tx
            } else {
                processed_tx
            };
            scope.spawn(move || {
                for batch in input_rx {
                    let result = batch.and_then(|b| {
                        process_batch(b, pool, config, before, after, hasher, recycle_limit)
                    });
                    let failed = result.is_err();
                    if processed_tx.send(result).is_err() || failed {
                        break;
                    }
                }
            });
            Some(processed_rx)
        } else {
            None
        };
        loop {
            let mut chunks = if let Some(receiver) = &processed_rx {
                match receiver.recv() {
                    Ok(result) => result?,
                    Err(_) => break,
                }
            } else {
                match serial_reader.as_mut().unwrap().next_batch()? {
                    Some(batch) => process_batch(
                        batch,
                        &pool,
                        &config,
                        &before_templates,
                        &after_templates,
                        hasher.as_ref(),
                        recycle_limit,
                    )?,
                    None => break,
                }
            };
            if processed_rx.is_none() {
                mark_duplicates(&mut chunks, duplication.as_mut());
            }
            let mut read_pool = ReadPool {
                reads: Vec::with_capacity(chunks.iter().map(|c| c.processed.len()).sum()),
                limit: max_bytes,
                ..ReadPool::default()
            };
            let mut read_pools = Vec::with_capacity(chunks.len() + 1);
            let profile_after_in_worker = !a.dedup && !a.overrepresentation_analysis;
            // Bloom decisions remain ordered. Determine global after-filter indices
            // before parallel profiling so sampling is independent of chunk boundaries.
            let mut offsets = std::array::from_fn(|i| after[i].simple.total_reads);
            for chunk in &mut chunks {
                chunk.after_offset = offsets;
                for p in &mut chunk.processed {
                    duplicates += usize::from(p.duplicate);
                    if !profile_after_in_worker {
                        for (i, _) in p.passing_reads(a).into_iter().flatten() {
                            offsets[i] += 1;
                        }
                    }
                }
            }
            if !profile_after_in_worker {
                pool.install(|| {
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
                for (dst, src) in before.iter_mut().zip(chunk.before) {
                    dst.merge(src);
                }
                for (dst, src) in after.iter_mut().zip(chunk.after) {
                    dst.merge(src);
                }
                if let Some(packed) = chunk.packed {
                    out.write_packed(main1, &packed.output[0])?;
                    out.write_packed(main2, &packed.output[1])?;
                    out.write_packed(merged_out, &packed.output[2])?;
                    count += packed.count;
                    merged_pairs += packed.merged_pairs;
                    duplicates += packed.duplicates;
                    adapter_reads += packed.adapter_reads;
                    adapter_bases += packed.adapter_bases;
                    corrected_reads += packed.corrected_reads;
                    corrected_bases += packed.corrected_bases;
                    unknown += packed.unknown;
                    for (dst, src) in failures.iter_mut().zip(packed.failures) {
                        *dst += src;
                    }
                    for (dst, src) in histogram.iter_mut().zip(packed.histogram) {
                        *dst += src;
                    }
                    for (dst, src) in poly_reads.iter_mut().zip(packed.poly_reads) {
                        *dst += src;
                    }
                    for (dst, src) in poly_bases.iter_mut().zip(packed.poly_bases) {
                        *dst += src;
                    }
                    read_pools.push(packed.recycled);
                    continue;
                }
                for p in chunk.processed {
                    let p = RecycleRead {
                        processed: p,
                        pool: &mut read_pool,
                    };
                    count += 1;
                    let duplicate = p.duplicate;
                    if p.skipped {
                        continue;
                    }
                    corrected_reads += p.corrected_reads;
                    corrected_bases += p.corrected_bases;
                    adapter_reads += p.adapter_reads;
                    adapter_bases += p.adapter_bases;
                    for &(base, n) in &p.poly_x {
                        poly_reads[base] += 1;
                        poly_bases[base] += n;
                    }
                    if pe {
                        if let Some(size) = p.insert {
                            histogram[size.min(512)] += 1;
                        } else {
                            unknown += 1;
                        }
                    }
                    if let Some(r) = &p.overlapped {
                        out.write(a.overlapped_out.as_deref(), r, None)?;
                    }
                    if let Some(r) = &p.merged {
                        failures[p.failures[0] as usize] += 2;
                        if p.failures[0] == Failure::Pass && !(a.dedup && duplicate) {
                            out.write(merged_out, r, None)?;
                            merged_pairs += 1;
                        }
                        continue;
                    }
                    if a.merge && a.include_unmerged {
                        for i in 0..2 {
                            if let Some(r) = &p.reads[i] {
                                failures[p.failures[i] as usize] += 1;
                                if p.failures[i] == Failure::Pass && !(a.dedup && duplicate) {
                                    out.write(merged_out, r, None)?;
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
                    failures[f as usize] += if pe { 2 } else { 1 };
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
                    "Processed {count} {}",
                    if pe { "read pairs" } else { "reads" }
                );
            }
            read_pools.push(read_pool);
            if processed_rx.is_some() {
                let _ = recycle_tx.try_send(read_pools);
            } else {
                for read_pool in read_pools {
                    serial_reader.as_mut().unwrap().pool.absorb(read_pool);
                }
            }
        }
        Ok(())
    })?;
    out.finish()?;
    for sample in before.iter_mut().chain(&mut after) {
        sample.finish_summary();
    }
    let mut report = make_report(&before, &after, pe, a.merge);
    report["filtering_result"] = json!({"passed_filter_reads":failures[0],"low_quality_reads":failures[4],"too_many_N_reads":failures[1],"too_short_reads":failures[2],"too_long_reads":failures[3],"low_complexity_reads":failures[5],"adapter_dimer_reads":failures[6]});
    if a.correction {
        report["filtering_result"]["corrected_reads"] = json!(corrected_reads);
        report["filtering_result"]["corrected_bases"] = json!(corrected_bases);
    }
    if !a.dont_eval_duplication || a.dedup {
        report["duplication"] =
            json!({"rate":if count>0 {duplicates as f64/count as f64} else {0.0}});
    }
    if !a.disable_adapter_trimming {
        report["adapter_cutting"] =
            json!({"adapter_trimmed_reads":adapter_reads,"adapter_trimmed_bases":adapter_bases});
    }
    if a.trim_poly_x {
        report["polyx_trimming"] = json!({"total_polyx_trimmed_reads":poly_reads.iter().sum::<usize>(),"total_polyx_trimmed_bases":poly_bases.iter().sum::<usize>(),"polyx_trimmed_reads":poly_reads,"polyx_trimmed_bases":poly_bases});
    }
    if pe {
        let peak = (0..512).max_by_key(|&i| (histogram[i], 512 - i)).unwrap();
        report["insert_size"] =
            json!({"peak":peak,"unknown":unknown+histogram[512],"histogram":&histogram[..512]});
    }
    if a.merge {
        report["merged_and_filtered"] = after[0].details();
        report["merged_pairs"] = json!(merged_pairs);
    }
    write_reports(a, &report, start)?;
    eprintln!(
        "Processed {count} {}; passed {} reads. Compatibility target: fastp {FASTP_VERSION}",
        if pe { "pairs" } else { "reads" },
        failures[0]
    );
    Ok(())
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
        let mut writer = BufWriter::new(std::fs::File::create(&a.json)?);
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
    let mut text = std::fs::read_to_string(&a.html)?;
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
    std::fs::write(&a.html, text.replace("</body>", &format!("{extra}</body>")))?;
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
    fn fused_statistics_match_separate_quality_and_kmer_passes() {
        let mut sample = Sample::default();
        let mut simple = stats::SimpleStats::default();
        let mut kmers = [0; 1024];
        let mut qualities = [0; 94];
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
            kmer::count_k5_2bit(&read.seq, &mut kmers);
            for &q in read.qual.iter() {
                qualities[usize::from(q - 33)] += 1;
            }
        }
        sample.finish_summary();
        assert_eq!(
            serde_json::to_value(sample.simple.to_read_stats()).unwrap(),
            serde_json::to_value(simple.to_read_stats()).unwrap()
        );
        assert_eq!(sample.kmers, kmers);
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
                    .map(|i| alphabet[(i * 17 + length) % alphabet.len()])
                    .collect();
                for split in [
                    length,
                    length / 3,
                    length / 2,
                    length.min(15),
                    length.min(17),
                ] {
                    let first = Read {
                        seq: seq[..split].to_vec().into(),
                        ..Read::default()
                    };
                    let second = Read {
                        seq: seq[split..].to_vec().into(),
                        ..Read::default()
                    };
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
                        hash.positions(&first, Some(&second)),
                        expected,
                        "level={level}, length={length}, split={split}"
                    );
                    if split == length {
                        assert_eq!(hash.positions(&first, None), expected);
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

    #[test]
    fn adaptive_bloom_page_preserves_every_bit_through_promotion() {
        let mut page = BloomPage::default();
        let mut reference = [0u8; 4096];
        let mut state = 140u64;
        for index in 0..100_000 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let bit = if index < 33 {
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
            if index == 31 {
                assert!(matches!(page, BloomPage::Small { len: 32, .. }));
            }
            if index == 32 {
                assert!(matches!(page, BloomPage::Bitmap(_)));
            }
        }
        for bit in 0..32768 {
            assert_eq!(
                page.check_and_set(bit as u16),
                reference[bit / 8] & (1 << (bit % 8)) != 0
            );
        }
    }

    #[test]
    fn bloom_pages_preserve_bit_addresses_and_order() {
        for level in [1, 3, 6] {
            let mut duplication = Duplication::new(level);
            for bit in [0, 1, 7, 8, 32767, 32768, duplication.bits - 1] {
                let positions = [bit; 8];
                assert!(
                    !duplication.check_positions(positions),
                    "level={level}, bit={bit}"
                );
                assert!(duplication.check_positions(positions));
            }
            // Addresses wrap modulo the same upstream bit count.
            assert!(duplication.check_positions([duplication.bits; 8]));
        }
    }
}
