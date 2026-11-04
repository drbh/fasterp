//! A FASTQ-aware deflate encoder.
//!
//! Standard deflate output (RFC 1951) that any inflater reads, produced by
//! exploiting FASTQ structure instead of a general-purpose match finder:
//!
//! * sequence lines are Huffman literals: on a four-letter alphabet short
//!   matches cost more bits than they save, and searching for them is where
//!   general compressors spend most of their time. One probe per read finds
//!   a duplicate of an earlier read that starts at the same base;
//! * header lines are matched against the previous header (common prefix
//!   and suffix), with the hash matcher for the part that changes;
//! * `+` lines and quality strings go through a four-way hash matcher
//!   (greedy, minimum match 4), which captures runs and repeated patterns of
//!   binned qualities.
//!
//! Line classification only affects compression efficiency, never validity:
//! every match is verified byte for byte, so any input compresses correctly.

const WINDOW: usize = 32768;
const MIN_MATCH: usize = 4;
const MAX_MATCH: usize = 258;
const NICE_MATCH: usize = 64;
const INSERT_IN_MATCH: usize = 16;

const HASH_BITS: u32 = 12;

/// Quality strings with at most this many distinct characters are binned,
/// as NovaSeq and NextSeq 2000 write them (four values). There this encoder
/// matches libdeflate level 4 in size at about 2.4 times its speed. With eight
/// bins or full-range qualities level 4 finds matches in the qualities that
/// greedy hashing misses, and is 2-8% smaller, so those blocks stay with it.
const MAX_BINS: u32 = 6;
const SAMPLE_RECORDS: usize = 64;

/// Whether the FASTQ-aware encoder should compress `input` instead of
/// libdeflate: the block's first records have binned qualities. The choice
/// depends only on the block, so output is the same on every run and thread
/// count.
pub(crate) fn suits(input: &[u8]) -> bool {
    binned(input)
}

fn binned(data: &[u8]) -> bool {
    let n = data.len();
    let mut pos = 0;
    while !record_starts(data, pos) {
        match memchr::memchr(b'\n', &data[pos.min(n)..]) {
            Some(i) => pos += i + 1,
            None => return false,
        }
    }
    let mut seen = [0u64; 4];
    let mut records = 0;
    while records < SAMPLE_RECORDS && pos < n && data[pos] == b'@' {
        let seq = line_end(data, pos) + 1;
        let plus = line_end(data, seq.min(n)) + 1;
        let qual = line_end(data, plus.min(n)) + 1;
        if qual >= n {
            break;
        }
        let end = line_end(data, qual);
        for &b in &data[qual..end] {
            seen[(b >> 6) as usize] |= 1 << (b & 63);
        }
        records += 1;
        pos = end + 1;
    }
    records > 0 && seen.iter().map(|w| w.count_ones()).sum::<u32>() <= MAX_BINS
}

// ---- static tables --------------------------------------------------------

const LEN_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LEN_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const CL_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

struct Tables {
    /// length (3..=258) -> length symbol (257..=285)
    len_sym: [u16; 259],
    /// distance (1..=256) -> distance symbol
    dist_sym_small: [u8; 257],
    /// (distance - 1) >> 7 for distance > 256 -> distance symbol
    dist_sym_large: [u8; 256],
}

fn tables() -> &'static Tables {
    use std::sync::OnceLock;
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let mut len_sym = [0u16; 259];
        for (i, &base) in LEN_BASE.iter().enumerate() {
            let next = if i + 1 < LEN_BASE.len() {
                LEN_BASE[i + 1] as usize
            } else {
                259
            };
            for sym in &mut len_sym[base as usize..next.min(259)] {
                *sym = 257 + i as u16;
            }
        }
        len_sym[258] = 285;
        let mut dist_sym_small = [0u8; 257];
        let mut dist_sym_large = [0u8; 256];
        for (i, &base) in DIST_BASE.iter().enumerate() {
            let next = if i + 1 < DIST_BASE.len() {
                DIST_BASE[i + 1] as usize
            } else {
                32769
            };
            for d in base as usize..next {
                if d <= 256 {
                    dist_sym_small[d] = i as u8;
                } else {
                    dist_sym_large[(d - 1) >> 7] = i as u8;
                }
            }
        }
        Tables {
            len_sym,
            dist_sym_small,
            dist_sym_large,
        }
    })
}

#[inline(always)]
fn dist_symbol(t: &Tables, d: usize) -> usize {
    if d <= 256 {
        t.dist_sym_small[d] as usize
    } else {
        t.dist_sym_large[(d - 1) >> 7] as usize
    }
}

// ---- tokens ---------------------------------------------------------------

/// A literal run (its length, bit 31 clear) or a match: bit 31 set, length
/// in bits 16..25, distance - 1 in bits 0..15. Tokens consume the input in
/// order, so the writer walks the input alongside them.
type Token = u32;
const MATCH: u32 = 1 << 31;

#[inline(always)]
fn token_match(len: usize, dist: usize) -> Token {
    MATCH | ((len as u32) << 16) | (dist as u32 - 1)
}

struct Tokens {
    buf: Vec<Token>,
    /// Sized to a power of two so masked indexes need no bounds check.
    lit_freq: [u32; 512],
    /// Literal counts from `literals`, added into `lit_freq` at the end.
    lit4: [[u32; 256]; 4],
    dist_freq: [u32; 32],
    /// Length of the literal run being accumulated.
    run: u32,
}

impl Tokens {
    #[inline(always)]
    fn literal(&mut self, b: u8) {
        self.lit_freq[b as usize] += 1;
        self.run += 1;
    }
    #[inline(always)]
    fn literals(&mut self, bytes: &[u8]) {
        // Four tables: sequence lines repeat four bytes, and counting into
        // one table would make each increment wait for the last.
        let mut quads = bytes.chunks_exact(4);
        for q in &mut quads {
            self.lit4[0][q[0] as usize] += 1;
            self.lit4[1][q[1] as usize] += 1;
            self.lit4[2][q[2] as usize] += 1;
            self.lit4[3][q[3] as usize] += 1;
        }
        for &b in quads.remainder() {
            self.lit_freq[b as usize] += 1;
        }
        self.run += bytes.len() as u32;
    }
    #[inline(always)]
    fn end_run(&mut self) {
        if self.run > 0 {
            self.buf.push(self.run);
            self.run = 0;
        }
    }
    /// Emits a match of any length >= 3, split into deflate matches.
    #[inline(always)]
    fn matched(&mut self, t: &Tables, mut len: usize, dist: usize) {
        debug_assert!((1..=WINDOW).contains(&dist));
        self.end_run();
        while len > 0 {
            let take = if len > MAX_MATCH {
                // Leave at least 3 for the remainder.
                if len - MAX_MATCH < 3 {
                    MAX_MATCH - 3
                } else {
                    MAX_MATCH
                }
            } else {
                len
            };
            self.buf.push(token_match(take, dist));
            self.lit_freq[t.len_sym[take] as usize & 511] += 1;
            self.dist_freq[dist_symbol(t, dist) & 31] += 1;
            len -= take;
        }
    }
}

// ---- hash matcher ---------------------------------------------------------

pub struct Matcher {
    /// Four candidates per bucket, most recent in the low 32 bits, so an
    /// insert is one shift and store.
    table: Box<[u128; 1 << HASH_BITS]>,
    /// Where each read starting with a given 16 bases was last seen.
    reads: Box<[u32; 1 << READ_BITS]>,
    /// Positions are stored as `base + pos + 1`. Each block raises `base`
    /// past everything stored before, which retires old entries without
    /// clearing the tables.
    base: u32,
    /// Token buffer, kept between blocks.
    tokens: Vec<Token>,
}

const READ_KEY: usize = 16;
const READ_BITS: u32 = 12;

impl Matcher {
    fn new() -> Self {
        Self {
            table: vec![0u128; 1 << HASH_BITS]
                .into_boxed_slice()
                .try_into()
                .unwrap(),
            reads: vec![0u32; 1 << READ_BITS]
                .into_boxed_slice()
                .try_into()
                .unwrap(),
            base: 0,
            tokens: Vec::new(),
        }
    }

    /// Start a block of `len` bytes.
    fn reset(&mut self, len: usize) {
        if self.base as usize + len + 1 >= u32::MAX as usize {
            self.table.fill(0);
            self.reads.fill(0);
            self.base = 0;
        }
    }

    /// End a block of `len` bytes: nothing stored so far matches again.
    fn retire(&mut self, len: usize) {
        self.base += len as u32 + 1;
    }

    #[inline(always)]
    fn stored(&self, pos: usize) -> u32 {
        self.base + pos as u32 + 1
    }

    /// The position in this block that a stored value refers to, if it was
    /// stored in this block before `pos`.
    #[inline(always)]
    fn position(&self, stored: u32, pos: usize) -> Option<usize> {
        let cand = stored.wrapping_sub(self.base).wrapping_sub(1) as usize;
        (cand < pos).then_some(cand)
    }

    /// Length and distance of an earlier read that this sequence line
    /// (`data[start..end]`) repeats from its first base, as PCR and optical
    /// duplicates do. One probe per read, so it costs next to nothing.
    #[inline(always)]
    fn duplicate(&mut self, data: &[u8], start: usize, end: usize) -> (usize, usize) {
        if end - start < READ_KEY {
            return (0, 0);
        }
        let a = u64::from_le_bytes(data[start..start + 8].try_into().unwrap());
        let b = u64::from_le_bytes(data[start + 8..start + 16].try_into().unwrap());
        let h = ((a ^ b.rotate_left(29)).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - READ_BITS))
            as usize;
        let c = self.reads[h];
        self.reads[h] = self.stored(start);
        let Some(cand) = self.position(c, start) else {
            return (0, 0);
        };
        let dist = start - cand;
        if dist > WINDOW || data[cand..cand + READ_KEY] != data[start..start + READ_KEY] {
            return (0, 0);
        }
        (extend(data, cand, start, end - start), dist)
    }

    #[inline(always)]
    fn hash(bytes: u32) -> usize {
        (bytes.wrapping_mul(0x9E37_79B1) >> (32 - HASH_BITS)) as usize
    }
    #[inline(always)]
    fn load4(data: &[u8], pos: usize) -> u32 {
        u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap())
    }
    #[inline(always)]
    fn insert(&mut self, data: &[u8], pos: usize) {
        let h = Self::hash(Self::load4(data, pos));
        self.table[h] = self.table[h] << 32 | u128::from(self.stored(pos));
    }

    /// Longest match at `pos` among the bucket candidates (and inserts `pos`).
    #[inline(always)]
    fn find(&mut self, data: &[u8], pos: usize, end: usize) -> (usize, usize) {
        let cur = Self::load4(data, pos);
        let h = Self::hash(cur);
        let bucket = self.table[h];
        self.table[h] = bucket << 32 | u128::from(self.stored(pos));
        let max = (end - pos).min(MAX_MATCH * 4);
        let mut best = (0usize, 0usize);
        let ways = [
            bucket as u32,
            (bucket >> 32) as u32,
            (bucket >> 64) as u32,
            (bucket >> 96) as u32,
        ];
        for stored in ways {
            let Some(cand) = self.position(stored, pos) else {
                break;
            };
            let dist = pos - cand;
            if dist > WINDOW || Self::load4(data, cand) != cur {
                continue;
            }
            let len = extend(data, cand, pos, max);
            if len > best.0 {
                best = (len, dist);
                if len >= NICE_MATCH {
                    break;
                }
            }
        }
        best
    }

    /// Greedy LZ over `data[start..end]`, with history anywhere before in `data`.
    fn encode(&mut self, t: &Tables, data: &[u8], start: usize, end: usize, tok: &mut Tokens) {
        let mut pos = start;
        // A match needs MIN_MATCH bytes up to `end`.
        let limit = (end + 1).saturating_sub(MIN_MATCH);
        while pos < limit {
            let (len, dist) = self.find(data, pos, end);
            if len == 0 {
                tok.literal(data[pos]);
                pos += 1;
                continue;
            }
            tok.matched(t, len, dist);
            // Insert the first positions inside the match so later data can
            // refer to them, then skip it. Inserting all of a long match
            // costs 10% of the speed for 0.5% of the size.
            let stop = (pos + len).min(limit);
            let mut p = pos + 1;
            let insert_until = (pos + 1 + INSERT_IN_MATCH).min(stop);
            while p < insert_until {
                self.insert(data, p);
                p += 1;
            }
            pos += len;
        }
        tok.literals(&data[pos.max(start)..end]);
    }
}

/// Match length between `data[cand..]` and `data[pos..]`, at least 4, at most `max`.
#[inline(always)]
fn extend(data: &[u8], cand: usize, pos: usize, max: usize) -> usize {
    // Slice once so the loops need no bounds checks.
    let (a, b) = (&data[cand..cand + max], &data[pos..pos + max]);
    let mut len = 4;
    while len + 8 <= max {
        let x = u64::from_le_bytes(a[len..len + 8].try_into().unwrap())
            ^ u64::from_le_bytes(b[len..len + 8].try_into().unwrap());
        if x != 0 {
            return len + (x.trailing_zeros() / 8) as usize;
        }
        len += 8;
    }
    while len < max && a[len] == b[len] {
        len += 1;
    }
    len
}

// ---- FASTQ tokenizer ------------------------------------------------------

#[inline(always)]
fn is_base(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
    )
}

fn line_end(data: &[u8], from: usize) -> usize {
    memchr::memchr(b'\n', &data[from..]).map_or(data.len(), |i| from + i)
}

/// Whether a record (header, sequence, plus) appears to start at `pos`.
fn record_starts(data: &[u8], pos: usize) -> bool {
    if pos >= data.len() || data[pos] != b'@' {
        return false;
    }
    let h_end = line_end(data, pos);
    if h_end + 1 >= data.len() {
        return true; // header only: harmless to treat as a record
    }
    let s = h_end + 1;
    let s_end = line_end(data, s);
    let probe = &data[s..(s + 8).min(s_end)];
    if !probe.iter().all(|&b| is_base(b)) {
        return false;
    }
    let p = s_end + 1;
    p >= data.len() || data[p] == b'+'
}

struct Tokenizer<'a> {
    data: &'a [u8],
    t: &'static Tables,
    matcher: &'a mut Matcher,
    tok: Tokens,
    /// Start of the previous header unit ("\n@..." or "@..."), with its length.
    prev_header: Option<(usize, usize)>,
}

impl<'a> Tokenizer<'a> {
    /// Header unit: data[start..end] is "\n@name...\n" (the leading newline
    /// is optional), including its trailing newline when present.
    fn header(&mut self, start: usize, end: usize) {
        let data = self.data;
        let mut pos = start;
        let mut tail = end;
        if let Some((pstart, plen)) = self.prev_header {
            let prev = &data[pstart..pstart + plen];
            let cur = &data[start..end];
            let n = prev.len().min(cur.len());
            let mut lcp = 0;
            while lcp < n && prev[lcp] == cur[lcp] {
                lcp += 1;
            }
            let dist = start - pstart;
            if lcp >= MIN_MATCH && dist <= WINDOW {
                self.tok.matched(self.t, lcp, dist);
                let lim = (start + lcp + 1).saturating_sub(MIN_MATCH);
                let mut p = start;
                while p < lim {
                    self.matcher.insert(data, p);
                    p += 1;
                }
                pos = start + lcp;
            }
            // Common suffix (at least 3 bytes) that does not overlap the prefix.
            let mut lcs = 0;
            let room = n - lcp.min(n);
            while lcs < room && prev[prev.len() - 1 - lcs] == cur[cur.len() - 1 - lcs] {
                lcs += 1;
            }
            let sdist = end - (pstart + plen);
            if lcs >= 3 && sdist <= WINDOW {
                tail = end - lcs;
                self.matcher.encode(self.t, data, pos, tail, &mut self.tok);
                self.tok.matched(self.t, lcs, sdist);
                self.prev_header = Some((start, end - start));
                return;
            }
        }
        self.matcher.encode(self.t, data, pos, tail, &mut self.tok);
        self.prev_header = Some((start, end - start));
    }

    fn sequence(&mut self, start: usize, end: usize) {
        let (len, dist) = self.matcher.duplicate(self.data, start, end);
        // Below about 24 bases a match costs more bits than two-bit literals.
        if len >= 24 {
            self.tok.matched(self.t, len, dist);
            self.tok.literals(&self.data[start + len..end]);
            return;
        }
        self.tok.literals(&self.data[start..end]);
    }

    fn generic(&mut self, start: usize, end: usize) {
        self.matcher
            .encode(self.t, self.data, start, end, &mut self.tok);
    }

    fn run(&mut self) {
        let data = self.data;
        let n = data.len();
        let mut pos = 0;
        // Phase: find the first record start, encoding what precedes generically.
        let mut phase_start = 0;
        'outer: while pos < n {
            // Search for a record start at a line beginning.
            let mut probe = pos;
            loop {
                if record_starts(data, probe) {
                    break;
                }
                let e = line_end(data, probe);
                if e >= n {
                    // No record in the rest: generic.
                    self.generic(phase_start, n);
                    break 'outer;
                }
                probe = e + 1;
            }
            // Bytes before the record (including the newline right before it,
            // which belongs to the header unit if there is one).
            let hdr_start = if probe > phase_start && data[probe - 1] == b'\n' {
                probe - 1
            } else {
                probe
            };
            if hdr_start > phase_start {
                self.generic(phase_start, hdr_start);
            }
            pos = hdr_start;
            // Records loop.
            loop {
                // header unit: data[pos..h_end]
                let h_end = line_end(data, if data[pos] == b'\n' { pos + 1 } else { pos });
                if h_end >= n {
                    self.header(pos, n);
                    break 'outer;
                }
                self.header(pos, h_end + 1);
                // sequence unit: bases only
                let s_end = line_end(data, h_end + 1);
                self.sequence(h_end + 1, s_end);
                if s_end >= n {
                    break 'outer;
                }
                // plus + quality unit: "\n+...\n" + quals (no trailing newline)
                let p_end = line_end(data, s_end + 1);
                let q_end = if p_end >= n {
                    n
                } else {
                    line_end(data, p_end + 1)
                };
                self.generic(s_end, q_end);
                if q_end >= n {
                    break 'outer;
                }
                // Next record expected at q_end + 1.
                let next = q_end + 1;
                if next < n && data[next] == b'@' && {
                    // cheap shape check: a '+' two lines later
                    let e1 = line_end(data, next);
                    e1 + 1 >= n || {
                        let e2 = line_end(data, e1 + 1);
                        e2 + 1 >= n || data[e2 + 1] == b'+'
                    }
                } {
                    pos = q_end; // header unit includes the newline
                    continue;
                }
                // Lost the phase (or the block ends at a record boundary):
                // the newline joins the generic segment.
                if next >= n {
                    self.generic(q_end, n);
                    break 'outer;
                }
                phase_start = q_end;
                pos = next;
                continue 'outer;
            }
        }
    }
}

// ---- Huffman --------------------------------------------------------------

/// Code lengths (max `limit`) for `freq`, with at least two non-zero lengths
/// when any symbol is used so the code is decodable by every inflater.
fn code_lengths(freq: &[u32], limit: usize, lens: &mut [u8]) {
    lens.iter_mut().for_each(|l| *l = 0);
    let mut syms: Vec<(u32, usize)> = freq
        .iter()
        .enumerate()
        .filter(|(_, &f)| f > 0)
        .map(|(i, &f)| (f, i))
        .collect();
    match syms.len() {
        0 => return,
        1 => {
            let s = syms[0].1;
            lens[s] = 1;
            lens[if s == 0 { 1 } else { 0 }] = 1;
            return;
        }
        _ => {}
    }
    syms.sort_unstable();
    // Huffman depths by the two-queue method.
    let n = syms.len();
    let mut node_freq: Vec<u64> = syms.iter().map(|s| s.0 as u64).collect();
    let mut parent: Vec<usize> = vec![usize::MAX; 2 * n];
    let mut leaf_next = 0usize;
    let mut inner_next = n;
    let mut inner_tail = n;
    node_freq.resize(2 * n, 0);
    for _ in 0..n - 1 {
        let pick = |node_freq: &Vec<u64>, leaf_next: &mut usize, inner_next: &mut usize| {
            if *leaf_next < n
                && (*inner_next >= inner_tail || node_freq[*leaf_next] <= node_freq[*inner_next])
            {
                *leaf_next += 1;
                *leaf_next - 1
            } else {
                *inner_next += 1;
                *inner_next - 1
            }
        };
        let a = pick(&node_freq, &mut leaf_next, &mut inner_next);
        let b = pick(&node_freq, &mut leaf_next, &mut inner_next);
        node_freq[inner_tail] = node_freq[a] + node_freq[b];
        parent[a] = inner_tail;
        parent[b] = inner_tail;
        inner_tail += 1;
    }
    let root = inner_tail - 1;
    let mut depth = vec![0u32; 2 * n];
    for i in (n..root).rev() {
        depth[i] = depth[parent[i]] + 1;
    }
    let mut count = [0u32; 64];
    for i in 0..n {
        let d = if i == root { 1 } else { depth[parent[i]] + 1 } as usize;
        count[d.min(limit)] += 1;
    }
    // Enforce the limit (miniz method): keep the Kraft sum at exactly 1.
    let mut total: u64 = (1..=limit).map(|i| (count[i] as u64) << (limit - i)).sum();
    while total != 1 << limit {
        count[limit] -= 1;
        for i in (1..limit).rev() {
            if count[i] > 0 {
                count[i] -= 1;
                count[i + 1] += 2;
                break;
            }
        }
        total -= 1;
    }
    // Longest codes to the least frequent symbols.
    let mut idx = 0;
    for len in (1..=limit).rev() {
        for _ in 0..count[len] {
            lens[syms[idx].1] = len as u8;
            idx += 1;
        }
    }
}

/// Canonical codes (bit-reversed for LSB-first writing).
fn canonical_codes(lens: &[u8], codes: &mut [u16]) {
    let mut bl_count = [0u16; 16];
    for &l in lens {
        bl_count[l as usize] += 1;
    }
    bl_count[0] = 0;
    let mut next = [0u16; 16];
    let mut code = 0u16;
    for bits in 1..16 {
        code = (code + bl_count[bits - 1]) << 1;
        next[bits] = code;
    }
    for (i, &l) in lens.iter().enumerate() {
        if l > 0 {
            let c = next[l as usize];
            next[l as usize] += 1;
            codes[i] = (c.reverse_bits()) >> (16 - l);
        } else {
            codes[i] = 0;
        }
    }
}

/// Writes bits into space reserved up front: every put stores eight bytes
/// at the write position and advances by the whole bytes, with no length
/// checks or branches.
struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    /// Write position in `out`'s spare capacity.
    at: usize,
    /// Reserved end; `put` never writes past it.
    end: usize,
    buf: u64,
    n: u32,
}
impl<'a> BitWriter<'a> {
    /// Space for `bytes` of output past what `out` holds.
    fn new(out: &'a mut Vec<u8>, bytes: usize) -> Self {
        out.reserve(bytes + 8);
        let at = out.len();
        let end = out.capacity();
        Self {
            out,
            at,
            end,
            buf: 0,
            n: 0,
        }
    }
    /// Up to 56 bits.
    #[inline(always)]
    fn put_small(&mut self, bits: u64, count: u32) {
        debug_assert!(count <= 56 && self.n <= 7 && bits >> count == 0);
        self.buf |= bits << self.n;
        self.n += count;
        assert!(self.at + 8 <= self.end);
        // SAFETY: `at + 8 <= end`, the capacity reserved in `new`; bytes past
        // the write position are rewritten by later puts or cut by `finish`.
        unsafe {
            self.out
                .as_mut_ptr()
                .add(self.at)
                .cast::<u64>()
                .write_unaligned(self.buf.to_le());
        }
        // At most 7 + 56 bits are pending, so the shift is below 64.
        let whole = self.n & !7;
        self.at += (whole / 8) as usize;
        self.buf >>= whole & 63;
        self.n -= whole;
    }
    fn finish(mut self) {
        if self.n > 0 {
            self.put_small(0, 7);
        }
        // SAFETY: every byte below `at` was written by `put_small`.
        unsafe { self.out.set_len(self.at) };
    }
}

/// Compresses `input` into a raw deflate stream appended to `out`. One
/// dynamic Huffman block with BFINAL set.
pub fn deflate(input: &[u8], matcher: &mut Matcher, out: &mut Vec<u8>) {
    let t = tables();
    matcher.reset(input.len());
    let mut buf = std::mem::take(&mut matcher.tokens);
    buf.clear();
    // Every token covers at least one input byte.
    buf.reserve(input.len() + 1);
    let mut tz = Tokenizer {
        data: input,
        t,
        matcher,
        tok: Tokens {
            buf,
            lit_freq: [0; 512],
            lit4: [[0; 256]; 4],
            dist_freq: [0; 32],
            run: 0,
        },
        prev_header: None,
    };
    tz.run();
    let matcher = tz.matcher;
    matcher.retire(input.len());
    let mut tok = tz.tok;
    tok.end_run();
    for table in &tok.lit4 {
        for (f, &c) in tok.lit_freq.iter_mut().zip(table) {
            *f += c;
        }
    }
    tok.lit_freq[256] = 1;

    let mut lit_lens = [0u8; 286];
    let mut dist_lens = [0u8; 30];
    code_lengths(&tok.lit_freq[..286], 15, &mut lit_lens);
    code_lengths(&tok.dist_freq[..30], 15, &mut dist_lens);
    let mut lit_codes = [0u16; 286];
    let mut dist_codes = [0u16; 30];
    canonical_codes(&lit_lens, &mut lit_codes);
    canonical_codes(&dist_lens, &mut dist_codes);
    let hlit = (257..=286)
        .rev()
        .find(|&i| lit_lens[i - 1] != 0)
        .unwrap_or(257)
        .max(257);
    let hdist = (1..=30)
        .rev()
        .find(|&i| dist_lens[i - 1] != 0)
        .unwrap_or(1)
        .max(1);

    // Run-length encode the code lengths.
    let mut all: Vec<u8> = Vec::with_capacity(hlit + hdist);
    all.extend_from_slice(&lit_lens[..hlit]);
    all.extend_from_slice(&dist_lens[..hdist]);
    let mut cl_syms: Vec<(u8, u8)> = Vec::new(); // (symbol, extra value)
    let mut i = 0;
    while i < all.len() {
        let v = all[i];
        let mut run = 1;
        while i + run < all.len() && all[i + run] == v {
            run += 1;
        }
        if v == 0 && run >= 3 {
            let mut r = run;
            while r >= 3 {
                if r >= 11 {
                    let take = r.min(138);
                    cl_syms.push((18, (take - 11) as u8));
                    r -= take;
                } else {
                    cl_syms.push((17, (r - 3) as u8));
                    r = 0;
                }
            }
            for _ in 0..r {
                cl_syms.push((0, 0));
            }
            i += run;
        } else if v != 0 && run >= 4 {
            cl_syms.push((v, 0));
            let mut r = run - 1;
            while r >= 3 {
                let take = r.min(6);
                cl_syms.push((16, (take - 3) as u8));
                r -= take;
            }
            for _ in 0..r {
                cl_syms.push((v, 0));
            }
            i += run;
        } else {
            cl_syms.push((v, 0));
            i += 1;
        }
    }
    let mut cl_freq = [0u32; 19];
    for &(s, _) in &cl_syms {
        cl_freq[s as usize] += 1;
    }
    let mut cl_lens = [0u8; 19];
    code_lengths(&cl_freq, 7, &mut cl_lens);
    let mut cl_codes = [0u16; 19];
    canonical_codes(&cl_lens, &mut cl_codes);
    let hclen = (4..=19)
        .rev()
        .find(|&i| cl_lens[CL_ORDER[i - 1]] != 0)
        .unwrap_or(4)
        .max(4);

    // Worst case: every byte a 15-bit literal, plus the block header.
    let mut w = BitWriter::new(out, input.len() * 2 + 1024);
    w.put_small(1, 1); // BFINAL
    w.put_small(2, 2); // dynamic
    w.put_small((hlit - 257) as u64, 5);
    w.put_small((hdist - 1) as u64, 5);
    w.put_small((hclen - 4) as u64, 4);
    for &o in &CL_ORDER[..hclen] {
        w.put_small(cl_lens[o] as u64, 3);
    }
    for &(s, extra) in &cl_syms {
        w.put_small(cl_codes[s as usize] as u64, cl_lens[s as usize] as u32);
        match s {
            16 => w.put_small(extra as u64, 2),
            17 => w.put_small(extra as u64, 3),
            18 => w.put_small(extra as u64, 7),
            _ => {}
        }
    }
    // Symbols.
    let four = lit_lens[..256].iter().max().is_none_or(|&m| m <= 14);
    // Code and length in one entry: one load per literal.
    let lit_table: [u32; 256] =
        std::array::from_fn(|b| u32::from(lit_codes[b]) | u32::from(lit_lens[b]) << 24);
    let mut ip = 0usize;
    for &tk in &tok.buf {
        if tk & MATCH == 0 {
            let run = &input[ip..ip + tk as usize];
            ip += tk as usize;
            // Several codes per put: up to 56 bits fit.
            let code = |b: u8| {
                let e = lit_table[b as usize];
                (u64::from(e & 0xff_ffff), e >> 24)
            };
            let mut rest = run;
            if four {
                let mut quads = rest.chunks_exact(4);
                for q in &mut quads {
                    let (c0, l0) = code(q[0]);
                    let (c1, l1) = code(q[1]);
                    let (c2, l2) = code(q[2]);
                    let (c3, l3) = code(q[3]);
                    let bits = c0 | c1 << l0 | c2 << (l0 + l1) | c3 << (l0 + l1 + l2);
                    w.put_small(bits, l0 + l1 + l2 + l3);
                }
                rest = quads.remainder();
            }
            let mut triples = rest.chunks_exact(3);
            for t3 in &mut triples {
                let (c0, l0) = code(t3[0]);
                let (c1, l1) = code(t3[1]);
                let (c2, l2) = code(t3[2]);
                w.put_small(c0 | c1 << l0 | c2 << (l0 + l1), l0 + l1 + l2);
            }
            for &b in triples.remainder() {
                let (c, l) = code(b);
                w.put_small(c, l);
            }
        } else {
            ip += ((tk >> 16) & 0x1ff) as usize;
            let len = ((tk >> 16) & 0x1ff) as usize;
            let dist = (tk & 0xffff) as usize + 1;
            let ls = t.len_sym[len] as usize;
            let li = ls - 257;
            let count = lit_lens[ls] as u32;
            let bits = lit_codes[ls] as u64 | (((len - LEN_BASE[li] as usize) as u64) << count);
            w.put_small(bits, count + LEN_EXTRA[li] as u32);
            let ds = dist_symbol(t, dist);
            let count = dist_lens[ds] as u32;
            let bits = dist_codes[ds] as u64 | (((dist - DIST_BASE[ds] as usize) as u64) << count);
            w.put_small(bits, count + DIST_EXTRA[ds] as u32);
        }
    }
    w.put_small(lit_codes[256] as u64, lit_lens[256] as u32);
    w.finish();
    matcher.tokens = tok.buf;
}

/// A complete gzip member for `input`, carrying the `FP` member-size extra
/// field that `io` uses to decompress members in parallel.
pub(crate) fn gzip_member(input: &[u8], out: &mut Vec<u8>) {
    thread_local! {
        static MATCHER: std::cell::RefCell<Matcher> = std::cell::RefCell::new(Matcher::new());
    }
    const FEXTRA: u8 = 4;
    // Header, the deflate worst case (see `deflate`) and trailer, reserved
    // once so nothing is copied as the member grows.
    out.reserve(20 + input.len() * 2 + 1024 + 8 + 8);
    let start = out.len();
    out.extend_from_slice(&[0x1f, 0x8b, 8, FEXTRA, 0, 0, 0, 0, 0, 3]);
    out.extend_from_slice(&8u16.to_le_bytes());
    out.extend_from_slice(&crate::gzip::MEMBER_SIZE_ID);
    out.extend_from_slice(&4u16.to_le_bytes());
    let size_at = out.len();
    out.extend_from_slice(&[0; 4]);
    MATCHER.with(|m| deflate(input, &mut m.borrow_mut(), out));
    out.extend_from_slice(&crc32fast::hash(input).to_le_bytes());
    out.extend_from_slice(&(input.len() as u32).to_le_bytes());
    let total = (out.len() - start) as u32;
    out[size_at..size_at + 4].copy_from_slice(&total.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn roundtrip(data: &[u8]) {
        let mut out = Vec::new();
        gzip_member(data, &mut out);
        let mut back = Vec::new();
        flate2::read::MultiGzDecoder::new(&out[..])
            .read_to_end(&mut back)
            .unwrap();
        assert_eq!(back, data, "roundtrip of {} bytes", data.len());
        let mut d = libdeflater::Decompressor::new();
        let mut buf = vec![0u8; data.len() + 1];
        let n = d.gzip_decompress(&out, &mut buf).unwrap();
        assert_eq!(&buf[..n], data);
    }

    #[test]
    fn positions_wrap_without_stale_matches() {
        let data: Vec<u8> = (0..200_000u32)
            .map(|i| b"ACGT\n@+F:"[(i.wrapping_mul(2654435761) >> 28) as usize % 9])
            .collect();
        let mut m = Matcher::new();
        // Blocks right below, across and after the point where the stored
        // positions would overflow and the tables are cleared.
        m.base = u32::MAX - 450_000;
        for _ in 0..4 {
            let mut out = Vec::new();
            deflate(&data, &mut m, &mut out);
            let mut back = Vec::new();
            flate2::read::DeflateDecoder::new(&out[..])
                .read_to_end(&mut back)
                .unwrap();
            assert!(back == data);
        }
        assert!(m.base < 1_000_000, "the tables were reset");
    }

    #[test]
    fn roundtrips_fastq_and_arbitrary_bytes() {
        let record = |i: usize| {
            format!(
                "@ERR3534523.{i} {i}/1\nACGT{}GGGGGG\n+\nFFFF,F:FF{}FF\n",
                "TGCA".repeat((i % 37) + 1),
                "F:".repeat(i % 5)
            )
        };
        let fastq: String = (1..3000).map(record).collect();
        let bytes = fastq.as_bytes();
        for len in [
            0,
            1,
            2,
            3,
            4,
            5,
            17,
            100,
            1000,
            bytes.len() / 2,
            bytes.len(),
        ] {
            roundtrip(&bytes[..len]);
        }
        for start in [1, 7, 23, 60, 61, 62, 63, 64, 65, 66] {
            roundtrip(&bytes[start..start + 2000]);
        }
        // Not FASTQ, runs, random-ish bytes, single symbol.
        roundtrip(&[b'A'; 100000]);
        roundtrip(&[b'A'; 1]);
        roundtrip(b"+\n+\n+\n@\n@@\n");
        let noise: Vec<u8> = (0..300000u32)
            .map(|i| (i.wrapping_mul(2654435761) >> 13) as u8)
            .collect();
        roundtrip(&noise);
        let mixed: Vec<u8> = bytes
            .iter()
            .chain(noise.iter())
            .chain(bytes.iter())
            .copied()
            .collect();
        roundtrip(&mixed);
    }
}

#[cfg(test)]
mod fuzz {
    use std::io::Read;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n.max(1) as u64) as usize
        }
        fn pick(&mut self, xs: &[u8]) -> u8 {
            xs[self.below(xs.len())]
        }
    }

    /// FASTQ-shaped text with the irregularities the tokenizer must survive:
    /// duplicate reads, long and empty lines, CRLF, stray '@' and '+',
    /// truncated records and random bytes.
    fn sample(r: &mut Rng) -> Vec<u8> {
        let mut out = Vec::new();
        let target = match r.below(4) {
            0 => r.below(64),
            1 => r.below(4096),
            _ => r.below(300_000),
        };
        let quals: Vec<u8> = match r.below(3) {
            0 => b"F:,#".to_vec(),
            1 => (b'!'..=b'J').collect(),
            _ => (0..=255).collect(),
        };
        let mut last_seq = Vec::new();
        let mut i = 0u64;
        while out.len() < target {
            i += 1;
            match r.below(40) {
                0 => out.extend((0..r.below(500)).map(|_| r.next() as u8)),
                1 => out.extend_from_slice(b"@\n+\n@@\n\n+\r\n"),
                2 => out.extend(std::iter::repeat_n(r.pick(b"AF@+\n"), r.below(70_000))),
                _ => {
                    let crlf = r.below(50) == 0;
                    let nl: &[u8] = if crlf { b"\r\n" } else { b"\n" };
                    out.extend_from_slice(
                        format!("@READ.{i} {}:{}/1", r.below(3), r.below(100_000)).as_bytes(),
                    );
                    if r.below(100) == 0 {
                        out.extend(std::iter::repeat_n(b'x', r.below(40_000)));
                    }
                    out.extend_from_slice(nl);
                    let len = if r.below(20) == 0 {
                        r.below(2000)
                    } else {
                        100 + r.below(51)
                    };
                    let seq: Vec<u8> = if r.below(5) == 0 && !last_seq.is_empty() {
                        last_seq.clone()
                    } else {
                        (0..len).map(|_| r.pick(b"ACGTN")).collect()
                    };
                    out.extend_from_slice(&seq);
                    out.extend_from_slice(nl);
                    out.push(b'+');
                    out.extend_from_slice(nl);
                    let n = if r.below(30) == 0 {
                        r.below(seq.len() + 1)
                    } else {
                        seq.len()
                    };
                    out.extend((0..n).map(|_| r.pick(&quals)));
                    out.extend_from_slice(nl);
                    last_seq = seq;
                }
            }
        }
        // Start and end anywhere, as 256 KiB members of a stream do.
        let a = r.below(out.len().min(200) + 1);
        let b = out.len() - r.below((out.len() - a).min(200) + 1);
        out[a..b].to_vec()
    }

    #[test]
    fn random_inputs_roundtrip() {
        let iters: usize = std::env::var("FQ_FUZZ_ITERS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(300);
        let seed: u64 = std::env::var("FQ_FUZZ_SEED")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0x5eed);
        let mut r = Rng(seed | 1);
        let mut d = libdeflater::Decompressor::new();
        let mut gz = Vec::new();
        for it in 0..iters {
            let data = sample(&mut r);
            gz.clear();
            super::gzip_member(&data, &mut gz);
            let mut back = Vec::with_capacity(data.len());
            flate2::read::GzDecoder::new(&gz[..])
                .read_to_end(&mut back)
                .unwrap_or_else(|e| panic!("iteration {it}: {e}"));
            assert!(
                back == data,
                "iteration {it}: flate2 mismatch on {} bytes",
                data.len()
            );
            let mut buf = vec![0; data.len() + 1];
            let n = d
                .gzip_decompress(&gz, &mut buf)
                .unwrap_or_else(|e| panic!("iteration {it}: {e:?}"));
            assert!(buf[..n] == data[..], "iteration {it}: libdeflate mismatch");
        }
    }
}

#[cfg(test)]
mod corpus {
    /// Size and speed against libdeflate on a real file, in 256 KiB members:
    /// FQ_CORPUS=path cargo test --release corpus -- --ignored --nocapture
    #[test]
    #[ignore]
    fn corpus() {
        let path = std::env::var("FQ_CORPUS").expect("FQ_CORPUS");
        let limit: usize = std::env::var("FQ_LIMIT")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(64 << 20);
        let mut data = std::fs::read(path).unwrap();
        data.truncate(limit);
        let chunks: Vec<&[u8]> = data.chunks(256 << 10).collect();
        let mut out = Vec::with_capacity(data.len());
        let t = std::time::Instant::now();
        let mut fq = 0;
        for c in &chunks {
            out.clear();
            super::gzip_member(c, &mut out);
            fq += out.len();
        }
        let fq_t = t.elapsed().as_secs_f64();
        let mb = data.len() as f64 / 1e6;
        println!(
            "fq     {:>11} ratio {:.4} {:>6.0} MB/s",
            fq,
            fq as f64 / data.len() as f64,
            mb / fq_t
        );
        let mut l4 = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(4).unwrap());
        let mut buf = vec![0; l4.gzip_compress_bound(256 << 10)];
        let t = std::time::Instant::now();
        let (mut gated, mut picked) = (0, 0);
        for c in &chunks {
            if super::suits(c) {
                out.clear();
                super::gzip_member(c, &mut out);
                gated += out.len();
                picked += 1;
            } else {
                gated += l4.gzip_compress(c, &mut buf).unwrap();
            }
        }
        let s = t.elapsed().as_secs_f64();
        println!(
            "gated  {:>11} ratio {:.4} {:>6.0} MB/s  fq blocks {picked}/{}",
            gated,
            gated as f64 / data.len() as f64,
            mb / s,
            chunks.len()
        );
        for level in [1, 2, 3, 4, 6] {
            let mut c =
                libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
            let mut buf = vec![0; c.gzip_compress_bound(256 << 10)];
            let t = std::time::Instant::now();
            let mut n = 0;
            for ch in &chunks {
                n += c.gzip_compress(ch, &mut buf).unwrap();
            }
            let s = t.elapsed().as_secs_f64();
            println!(
                "ld L{level}  {:>11} ratio {:.4} {:>6.0} MB/s  fq/ld {:+.2}%",
                n,
                n as f64 / data.len() as f64,
                mb / s,
                (fq as f64 / n as f64 - 1.0) * 100.0
            );
        }
    }
}
