//! Parallel decompression of a single gzip member, the approach of pugz and
//! rapidgzip. A deflate stream can only be decoded from its start, so:
//!
//! 1. The compressed file is cut into regions. A worker finds a dynamic
//!    Huffman block that starts in its region by trying every bit offset and
//!    keeping the first whose header is valid and whose whole block decodes.
//! 2. It decodes from there, without the 32 KiB of history before the region:
//!    a back-reference into that unknown window becomes a placeholder symbol.
//!    Once 32 KiB in a row have no placeholder no later copy can reach one, so
//!    from there it decodes straight to bytes.
//! 3. It stops at the first block boundary at or after the next region's
//!    start, which is where the next worker should have started.
//! 4. A coordinator takes regions in order, fills in placeholders from the
//!    previous region's last 32 KiB, checks the CRC and length, and sends the
//!    bytes on. A region whose start does not match the previous region's end
//!    (no block found, or a false start) is decoded again from the known end.
//!
//! A wrong guess costs time, never output: only regions that continue exactly
//! from verified output are used, and the trailer CRC covers everything.
use std::fs::File;
use std::io::{self, Read};

const WINDOW: usize = 32768;
/// Compressed bytes per region.
const REGION: usize = 2 << 20;
/// Bytes read past a region for the block that crosses its end.
const MARGIN: usize = 1 << 20;

// bits

struct Bits<'a> {
    data: &'a [u8],
    /// Next byte to load.
    pos: usize,
    buf: u64,
    n: u32,
}

impl<'a> Bits<'a> {
    fn at(data: &'a [u8], bit: usize) -> Self {
        let mut bits = Self {
            data,
            pos: bit / 8,
            buf: 0,
            n: 0,
        };
        bits.refill();
        bits.consume((bit % 8) as u32);
        bits
    }

    /// Ensure at least 56 bits. Past the end it reads zeros; `overrun` tells.
    #[inline(always)]
    fn refill(&mut self) {
        if self.pos + 8 <= self.data.len() {
            let v = u64::from_le_bytes(self.data[self.pos..self.pos + 8].try_into().unwrap());
            // Bytes beyond those counted are loaded again at the same bit
            // positions next time, so or-ing them in early is harmless.
            self.buf |= v << self.n;
            let take = (63 - self.n) / 8;
            self.pos += take as usize;
            self.n += take * 8;
        } else {
            while self.n <= 56 {
                let b = self.data.get(self.pos).copied().unwrap_or(0);
                self.buf |= u64::from(b) << self.n;
                self.pos += 1;
                self.n += 8;
            }
        }
    }

    #[inline(always)]
    fn consume(&mut self, k: u32) {
        self.buf >>= k;
        self.n -= k;
    }

    #[inline(always)]
    fn take(&mut self, k: u32) -> u32 {
        if self.n < k {
            self.refill();
        }
        let v = (self.buf & ((1u64 << k) - 1)) as u32;
        self.consume(k);
        v
    }

    fn bit_pos(&self) -> usize {
        self.pos * 8 - self.n as usize
    }

    fn overrun(&self) -> bool {
        self.bit_pos() > self.data.len() * 8
    }
}

// huffman

/// A two-level decode table for a canonical Huffman code, indexed by the
/// next `P` bits. Its entries carry everything a symbol needs, so decoding
/// one is a lookup, a shift and at most one add:
///
/// * bits 0..8: bits to consume, the code length plus the symbol's extra bits;
/// * bits 8..12: the code length, to shift the extra bits down;
/// * bits 16..32: the symbol, a literal byte, a length or distance base, or
///   the start of a subtable for longer codes (then bits 8..12 hold the
///   subtable's index width and the entry consumes the `P` bits that led to
///   it; subtable entries count their code length from there).
///
/// The buffer is kept between blocks, which each bring their own codes.
#[derive(Default)]
struct Table<const P: u32> {
    entries: Vec<u32>,
}

/// Set on literal entries.
const LITERAL: u32 = 1 << 15;
/// Set on every entry that is not a symbol, literal, length or distance.
const SPECIAL: u32 = 1 << 14;
/// With `SPECIAL`: the entry points at a subtable.
const SUBTABLE: u32 = 1 << 13;
/// With `SPECIAL`: the end-of-block symbol.
const END: u32 = 1 << 12;
/// An unused code, or a symbol deflate does not define.
const INVALID: u32 = SPECIAL;

/// Index bits for literal/length codes, distance codes and the code that
/// compresses their lengths, which is at most seven bits long.
const LITLEN: u32 = 10;
const DISTANCE: u32 = 8;
const LENGTHS: u32 = 7;

impl<const P: u32> Table<P> {
    /// Fill the table for the canonical code `lens`, with `entry` giving a
    /// symbol's entry from its code length. Returns false if the code is
    /// over-subscribed, or incomplete with more than a single one-bit code
    /// (zlib's rule).
    fn build(&mut self, lens: &[u8], entry: impl Fn(usize, u32) -> u32) -> bool {
        let mut count = [0u16; 16];
        for &l in lens {
            count[usize::from(l)] += 1;
        }
        count[0] = 0;
        let mut left = 1i32;
        let mut max = 0;
        for (len, &c) in count.iter().enumerate().skip(1) {
            left = (left << 1) - i32::from(c);
            if left < 0 {
                return false;
            }
            if c > 0 {
                max = len;
            }
        }
        if left > 0 && max > 1 {
            return false;
        }
        let mut next = [0u32; 16];
        let mut code = 0u32;
        for len in 1..16 {
            code = (code + u32::from(count[len - 1])) << 1;
            next[len] = code;
        }
        let primary = 1usize << P;
        let table = &mut self.entries;
        table.clear();
        table.resize(primary, INVALID);
        let reverse = |code: u32, len: usize| (code.reverse_bits() >> (32 - len)) as usize;
        if max <= P as usize {
            for (sym, &l) in lens.iter().enumerate() {
                let len = usize::from(l);
                if len == 0 {
                    continue;
                }
                let rev = reverse(next[len], len);
                next[len] += 1;
                let e = entry(sym, len as u32);
                for j in (rev..primary).step_by(1 << len) {
                    table[j] = e;
                }
            }
            return true;
        }
        // Widths of the subtables that long codes need, by primary index.
        let mut starts = next;
        let mut width = [0u8; 1 << 10];
        for &l in lens {
            let len = usize::from(l);
            if len == 0 {
                continue;
            }
            let rev = reverse(starts[len], len);
            starts[len] += 1;
            if len > P as usize {
                let slot = &mut width[rev & (primary - 1)];
                *slot = (*slot).max((len - P as usize) as u8);
            }
        }
        for (p, &bits) in width[..primary].iter().enumerate() {
            if bits > 0 {
                table[p] =
                    (table.len() as u32) << 16 | SPECIAL | SUBTABLE | u32::from(bits) << 8 | P;
                table.resize(table.len() + (1 << bits), INVALID);
            }
        }
        for (sym, &l) in lens.iter().enumerate() {
            let len = usize::from(l);
            if len == 0 {
                continue;
            }
            let rev = reverse(next[len], len);
            next[len] += 1;
            if len <= P as usize {
                let e = entry(sym, len as u32);
                for j in (rev..primary).step_by(1 << len) {
                    table[j] = e;
                }
            } else {
                let p = rev & (primary - 1);
                let offset = (table[p] >> 16) as usize;
                let extra = len - P as usize;
                let e = entry(sym, extra as u32);
                for j in ((rev >> P)..(1 << width[p])).step_by(1 << extra) {
                    table[offset + j] = e;
                }
            }
        }
        true
    }

    /// The entry for the next code in `bits`, following a subtable pointer
    /// and consuming its primary bits. Needs at least 15 bits in `bits`.
    #[inline(always)]
    fn entry(&self, bits: &mut Bits) -> u32 {
        // SAFETY: the primary table has `1 << P` entries, and a subtable
        // pointer addresses `1 << width` entries from its start.
        unsafe {
            let e = *self
                .entries
                .get_unchecked((bits.buf & ((1 << P) - 1)) as usize);
            if e & (SPECIAL | SUBTABLE) != SPECIAL | SUBTABLE {
                return e;
            }
            bits.consume(P);
            let width = (e >> 8) & 15;
            let index = (e >> 16) as usize + (bits.buf & ((1 << width) - 1)) as usize;
            *self.entries.get_unchecked(index)
        }
    }
}

impl Table<LITLEN> {
    fn litlen(&mut self, lens: &[u8]) -> bool {
        self.build(lens, |sym, len| match sym {
            0..=255 => (sym as u32) << 16 | LITERAL | len << 8 | len,
            256 => SPECIAL | END | len << 8 | len,
            257..=285 => {
                let extra = u32::from(LEN_EXTRA[sym - 257]);
                u32::from(LEN_BASE[sym - 257]) << 16 | len << 8 | (len + extra)
            }
            _ => INVALID,
        })
    }
}

impl Table<DISTANCE> {
    fn distance(&mut self, lens: &[u8]) -> bool {
        self.build(lens, |sym, len| match sym {
            0..=29 => {
                let extra = u32::from(DIST_EXTRA[sym]);
                u32::from(DIST_BASE[sym]) << 16 | len << 8 | (len + extra)
            }
            _ => INVALID,
        })
    }
}

impl Table<LENGTHS> {
    fn lengths(&mut self, lens: &[u8]) -> bool {
        self.build(lens, |sym, len| (sym as u32) << 16 | len << 8 | len)
    }

    /// The next code-length symbol. Needs at least 7 bits in `bits`.
    #[inline(always)]
    fn decode(&self, bits: &mut Bits) -> Option<u32> {
        let e = self.entries[(bits.buf & ((1 << LENGTHS) - 1)) as usize];
        if e & SPECIAL != 0 {
            return None;
        }
        bits.consume(e & 0xff);
        Some(e >> 16)
    }
}

/// The base of a length or distance entry plus its extra bits, taken from
/// `saved`, the bit buffer before the entry was consumed.
#[inline(always)]
fn with_extra(entry: u32, saved: u64) -> usize {
    let extra = (saved & ((1u64 << (entry & 0xff)) - 1)) >> ((entry >> 8) & 15);
    (entry >> 16) as usize + extra as usize
}

/// The literal/length and distance tables of a block.
#[derive(Default)]
struct Tables {
    litlen: Table<LITLEN>,
    distance: Table<DISTANCE>,
    lengths: Table<LENGTHS>,
}

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
const CODE_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

fn fixed_tables() -> Tables {
    let mut lens = [0u8; 288];
    lens[..144].fill(8);
    lens[144..256].fill(9);
    lens[256..280].fill(7);
    lens[280..].fill(8);
    let mut tables = Tables::default();
    assert!(tables.litlen.litlen(&lens) && tables.distance.distance(&[5u8; 32]));
    tables
}

/// Read the literal/length and distance codes of a dynamic block, after its
/// 3 header bits, into `tables`. Returns false if they are invalid.
fn dynamic_tables(bits: &mut Bits, tables: &mut Tables) -> bool {
    let hlit = bits.take(5) as usize + 257;
    let hdist = bits.take(5) as usize + 1;
    let hclen = bits.take(4) as usize + 4;
    if hlit > 286 || hdist > 30 {
        return false;
    }
    let mut code_lens = [0u8; 19];
    for &i in &CODE_ORDER[..hclen] {
        code_lens[i] = bits.take(3) as u8;
    }
    if !tables.lengths.lengths(&code_lens) {
        return false;
    }
    let mut lens = [0u8; 286 + 30];
    let total = hlit + hdist;
    let mut i = 0;
    while i < total {
        bits.refill();
        let Some(sym) = tables.lengths.decode(bits) else {
            return false;
        };
        let (value, repeat) = match sym {
            0..=15 => (sym as u8, 1),
            16 => {
                let Some(&previous) = lens[..i].last() else {
                    return false;
                };
                (previous, 3 + bits.take(2) as usize)
            }
            17 => (0, 3 + bits.take(3) as usize),
            _ => (0, 11 + bits.take(7) as usize),
        };
        if i + repeat > total {
            return false;
        }
        lens[i..i + repeat].fill(value);
        i += repeat;
    }
    lens[256] != 0
        && tables.litlen.litlen(&lens[..hlit])
        && tables.distance.distance(&lens[hlit..total])
}

// inflate

/// Output symbols: bytes, or (in `u16`) placeholders for the unknown window.
trait Symbol: Copy {
    fn byte(b: u8) -> Self;
    /// The symbol for window byte `index` (0 is 32 KiB back), if representable.
    fn window(index: usize) -> Option<Self>;
}

impl Symbol for u8 {
    fn byte(b: u8) -> Self {
        b
    }
    fn window(_: usize) -> Option<Self> {
        None
    }
}

impl Symbol for u16 {
    fn byte(b: u8) -> Self {
        u16::from(b)
    }
    fn window(index: usize) -> Option<Self> {
        Some(256 + index as u16)
    }
}

#[derive(Debug)]
enum Error {
    /// Invalid deflate data.
    Bad,
    /// Ran past the end of the bytes given.
    More,
}

struct Inflater<'a> {
    bits: Bits<'a>,
    fixed: Option<Tables>,
    /// The tables of the dynamic block being decoded.
    dynamic: Tables,
}

impl<'a> Inflater<'a> {
    /// Decode one block into `out`; returns whether it was the final block.
    /// A copy reaching before `out[0]` becomes window placeholders if `T`
    /// has them, and is an error otherwise. `last_window` is set to one past
    /// the last placeholder written.
    /// Decode one block into `out[*n..]`, advancing `*n`; returns whether it
    /// was the final block. `out` is used as a buffer with spare length past
    /// `*n` and grows as needed; callers truncate it to `*n` at the end.
    /// A copy reaching before `out[0]` becomes window placeholders if `T`
    /// has them, and is an error otherwise. `last_window` is kept at or past
    /// the end of the last placeholder written.
    fn block<T: Symbol>(
        &mut self,
        out: &mut Vec<T>,
        n: &mut usize,
        last_window: &mut usize,
    ) -> Result<bool, Error> {
        let bits = &mut self.bits;
        let last = bits.take(1) == 1;
        let (litlen, dist) = match bits.take(2) {
            0 => {
                let skip = bits.n % 8;
                bits.consume(skip);
                let len = bits.take(16) as usize;
                if len as u32 != !bits.take(16) & 0xffff {
                    return Err(Error::Bad);
                }
                grow(out, *n + len);
                for slot in &mut out[*n..*n + len] {
                    *slot = T::byte(bits.take(8) as u8);
                }
                *n += len;
                return if bits.overrun() {
                    Err(Error::More)
                } else {
                    Ok(last)
                };
            }
            1 => {
                let fixed = self.fixed.get_or_insert_with(fixed_tables);
                (&fixed.litlen, &fixed.distance)
            }
            2 => {
                if !dynamic_tables(bits, &mut self.dynamic) {
                    return Err(Error::Bad);
                }
                (&self.dynamic.litlen, &self.dynamic.distance)
            }
            _ => return Err(Error::Bad),
        };
        let mut at = *n;
        let mut window_end = *last_window;
        // The bit reader lives in locals while symbols are decoded: output is
        // written through a raw pointer, and with the reader behind a
        // reference every such write would force its fields to be reloaded.
        let mut local = Bits {
            data: bits.data,
            pos: bits.pos,
            buf: bits.buf,
            n: bits.n,
        };
        let outcome = 'symbols: loop {
            // Room for three literals and the longest copy, which may write up
            // to a chunk past its end.
            if at + 3 + 258 + CHUNK > out.len() {
                grow(out, at + 3 + 258 + CHUNK);
            }
            let base = out.as_mut_ptr();
            // SAFETY (for every write through `base` in this iteration): the
            // check above leaves room for three literals at `at` and then a
            // copy of at most 258 elements rounded up to a whole chunk.
            let room = out.len();
            let put = |slot: usize, b: u8| {
                debug_assert!(slot < room);
                unsafe { *base.add(slot) = T::byte(b) }
            };
            // A refill leaves at least 56 bits: three literal codes of up to
            // 15 bits each, or one length, distance and their extra local.
            local.refill();
            let mut entry = litlen.entry(&mut local);
            if entry & LITERAL != 0 {
                local.consume(entry & 0xff);
                put(at, (entry >> 16) as u8);
                at += 1;
                entry = litlen.entry(&mut local);
                if entry & LITERAL != 0 {
                    local.consume(entry & 0xff);
                    put(at, (entry >> 16) as u8);
                    at += 1;
                    entry = litlen.entry(&mut local);
                    if entry & LITERAL != 0 {
                        local.consume(entry & 0xff);
                        put(at, (entry >> 16) as u8);
                        at += 1;
                        // Past the end the bits read as zeros, which can
                        // decode as literals forever: stop there as copies do.
                        if local.overrun() {
                            break 'symbols Err(Error::More);
                        }
                        continue;
                    }
                }
                // The entry read before this refill is still the next code:
                // a refill only adds bits above those already present.
                if local.n < 48 {
                    local.refill();
                }
            }
            if entry & SPECIAL != 0 {
                if entry & END == 0 {
                    break 'symbols Err(Error::Bad);
                }
                local.consume(entry & 0xff);
                break 'symbols Ok(());
            }
            let saved = local.buf;
            local.consume(entry & 0xff);
            let len = with_extra(entry, saved);
            let entry = dist.entry(&mut local);
            if entry & SPECIAL != 0 {
                break 'symbols Err(Error::Bad);
            }
            let saved = local.buf;
            local.consume(entry & 0xff);
            let distance = with_extra(entry, saved);
            let p = at;
            // Copying from before the last placeholder may copy placeholders.
            if p.saturating_sub(distance) < window_end {
                window_end = p + len;
            }
            if distance <= p {
                let from = p - distance;
                debug_assert!(p + len.next_multiple_of(CHUNK) <= room);
                // SAFETY: `from + len <= p + len`, within the room checked
                // above, and everything before `p` has been written.
                unsafe {
                    if distance >= CHUNK {
                        // Copy whole chunks, which compile to vector loads
                        // and stores, and may write past `len` into spare
                        // buffer. Chunks never overlap the bytes they read.
                        let mut i = 0;
                        while i < len {
                            std::ptr::copy_nonoverlapping(
                                base.add(from + i),
                                base.add(p + i),
                                CHUNK,
                            );
                            i += CHUNK;
                        }
                    } else if distance == 1 {
                        let v = *base.add(from);
                        for i in 0..len {
                            *base.add(p + i) = v;
                        }
                    } else {
                        for i in 0..len {
                            *base.add(p + i) = *base.add(from + i);
                        }
                    }
                }
            } else {
                for i in 0..len {
                    let back = (p + i) as isize - distance as isize;
                    out[p + i] = if back >= 0 {
                        out[back as usize]
                    } else {
                        match T::window((WINDOW as isize + back) as usize) {
                            Some(symbol) => symbol,
                            None => break 'symbols Err(Error::Bad),
                        }
                    };
                }
                window_end = p + len;
            }
            at += len;
            if local.overrun() {
                break 'symbols Err(Error::More);
            }
        };
        *bits = local;
        *n = at;
        *last_window = window_end;
        outcome?;
        if bits.overrun() {
            return Err(Error::More);
        }
        Ok(last)
    }
}

/// Elements per copy step for back-references at least this far back.
const CHUNK: usize = 16;

/// Make `out` at least `need` long. Capacity doubles so growth is amortized,
/// but only a little past `need` is initialized: zero-filling the whole
/// doubled buffer would make all of it resident and write it twice.
fn grow<T: Symbol>(out: &mut Vec<T>, need: usize) {
    if out.len() < need {
        let len = need + (1 << 18);
        if out.capacity() < len {
            out.reserve(len.max(out.capacity() * 2) - out.len());
        }
        out.resize(len, T::byte(0));
    }
}

/// Buffers kept for reuse between regions, so that the decoder does not
/// allocate and free several megabytes per region. Freed large blocks stay
/// resident in most allocators (macOS kept 0.5 to 1 GB of them), and fresh
/// ones fault in new pages. A buffer is allocated only when no idle one is
/// available, so the buffers that exist never outnumber the most that were
/// in use at once; each channel's capacity caps the idle ones.
///
/// Reused buffers keep their length and old contents: the decoder tracks how
/// much it has written and grows only past that, so pages are not zeroed
/// again for every region.
struct Pool {
    /// Compressed input, one per decoding thread.
    input: (
        crossbeam_channel::Sender<Vec<u8>>,
        crossbeam_channel::Receiver<Vec<u8>>,
    ),
    /// Output that may hold window placeholders.
    words: (
        crossbeam_channel::Sender<Vec<u16>>,
        crossbeam_channel::Receiver<Vec<u16>>,
    ),
    /// Decoded output: buffers the reader is done with come back here.
    output: (
        crossbeam_channel::Sender<Vec<u8>>,
        crossbeam_channel::Receiver<Vec<u8>>,
    ),
    /// Capacity for new output and placeholder buffers.
    output_size: Estimate,
    words_size: Estimate,
}

/// How large a buffer must be: the largest size seen so far, or a guess
/// before the first.
struct Estimate {
    largest: std::sync::atomic::AtomicUsize,
    guess: usize,
}

impl Estimate {
    fn new(guess: usize) -> Self {
        Self {
            largest: std::sync::atomic::AtomicUsize::new(0),
            guess,
        }
    }
    fn get(&self) -> usize {
        match self.largest.load(std::sync::atomic::Ordering::Relaxed) {
            0 => self.guess,
            n => n,
        }
    }
    fn saw(&self, n: usize) {
        self.largest
            .fetch_max(n, std::sync::atomic::Ordering::Relaxed);
    }
}

/// A buffer from `idle` with capacity for at least `need` elements. A new
/// one gets 1/8 more, so that later regions rarely outgrow it, and it is
/// reserved exactly rather than grown by doubling, which copies and frees
/// every smaller block on the way (freed blocks stay resident).
fn reuse<T>(idle: &crossbeam_channel::Receiver<Vec<T>>, need: usize) -> Vec<T> {
    let v = idle.try_recv().unwrap_or_default();
    if v.capacity() >= need {
        return v;
    }
    // Not `reserve`: reallocating would copy the old contents.
    drop(v);
    Vec::with_capacity(need + need / 8)
}

impl Pool {
    fn new(threads: usize, region: usize, output: Recycle) -> Self {
        Self {
            input: crossbeam_channel::bounded(threads + 2),
            words: crossbeam_channel::bounded(threads + 2),
            output,
            // About a FASTQ region's output (ratio 4 to 5); for regions with
            // less output the rest is only address space, never touched.
            output_size: Estimate::new(region.saturating_mul(5)),
            // Placeholder output is often the whole region: copies of copies
            // keep reaching the unknown window (FASTQ names and qualities).
            words_size: Estimate::new(region.saturating_mul(5)),
        }
    }
    fn input(&self) -> Vec<u8> {
        self.input.1.try_recv().unwrap_or_default()
    }
    /// A placeholder buffer with room for a region's output.
    fn words(&self) -> Vec<u16> {
        reuse(&self.words.1, self.words_size.get())
    }
    /// An output buffer with room for a region's output.
    fn output(&self) -> Vec<u8> {
        reuse(&self.output.1, self.output_size.get())
    }
    fn give_input(&self, v: Vec<u8>) {
        if v.capacity() > 0 {
            let _ = self.input.0.try_send(v);
        }
    }
    fn give_words(&self, v: Vec<u16>) {
        if v.capacity() > 0 {
            let _ = self.words.0.try_send(v);
        }
    }
    fn give_output(&self, v: Vec<u8>) {
        if v.capacity() > 0 {
            let _ = self.output.0.try_send(v);
        }
    }
}

/// Output buffers handed back for reuse: the reader sends on one end, the
/// decoder takes from the other.
pub(crate) type Recycle = (
    crossbeam_channel::Sender<Vec<u8>>,
    crossbeam_channel::Receiver<Vec<u8>>,
);

/// Output buffers that can be in use at once with `ahead` regions in
/// flight: one per region in flight, the output queue's [`OUTPUT_DEPTH`],
/// and the reader's own, plus one to spare. The reader's recycling channel
/// holds this many, so a returned buffer is never freed.
pub(crate) fn output_buffers(ahead: usize) -> usize {
    ahead + OUTPUT_DEPTH + 2
}

/// Regions queued between the region decoder and the reader.
pub(crate) const OUTPUT_DEPTH: usize = 2;

/// A region decoded from a block start up to a block boundary.
struct Decoded {
    start: usize,
    end: usize,
    /// The member ended here (the final block was decoded).
    last: bool,
    /// Output that may hold window placeholders.
    head: Vec<u16>,
    /// The region's whole output: `head.len()` bytes for the head, filled in
    /// by [`resolve`], then the plain bytes decoded after it (the tail).
    /// Empty until needed if there is no tail.
    out: Vec<u8>,
    /// CRC32 of the tail.
    tail_crc: u32,
}

impl Decoded {
    fn tail(&self) -> &[u8] {
        self.out.get(self.head.len()..).unwrap_or_default()
    }
    fn len(&self) -> usize {
        self.head.len() + self.tail().len()
    }
}

/// Decode from bit `start` of `data` until a block boundary at or after bit
/// `stop`, or the final block. `first` means the stream starts at `start`, so
/// there is no window and copies may not reach before it.
///
/// The tail is decoded straight into the output buffer after room for the
/// head, so the region is delivered as one buffer without a copy.
fn decode(
    data: &[u8],
    start: usize,
    stop: usize,
    first: bool,
    pool: &Pool,
) -> Result<Decoded, Error> {
    let mut inflater = Inflater {
        bits: Bits::at(data, start),
        fixed: None,
        dynamic: Tables::default(),
    };
    let mut head = Vec::new();
    // Taken when the tail starts: a region that is all head needs its
    // output buffer only once it can be resolved.
    let mut out = Vec::new();
    if first {
        out = pool.output();
    } else {
        head = pool.words();
    }
    let (mut head_n, mut out_n) = (0, 0);
    let mut in_tail = first;
    let mut last_window = 0;
    let mut last = false;
    while inflater.bits.bit_pos() < stop || inflater.bits.bit_pos() == start {
        let result = if in_tail {
            inflater.block(&mut out, &mut out_n, &mut last_window)
        } else {
            inflater.block(&mut head, &mut head_n, &mut last_window)
        };
        last = match result {
            Ok(last) => last,
            Err(e) => {
                pool.give_words(head);
                pool.give_output(out);
                return Err(e);
            }
        };
        if !in_tail && head_n - last_window.min(head_n) >= WINDOW {
            // No later copy can reach a placeholder: continue in bytes after
            // the head, whose last 32 KiB are plain bytes and the history.
            out = pool.output();
            grow(&mut out, head_n);
            for (o, &s) in out[head_n - WINDOW..head_n]
                .iter_mut()
                .zip(&head[head_n - WINDOW..head_n])
            {
                *o = s as u8;
            }
            out_n = head_n;
            in_tail = true;
        }
        if last {
            break;
        }
    }
    head.truncate(head_n);
    out.truncate(out_n);
    pool.words_size.saw(head_n);
    pool.output_size.saw(out_n.max(head_n));
    let tail_crc = crc32fast::hash(out.get(head_n..).unwrap_or_default());
    Ok(Decoded {
        start,
        end: inflater.bits.bit_pos(),
        last,
        head,
        out,
        tail_crc,
    })
}

/// The first bit at or after `from` and before `to` where a non-final
/// dynamic block starts whose header is valid and whose whole block decodes.
fn find_block(data: &[u8], from: usize, to: usize, scratch: &mut Vec<u16>) -> Option<usize> {
    let to = to.min(data.len().saturating_sub(8) * 8);
    // Table buffers, reused for every candidate.
    let mut tables = Tables::default();
    for bit in from..to {
        let byte = bit / 8;
        let w = u32::from_le_bytes(data[byte..byte + 4].try_into().unwrap()) >> (bit % 8);
        // BFINAL = 0, BTYPE = 2, HLIT <= 29, HDIST <= 29.
        if w & 7 != 4 || (w >> 3) & 31 > 29 || (w >> 8) & 31 > 29 {
            continue;
        }
        // The code-length code must be complete or a single one-bit code.
        let hclen = ((w >> 13) & 15) as usize + 4;
        let mut bits = Bits::at(data, bit + 17);
        let mut lens = [0u8; 19];
        for &i in &CODE_ORDER[..hclen] {
            lens[i] = bits.take(3) as u8;
        }
        if !tables.lengths.lengths(&lens) {
            continue;
        }
        let mut inflater = Inflater {
            bits: Bits::at(data, bit),
            fixed: None,
            dynamic: std::mem::take(&mut tables),
        };
        let (mut n, mut last_window) = (0, 0);
        if inflater.block(scratch, &mut n, &mut last_window).is_ok() {
            return Some(bit);
        }
        tables = inflater.dynamic;
    }
    None
}

// driver

/// Compressed input: a file, or bytes in memory for tests and fuzzing.
pub(crate) trait Source: Sync {
    fn len(&self) -> usize;
    /// Fill `buf` with the bytes at `offset`.
    fn read_at(&self, buf: &mut [u8], offset: usize) -> io::Result<()>;
    /// The bytes from `offset` to the end, for the members after the first.
    fn rest(&self, offset: usize) -> io::Result<Box<dyn Read + Send + '_>>;
}

pub(crate) struct FileSource {
    file: File,
    len: usize,
}

impl FileSource {
    pub(crate) fn open(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let len = file.metadata()?.len() as usize;
        Ok(Self { file, len })
    }
}

impl Source for FileSource {
    fn len(&self) -> usize {
        self.len
    }

    #[cfg(unix)]
    fn read_at(&self, buf: &mut [u8], offset: usize) -> io::Result<()> {
        std::os::unix::fs::FileExt::read_exact_at(&self.file, buf, offset as u64)
    }

    #[cfg(windows)]
    fn read_at(&self, mut buf: &mut [u8], offset: usize) -> io::Result<()> {
        let mut offset = offset as u64;
        while !buf.is_empty() {
            match std::os::windows::fs::FileExt::seek_read(&self.file, buf, offset)? {
                0 => return Err(io::ErrorKind::UnexpectedEof.into()),
                n => {
                    buf = &mut buf[n..];
                    offset += n as u64;
                }
            }
        }
        Ok(())
    }

    fn rest(&self, offset: usize) -> io::Result<Box<dyn Read + Send + '_>> {
        let mut file = self.file.try_clone()?;
        io::Seek::seek(&mut file, io::SeekFrom::Start(offset as u64))?;
        Ok(Box::new(file))
    }
}

impl Source for &[u8] {
    fn len(&self) -> usize {
        <[u8]>::len(self)
    }

    fn read_at(&self, buf: &mut [u8], offset: usize) -> io::Result<()> {
        let bytes = self
            .get(offset..offset + buf.len())
            .ok_or(io::ErrorKind::UnexpectedEof)?;
        buf.copy_from_slice(bytes);
        Ok(())
    }

    fn rest(&self, offset: usize) -> io::Result<Box<dyn Read + Send + '_>> {
        Ok(Box::new(self.get(offset..).unwrap_or_default()))
    }
}

/// Bytes `from..to` of the input.
fn read_range(src: &dyn Source, from: usize, to: usize) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    read_range_into(src, from, to, &mut buf)?;
    Ok(buf)
}

/// Bytes `from..to` of the input into `buf`, reusing its capacity.
fn read_range_into(src: &dyn Source, from: usize, to: usize, buf: &mut Vec<u8>) -> io::Result<()> {
    buf.resize(to - from, 0);
    src.read_at(buf, from)
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

/// Bit offset of the deflate stream after a gzip member header at `data[0]`.
fn header_len(data: &[u8]) -> io::Result<usize> {
    let bad = || invalid("Invalid gzip header");
    // Reserved flag bits must be zero.
    if data.len() < 10 || data[..3] != [0x1f, 0x8b, 8] || data[3] & 0xe0 != 0 {
        return Err(bad());
    }
    let flags = data[3];
    let mut at = 10;
    if flags & 4 != 0 {
        let xlen = data.get(at..at + 2).ok_or_else(bad)?;
        at += 2 + usize::from(u16::from_le_bytes([xlen[0], xlen[1]]));
    }
    for flag in [8, 16] {
        if flags & flag != 0 {
            at += 1 + data
                .get(at..)
                .ok_or_else(bad)?
                .iter()
                .position(|&b| b == 0)
                .ok_or_else(bad)?;
        }
    }
    if flags & 2 != 0 {
        // The header CRC covers every header byte before it.
        let crc = data.get(at..at + 2).ok_or_else(bad)?;
        if u16::from_le_bytes([crc[0], crc[1]]) != crc32fast::hash(&data[..at]) as u16 {
            return Err(invalid("Invalid gzip header checksum"));
        }
        at += 2;
    }
    if at > data.len() {
        return Err(bad());
    }
    Ok(at * 8)
}

/// Region `k`: decode from its first block, or from `start` when given,
/// reading more past the region if the block crossing its end needs it.
/// `Ok(None)` means no block start was found in the region.
fn region(
    src: &dyn Source,
    size: usize,
    k: usize,
    deflate: usize,
    start: Option<usize>,
    pool: &Pool,
) -> io::Result<Option<Decoded>> {
    let len = src.len();
    let begin = deflate / 8 + k * size;
    let end = (begin + size).min(len);
    let from = start.map_or(begin, |s| s / 8);
    let mut margin = MARGIN;
    let mut data = pool.input();
    let result = loop {
        let read_to = (end + margin).min(len);
        if let Err(e) = read_range_into(src, from, read_to, &mut data) {
            break Err(e);
        }
        let base = from * 8;
        let local_start = match start {
            Some(s) => s - base,
            None if k == 0 => deflate - base,
            None => {
                let mut scratch = pool.words();
                let found = find_block(&data, 0, (end - from) * 8, &mut scratch);
                pool.give_words(scratch);
                match found {
                    Some(bit) => bit,
                    None => break Ok(None),
                }
            }
        };
        let first = k == 0 && start.is_none_or(|s| s == deflate);
        match decode(&data, local_start, (end - from) * 8, first, pool) {
            Ok(mut decoded) => {
                decoded.start += base;
                decoded.end += base;
                break Ok(Some(decoded));
            }
            Err(Error::More) if read_to < len => margin *= 4,
            Err(Error::More) => break Err(invalid("Truncated gzip input")),
            Err(Error::Bad) => break Err(invalid("Invalid deflate data")),
        }
    };
    pool.give_input(data);
    result
}

/// A resolved region: its dispatch number, the head's CRC32 (or the error
/// resolving it), the head's length, and the region.
type Resolved = (usize, io::Result<u32>, usize, Decoded);

enum Job {
    Decode(usize),
    /// Fill a decoded region's placeholders from the 32 KiB before it.
    Resolve(usize, Decoded, std::sync::Arc<Vec<u8>>),
}

/// Write the bytes for `head` into `out` (as long as `head`), filling
/// placeholders from `window`, the bytes before it.
fn resolve(head: &[u16], window: &[u8], out: &mut [u8]) -> io::Result<()> {
    let mut lut = [0u8; 256 + WINDOW];
    for (i, b) in lut[..256].iter_mut().enumerate() {
        *b = i as u8;
    }
    let known = window.len().min(WINDOW);
    lut[256 + WINDOW - known..].copy_from_slice(&window[window.len() - known..]);
    // Near the start of the stream the window is short: a placeholder for a
    // byte before the stream is invalid data.
    if known < WINDOW
        && head
            .iter()
            .any(|&s| s >= 256 && usize::from(s - 256) < WINDOW - known)
    {
        return Err(invalid("Invalid deflate distance"));
    }
    // Most groups of symbols hold no placeholder and only need narrowing,
    // which compiles to vector code; the rest go through the table.
    let mut bytes = out.chunks_exact_mut(32);
    let mut symbols = head.chunks_exact(32);
    for (bytes, symbols) in (&mut bytes).zip(&mut symbols) {
        if symbols.iter().fold(0, |all, &s| all | s) < 256 {
            for (o, &s) in bytes.iter_mut().zip(symbols) {
                *o = s as u8;
            }
        } else {
            for (o, &s) in bytes.iter_mut().zip(symbols) {
                *o = lut[usize::from(s)];
            }
        }
    }
    for (o, &s) in bytes.into_remainder().iter_mut().zip(symbols.remainder()) {
        *o = lut[usize::from(s)];
    }
    Ok(())
}

/// Whether a gzip file of `len` bytes spans enough regions to decode in parallel.
pub(crate) fn worth_it(len: u64) -> bool {
    len >= 3 * REGION as u64
}

/// Decompress the gzip file at `path` on `threads` workers, sending output in
/// order on `sender`, then an empty buffer at the end. Members after the first
/// are decoded sequentially.
/// The reader hands buffers it is done with back on `recycle`; its channel
/// should hold [`output_buffers`]`(threads + 2)` so that none is freed.
pub(crate) fn decompress(
    path: &str,
    threads: usize,
    sender: crossbeam_channel::Sender<io::Result<Vec<u8>>>,
    recycle: Recycle,
) {
    let result =
        FileSource::open(path).and_then(|src| run(&src, region_size(), threads, &sender, recycle));
    if let Err(e) = result {
        let _ = sender.send(Err(e));
    }
}

/// `decompress` from any source in regions of `size` compressed bytes.
pub(crate) fn decompress_source(
    src: &dyn Source,
    size: usize,
    threads: usize,
    sender: crossbeam_channel::Sender<io::Result<Vec<u8>>>,
    recycle: Recycle,
) {
    if let Err(e) = run(src, size, threads, &sender, recycle) {
        let _ = sender.send(Err(e));
    }
}

/// Region size for files: `REGION` shared between the inputs being decoded,
/// unless tests asked for smaller regions so that small inputs cross many of
/// them. Two inputs with full regions each ran no faster and held twice the
/// decoded output in memory.
fn region_size() -> usize {
    match crate::io::forced_split_region() {
        0 => REGION / crate::io::split_region_shares().max(1),
        size => size,
    }
}

fn run(
    src: &dyn Source,
    size: usize,
    threads: usize,
    sender: &crossbeam_channel::Sender<io::Result<Vec<u8>>>,
    recycle: Recycle,
) -> io::Result<()> {
    let len = src.len();
    let deflate = header_len(&read_range(src, 0, len.min(64 << 10))?)?;
    let regions = (len - deflate / 8).div_ceil(size).max(1);
    // Stops the workers early if the reader goes away.
    let stopped = std::sync::atomic::AtomicBool::new(false);
    let stopped = &stopped;
    let pool = &Pool::new(threads, size, recycle);
    let result = std::thread::scope(|scope| -> io::Result<()> {
        // Regions decode to about 3x their size; keep a few in flight per worker.
        let ahead = threads + 2;
        let (jobs, queue) = crossbeam_channel::bounded::<Job>(2 * ahead);
        // Resolve jobs go first: each holds a region's placeholder output,
        // the largest buffer, until it is done.
        let (resolves, urgent) = crossbeam_channel::unbounded::<Job>();
        let (decoded_tx, decoded_rx) = crossbeam_channel::unbounded();
        let (resolved_tx, resolved_rx) = crossbeam_channel::unbounded();
        for _ in 0..threads {
            let (queue, urgent, decoded_tx, resolved_tx) = (
                queue.clone(),
                urgent.clone(),
                decoded_tx.clone(),
                resolved_tx.clone(),
            );
            scope.spawn(move || {
                let next = || match urgent.try_recv() {
                    Ok(job) => Ok(job),
                    Err(_) => {
                        let mut select = crossbeam_channel::Select::new();
                        select.recv(&urgent);
                        select.recv(&queue);
                        let op = select.select();
                        match op.index() {
                            0 => op.recv(&urgent),
                            _ => op.recv(&queue),
                        }
                    }
                };
                // Both channels close together, at the end of the member.
                while let Ok(job) = next() {
                    if stopped.load(std::sync::atomic::Ordering::Relaxed) {
                        break;
                    }
                    match job {
                        // A failure here may only be a false block start, so
                        // the coordinator decides; it re-decodes from known output.
                        Job::Decode(k) => {
                            let _ = decoded_tx.send((k, region(src, size, k, deflate, None, pool)));
                        }
                        Job::Resolve(k, mut decoded, window) => {
                            let n = decoded.head.len();
                            if decoded.out.len() < n {
                                decoded.out = pool.output();
                                decoded.out.resize(n, 0);
                                decoded.out.truncate(n);
                            }
                            let head = &mut decoded.out[..n];
                            let crc = resolve(&decoded.head, &window, head)
                                .map(|()| crc32fast::hash(head));
                            // Keep only the length; the buffer is free now.
                            let head = std::mem::take(&mut decoded.head);
                            let _ = resolved_tx.send((k, crc, n, decoded));
                            pool.give_words(head);
                        }
                    }
                }
            });
        }
        drop((queue, decoded_tx, resolved_tx));
        let stop = || stopped.store(true, std::sync::atomic::Ordering::Relaxed);
        let workers_stopped = || invalid("Decompression workers stopped");
        let mut queued = 0;
        let mut decoded_pending = std::collections::BTreeMap::new();
        let mut resolved_pending = std::collections::BTreeMap::new();
        let mut window: Vec<u8> = Vec::with_capacity(2 * WINDOW);
        let mut crc = crc32fast::Hasher::new();
        let mut total = 0u64;
        let mut pos = deflate;
        // Resolve jobs are numbered in dispatch order: skipped regions have
        // none; `next_send` is the next of them to send.
        let (mut k, mut dispatched, mut next_send) = (0, 0, 0);
        let mut finished = None;
        // Take a resolved region and send what is now in order; `false` if
        // the reader went away.
        let deliver = |item: Resolved,
                       next_send: &mut usize,
                       resolved_pending: &mut std::collections::BTreeMap<_, _>,
                       crc: &mut crc32fast::Hasher|
         -> io::Result<bool> {
            let (i, head_crc, head_len, decoded) = item;
            resolved_pending.insert(i, (head_crc, head_len, decoded));
            while let Some((head_crc, head_len, decoded)) = resolved_pending.remove(next_send) {
                let head_crc: u32 = head_crc?;
                let decoded: Decoded = decoded;
                crc.combine(&crc32fast::Hasher::new_with_initial_len(
                    head_crc,
                    head_len as u64,
                ));
                crc.combine(&crc32fast::Hasher::new_with_initial_len(
                    decoded.tail_crc,
                    (decoded.out.len() - head_len) as u64,
                ));
                if decoded.out.is_empty() {
                    pool.give_output(decoded.out);
                } else if sender.send(Ok(decoded.out)).is_err() {
                    return Ok(false);
                }
                *next_send += 1;
            }
            Ok(true)
        };
        macro_rules! deliver {
            ($item:expr) => {
                if !deliver($item, &mut next_send, &mut resolved_pending, &mut crc)? {
                    stop();
                    return Ok(());
                }
            };
        }
        // At most `ahead` regions are in flight, from queued to sent: queued
        // from `k` on, and taken but not yet sent. Regions decoded out of
        // order wait in `decoded_pending` within that bound.
        macro_rules! top_up {
            () => {
                while queued < regions && queued - k + dispatched - next_send < ahead {
                    jobs.send(Job::Decode(queued)).ok();
                    queued += 1;
                }
            };
        }
        while k < regions && finished.is_none() {
            top_up!();
            if queued == k {
                // Every slot is taken by regions being resolved.
                deliver!(resolved_rx.recv().map_err(|_| workers_stopped())?);
                continue;
            }
            let result = loop {
                if let Some(r) = decoded_pending.remove(&k) {
                    break r;
                }
                crossbeam_channel::select! {
                    recv(decoded_rx) -> r => {
                        let (i, r) = r.map_err(|_| workers_stopped())?;
                        decoded_pending.insert(i, r);
                    }
                    recv(resolved_rx) -> r => {
                        deliver!(r.map_err(|_| workers_stopped())?);
                        top_up!();
                    }
                }
            };
            let region_end = (deflate / 8 + (k + 1) * size).min(len) * 8;
            if pos >= region_end {
                // An earlier region's last block covered this one.
                if let Ok(Some(d)) = result {
                    pool.give_words(d.head);
                    pool.give_output(d.out);
                }
                k += 1;
                continue;
            }
            let decoded = match result {
                Ok(Some(d)) if d.start == pos => d,
                other => {
                    // A false start: its buffers are reused for the redo.
                    if let Ok(Some(d)) = other {
                        pool.give_words(d.head);
                        pool.give_output(d.out);
                    }
                    region(src, size, k, deflate, Some(pos), pool)?
                        .expect("a region decoded from a known start always has a block")
                }
            };
            // The next region's window is this region's last 32 KiB. Only that
            // much of the head is resolved here; workers resolve the rest.
            let need = WINDOW.saturating_sub(decoded.tail().len());
            let head_end = &decoded.head[decoded.head.len().saturating_sub(need)..];
            let mut head_bytes = vec![0; head_end.len()];
            resolve(head_end, &window, &mut head_bytes)?;
            let this_window = std::sync::Arc::new(window.clone());
            for part in [&head_bytes[..], decoded.tail()] {
                window.extend_from_slice(&part[part.len().saturating_sub(WINDOW)..]);
                if window.len() > WINDOW {
                    window.drain(..window.len() - WINDOW);
                }
            }
            total += decoded.len() as u64;
            pos = decoded.end;
            if decoded.last {
                finished = Some(pos);
            }
            resolves
                .send(Job::Resolve(dispatched, decoded, this_window))
                .ok();
            dispatched += 1;
            k += 1;
            while let Ok(item) = resolved_rx.try_recv() {
                deliver!(item);
            }
        }
        while next_send < dispatched {
            deliver!(resolved_rx.recv().map_err(|_| workers_stopped())?);
        }
        stop();
        drop((jobs, resolves));
        let end = finished.ok_or_else(|| invalid("Truncated gzip input"))?;
        let trailer_at = end.div_ceil(8);
        if trailer_at + 8 > len {
            return Err(invalid("Truncated gzip input"));
        }
        let trailer = read_range(src, trailer_at, trailer_at + 8)?;
        if u32::from_le_bytes(trailer[..4].try_into().unwrap()) != crc.finalize()
            || u32::from_le_bytes(trailer[4..].try_into().unwrap()) != total as u32
        {
            return Err(invalid("Corrupt gzip member: CRC or length mismatch"));
        }
        // Any further members are decoded sequentially, as MultiGzDecoder would.
        if trailer_at + 8 < len {
            let rest = src.rest(trailer_at + 8)?;
            let mut decoder =
                flate2::bufread::MultiGzDecoder::new(io::BufReader::with_capacity(1 << 20, rest));
            loop {
                let mut chunk = pool.output.1.try_recv().unwrap_or_default();
                chunk.resize(1 << 20, 0);
                let mut filled = 0;
                while filled < chunk.len() {
                    match decoder.read(&mut chunk[filled..]) {
                        Ok(0) => break,
                        Ok(n) => filled += n,
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                        Err(e) => return Err(e),
                    }
                }
                chunk.truncate(filled);
                if filled == 0 {
                    break;
                }
                if sender.send(Ok(chunk)).is_err() {
                    return Ok(());
                }
            }
        }
        Ok(())
    });
    result?;
    let _ = sender.send(Ok(Vec::new()));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Decode a raw deflate stream block by block, as bytes or as symbols
    /// with window placeholders.
    fn inflate_all<T: Symbol>(deflated: &[u8]) -> Result<Vec<T>, Error> {
        let mut inflater = Inflater {
            bits: Bits::at(deflated, 0),
            fixed: None,
            dynamic: Tables::default(),
        };
        let mut out = Vec::new();
        let (mut n, mut last_window) = (0, 0);
        while !inflater.block(&mut out, &mut n, &mut last_window)? {}
        out.truncate(n);
        Ok(out)
    }

    #[test]
    fn inflater_matches_the_input_for_every_block_type_and_survives_damage() {
        use std::io::Write;
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut random = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut inputs: Vec<Vec<u8>> = vec![
            Vec::new(),
            b"a".to_vec(),
            b"abcabcabcabcabcabc".to_vec(),
            vec![b'A'; 70_000],
            fastq(3000, 5),
        ];
        // Incompressible bytes are stored; skewed ones need long codes and
        // subtables; short periods give overlapping copies at every distance.
        inputs.push((0..40_000).map(|_| random() as u8).collect());
        inputs.push(
            (0..120_000)
                .map(|_| (random() % 97).min(random() % 251) as u8)
                .collect(),
        );
        for period in [1usize, 2, 3, 7, 15, 16, 17, 31, 255, 300] {
            let unit: Vec<u8> = (0..period).map(|_| random() as u8).collect();
            inputs.push(unit.iter().cycle().take(5000 + period).copied().collect());
        }
        for input in &inputs {
            for level in [0, 1, 4, 6, 9] {
                let mut encoder =
                    flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::new(level));
                encoder.write_all(input).unwrap();
                let deflated = encoder.finish().unwrap();
                let bytes: Vec<u8> = inflate_all(&deflated).ok().unwrap();
                assert!(
                    bytes == *input,
                    "bytes, level {level}, {} long",
                    input.len()
                );
                let symbols: Vec<u16> = inflate_all(&deflated).ok().unwrap();
                assert!(
                    symbols.iter().map(|&s| s as u8).eq(input.iter().copied())
                        && symbols.iter().all(|&s| s < 256),
                    "symbols, level {level}, {} long",
                    input.len()
                );
                // Damaged and truncated streams may decode to anything or
                // fail, but must not read or write out of bounds.
                for _ in 0..40 {
                    let mut damaged = deflated.clone();
                    if damaged.is_empty() {
                        break;
                    }
                    for _ in 0..1 + random() % 3 {
                        let at = random() as usize % damaged.len();
                        damaged[at] ^= 1 << (random() % 8);
                    }
                    damaged.truncate(1 + random() as usize % damaged.len());
                    let _ = inflate_all::<u8>(&damaged);
                    let _ = inflate_all::<u16>(&damaged);
                }
            }
        }
    }
    use std::io::Write;

    fn fastq(n: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        let mut next = |m: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % m
        };
        let mut out = Vec::new();
        for i in 0..n {
            let len = 50 + next(100) as usize;
            write!(out, "@read.{i} length={len}\n").unwrap();
            out.extend((0..len).map(|_| b"ACGTN"[next(if i % 50 == 0 { 5 } else { 4 }) as usize]));
            out.extend_from_slice(b"\n+\n");
            out.extend((0..len).map(|_| b"#,:FFFF"[next(7) as usize]));
            out.push(b'\n');
        }
        out
    }

    fn decompress_bytes(gz: &[u8], threads: usize) -> io::Result<Vec<u8>> {
        let file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(file.path(), gz).unwrap();
        let (sender, receiver) = crossbeam_channel::unbounded();
        decompress(
            file.path().to_str().unwrap(),
            threads,
            sender,
            crossbeam_channel::bounded(1),
        );
        let mut out = Vec::new();
        for part in receiver {
            let part = part?;
            if part.is_empty() {
                return Ok(out);
            }
            out.extend(part);
        }
        Err(invalid("no end"))
    }

    fn zlib_gzip(data: &[u8], level: u32) -> Vec<u8> {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn regions_reassemble_exactly() {
        // Several regions at every level, including stored (0) and fast (1).
        let data = fastq(160_000, 3);
        for level in [0, 1, 6, 9] {
            let gz = zlib_gzip(&data, level);
            assert!(
                gz.len() > 2 * REGION || level == 9 || level == 6,
                "level {level} spans regions"
            );
            assert_eq!(decompress_bytes(&gz, 4).unwrap(), data, "level {level}");
        }
        let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(4).unwrap());
        let mut gz = vec![0; c.gzip_compress_bound(data.len())];
        let n = c.gzip_compress(&data, &mut gz).unwrap();
        gz.truncate(n);
        assert_eq!(decompress_bytes(&gz, 3).unwrap(), data, "libdeflate");
    }

    /// `PGZIP_BENCH=file.gz cargo test --release pgzip_throughput -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn pgzip_throughput() {
        let path = std::env::var("PGZIP_BENCH").unwrap();
        let counts: Vec<usize> = std::env::var("PGZIP_THREADS")
            .map(|v| v.split(',').map(|x| x.parse().unwrap()).collect())
            .unwrap_or(vec![1, 2, 4, 8]);
        for threads in counts {
            let t = std::time::Instant::now();
            let (sender, receiver) = crossbeam_channel::bounded(4);
            let p = path.clone();
            let worker = std::thread::spawn(move || {
                decompress(&p, threads, sender, crossbeam_channel::bounded(1))
            });
            let mut total = 0;
            for part in receiver {
                let part = part.unwrap();
                if part.is_empty() {
                    break;
                }
                total += part.len();
            }
            worker.join().unwrap();
            eprintln!(
                "threads {threads}: {:.3}s for {} MB",
                t.elapsed().as_secs_f64(),
                total >> 20
            );
        }
    }

    #[test]
    fn damaged_stream_fails_instead_of_decoding_zeros() {
        // Found by fuzzing: a 40-byte damaged member whose block runs past
        // the end of the input. The bits read there are zeros, which decoded
        // as literals forever, until the decoder ran out of memory.
        let gz = [
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0x03, 0x08, 0, 0x46, 0x50, 0x04, 0, 0x28, 0, 0,
            0, 0x05, 0xc1, 0x81, 0x08, 0, 0, 0, 0, 0xa0, 0xfd, 0xa9, 0x2f, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        assert!(crate::testing::pgzip(&gz, 64, 1).is_err());
    }

    #[test]
    fn later_members_and_errors() {
        let data = fastq(160_000, 5);
        let mut gz = zlib_gzip(&data, 1);
        gz.extend(zlib_gzip(b"second member\n", 6));
        assert_eq!(
            decompress_bytes(&gz, 4).unwrap(),
            [&data[..], b"second member\n"].concat()
        );
        let gz = zlib_gzip(&data, 1);
        let mut corrupt = gz.clone();
        let at = corrupt.len() / 2;
        corrupt[at] ^= 0x10;
        assert!(decompress_bytes(&corrupt, 4).is_err(), "corrupt");
        assert!(
            decompress_bytes(&gz[..gz.len() - 3], 4).is_err(),
            "truncated trailer"
        );
        assert!(
            decompress_bytes(&gz[..gz.len() / 2], 4).is_err(),
            "truncated stream"
        );
    }
}
