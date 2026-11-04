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

const PRIMARY: u32 = 10;

/// Two-level decode table. A primary entry is `symbol << 8 | length`, or
/// `offset << 8 | 0x80 | sub_bits` pointing at a subtable for longer codes,
/// whose entries are `symbol << 8 | (length - PRIMARY)`. Length 0 is invalid.
struct Huffman {
    table: Vec<u32>,
}

impl Huffman {
    /// Canonical code for `lens`, or `None` if it is over-subscribed, or
    /// incomplete with more than a single one-bit code (zlib's rule).
    fn new(lens: &[u8]) -> Option<Self> {
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
                return None;
            }
            if c > 0 {
                max = len;
            }
        }
        if left > 0 && max > 1 {
            return None;
        }
        let mut next = [0u32; 16];
        let mut code = 0u32;
        for len in 1..16 {
            code = (code + u32::from(count[len - 1])) << 1;
            next[len] = code;
        }
        let primary = 1usize << PRIMARY;
        let mut table = vec![0u32; primary];
        let reverse = |code: u32, len: usize| code.reverse_bits() >> (32 - len);
        let mut long = Vec::new();
        for (sym, &l) in lens.iter().enumerate() {
            let len = usize::from(l);
            if len == 0 {
                continue;
            }
            let rev = reverse(next[len], len) as usize;
            next[len] += 1;
            if len <= PRIMARY as usize {
                for j in (rev..primary).step_by(1 << len) {
                    table[j] = (sym as u32) << 8 | len as u32;
                }
            } else {
                long.push((sym, len, rev));
            }
        }
        if !long.is_empty() {
            let mut sub_len = vec![0usize; primary];
            for &(_, len, rev) in &long {
                let p = rev & (primary - 1);
                sub_len[p] = sub_len[p].max(len - PRIMARY as usize);
            }
            for (p, &bits) in sub_len.iter().enumerate() {
                if bits > 0 {
                    table[p] = (table.len() as u32) << 8 | 0x80 | bits as u32;
                    table.resize(table.len() + (1 << bits), 0);
                }
            }
            for (sym, len, rev) in long {
                let p = rev & (primary - 1);
                let bits = sub_len[p];
                let offset = (table[p] >> 8) as usize;
                let extra = len - PRIMARY as usize;
                for j in ((rev >> PRIMARY)..(1 << bits)).step_by(1 << extra) {
                    table[offset + j] = (sym as u32) << 8 | extra as u32;
                }
            }
        }
        Some(Self { table })
    }

    /// Needs at least 15 bits in `bits`.
    #[inline(always)]
    fn decode(&self, bits: &mut Bits) -> Option<u32> {
        let e = self.table[(bits.buf & ((1 << PRIMARY) - 1)) as usize];
        if e & 0x80 == 0 {
            let len = e & 15;
            if len == 0 {
                return None;
            }
            bits.consume(len);
            Some(e >> 8)
        } else {
            let sub = ((bits.buf >> PRIMARY) & ((1 << (e & 15)) - 1)) as usize;
            let e2 = self.table[(e >> 8) as usize + sub];
            let len = e2 & 15;
            if len == 0 {
                return None;
            }
            bits.consume(PRIMARY + len);
            Some(e2 >> 8)
        }
    }
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

fn fixed_tables() -> (Huffman, Huffman) {
    let mut lens = [0u8; 288];
    lens[..144].fill(8);
    lens[144..256].fill(9);
    lens[256..280].fill(7);
    lens[280..].fill(8);
    (
        Huffman::new(&lens).unwrap(),
        Huffman::new(&[5u8; 32]).unwrap(),
    )
}

/// The litlen and distance codes of a dynamic block, after its 3 header bits.
fn dynamic_tables(bits: &mut Bits) -> Option<(Huffman, Huffman)> {
    let hlit = bits.take(5) as usize + 257;
    let hdist = bits.take(5) as usize + 1;
    let hclen = bits.take(4) as usize + 4;
    if hlit > 286 || hdist > 30 {
        return None;
    }
    let mut code_lens = [0u8; 19];
    for &i in &CODE_ORDER[..hclen] {
        code_lens[i] = bits.take(3) as u8;
    }
    let codes = Huffman::new(&code_lens)?;
    let mut lens = [0u8; 286 + 30];
    let total = hlit + hdist;
    let mut i = 0;
    while i < total {
        bits.refill();
        let sym = codes.decode(bits)?;
        let (value, repeat) = match sym {
            0..=15 => (sym as u8, 1),
            16 => (*lens[..i].last()?, 3 + bits.take(2) as usize),
            17 => (0, 3 + bits.take(3) as usize),
            _ => (0, 11 + bits.take(7) as usize),
        };
        if i + repeat > total {
            return None;
        }
        lens[i..i + repeat].fill(value);
        i += repeat;
    }
    if lens[256] == 0 {
        return None;
    }
    Some((
        Huffman::new(&lens[..hlit])?,
        Huffman::new(&lens[hlit..total])?,
    ))
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
    fixed: Option<(Huffman, Huffman)>,
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
        let dynamic;
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
                (&fixed.0, &fixed.1)
            }
            2 => {
                dynamic = dynamic_tables(bits).ok_or(Error::Bad)?;
                (&dynamic.0, &dynamic.1)
            }
            _ => return Err(Error::Bad),
        };
        loop {
            // Room for three literals and the longest copy, which may write up
            // to a chunk past its end.
            if *n + 3 + 258 + CHUNK > out.len() {
                grow(out, *n + 3 + 258 + CHUNK);
            }
            // A refill leaves at least 56 bits: three literal codes of up to
            // 15 bits each, or one length, distance and their extra bits.
            bits.refill();
            let mut sym = litlen.decode(bits).ok_or(Error::Bad)?;
            if sym < 256 {
                out[*n] = T::byte(sym as u8);
                *n += 1;
                sym = litlen.decode(bits).ok_or(Error::Bad)?;
                if sym < 256 {
                    out[*n] = T::byte(sym as u8);
                    *n += 1;
                    sym = litlen.decode(bits).ok_or(Error::Bad)?;
                    if sym < 256 {
                        out[*n] = T::byte(sym as u8);
                        *n += 1;
                        continue;
                    }
                }
                if bits.n < 48 {
                    bits.refill();
                }
            }
            if sym == 256 {
                break;
            }
            let s = (sym - 257) as usize;
            if s >= 29 {
                return Err(Error::Bad);
            }
            let len = usize::from(LEN_BASE[s]) + bits.take(u32::from(LEN_EXTRA[s])) as usize;
            let d = dist.decode(bits).ok_or(Error::Bad)? as usize;
            if d >= 30 {
                return Err(Error::Bad);
            }
            let distance = usize::from(DIST_BASE[d]) + bits.take(u32::from(DIST_EXTRA[d])) as usize;
            let p = *n;
            // Copying from before the last placeholder may copy placeholders.
            if p.saturating_sub(distance) < *last_window {
                *last_window = p + len;
            }
            if distance <= p {
                let from = p - distance;
                if distance >= CHUNK {
                    // Copy whole chunks, which compile to vector loads and
                    // stores, and may write past `len` into spare buffer.
                    // Chunks never overlap the bytes they read.
                    let mut i = 0;
                    while i < len {
                        let chunk: [T; CHUNK] = out[from + i..from + i + CHUNK].try_into().unwrap();
                        out[p + i..p + i + CHUNK].copy_from_slice(&chunk);
                        i += CHUNK;
                    }
                } else if distance >= len {
                    out.copy_within(from..from + len, p);
                } else if distance == 1 {
                    let v = out[from];
                    out[p..p + len].fill(v);
                } else {
                    for i in 0..len {
                        out[p + i] = out[from + i];
                    }
                }
            } else {
                for i in 0..len {
                    let at = (p + i) as isize - distance as isize;
                    out[p + i] = if at >= 0 {
                        out[at as usize]
                    } else {
                        T::window((WINDOW as isize + at) as usize).ok_or(Error::Bad)?
                    };
                }
                *last_window = p + len;
            }
            *n += len;
            if bits.overrun() {
                return Err(Error::More);
            }
        }
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

/// Buffers kept for reuse between regions. Every region needs several
/// megabytes of buffers; allocating them fresh faults in new zeroed pages, and
/// freeing them unmaps pages, which interrupts every core.
struct Pool {
    bytes: (
        crossbeam_channel::Sender<Vec<u8>>,
        crossbeam_channel::Receiver<Vec<u8>>,
    ),
    words: (
        crossbeam_channel::Sender<Vec<u16>>,
        crossbeam_channel::Receiver<Vec<u16>>,
    ),
    /// Output buffers the reader is done with.
    returned: crossbeam_channel::Receiver<Vec<u8>>,
}

impl Pool {
    fn new(size: usize, returned: crossbeam_channel::Receiver<Vec<u8>>) -> Self {
        Self {
            bytes: crossbeam_channel::bounded(size),
            words: crossbeam_channel::bounded(size),
            returned,
        }
    }
    fn bytes(&self) -> Vec<u8> {
        let mut v = (self
            .bytes
            .1
            .try_recv()
            .or_else(|_| self.returned.try_recv()))
        .unwrap_or_default();
        v.clear();
        v
    }
    fn words(&self) -> Vec<u16> {
        let mut v = self.words.1.try_recv().unwrap_or_default();
        v.clear();
        v
    }
    fn give_bytes(&self, v: Vec<u8>) {
        let _ = self.bytes.0.try_send(v);
    }
    fn give_words(&self, v: Vec<u16>) {
        let _ = self.words.0.try_send(v);
    }
}

/// A region decoded from a block start up to a block boundary.
struct Decoded {
    start: usize,
    end: usize,
    /// The member ended here (the final block was decoded).
    last: bool,
    /// Output that may hold window placeholders.
    head: Vec<u16>,
    /// Output after `head`, plain bytes, and its CRC32.
    tail: Vec<u8>,
    tail_crc: u32,
}

/// Decode from bit `start` of `data` until a block boundary at or after bit
/// `stop`, or the final block. `first` means the stream starts at `start`, so
/// there is no window and copies may not reach before it.
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
    };
    let mut head = pool.words();
    let mut tail = pool.bytes();
    let (mut head_n, mut tail_n) = (0, 0);
    let mut in_tail = first;
    let mut last_window = 0;
    let mut last = false;
    while inflater.bits.bit_pos() < stop || inflater.bits.bit_pos() == start {
        if in_tail {
            last = inflater.block(&mut tail, &mut tail_n, &mut last_window)?;
        } else {
            last = inflater.block(&mut head, &mut head_n, &mut last_window)?;
            if head_n - last_window.min(head_n) >= WINDOW {
                // No later copy can reach a placeholder: continue in bytes,
                // seeded with the last 32 KiB as history.
                tail.extend(head[head_n - WINDOW..head_n].iter().map(|&s| s as u8));
                tail_n = WINDOW;
                in_tail = true;
            }
        }
        if last {
            break;
        }
    }
    head.truncate(head_n);
    tail.truncate(tail_n);
    let seeded = if in_tail && !first { WINDOW } else { 0 };
    tail.drain(..seeded);
    let tail_crc = crc32fast::hash(&tail);
    Ok(Decoded {
        start,
        end: inflater.bits.bit_pos(),
        last,
        head,
        tail,
        tail_crc,
    })
}

/// The first bit at or after `from` and before `to` where a non-final
/// dynamic block starts whose header is valid and whose whole block decodes.
fn find_block(data: &[u8], from: usize, to: usize) -> Option<usize> {
    let to = to.min(data.len().saturating_sub(8) * 8);
    let mut scratch: Vec<u16> = Vec::new();
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
        if Huffman::new(&lens).is_none() {
            continue;
        }
        let mut inflater = Inflater {
            bits: Bits::at(data, bit),
            fixed: None,
        };
        let (mut n, mut last_window) = (0, 0);
        if inflater
            .block(&mut scratch, &mut n, &mut last_window)
            .is_ok()
        {
            return Some(bit);
        }
    }
    None
}

// driver

#[cfg(unix)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    std::os::unix::fs::FileExt::read_exact_at(file, buf, offset)
}

#[cfg(windows)]
fn read_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> io::Result<()> {
    while !buf.is_empty() {
        match std::os::windows::fs::FileExt::seek_read(file, buf, offset)? {
            0 => return Err(io::ErrorKind::UnexpectedEof.into()),
            n => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
        }
    }
    Ok(())
}

/// Bytes `from..to` of the file.
fn read_range(file: &File, from: usize, to: usize) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    read_range_into(file, from, to, &mut buf)?;
    Ok(buf)
}

/// Bytes `from..to` of the file into `buf`, reusing its capacity.
fn read_range_into(file: &File, from: usize, to: usize, buf: &mut Vec<u8>) -> io::Result<()> {
    buf.resize(to - from, 0);
    read_at(file, buf, from as u64)
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

/// Bit offset of the deflate stream after a gzip member header at `data[0]`.
fn header_len(data: &[u8]) -> io::Result<usize> {
    let bad = || invalid("Invalid gzip header");
    if data.len() < 10 || data[..3] != [0x1f, 0x8b, 8] {
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
    file: &File,
    len: usize,
    k: usize,
    deflate: usize,
    start: Option<usize>,
    pool: &Pool,
) -> io::Result<Option<Decoded>> {
    let begin = deflate / 8 + k * REGION;
    let end = (begin + REGION).min(len);
    let from = start.map_or(begin, |s| s / 8);
    let mut margin = MARGIN;
    let mut data = pool.bytes();
    let result = loop {
        let read_to = (end + margin).min(len);
        if let Err(e) = read_range_into(file, from, read_to, &mut data) {
            break Err(e);
        }
        let base = from * 8;
        let local_start = match start {
            Some(s) => s - base,
            None if k == 0 => deflate - base,
            None => match find_block(&data, 0, (end - from) * 8) {
                Some(bit) => bit,
                None => break Ok(None),
            },
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
    pool.give_bytes(data);
    result
}

enum Job {
    Decode(usize),
    /// Fill a decoded region's placeholders from the 32 KiB before it.
    Resolve(usize, Decoded, std::sync::Arc<Vec<u8>>),
}

/// Bytes for `head`, filling placeholders from `window`, the bytes before it.
fn resolve(head: &[u16], window: &[u8], mut out: Vec<u8>) -> io::Result<Vec<u8>> {
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
    out.clear();
    out.extend(head.iter().map(|&s| lut[usize::from(s)]));
    Ok(out)
}

/// Whether a gzip file of `len` bytes spans enough regions to decode in parallel.
pub(crate) fn worth_it(len: u64) -> bool {
    len >= 3 * REGION as u64
}

/// Decompress the gzip file at `path` on `threads` workers, sending output in
/// order on `sender`, then an empty buffer at the end. Members after the first
/// are decoded sequentially.
/// The reader hands buffers it is done with back on `returned`.
pub(crate) fn decompress(
    path: &str,
    threads: usize,
    sender: crossbeam_channel::Sender<io::Result<Vec<u8>>>,
    returned: crossbeam_channel::Receiver<Vec<u8>>,
) {
    if let Err(e) = run(path, threads, &sender, returned) {
        let _ = sender.send(Err(e));
    }
}

fn run(
    path: &str,
    threads: usize,
    sender: &crossbeam_channel::Sender<io::Result<Vec<u8>>>,
    returned: crossbeam_channel::Receiver<Vec<u8>>,
) -> io::Result<()> {
    let file = File::open(path)?;
    let len = file.metadata()?.len() as usize;
    let deflate = header_len(&read_range(&file, 0, len.min(64 << 10))?)?;
    let regions = (len - deflate / 8).div_ceil(REGION).max(1);
    let file = &file;
    // Stops the workers early if the reader goes away.
    let stopped = std::sync::atomic::AtomicBool::new(false);
    let stopped = &stopped;
    // Enough for the regions in flight and the reader's buffers.
    let pool = &Pool::new(4 * (threads + 2), returned);
    let result = std::thread::scope(|scope| -> io::Result<()> {
        // Regions decode to about 3x their size; keep a few in flight per worker.
        let ahead = threads + 2;
        let (jobs, queue) = crossbeam_channel::bounded::<Job>(2 * ahead);
        let (decoded_tx, decoded_rx) = crossbeam_channel::unbounded();
        let (resolved_tx, resolved_rx) = crossbeam_channel::unbounded();
        for _ in 0..threads {
            let (queue, decoded_tx, resolved_tx) =
                (queue.clone(), decoded_tx.clone(), resolved_tx.clone());
            scope.spawn(move || {
                for job in queue {
                    if stopped.load(std::sync::atomic::Ordering::Relaxed) {
                        break;
                    }
                    match job {
                        // A failure here may only be a false block start, so
                        // the coordinator decides; it re-decodes from known output.
                        Job::Decode(k) => {
                            let _ = decoded_tx.send((k, region(file, len, k, deflate, None, pool)));
                        }
                        Job::Resolve(k, mut decoded, window) => {
                            let head = resolve(&decoded.head, &window, pool.bytes());
                            pool.give_words(std::mem::take(&mut decoded.head));
                            let crc = head.as_ref().map(|h| crc32fast::hash(h)).unwrap_or(0);
                            let _ = resolved_tx.send((k, head, crc, decoded));
                        }
                    }
                }
            });
        }
        drop((queue, decoded_tx, resolved_tx));
        let mut queued = 0;
        while queued < regions.min(ahead) {
            jobs.send(Job::Decode(queued)).ok();
            queued += 1;
        }
        let mut decoded_pending = std::collections::BTreeMap::new();
        let mut resolved_pending = std::collections::BTreeMap::new();
        let mut window: Vec<u8> = Vec::with_capacity(2 * WINDOW);
        let mut crc = crc32fast::Hasher::new();
        let mut total = 0u64;
        let mut pos = deflate;
        // Resolve jobs are numbered in dispatch order: skipped regions have none.
        let (mut k, mut resolving, mut next_send, mut dispatched) = (0, 0, 0, 0);
        let mut finished = None;
        // Send resolved regions in order; with `wait`, block for at least one.
        let mut flush = |wait: bool,
                         resolving: &mut usize,
                         resolved_pending: &mut std::collections::BTreeMap<_, _>,
                         crc: &mut crc32fast::Hasher|
         -> io::Result<bool> {
            let mut got = !wait;
            loop {
                let next = if got {
                    resolved_rx.try_recv().ok()
                } else {
                    Some(
                        resolved_rx
                            .recv()
                            .map_err(|_| invalid("Decompression workers stopped"))?,
                    )
                };
                let Some((i, head, head_crc, decoded)) = next else {
                    return Ok(true);
                };
                got = true;
                resolved_pending.insert(i, (head, head_crc, decoded));
                while let Some((head, head_crc, decoded)) = resolved_pending.remove(&next_send) {
                    let head: Vec<u8> = head?;
                    let decoded: Decoded = decoded;
                    crc.combine(&crc32fast::Hasher::new_with_initial_len(
                        head_crc,
                        head.len() as u64,
                    ));
                    crc.combine(&crc32fast::Hasher::new_with_initial_len(
                        decoded.tail_crc,
                        decoded.tail.len() as u64,
                    ));
                    for part in [head, decoded.tail] {
                        if !part.is_empty() && sender.send(Ok(part)).is_err() {
                            return Ok(false);
                        }
                    }
                    next_send += 1;
                    *resolving -= 1;
                }
            }
        };
        while k < regions && finished.is_none() {
            let result = loop {
                if let Some(r) = decoded_pending.remove(&k) {
                    break r;
                }
                let (i, r) = decoded_rx
                    .recv()
                    .map_err(|_| invalid("Decompression workers stopped"))?;
                decoded_pending.insert(i, r);
            };
            if queued < regions {
                jobs.send(Job::Decode(queued)).ok();
                queued += 1;
            }
            let region_end = (deflate / 8 + (k + 1) * REGION).min(len) * 8;
            if pos >= region_end {
                // An earlier region's last block covered this one.
                k += 1;
                continue;
            }
            let decoded = match result {
                Ok(Some(d)) if d.start == pos => d,
                _ => region(file, len, k, deflate, Some(pos), pool)?.unwrap(),
            };
            // The next region's window is this region's last 32 KiB. Only that
            // much of the head is resolved here; workers resolve the rest.
            let need = WINDOW.saturating_sub(decoded.tail.len());
            let head_end = resolve(
                &decoded.head[decoded.head.len().saturating_sub(need)..],
                &window,
                Vec::with_capacity(need),
            )?;
            let this_window = std::sync::Arc::new(window.clone());
            for part in [&head_end[..], &decoded.tail[..]] {
                window.extend_from_slice(&part[part.len().saturating_sub(WINDOW)..]);
                if window.len() > WINDOW {
                    window.drain(..window.len() - WINDOW);
                }
            }
            total += (decoded.head.len() + decoded.tail.len()) as u64;
            pos = decoded.end;
            if decoded.last {
                finished = Some(pos);
            }
            jobs.send(Job::Resolve(dispatched, decoded, this_window))
                .ok();
            dispatched += 1;
            resolving += 1;
            if !flush(
                resolving > ahead,
                &mut resolving,
                &mut resolved_pending,
                &mut crc,
            )? {
                stopped.store(true, std::sync::atomic::Ordering::Relaxed);
                return Ok(());
            }
            k += 1;
        }
        while resolving > 0 {
            if !flush(true, &mut resolving, &mut resolved_pending, &mut crc)? {
                stopped.store(true, std::sync::atomic::Ordering::Relaxed);
                return Ok(());
            }
        }
        stopped.store(true, std::sync::atomic::Ordering::Relaxed);
        drop(jobs);
        let end = finished.ok_or_else(|| invalid("Truncated gzip input"))?;
        let trailer_at = end.div_ceil(8);
        if trailer_at + 8 > len {
            return Err(invalid("Truncated gzip input"));
        }
        let trailer = read_range(file, trailer_at, trailer_at + 8)?;
        if u32::from_le_bytes(trailer[..4].try_into().unwrap()) != crc.finalize()
            || u32::from_le_bytes(trailer[4..].try_into().unwrap()) != total as u32
        {
            return Err(invalid("Corrupt gzip member: CRC or length mismatch"));
        }
        // Any further members are decoded sequentially, as MultiGzDecoder would.
        if trailer_at + 8 < len {
            let mut rest = File::open(path)?;
            io::Seek::seek(&mut rest, io::SeekFrom::Start((trailer_at + 8) as u64))?;
            let mut decoder =
                flate2::bufread::MultiGzDecoder::new(io::BufReader::with_capacity(1 << 20, rest));
            loop {
                let mut chunk = vec![0; 1 << 20];
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
            crossbeam_channel::never(),
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
                decompress(&p, threads, sender, crossbeam_channel::never())
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
