//! Adapter evaluation ported from OpenGene/fastp v1.4.0 evaluator.cpp and
//! nucleotidetree.cpp. Upstream copyright and MIT terms: LICENSE.fastp.
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::sync::Mutex;

#[derive(Default)]
struct SeedCandidates {
    positions: Vec<usize>,
    dense: bool,
}

struct KnownSeeds {
    index: rustc_hash::FxHashMap<u64, Vec<(usize, usize)>>,
    occupied: Box<[bool; 8192]>,
    candidates: Vec<SeedCandidates>,
    short: Vec<usize>,
    active: Vec<usize>,
}

#[inline]
fn seed_bucket<const BITS: u32>(key: u64) -> usize {
    (key.wrapping_mul(0x9e3779b97f4a7c15) >> (64 - BITS)) as usize
}

impl KnownSeeds {
    fn new(adapters: &[&[u8]]) -> Self {
        let mut index: rustc_hash::FxHashMap<u64, Vec<(usize, usize)>> = Default::default();
        for (adapter, seq) in adapters.iter().enumerate() {
            if seq.len() < 8 {
                continue;
            }
            // With floor(L/16) mismatches, one of the first floor(L/16)+1
            // disjoint 8-byte seeds must be exact. Full-adapter seeds are a
            // superset of the required seeds for every valid partial tail L>=9.
            for seed in 0..=seq.len() / 16 {
                let offset = seed * 8;
                let key = u64::from_le_bytes(seq[offset..offset + 8].try_into().unwrap());
                index.entry(key).or_default().push((adapter, offset));
            }
        }
        let mut occupied = Box::new([false; 8192]);
        for &key in index.keys() {
            occupied[seed_bucket::<13>(key)] = true;
        }
        Self {
            index,
            occupied,
            candidates: (0..adapters.len())
                .map(|_| SeedCandidates::default())
                .collect(),
            short: (0..adapters.len())
                .filter(|&i| adapters[i].len() < 8)
                .collect(),
            active: Vec::new(),
        }
    }

    fn collect(&mut self, read: &[u8]) {
        for i in self.active.drain(..) {
            let c = &mut self.candidates[i];
            c.positions.clear();
            c.dense = false;
        }
        for &i in &self.short {
            self.candidates[i].dense = true;
        }
        self.active.extend_from_slice(&self.short);
        // Scan each read once for all known adapters. Exact 8-byte keys retain
        // case and ambiguous bytes; no DNA normalization or lossy hash is used.
        for (found, bytes) in read.windows(8).enumerate() {
            let key = u64::from_le_bytes(bytes.try_into().unwrap());
            // A small presence table rejects most windows before probing the
            // full catalog. Bucket collisions still get the exact-key lookup.
            if !self.occupied[seed_bucket::<13>(key)] {
                continue;
            }
            if let Some(seeds) = self.index.get(&key) {
                for &(adapter, offset) in seeds {
                    let c = &mut self.candidates[adapter];
                    if !c.dense && found >= offset && found - offset < read.len().saturating_sub(8)
                    {
                        if c.positions.is_empty() {
                            self.active.push(adapter);
                        }
                        c.positions.push(found - offset);
                        // Avoid unbounded candidate lists on repeated sequence.
                        // The original ordered search is cheap for dense hits.
                        if c.positions.len() > 64 {
                            c.dense = true;
                        }
                    }
                }
            }
        }
        // The leader count can change during the scan of a single read.
        // Visit only possible matches, in the same catalog order as fastp.
        self.active.sort_unstable();
    }

    fn first_match(&mut self, index: usize, read: &[u8], adapter: &[u8]) -> Option<usize> {
        let at = |pos: usize| {
            let len = (read.len() - pos).min(adapter.len());
            let limit = len / 16;
            let (diff, _) =
                crate::simd::count_mismatches(&read[pos..pos + len], &adapter[..len], limit, len);
            (diff <= limit).then_some(diff)
        };
        let c = &mut self.candidates[index];
        if c.dense {
            return (0..read.len().saturating_sub(8)).find_map(at);
        }
        // Preserve the first accepted offset and its mismatch count, including
        // duplicate seeds and matches discovered at different seed offsets.
        c.positions.sort_unstable();
        c.positions.dedup();
        c.positions.iter().copied().find_map(at)
    }
}

fn known(reads: &[&[u8]]) -> Option<Vec<u8>> {
    // Whether a read matches an adapter does not depend on the other reads, so
    // find every match in parallel, in the order the serial loop visits them,
    // then replay fastp's loop (which skips adapters that fall behind and
    // stops early) over those results. Only the first 100000 reads can count.
    let matches: Vec<Vec<(usize, usize)>> = reads[..reads.len().min(100000)]
        .par_iter()
        .map_init(
            || KnownSeeds::new(KNOWN),
            |seeds, r| {
                seeds.collect(r);
                let mut found = Vec::new();
                for j in 0..seeds.active.len() {
                    let i = seeds.active[j];
                    let adapter = KNOWN[i];
                    if adapter.len() < r.len() {
                        if let Some(diff) = seeds.first_match(i, r, adapter) {
                            found.push((i, diff));
                        }
                    }
                }
                found
            },
        )
        .collect();
    let mut counts = vec![0usize; KNOWN.len()];
    let mut mismatches = vec![0usize; KNOWN.len()];
    let mut checked = 0;
    let mut bases = 0;
    let mut max = 0;
    for (r, found) in reads.iter().zip(&matches) {
        checked += 1;
        bases += r.len();
        if checked > 100000 || bases > 100000000 || max > 1000 {
            break;
        }
        for &(i, diff) in found {
            if max > 20 && counts[i] < max / 10 {
                continue;
            }
            counts[i] += 1;
            max = max.max(counts[i]);
            mismatches[i] += diff;
        }
    }
    let mut best = 0;
    for i in 1..counts.len() {
        if counts[i] > counts[best] {
            best = i;
        }
    }
    (counts[best] > checked / 50 || (counts[best] > checked / 200 && mismatches[best] < checked))
        .then(|| KNOWN[best].to_vec())
}

fn count_seeds<T: From<u8> + std::ops::AddAssign>(read: &[u8], shift: usize, counts: &mut [T]) {
    each_seed(read, shift, |key| counts[key] += T::from(1));
}

/// Call `seen` with the key of every ten-base seed that discovery counts.
#[inline(always)]
fn each_seed(read: &[u8], shift: usize, mut seen: impl FnMut(usize)) {
    let end = read.len().saturating_sub(10 + shift).saturating_add(1);
    if end <= 20 {
        return;
    }
    // A table lookup: branching on each base mispredicts on most of them.
    const CODES: [u8; 256] = {
        let mut codes = [4; 256];
        codes[b'A' as usize] = 0;
        codes[b'T' as usize] = 1;
        codes[b'C' as usize] = 2;
        codes[b'G' as usize] = 3;
        codes
    };
    let mut key = 0;
    let mut valid = 0;
    for &b in &read[20..end + 9] {
        let base = CODES[b as usize] as usize;
        if base > 3 {
            key = 0;
            valid = 0;
            continue;
        }
        key = ((key << 2) | base) & ((1 << 20) - 1);
        valid += 1;
        if valid >= 10 {
            seen(key);
        }
    }
}

/// Seed counts for one worker. A megabyte of one-byte counters stays in
/// cache where a table of full counters per worker does not; the few seeds
/// seen more than 255 times record each wrap, so counts remain exact.
#[derive(Default)]
struct SeedCounts {
    low: Vec<u8>,
    wraps: FxHashMap<usize, u32>,
}

impl SeedCounts {
    #[inline(always)]
    fn add(&mut self, key: usize) {
        let count = &mut self.low[key];
        *count = count.wrapping_add(1);
        if *count == 0 {
            *self.wraps.entry(key).or_default() += 1;
        }
    }
}

// A radix-tree edge borrows the first observed spelling of an unbranched
// path. Split at both forks and sequence ends so every base on an edge has
// the same support. Index zero is the root and cannot be a child.
#[derive(Default)]
struct Node<'a> {
    children: [usize; 8],
    count: usize,
    label: &'a [u8],
}

struct Tree<'a, const REVERSE: bool> {
    nodes: Vec<Node<'a>>,
}

impl<'a, const REVERSE: bool> Tree<'a, REVERSE> {
    fn new() -> Self {
        Self {
            nodes: vec![Node::default()],
        }
    }

    fn slot(seq: &[u8]) -> usize {
        let b = if REVERSE { seq[seq.len() - 1] } else { seq[0] };
        (b & 7) as usize
    }

    fn split(seq: &'a [u8], n: usize) -> (&'a [u8], &'a [u8]) {
        if REVERSE {
            let (rest, prefix) = seq.split_at(seq.len() - n);
            (prefix, rest)
        } else {
            seq.split_at(n)
        }
    }

    fn add(&mut self, seq: &'a [u8]) {
        let mut seq = if REVERSE {
            memchr::memrchr(b'N', seq).map_or(seq, |n| &seq[n + 1..])
        } else {
            &seq[..memchr::memchr(b'N', seq).unwrap_or(seq.len())]
        };
        let mut parent = 0;
        while !seq.is_empty() {
            let slot = Self::slot(seq);
            let child = self.nodes[parent].children[slot];
            if child == 0 {
                self.nodes[parent].children[slot] = self.nodes.len();
                self.nodes.push(Node {
                    label: seq,
                    count: 1,
                    ..Node::default()
                });
                return;
            }
            let label = self.nodes[child].label;
            // The upstream tree groups bases by their low three bits.
            let same = |(a, b): &(&u8, &u8)| (*a ^ *b) & 7 == 0;
            let shared = if REVERSE {
                label
                    .iter()
                    .rev()
                    .zip(seq.iter().rev())
                    .take_while(same)
                    .count()
            } else {
                label.iter().zip(seq).take_while(same).count()
            };
            if shared < label.len() {
                let (prefix, rest) = Self::split(label, shared);
                let suffix = self.nodes.len();
                let old = std::mem::take(&mut self.nodes[child]);
                self.nodes[child].label = prefix;
                self.nodes[child].count = old.count;
                self.nodes[child].children[Self::slot(rest)] = suffix;
                self.nodes.push(Node { label: rest, ..old });
            }
            self.nodes[child].count += 1;
            parent = child;
            seq = Self::split(seq, shared).1;
        }
    }

    fn path(&self, leaf: &mut bool) -> Vec<u8> {
        let mut node = &self.nodes[0];
        let mut path = Vec::new();
        loop {
            let mut children = node
                .children
                .iter()
                .filter(|&&i| i != 0)
                .map(|&i| &self.nodes[i]);
            let total: usize = children.clone().map(|child| child.count).sum();
            if total < 50 {
                break;
            }
            let Some(next) = children.find(|c| c.count as f64 / total as f64 >= 0.95) else {
                *leaf = false;
                break;
            };
            // The first base qualifies against the parent's total. Continuing
            // down a unary edge requires at least 50 observations of that edge.
            let label = if next.count < 50 {
                Self::split(next.label, 1).0
            } else {
                next.label
            };
            if REVERSE {
                path.extend(label.iter().rev());
            } else {
                path.extend_from_slice(label);
            }
            if next.count < 50 {
                break;
            }
            node = next;
        }
        path
    }
}

// Index sparse hits in input order. The frequencies from the first pass are
// upper bounds since extension only inspects the first 500 start positions.
// Dense seeds stay on the streaming search, which stores no hit list.
fn index_seed_hits(
    reads: &[&[u8]],
    shift: usize,
    seeds: &[(Vec<u8>, usize)],
) -> Option<Vec<Vec<(usize, usize)>>> {
    if seeds.len() < 2 || seeds.iter().map(|(_, n)| n).sum::<usize>() > 256 * 1024 {
        return None;
    }
    fn bucket(window: &[u8]) -> usize {
        let prefix = u64::from_le_bytes(window[..8].try_into().unwrap());
        seed_bucket::<10>(prefix)
    }
    let mut buckets = [0u16; 1024];
    for (rank, (seed, _)) in seeds.iter().enumerate() {
        buckets[bucket(seed)] |= 1 << rank;
    }
    // Each thread scans a contiguous share of the sample; concatenating the
    // shares in order keeps every hit list in input order.
    let share = reads.len().div_ceil(rayon::current_num_threads()).max(1);
    let parts: Vec<Vec<Vec<(usize, usize)>>> = reads
        .par_chunks(share)
        .enumerate()
        .map(|(c, chunk)| {
            let mut hits = vec![Vec::new(); seeds.len()];
            for (k, &read) in chunk.iter().enumerate() {
                let end = 500.min(read.len().saturating_sub(10 + shift).saturating_add(1));
                if end <= 20 {
                    continue;
                }
                for (offset, window) in read[20..end + 9].windows(10).enumerate() {
                    let mut candidates = buckets[bucket(window)];
                    while candidates != 0 {
                        let rank = candidates.trailing_zeros() as usize;
                        if window == seeds[rank].0 {
                            hits[rank].push((c * share + k, offset + 20));
                        }
                        candidates &= candidates - 1;
                    }
                }
            }
            hits
        })
        .collect();
    let mut hits: Vec<Vec<_>> = seeds.iter().map(|(_, n)| Vec::with_capacity(*n)).collect();
    for part in parts {
        for (all, some) in hits.iter_mut().zip(part) {
            all.extend(some);
        }
    }
    Some(hits)
}

fn seed_histogram<T: From<u8> + std::ops::AddAssign + Copy>(
    reads: &[&[u8]],
    shift: usize,
) -> Vec<T> {
    let mut counts = vec![T::from(0); 1 << 20];
    for r in reads {
        count_seeds(r, shift, &mut counts);
    }
    counts
}

/// Candidate seeds among `counts`, whose first key is `first`: the total of
/// the eligible counts and the ten highest (count, key) pairs.
fn seed_candidates(
    first: usize,
    counts: impl Iterator<Item = usize>,
) -> (usize, Vec<(usize, usize)>) {
    let mut candidates = BinaryHeap::with_capacity(10);
    let mut total = 0;
    for (key, n) in counts.enumerate().map(|(i, n)| (first + i, n)) {
        if n == 0 || key == 0 {
            continue;
        }
        let mut composition = [0; 4];
        for i in 0..10 {
            composition[(key >> (i * 2)) & 3] += 1;
        }
        if composition.iter().any(|&n| n >= 6)
            || composition[2] + composition[3] >= 8
            || key >> 12 == 255
        {
            continue;
        }
        total += n;
        if candidates.len() < 10 {
            candidates.push(Reverse((n, key)));
        } else if (n, key) > candidates.peek().unwrap().0 {
            candidates.pop();
            candidates.push(Reverse((n, key)));
        }
    }
    let candidates = candidates
        .into_iter()
        .map(|Reverse((n, key))| (key, n))
        .collect();
    (total, candidates)
}

fn rank_seeds(counts: impl Iterator<Item = usize>) -> Vec<(Vec<u8>, usize)> {
    let (total, candidates) = seed_candidates(0, counts);
    seeds_from(total, candidates)
}

/// Only the ten highest (count, key) pairs are used by discovery. Keys are
/// distinct, so the ten highest overall are the ten highest of any split's
/// per-part candidates.
fn seeds_from(total: usize, mut candidates: Vec<(usize, usize)>) -> Vec<(Vec<u8>, usize)> {
    candidates.sort_unstable_by_key(|&(key, n)| (std::cmp::Reverse(n), std::cmp::Reverse(key)));
    let mut seeds = Vec::new();
    for (key, n) in candidates.into_iter().take(10) {
        if n < 10 || n * (1 << 20) < total * 20 {
            break;
        }
        let seed: Vec<_> = (0..10)
            .rev()
            .map(|i| b"ATCG"[(key >> (i * 2)) & 3])
            .collect();
        if seed.windows(2).filter(|w| w[0] != w[1]).count() < 3 {
            continue;
        }
        seeds.push((seed, n));
    }
    seeds
}

/// `rank_seeds` over the sum of `tables`, summing and ranking ranges of keys
/// in parallel.
fn rank_tables(tables: &[SeedCounts]) -> Vec<(Vec<u8>, usize)> {
    const RANGE: usize = 1 << 14;
    let mut wraps: Vec<(usize, u32)> = tables
        .iter()
        .flat_map(|t| t.wraps.iter().map(|(&key, &n)| (key, n * 256)))
        .collect();
    wraps.sort_unstable();
    let parts: Vec<(usize, Vec<(usize, usize)>)> = (0..(1 << 20) / RANGE)
        .into_par_iter()
        .map(|part| {
            let first = part * RANGE;
            let mut sum = [0u32; RANGE];
            for table in tables {
                let counts = &table.low[first..][..RANGE];
                sum.iter_mut()
                    .zip(counts)
                    .for_each(|(a, &b)| *a += u32::from(b));
            }
            let from = wraps.partition_point(|&(key, _)| key < first);
            for &(key, n) in wraps[from..].iter().take_while(|w| w.0 < first + RANGE) {
                sum[key - first] += n;
            }
            seed_candidates(first, sum.iter().map(|&n| n as usize))
        })
        .collect();
    let total = parts.iter().map(|(n, _)| n).sum();
    seeds_from(total, parts.into_iter().flat_map(|(_, c)| c).collect())
}

pub fn detect<'a>(reads: &[&'a [u8]], tail: usize) -> Option<Vec<u8>> {
    if reads.len() < 10000 {
        return None;
    }
    if let Some(known) = known(reads) {
        return Some(known);
    }
    let shift = tail.max(1);
    // Each base contributes at most one count. Two independent u32 tables
    // let the existing workers count in parallel without growing the usual
    // 64-bit histogram. A sample containing an exceptionally long read
    // retains the full-width counters.
    let bases = reads
        .iter()
        .try_fold(0u32, |n, r| n.checked_add(u32::try_from(r.len()).ok()?));
    let threads = rayon::current_num_threads();
    let seeds = if bases.is_some() && threads > 1 && usize::BITS >= 64 {
        // Every worker counts into its own table, then sums and ranks a
        // range of keys across the tables. Workers take small pieces of the
        // sample rather than an equal share each, which left the step
        // waiting for the slowest core.
        let tables: Vec<_> = (0..threads)
            .map(|_| Mutex::new(SeedCounts::default()))
            .collect();
        reads.par_chunks(512).for_each(|chunk| {
            let worker = rayon::current_thread_index().unwrap_or(0) % threads;
            let mut counts = tables[worker].lock().unwrap_or_else(|e| e.into_inner());
            if counts.low.is_empty() {
                counts.low.resize(1 << 20, 0);
            }
            for r in chunk {
                each_seed(r, shift, |key| counts.add(key));
            }
        });
        let tables: Vec<SeedCounts> = tables
            .into_iter()
            .map(|t| t.into_inner().unwrap_or_else(|e| e.into_inner()))
            .filter(|t| !t.low.is_empty())
            .collect();
        rank_tables(&tables)
    } else {
        rank_seeds(seed_histogram::<usize>(reads, shift).into_iter())
    };

    // The adapter a seed leads to, if the reads around its hits agree on one.
    // `indexed` holds the seed's hits when a shared scan already found them.
    let try_seed = |seed: &[u8], indexed: Option<&[(usize, usize)]>| -> Option<Vec<u8>> {
        let found;
        let hits: &[(usize, usize)] = if let Some(hits) = indexed {
            hits
        } else {
            // Find every hit in parallel; the trees take them in read order,
            // as a single pass would.
            let finder = memchr::memmem::Finder::new(seed);
            found = reads
                .par_iter()
                .enumerate()
                .flat_map_iter(|(i, &r)| {
                    let end = 500.min(r.len().saturating_sub(10 + shift).saturating_add(1));
                    let mut found = Vec::new();
                    let mut start = 20;
                    while start < end {
                        let Some(offset) = finder.find(&r[start..end + 9]) else {
                            break;
                        };
                        let pos = start + offset;
                        found.push((i, pos));
                        start = pos + 1;
                    }
                    found
                })
                .collect::<Vec<_>>();
            &found
        };
        // The two sides of the seed are independent trees.
        let ((right, right_leaf), (left, left_leaf)) = rayon::join(
            || {
                let mut forward = Tree::<false>::new();
                for &(i, pos) in hits {
                    let r: &'a [u8] = reads[i];
                    forward.add(&r[pos + 10..r.len() - shift]);
                }
                let mut leaf = true;
                (forward.path(&mut leaf), leaf)
            },
            || {
                let mut backward = Tree::<true>::new();
                for &(i, pos) in hits {
                    let r: &'a [u8] = reads[i];
                    backward.add(&r[..pos]);
                }
                let mut leaf = true;
                (backward.path(&mut leaf), leaf)
            },
        );
        let adapter: Vec<_> = left
            .into_iter()
            .rev()
            .chain(seed.iter().copied())
            .chain(right)
            .take(60)
            .collect();
        if let Some(known) = KNOWN.iter().find(|seq| adapter.starts_with(seq)) {
            return Some(known.to_vec());
        }
        (right_leaf && left_leaf).then_some(adapter)
    };
    // Seeds are tried at the same time, and the first in rank order that
    // leads to an adapter wins, as when they are tried in turn. A single scan
    // locates the hits of every seed when they are sparse enough to index.
    if let Some(indexed) = index_seed_hits(reads, shift, &seeds) {
        return seeds
            .par_iter()
            .zip(&indexed)
            .find_map_first(|((seed, _), hits)| try_seed(seed, Some(hits)));
    }
    // Otherwise try the first seed with a scan of its own; when it fails,
    // the remaining seeds may still be sparse enough to share one.
    let ((first, _), rest) = seeds.split_first()?;
    if let Some(adapter) = try_seed(first, None) {
        return Some(adapter);
    }
    let indexed = index_seed_hits(reads, shift, rest);
    rest.par_iter()
        .enumerate()
        .find_map_first(|(k, (seed, _))| {
            try_seed(seed, indexed.as_ref().map(|hits| hits[k].as_slice()))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranking_summed_ranges_matches_serial_ranking() {
        let mut state = 7u64;
        let mut random = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 33) as usize
        };
        for trial in 0..4 {
            let mut tables: Vec<_> = (0..3)
                .map(|_| SeedCounts {
                    low: vec![0; 1 << 20],
                    ..SeedCounts::default()
                })
                .collect();
            let mut summed = vec![0usize; 1 << 20];
            let mut add = |table: usize, key: usize, n: usize| {
                summed[key] += n;
                (0..n).for_each(|_| tables[table].add(key));
            };
            for _ in 0..20000 {
                add(random() % 3, random() % (1 << 20), random() % 50);
            }
            // Ties on count across range boundaries are broken by key, and
            // these counts wrap the one-byte counters several times.
            for key in [3usize, 1 << 14, (1 << 14) + 7, 900_001] {
                add(0, key, 5000 + trial);
                add(2, key, 255);
            }
            assert_eq!(rank_tables(&tables), rank_seeds(summed.into_iter()));
        }
    }

    #[test]
    fn seed_hit_index_preserves_order_overlaps_and_search_bounds() {
        let mut long = b"AC".repeat(300);
        long[22] = b'N';
        long[34] = b'c';
        let reads = [long.as_slice(), b"short", &long[1..49]];
        // The first eight bytes deliberately collide in the bucket index.
        let seeds = [(b"ACACACACAC".to_vec(), 1000), (b"ACACACACCA".to_vec(), 10)];
        for shift in [1, 16, 600] {
            let hits = index_seed_hits(&reads, shift, &seeds).unwrap();
            for ((seed, _), positions) in seeds.iter().zip(hits) {
                let mut expected = Vec::new();
                for (i, r) in reads.iter().enumerate() {
                    let end = 500.min(r.len().saturating_sub(10 + shift).saturating_add(1));
                    for p in 20..end {
                        if r[p..p + 10] == **seed {
                            expected.push((i, p));
                        }
                    }
                }
                assert_eq!(positions, expected);
            }
        }
    }

    #[test]
    fn compressed_edges_match_base_by_base_consensus() {
        fn check<const REVERSE: bool>(reads: &[Vec<u8>]) {
            let mut tree = Tree::<REVERSE>::new();
            let mut paths = Vec::new();
            for read in reads {
                tree.add(read);
                let mut path = read.clone();
                if REVERSE {
                    path.reverse();
                }
                path.truncate(path.iter().position(|&b| b == b'N').unwrap_or(path.len()));
                paths.push(path);
            }
            let mut expected = Vec::new();
            let mut complete = true;
            loop {
                let mut counts = [0usize; 8];
                let mut spelling = [None; 8];
                let depth = expected.len();
                for b in paths.iter().filter_map(|path| path.get(depth)) {
                    let slot = (b & 7) as usize;
                    counts[slot] += 1;
                    spelling[slot].get_or_insert(*b);
                }
                let total = counts.iter().sum::<usize>();
                if total < 50 {
                    break;
                }
                let Some(slot) = counts.iter().position(|&n| 100 * n >= 95 * total) else {
                    complete = false;
                    break;
                };
                expected.push(spelling[slot].unwrap());
                paths.retain(|path| path.get(depth).is_some_and(|b| (b & 7) as usize == slot));
            }
            let mut leaf = true;
            assert_eq!(tree.path(&mut leaf), expected);
            assert_eq!(leaf, complete);
        }

        // Upper/lowercase aliases, paths ending inside another path, N
        // terminators, forks, and the 50-observation / 95% boundary.
        let mut reads = vec![b"acgtACGT".to_vec()];
        reads.extend(vec![b"ACGTACGT".to_vec(); 47]);
        reads.extend(vec![b"TCCCTCCT".to_vec(); 2]);
        for seq in [
            b"ACGT".as_slice(),
            b"ACGNACGT",
            b"ACGTGCGT",
            b"AGGTACGT",
            b"",
        ] {
            reads.push(seq.to_vec());
        }
        reads.extend(vec![b"ACGTACGT".to_vec(); 75]);
        for n in 0..=reads.len() {
            check::<false>(&reads[..n]);
            check::<true>(&reads[..n]);
        }
    }

    #[test]
    fn known_seed_filter_matches_exhaustive_search() {
        let mut state = 140u64;
        let mut random = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as usize
        };
        let adapters: Vec<_> = KNOWN
            .iter()
            .copied()
            .chain([&b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"[..], &b"ACGT"[..]])
            .collect();
        let mut seeds = KnownSeeds::new(&adapters);
        for (adapter_index, &adapter) in adapters.iter().enumerate() {
            for case in 0..120 {
                let len = (adapter.len() + 1 + random() % 160).max(9);
                let mut read: Vec<_> = (0..len).map(|_| b"ATCGNnatcg"[random() % 10]).collect();
                if case % 4 != 0 {
                    let pos = random() % (len - 8);
                    let overlap = (len - pos).min(adapter.len());
                    read[pos..pos + overlap].copy_from_slice(&adapter[..overlap]);
                    for _ in 0..case % (overlap / 16 + 2) {
                        let i = pos + random() % overlap;
                        read[i] = b"ATCGN"[random() % 5];
                    }
                }
                if case == 0 {
                    read.fill(b'A');
                }
                let expected = (0..len.saturating_sub(8)).find_map(|pos| {
                    let overlap = (len - pos).min(adapter.len());
                    let diff = read[pos..pos + overlap]
                        .iter()
                        .zip(adapter)
                        .filter(|(a, b)| a != b)
                        .count();
                    (diff <= overlap / 16).then_some(diff)
                });
                seeds.collect(&read);
                if expected.is_some() {
                    assert!(seeds.active.binary_search(&adapter_index).is_ok());
                }
                assert_eq!(
                    seeds.first_match(adapter_index, &read, adapter),
                    expected,
                    "adapter={adapter:?}, case={case}, read={read:?}"
                );
            }
        }
    }

    #[test]
    fn rolling_seeds_match_independent_window_encoding() {
        let mut seed = 140u64;
        let mut read = Vec::new();
        for _ in 0..300 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            read.push(b"ATCGNaatcg"[(seed >> 32) as usize % 10]);
        }
        read[40..100].copy_from_slice(&b"ATCG".repeat(15));
        for len in [0, 20, 30, 31, 32, 60, 150, 300] {
            for shift in [1, 4, 30, 400] {
                let mut actual = vec![0; 1 << 20];
                let mut expected = vec![0; 1 << 20];
                count_seeds(&read[..len], shift, &mut actual);
                for pos in 20..len.saturating_sub(10 + shift).saturating_add(1) {
                    let encoded = read[pos..pos + 10].iter().try_fold(0usize, |key, b| {
                        b"ATCG"
                            .iter()
                            .position(|v| v == b)
                            .map(|base| key * 4 + base)
                    });
                    if let Some(key) = encoded {
                        expected[key] += 1;
                    }
                }
                assert_eq!(actual, expected, "length={len}, shift={shift}");
            }
        }
    }
}

// The v1.4.0 known adapter table, in the upstream std::map iteration order.
const KNOWN: &[&[u8]] = &[
    b"AACTGTAGGCACCATCAAT",
    b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    b"AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
    b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    b"AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
    b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
    b"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    b"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    b"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA",
    b"AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    b"CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    b"CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC",
    b"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    b"CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    b"CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT",
    b"CCGACAGGTTCAGAGTTCTACAGTCCGACATG",
    b"CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG",
    b"CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG",
    b"CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
    b"CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    b"CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT",
    b"CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    b"CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG",
    b"CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
    b"CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC",
    b"CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC",
    b"CTGTCTCTTATACACATCTGACGCTGCCGACGA",
    b"GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    b"GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    b"GATCGTCGGACTGTAGAACTCTGAACGTGTAGA",
    b"GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    b"TACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    b"TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT",
    b"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
    b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
    b"TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC",
    b"TTTTTTTTTTCAAGCAGAAGACGGCATACGA",
];
