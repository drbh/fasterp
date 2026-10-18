//! Adapter evaluation ported from OpenGene/fastp v1.4.0 evaluator.cpp and
//! nucleotidetree.cpp. Upstream copyright and MIT terms: LICENSE.fastp.
use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(Default)]
struct SeedCandidates {
    positions: Vec<usize>,
    dense: bool,
}

struct KnownSeeds {
    index: rustc_hash::FxHashMap<u64, Vec<(usize, usize)>>,
    candidates: Vec<SeedCandidates>,
    short: Vec<bool>,
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
        Self {
            index,
            candidates: (0..adapters.len())
                .map(|_| SeedCandidates::default())
                .collect(),
            short: adapters.iter().map(|a| a.len() < 8).collect(),
        }
    }

    fn collect(&mut self, read: &[u8]) {
        for (c, &short) in self.candidates.iter_mut().zip(&self.short) {
            c.positions.clear();
            c.dense = short;
        }
        // Scan each read once for all known adapters. Exact 8-byte keys retain
        // case and ambiguous bytes; no DNA normalization or lossy hash is used.
        for (found, bytes) in read.windows(8).enumerate() {
            let key = u64::from_le_bytes(bytes.try_into().unwrap());
            if let Some(seeds) = self.index.get(&key) {
                for &(adapter, offset) in seeds {
                    let c = &mut self.candidates[adapter];
                    if !c.dense && found >= offset && found - offset < read.len().saturating_sub(8)
                    {
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
    let mut seeds = KnownSeeds::new(KNOWN);
    let mut counts = vec![0usize; KNOWN.len()];
    let mut mismatches = vec![0usize; KNOWN.len()];
    let mut checked = 0;
    let mut bases = 0;
    let mut max = 0;
    for r in reads {
        checked += 1;
        bases += r.len();
        if checked > 100000 || bases > 100000000 || max > 1000 {
            break;
        }
        seeds.collect(r);
        for (i, adapter) in KNOWN.iter().enumerate() {
            if adapter.len() >= r.len() || (max > 20 && counts[i] < max / 10) {
                continue;
            }
            if let Some(diff) = seeds.first_match(i, r, adapter) {
                counts[i] += 1;
                max = max.max(counts[i]);
                mismatches[i] += diff;
            }
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

fn count_seeds(read: &[u8], shift: usize, counts: &mut [usize]) {
    let end = read.len().saturating_sub(10 + shift).saturating_add(1);
    if end <= 20 {
        return;
    }
    let mut key = 0;
    let mut valid = 0;
    for b in &read[20..end + 9] {
        let base = match b {
            b'A' => 0,
            b'T' => 1,
            b'C' => 2,
            b'G' => 3,
            _ => {
                key = 0;
                valid = 0;
                continue;
            }
        };
        key = ((key << 2) | base) & ((1 << 20) - 1);
        valid += 1;
        if valid >= 10 {
            counts[key] += 1;
        }
    }
}

#[derive(Default)]
struct Node {
    children: [Option<usize>; 8],
    count: usize,
    base: u8,
}

struct Tree {
    nodes: Vec<Node>,
}
impl Tree {
    fn new() -> Self {
        Self {
            nodes: vec![Node::default()],
        }
    }
    fn add(&mut self, seq: impl Iterator<Item = u8>) {
        let mut node = 0;
        for b in seq {
            if b == b'N' {
                break;
            }
            let idx = (b & 7) as usize;
            let child = if let Some(child) = self.nodes[node].children[idx] {
                child
            } else {
                let child = self.nodes.len();
                self.nodes.push(Node {
                    base: b,
                    ..Node::default()
                });
                self.nodes[node].children[idx] = Some(child);
                child
            };
            self.nodes[child].count += 1;
            node = child;
        }
    }
    fn path(&self, leaf: &mut bool) -> Vec<u8> {
        let mut node = 0;
        let mut path = Vec::new();
        loop {
            let total: usize = self.nodes[node]
                .children
                .iter()
                .flatten()
                .map(|&n| self.nodes[n].count)
                .sum();
            if total < 50 {
                break;
            }
            let next = self.nodes[node]
                .children
                .iter()
                .flatten()
                .find(|&&n| self.nodes[n].count as f64 / total as f64 >= 0.95);
            if let Some(&next) = next {
                path.push(self.nodes[next].base);
                node = next;
            } else {
                *leaf = false;
                break;
            }
        }
        path
    }
}

pub fn detect(reads: &[&[u8]], tail: usize) -> Option<Vec<u8>> {
    if reads.len() < 10000 {
        return None;
    }
    if let Some(known) = known(reads) {
        return Some(known);
    }
    let shift = tail.max(1);
    let mut counts = vec![0usize; 1 << 20];
    for r in reads {
        count_seeds(r, shift, &mut counts);
    }
    counts[0] = 0;
    // Only the ten highest (count, key) pairs are used by discovery.
    let mut candidates = BinaryHeap::with_capacity(10);
    let mut total = 0;
    for (key, &n) in counts.iter().enumerate() {
        if n == 0 {
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
    let mut candidates: Vec<_> = candidates
        .into_iter()
        .map(|Reverse((n, key))| (key, n))
        .collect();
    candidates.sort_unstable_by_key(|&(key, n)| (std::cmp::Reverse(n), std::cmp::Reverse(key)));
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
        let mut forward = Tree::new();
        let mut backward = Tree::new();
        for r in reads {
            for pos in 20..500.min(r.len().saturating_sub(10 + shift).saturating_add(1)) {
                if r[pos..pos + 10] == seed {
                    forward.add(r[pos + 10..r.len() - shift].iter().copied());
                    backward.add(r[..pos].iter().rev().copied());
                }
            }
        }
        let mut leaf = true;
        let right = forward.path(&mut leaf);
        let left = backward.path(&mut leaf);
        let mut adapter: Vec<_> = left
            .into_iter()
            .rev()
            .chain(seed)
            .chain(right)
            .take(60)
            .collect();
        if let Some(known) = KNOWN.iter().find(|seq| adapter.starts_with(seq)) {
            adapter = known.to_vec();
            return Some(adapter);
        }
        if leaf {
            return Some(adapter);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

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
