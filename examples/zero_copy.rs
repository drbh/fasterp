//! fasterp as a library, with no FASTQ text anywhere.
//!
//! A decoder fills readkit batch columns, fasterp trims and filters them,
//! and a k-mer counter reads fasterp's output columns in place. Buffers flow
//! back to the pool they came from, so the steady state allocates nothing.
//!
//! ```text
//! cargo run --release --example zero_copy
//! ```
use fasterp::readkit_batch::{Batch, BatchBuilder, Pool};
use fasterp::readkit_io::BatchSource;
use fasterp::{Options, Paired, Run};

const ADAPTER: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

/// A stand-in for a real decoder (of SRA, say): it writes bases and
/// qualities straight into batch columns, 10,000 pairs per batch. Each pair
/// reads a random insert from both ends; short inserts run into the adapter.
struct Decoder {
    pool: Pool,
    pairs: usize,
    rng: u64,
    insert: Vec<u8>,
}

impl Decoder {
    fn random(&mut self, n: u64) -> u64 {
        self.rng = self.rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        (self.rng >> 33) % n
    }

    /// Mate `segment` of the current insert: read 1 reads it forwards, read
    /// 2 backwards (reverse complement), then the adapter, then random
    /// bases. Qualities fall off after base 120.
    fn mate(&mut self, segment: usize, seq: &mut [u8], qual: &mut [u8]) {
        let n = self.insert.len();
        for (i, (s, q)) in seq.iter_mut().zip(qual.iter_mut()).enumerate() {
            *s = if i < n && segment == 0 {
                self.insert[i]
            } else if i < n {
                match self.insert[n - 1 - i] {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                }
            } else if i - n < ADAPTER.len() {
                ADAPTER[i - n]
            } else {
                b"ACGT"[self.random(4) as usize]
            };
            *q = b'!' + if i < 120 { 37 } else { self.random(38) as u8 };
        }
    }
}

impl BatchSource<Paired> for Decoder {
    type Error = std::convert::Infallible;

    fn next_batch(&mut self) -> Result<Option<Batch<Paired>>, Self::Error> {
        let n = self.pairs.min(10_000);
        if n == 0 {
            return Ok(None);
        }
        self.pairs -= n;
        let pool = self.pool.clone();
        let mut batch = BatchBuilder::<Paired>::with_pool(&pool, n);
        for _ in 0..n {
            let len = 10 + self.random(250) as usize;
            self.insert.clear();
            for _ in 0..len {
                let base = b"ACGT"[self.random(4) as usize];
                self.insert.push(base);
            }
            // No names: a quantifier does not need them, so none are made.
            // `push_with` hands out the column space to fill in place.
            batch
                .push_with(&[150, 150], None, |s, seq, qual| self.mate(s, seq, qual))
                .expect("both mates have 150 bases and qualities");
        }
        Ok(Some(batch.finish()))
    }
}

fn main() -> Result<(), fasterp::Error> {
    let decoder = Decoder {
        pool: Pool::new(),
        pairs: 200_000,
        rng: 7,
        insert: Vec::new(),
    };
    let decoder_pool = decoder.pool.clone();

    // fastp's defaults plus adapter and quality trimming.
    let options = Options::new()
        .threads(4)
        .adapter(std::str::from_utf8(ADAPTER).unwrap())
        .cut_tail(20);
    let mut run = Run::from_source(&options, decoder)?;

    // Downstream: count 8-mers of every trimmed read straight from fasterp's
    // sequence columns (one contiguous buffer and offsets per mate).
    const K: usize = 8;
    let mut kmers = vec![0u32; 1 << (2 * K)];
    while let Some(batch) = run.next_batch()? {
        for column in [batch.read1(), batch.read2()] {
            for seq in column.seqs().iter() {
                count_kmers::<K>(seq, &mut kmers);
            }
        }
    }
    let output_pool = run.pool().clone();
    let report = run.finish()?;

    let top = (0..kmers.len()).max_by_key(|&i| kmers[i]).unwrap_or(0);
    println!(
        "{} of {} reads passed; most common {K}-mer {} ({} times)",
        report.passed_reads(),
        report.total_reads(),
        decode::<K>(top),
        kmers[top]
    );
    println!("filtering: {}", report.json()["filtering_result"]);
    for (pool, label) in [(decoder_pool, "decoder"), (output_pool, "fasterp")] {
        let s = pool.stats();
        println!(
            "{label} pool: {} buffers allocated, {} reused",
            s.fresh, s.reused
        );
    }
    Ok(())
}

/// Adds every `K`-mer of `seq` without an N to `counts` (2 bits per base).
fn count_kmers<const K: usize>(seq: &[u8], counts: &mut [u32]) {
    let mask = (1usize << (2 * K)) - 1;
    let (mut code, mut valid) = (0usize, 0);
    for &base in seq {
        let bits = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                valid = 0;
                continue;
            }
        };
        code = (code << 2 | bits) & mask;
        valid += 1;
        if valid >= K {
            counts[code] += 1;
        }
    }
}

fn decode<const K: usize>(code: usize) -> String {
    (0..K)
        .rev()
        .map(|i| b"ACGT"[(code >> (2 * i)) & 3] as char)
        .collect()
}
