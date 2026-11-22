//! An `.sra` run through fasterp with no FASTQ text anywhere.
//!
//! readkit-sra's in-process VDB decoder (its `ngs` backend) writes bases and
//! qualities straight into batch columns, fasterp trims and filters them,
//! and a k-mer counter reads fasterp's output columns in place. Buffers flow
//! back to the pool they came from.
//!
//! ```text
//! cargo run --release --manifest-path examples/sra-pipeline/Cargo.toml -- RUN.sra [THREADS] [--fastq OUT]
//! ```
//!
//! Needs the SRA Toolkit's `libncbi-ngs` (set `NGS_PREFIX` if it is not
//! under `/opt/homebrew`). With `--fastq OUT` the decoder also names reads
//! as `fasterq-dump --split-spot` does and the passing reads are written to
//! `OUT` (with a bare `+` line), for comparing with the command line's `-o`
//! on `fasterq-dump --stdout --split-spot RUN.sra`.
use fasterp::{Options, Run, Single};
use readkit_sra::{DecodeOptions, Names, NgsBackend, SingleDecoder};
use std::io::Write;
use std::sync::Arc;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let mut path = None;
    let mut threads = 4;
    let mut fastq = None;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--fastq" => fastq = Some(args.next().ok_or("--fastq needs a path")?),
            _ if path.is_none() => path = Some(arg),
            _ => threads = arg.parse()?,
        }
    }
    let path = path.ok_or("usage: sra_pipeline RUN.sra [THREADS] [--fastq OUT]")?;

    // The default segment selection keeps biological reads only, as
    // fasterq-dump does: one 301-base read per spot for the yeast runs, with
    // the empty technical read left out. Names cost a formatted string per
    // read and nothing downstream needs them, so they are made only for the
    // FASTQ check.
    let decode_options = DecodeOptions {
        names: if fastq.is_some() {
            Names::FasterqSplitSpot
        } else {
            Names::None
        },
        threads: 2,
        // Smaller than the default 16 384 spots, so that a small run still
        // makes enough batches for fasterp's buffers to come back and be
        // reused.
        batch_spots: 4096,
        ..DecodeOptions::default()
    };
    let decoder = SingleDecoder::open_with(Arc::new(NgsBackend::new()), &path, decode_options)?;
    println!(
        "{path}: {} spots, read types {:?}, decoded by {}",
        decoder.info().spots,
        decoder.read_types(),
        decoder.backend_name()
    );
    let decoder_pool = decoder.pool().clone();

    // fastp's defaults: adapter detection, quality and length filters.
    let options = Options::new().threads(threads);
    let started = Instant::now();
    let mut run: Run<Single> = Run::from_source(&options, decoder)?;

    // Downstream: count 8-mers of every passing read straight from
    // fasterp's sequence column (one contiguous buffer plus offsets).
    const K: usize = 8;
    let mut kmers = vec![0u32; 1 << (2 * K)];
    let (mut bases, mut gc, mut longest) = (0usize, 0usize, 0usize);
    let mut out = match &fastq {
        Some(p) => Some(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => None,
    };
    while let Some(batch) = run.next_batch()? {
        let seqs = batch.column().seqs();
        let data = seqs.data();
        bases += data.len();
        gc += data.iter().filter(|&&b| b == b'G' || b == b'C').count();
        for seq in seqs.iter() {
            longest = longest.max(seq.len());
            count_kmers::<K>(seq, &mut kmers);
        }
        if let Some(out) = &mut out {
            for read in batch.reads() {
                let name = read.name().unwrap_or_default();
                for part in [b"@", name, b"\n", read.seq(), b"\n+\n", read.qual(), b"\n"] {
                    out.write_all(part)?;
                }
            }
        }
    }
    if let Some(mut out) = out {
        out.flush()?;
    }
    let output_pool = run.pool().clone();
    let progress = run.progress();
    let report = run.finish()?;
    let elapsed = started.elapsed();

    let filtering = &report.json()["filtering_result"];
    println!(
        "{} of {} reads passed in {} batches ({:.2} s)",
        report.passed_reads(),
        report.total_reads(),
        progress.batches,
        elapsed.as_secs_f64()
    );
    println!("filtering: {filtering}");
    println!(
        "passing bases: {bases} (mean length {:.1}, longest {longest}, GC {:.1}%)",
        bases as f64 / report.passed_reads().max(1) as f64,
        100.0 * gc as f64 / bases.max(1) as f64
    );
    let top = (0..kmers.len()).max_by_key(|&i| kmers[i]).unwrap_or(0);
    println!(
        "most common {K}-mer: {} ({} times)",
        decode::<K>(top),
        kmers[top]
    );
    // fasterp holds the input it samples for adapter detection (as fastp,
    // the first 256Ki reads) until detection is done, so the decoder's
    // buffers start coming back only after that: a run smaller than the
    // sample allocates every batch. fasterp's output buffers come back as
    // soon as each batch is dropped.
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
