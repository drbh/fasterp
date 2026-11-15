//! Results must not depend on how the work is split: worker count, batch
//! size, how gzip input is decoded, and whether the command line or the
//! library runs it (on readers, or on readkit batches from a source) change
//! only speed. Each case runs the same input several ways and compares the
//! decompressed outputs and the JSON report byte for byte; library batches
//! are written as FASTQ for the comparison.
use clap::Parser;
use fasterp::readkit_batch::{Batch, BatchBuilder, Pool, Read as Record};
use fasterp::readkit_io::iter_source as source;
use fasterp::{Failure, Layout, Options, Paired, Run, Single};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Mutex;

/// Runs share the process-wide decoding override, so they take turns.
static SERIAL: Mutex<()> = Mutex::new(());

/// A small deterministic generator, so the input is the same on every run.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

const ADAPTER: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

/// Paired reads with what the processing steps act on: adapters after short
/// inserts, poly-G tails, Ns, low-quality tails, duplicates, varying lengths
/// and both binned and full-range qualities.
fn pairs(n: usize, binned: bool) -> (Vec<u8>, Vec<u8>) {
    let mut rng = Rng(7 + binned as u64);
    let quals: Vec<u8> = if binned {
        b"FF:F,#".to_vec()
    } else {
        (b'#'..=b'J').collect()
    };
    let (mut r1, mut r2) = (Vec::new(), Vec::new());
    let mut last: Option<(Vec<u8>, Vec<u8>)> = None;
    for i in 0..n {
        let (s1, s2) = match &last {
            Some(prev) if rng.below(20) == 0 => prev.clone(),
            _ => {
                let len = 60 + rng.below(91);
                let insert: Vec<u8> = (0..40 + rng.below(200))
                    .map(|_| b"ACGT"[rng.below(4)])
                    .collect();
                let mate = |fragment: &[u8], rng: &mut Rng| {
                    let mut s: Vec<u8> =
                        fragment.iter().chain(ADAPTER).copied().take(len).collect();
                    while s.len() < len {
                        s.push(b"ACGT"[rng.below(4)]);
                    }
                    if rng.below(10) == 0 {
                        let at = rng.below(len);
                        s[at..].fill(b'G');
                    }
                    if rng.below(15) == 0 {
                        let at = rng.below(len);
                        s[at] = b'N';
                    }
                    s
                };
                let reverse: Vec<u8> = insert
                    .iter()
                    .rev()
                    .map(|b| match b {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        _ => b'A',
                    })
                    .collect();
                (mate(&insert, &mut rng), mate(&reverse, &mut rng))
            }
        };
        for (out, seq, mate) in [(&mut r1, &s1, 1), (&mut r2, &s2, 2)] {
            writeln!(
                out,
                "@A00123:8:H7:1:{}:{}:{} {mate}:N:0:ACGT",
                i % 7,
                i,
                i * 3
            )
            .unwrap();
            out.extend_from_slice(seq);
            out.extend_from_slice(b"\n+\n");
            let tail = rng.below(seq.len());
            for k in 0..seq.len() {
                let q = if k > tail && rng.below(3) == 0 {
                    quals[quals.len() - 1]
                } else {
                    quals[rng.below(quals.len())]
                };
                out.push(q);
            }
            out.push(b'\n');
        }
        last = Some((s1, s2));
    }
    (r1, r2)
}

/// Both mates of each pair, one after the other.
fn interleave(r1: &[u8], r2: &[u8]) -> Vec<u8> {
    let (a, b) = (records(r1), records(r2));
    assert_eq!(a.len(), b.len());
    let mut out = Vec::new();
    for (a, b) in a.iter().zip(&b) {
        for (name, seq, qual) in [a, b] {
            write_fastq(&mut out, Record::named(name, seq, qual), None);
        }
    }
    out
}

/// One single-member gzip file, which the parallel decoder can split.
fn write_gz(path: &Path, data: &[u8]) {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
    e.write_all(data).unwrap();
    std::fs::write(path, e.finish().unwrap()).unwrap();
}

fn decompressed(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Who runs: the command line, or the library on readers or on readkit
/// batches of this many records.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Runner {
    Cli,
    Readers,
    Source(usize),
}

/// One way of running: worker count, batch size, parallel decoding and
/// the runner.
#[derive(Debug, Clone, Copy)]
struct Split {
    workers: usize,
    batch_bytes: Option<usize>,
    pgzip_region: usize,
    runner: Runner,
}

const SPLITS: [Split; 10] = [
    Split {
        workers: 1,
        batch_bytes: None,
        pgzip_region: 0,
        runner: Runner::Cli,
    },
    Split {
        workers: 2,
        batch_bytes: Some(64 << 10),
        pgzip_region: 0,
        runner: Runner::Cli,
    },
    Split {
        workers: 3,
        batch_bytes: Some(200 << 10),
        pgzip_region: 16 << 10,
        runner: Runner::Cli,
    },
    Split {
        workers: 8,
        batch_bytes: None,
        pgzip_region: 0,
        runner: Runner::Cli,
    },
    Split {
        workers: 8,
        batch_bytes: Some(64 << 10),
        pgzip_region: 16 << 10,
        runner: Runner::Cli,
    },
    Split {
        workers: 5,
        batch_bytes: None,
        pgzip_region: 3 << 10,
        runner: Runner::Cli,
    },
    Split {
        workers: 1,
        batch_bytes: None,
        pgzip_region: 0,
        runner: Runner::Readers,
    },
    Split {
        workers: 8,
        batch_bytes: Some(64 << 10),
        pgzip_region: 0,
        runner: Runner::Readers,
    },
    Split {
        workers: 3,
        batch_bytes: None,
        pgzip_region: 0,
        runner: Runner::Source(5000),
    },
    Split {
        workers: 8,
        batch_bytes: None,
        pgzip_region: 0,
        runner: Runner::Source(997),
    },
];

/// Run fasterp with `args`, writing to `outputs` in `dir`, and return the
/// decompressed outputs and the report.
fn run(dir: &Path, args: &[&str], outputs: &[&str], split: Split) -> (Vec<Vec<u8>>, String) {
    let _turn = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    for name in outputs {
        let _ = std::fs::remove_file(dir.join(name));
    }
    let report = dir.join("report.json");
    // Arguments that name files are relative to `dir`.
    let mut argv: Vec<String> = std::iter::once("fasterp".to_string())
        .chain(args.iter().map(|a| match a.ends_with(".gz") {
            true => dir.join(a).to_str().unwrap().to_string(),
            false => a.to_string(),
        }))
        .collect();
    argv.extend([
        "-w".into(),
        split.workers.to_string(),
        "-j".into(),
        report.to_str().unwrap().into(),
        "-h".into(),
        dir.join("report.html").to_str().unwrap().into(),
    ]);
    if let Some(bytes) = split.batch_bytes {
        argv.extend(["--batch-bytes".into(), bytes.to_string()]);
    }
    if split.runner != Runner::Cli {
        return run_library(argv, outputs, split.runner);
    }
    fasterp::testing::force_split_decoding(split.pgzip_region);
    let result = fasterp::compat::run(fasterp::Args::parse_from(&argv));
    fasterp::testing::force_split_decoding(0);
    result.unwrap_or_else(|e| panic!("{split:?}: {e:#}"));
    let files = outputs
        .iter()
        .map(|name| decompressed(&dir.join(name)))
        .collect();
    (files, std::fs::read_to_string(report).unwrap())
}

/// The fields of each FASTQ record: name without `@`, sequence, quality.
fn records(fastq: &[u8]) -> Vec<(&[u8], &[u8], &[u8])> {
    let lines: Vec<&[u8]> = fastq.split(|&b| b == b'\n').collect();
    lines
        .chunks_exact(4)
        .map(|l| (&l[0][1..], l[1], l[3]))
        .collect()
}

/// A record as FASTQ, with fastp's failure tag after the name if given.
fn write_fastq(out: &mut Vec<u8>, read: Record, tag: Option<&str>) {
    out.push(b'@');
    out.extend_from_slice(read.name().unwrap());
    if let Some(tag) = tag {
        out.push(b' ');
        out.extend_from_slice(tag.as_bytes());
    }
    for line in [read.seq(), b"+", read.qual()] {
        out.push(b'\n');
        out.extend_from_slice(line);
    }
    out.push(b'\n');
}

/// Readkit batches of `per_batch` records holding `mates` (one FASTQ text
/// per segment): what a decoder with no FASTQ text in between would yield.
fn batches<L: Layout>(mates: &[&[u8]], per_batch: usize) -> Vec<Result<Batch<L>, std::io::Error>> {
    let pool = Pool::new();
    let mates: Vec<_> = mates.iter().map(|m| records(m)).collect();
    let mut out = Vec::new();
    for start in (0..mates[0].len()).step_by(per_batch) {
        let end = (start + per_batch).min(mates[0].len());
        let mut b = BatchBuilder::with_segments(&pool, mates.len(), end - start);
        for i in start..end {
            let record: Vec<_> = mates
                .iter()
                .map(|m| Record::named(m[i].0, m[i].1, m[i].2))
                .collect();
            b.push(&record).unwrap();
        }
        let batch = b.finish();
        out.push(Ok(match L::SEGMENTS {
            Some(1) => relabel(batch.try_into_layout::<Single>().ok().unwrap()),
            _ => relabel(batch.try_into_layout::<Paired>().ok().unwrap()),
        }));
    }
    out
}

fn relabel<A: 'static, B: 'static>(batch: Batch<A>) -> Batch<B> {
    *(Box::new(batch) as Box<dyn std::any::Any>)
        .downcast::<Batch<B>>()
        .unwrap()
}

/// FASTQ texts of a library run, keyed by the command line option that
/// writes them.
#[derive(Default)]
struct Texts(BTreeMap<&'static str, Vec<u8>>);

impl Texts {
    fn drain<L: Layout>(&mut self, run: Run<L>) -> fasterp::Report {
        let mut run = run;
        while let Some(batch) = run
            .next_batch()
            .unwrap_or_else(|e| panic!("library: {e:#}"))
        {
            for (s, column) in batch.records().columns().enumerate() {
                let out = self.0.entry(["-o", "-O"][s]).or_default();
                for i in 0..column.len() {
                    write_fastq(out, column.read(i), None);
                }
            }
            for read in batch.merged().reads() {
                write_fastq(self.0.entry("--merged_out").or_default(), read, None);
            }
            for (read, reason, _) in batch.failures() {
                let out = self.0.entry("--failed_out").or_default();
                write_fastq(out, read, Some(reason.fastp_name()));
            }
            let tags = batch.unpaired().tags().unwrap_or(&[]);
            for (read, &tag) in batch.unpaired().reads().zip(tags) {
                let flag = ["--unpaired1", "--unpaired2"][Failure::segment_of_tag(tag)];
                write_fastq(self.0.entry(flag).or_default(), read, None);
            }
        }
        run.finish().unwrap_or_else(|e| panic!("library: {e:#}"))
    }
}

/// Like `run`, but through the library: inputs from readers or readkit
/// batches, outputs written from the returned batches. The report is
/// serialized the way `-j` writes it by default.
/// The caller holds the `SERIAL` turn.
fn run_library(argv: Vec<String>, outputs: &[&str], runner: Runner) -> (Vec<Vec<u8>>, String) {
    const PATHS: [&str; 11] = [
        "-i",
        "-I",
        "-o",
        "-O",
        "-j",
        "-h",
        "--merged_out",
        "--failed_out",
        "--unpaired1",
        "--unpaired2",
        "--adapter_fasta",
    ];
    let mut paths = BTreeMap::new();
    let mut rest = Vec::new();
    let mut interleaved = false;
    // Skip the program name; the library takes options only.
    let mut argv = argv.into_iter().skip(1);
    while let Some(arg) = argv.next() {
        match arg.as_str() {
            "--interleaved_in" => interleaved = true,
            flag if PATHS.contains(&flag) => {
                paths.insert(
                    PATHS.iter().find(|p| **p == flag).unwrap(),
                    argv.next().unwrap(),
                );
            }
            _ => rest.push(arg),
        }
    }
    let mut options = Options::from_fastp_options(&rest).unwrap();
    if paths.contains_key(&"--failed_out") {
        options = options.keep_failed();
    }
    if paths.contains_key(&"--unpaired1") {
        options = options.keep_unpaired();
    }
    if let Some(path) = paths.get(&"--adapter_fasta") {
        let fasta = std::fs::read_to_string(path).unwrap();
        options = options.adapters(fasta.lines().filter(|l| !l.starts_with('>')));
    }
    let input = |flag: &str| paths.get(&flag).map(|p| decompressed(Path::new(p)));
    let (in1, in2) = (input("-i").unwrap(), input("-I"));
    let cursor = std::io::Cursor::new;
    let mut texts = Texts::default();
    let report = match (runner, in2) {
        (Runner::Readers, Some(in2)) => {
            texts.drain(Run::paired(&options, cursor(in1), cursor(in2)).unwrap())
        }
        (Runner::Readers, None) if interleaved => {
            texts.drain(Run::interleaved(&options, cursor(in1)).unwrap())
        }
        (Runner::Readers, None) => texts.drain(Run::single(&options, cursor(in1)).unwrap()),
        (Runner::Source(n), in2) => {
            match in2 {
                Some(in2) => texts.drain(
                    Run::from_source(&options, source(batches::<Paired>(&[&in1, &in2], n)))
                        .unwrap(),
                ),
                None if interleaved => {
                    // An interleaved file holds pairs: give them as such.
                    let r: Vec<_> = records(&in1);
                    let mates: [Vec<u8>; 2] = std::array::from_fn(|m| {
                        let mut out = Vec::new();
                        for (name, seq, qual) in r.iter().skip(m).step_by(2) {
                            write_fastq(&mut out, Record::named(name, seq, qual), None);
                        }
                        out
                    });
                    let batches = batches::<Paired>(&[&mates[0], &mates[1]], n);
                    texts.drain(Run::from_source(&options, source(batches)).unwrap())
                }
                None => texts.drain(
                    Run::from_source(&options, source(batches::<Single>(&[&in1], n))).unwrap(),
                ),
            }
        }
        (Runner::Cli, _) => unreachable!(),
    };
    let files = outputs
        .iter()
        .map(|name| {
            let flag = paths
                .iter()
                .find(|(_, path)| Path::new(path).file_name() == Some(std::ffi::OsStr::new(name)))
                .map(|(flag, _)| **flag)
                .unwrap_or_else(|| panic!("{name} is not an output"));
            texts.0.get(flag).cloned().unwrap_or_default()
        })
        .collect();
    (files, serde_json::to_string(report.json()).unwrap())
}

/// Every split gives the first split's outputs and report.
fn check(name: &str, dir: &Path, args: &[&str], outputs: &[&str]) {
    let (files, report) = run(dir, args, outputs, SPLITS[0]);
    for (file, output) in files.iter().zip(outputs) {
        assert!(!file.is_empty(), "{name}: {output} is empty");
    }
    for split in &SPLITS[1..] {
        let (other_files, other_report) = run(dir, args, outputs, *split);
        for (k, (a, b)) in files.iter().zip(&other_files).enumerate() {
            assert!(
                a == b,
                "{name}: output {} differs with {split:?}",
                outputs[k]
            );
        }
        assert!(
            report == other_report,
            "{name}: report differs with {split:?}"
        );
    }
}

#[test]
fn results_do_not_depend_on_how_work_is_split() {
    for binned in [true, false] {
        let dir = tempfile::tempdir().unwrap();
        let d = dir.path();
        let (r1, r2) = pairs(30_000, binned);
        write_gz(&d.join("in1.fq.gz"), &r1);
        write_gz(&d.join("in2.fq.gz"), &r2);
        write_gz(&d.join("interleaved.fq.gz"), &interleave(&r1, &r2));
        let label = |w: &str| format!("{w} ({})", if binned { "binned" } else { "full-range" });

        check(
            &label("single-end"),
            d,
            &["-i", "in1.fq.gz", "-o", "out1.fq.gz"],
            &["out1.fq.gz"],
        );
        check(
            &label("single-end, poly-G, low complexity, overrepresentation"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-o",
                "out1.fq.gz",
                "-g",
                "-y",
                "-p",
                "-5",
                "-3",
            ],
            &["out1.fq.gz"],
        );
        check(
            &label("paired-end"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-I",
                "in2.fq.gz",
                "-o",
                "out1.fq.gz",
                "-O",
                "out2.fq.gz",
            ],
            &["out1.fq.gz", "out2.fq.gz"],
        );
        std::fs::write(
            d.join("adapters.fa"),
            format!(
                ">a\n{}\n>b\nCTGTCTCTTATACACATCT\n",
                String::from_utf8_lossy(&ADAPTER[..20])
            ),
        )
        .unwrap();
        check(
            &label("single-end adapter fasta"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-o",
                "out1.fq.gz",
                "--adapter_fasta",
                &d.join("adapters.fa").to_string_lossy(),
            ],
            &["out1.fq.gz"],
        );
        check(
            &label("interleaved paired-end"),
            d,
            &[
                "--interleaved_in",
                "-i",
                "interleaved.fq.gz",
                "-o",
                "out1.fq.gz",
                "-O",
                "out2.fq.gz",
            ],
            &["out1.fq.gz", "out2.fq.gz"],
        );
        check(
            &label("single-end failed reads"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-o",
                "out1.fq.gz",
                "--failed_out",
                "failed.fq.gz",
                "-l",
                "70",
                "--length_limit",
                "140",
                "-n",
                "0",
                "-y",
                "-e",
                "25",
            ],
            &["out1.fq.gz", "failed.fq.gz"],
        );
        check(
            &label("paired failed and unpaired reads"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-I",
                "in2.fq.gz",
                "-o",
                "out1.fq.gz",
                "-O",
                "out2.fq.gz",
                "--failed_out",
                "failed.fq.gz",
                "--unpaired1",
                "unpaired1.fq.gz",
                "--unpaired2",
                "unpaired2.fq.gz",
                "-l",
                "70",
                "--length_limit",
                "140",
                "-n",
                "0",
                "-y",
            ],
            &[
                "out1.fq.gz",
                "out2.fq.gz",
                "failed.fq.gz",
                "unpaired1.fq.gz",
                "unpaired2.fq.gz",
            ],
        );
        check(
            &label("paired failed reads, mates not kept"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-I",
                "in2.fq.gz",
                "-o",
                "out1.fq.gz",
                "-O",
                "out2.fq.gz",
                "--failed_out",
                "failed.fq.gz",
                "-l",
                "90",
                "-n",
                "1",
            ],
            &["out1.fq.gz", "out2.fq.gz", "failed.fq.gz"],
        );
        check(
            &label("paired correct, merge, dedup"),
            d,
            &[
                "-i",
                "in1.fq.gz",
                "-I",
                "in2.fq.gz",
                "-o",
                "out1.fq.gz",
                "-O",
                "out2.fq.gz",
                "-c",
                "-m",
                "--merged_out",
                "merged.fq.gz",
                "-D",
                "--detect_adapter_for_pe",
            ],
            &["out1.fq.gz", "out2.fq.gz", "merged.fq.gz"],
        );
    }
}
