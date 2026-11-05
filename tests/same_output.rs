//! Results must not depend on how the work is split: worker count, batch
//! size and how gzip input is decoded change only speed. Each case runs the
//! same input several ways and compares the decompressed outputs and the JSON
//! report byte for byte.
use clap::Parser;
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

/// One way of running: worker count, batch size and parallel decoding.
#[derive(Debug, Clone, Copy)]
struct Split {
    workers: usize,
    batch_bytes: Option<usize>,
    pgzip_region: usize,
}

const SPLITS: [Split; 6] = [
    Split {
        workers: 1,
        batch_bytes: None,
        pgzip_region: 0,
    },
    Split {
        workers: 2,
        batch_bytes: Some(64 << 10),
        pgzip_region: 0,
    },
    Split {
        workers: 3,
        batch_bytes: Some(200 << 10),
        pgzip_region: 16 << 10,
    },
    Split {
        workers: 8,
        batch_bytes: None,
        pgzip_region: 0,
    },
    Split {
        workers: 8,
        batch_bytes: Some(64 << 10),
        pgzip_region: 16 << 10,
    },
    Split {
        workers: 5,
        batch_bytes: None,
        pgzip_region: 3 << 10,
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

/// Every split gives the first split's outputs and report.
fn check(name: &str, dir: &Path, args: &[&str], outputs: &[&str]) {
    let (files, report) = run(dir, args, outputs, SPLITS[0]);
    assert!(files.iter().any(|f| !f.is_empty()), "{name}: no output");
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
