//! Self-contained regressions for the fastp 1.4.0 CLI port.
use flate2::{write::GzEncoder, Compression};
use serde_json::Value;
use std::{
    fs,
    io::Write,
    path::Path,
    process::{Command, Output, Stdio},
};
use tempfile::TempDir;

fn record(name: &str, seq: &str, qual: char) -> String {
    format!(
        "@{name}\n{seq}\n+\n{}\n",
        qual.to_string().repeat(seq.len())
    )
}
fn command(root: &Path) -> Command {
    let mut c = Command::new(assert_cmd::cargo::cargo_bin!("fasterp"));
    c.current_dir(root).args([
        "--dont_eval_duplication",
        "-A",
        "-G",
        "-j",
        "report.json",
        "-h",
        "report.html",
    ]);
    c
}
fn success(output: Output) -> Output {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
fn report(root: &Path) -> Value {
    serde_json::from_slice(&fs::read(root.join("report.json")).unwrap()).unwrap()
}
fn fixture() -> (TempDir, String) {
    let d = TempDir::new().unwrap();
    let data = (0..30)
        .map(|i| record(&format!("read{i}"), "ACGTGTCAGTACGATCGTAGCTAGCTACGATC", 'I'))
        .collect::<String>();
    fs::write(d.path().join("in.fq"), &data).unwrap();
    (d, data)
}

fn bounded_exit(mut child: std::process::Child) -> std::process::ExitStatus {
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(15);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            return status;
        }
        if std::time::Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            panic!("FASTQ pipeline did not terminate after an input/output error");
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
}

#[test]
fn pipeline_stops_when_output_pipe_closes() {
    let d = TempDir::new().unwrap();
    let mut file = std::io::BufWriter::new(fs::File::create(d.path().join("in.fq")).unwrap());
    let read = record("read", &"ACGT".repeat(38), 'I');
    // Exceed the stdout writer's buffer while the reader still has work queued.
    for _ in 0..100_000 {
        file.write_all(read.as_bytes()).unwrap();
    }
    drop(file);
    for output in [&["--stdout"][..], &["-o", "/dev/stdout"][..]] {
        let mut child = command(d.path())
            .args(["-i", "in.fq", "-w", "4", "--batch_bytes", "65536"])
            .args(output)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        drop(child.stdout.take());
        assert!(!bounded_exit(child).success());
    }
}

#[test]
fn pipeline_propagates_late_input_error() {
    let d = TempDir::new().unwrap();
    let good = record("read", &"ACGT".repeat(38), 'I');
    let mut data = good.repeat(1_000);
    data.push_str("@broken\nACGT\n+\nI\n");
    fs::write(d.path().join("in.fq"), data).unwrap();
    let child = command(d.path())
        .args([
            "-i",
            "in.fq",
            "-o",
            "out.fq",
            "-w",
            "4",
            "--batch_bytes",
            "4096",
        ])
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    assert!(!bounded_exit(child).success());
    assert!(!d.path().join("report.json").exists());
}
fn bgzf(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    let gz = encoder.finish().unwrap();
    let size = (gz.len() + 7) as u16;
    let mut block = vec![31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0];
    block.extend(size.to_le_bytes());
    block.extend_from_slice(&gz[10..]);
    block
}

#[test]
fn gzip_checksum_failure_is_not_reported_as_clean_eof() {
    let (d, data) = fixture();
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data.as_bytes()).unwrap();
    let mut bytes = encoder.finish().unwrap();
    let crc = bytes.len() - 8;
    bytes[crc] ^= 1;
    fs::write(d.path().join("corrupt.fq.gz"), bytes).unwrap();
    let output = command(d.path())
        .args(["-i", "corrupt.fq.gz", "-o", "out.fq"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(!d.path().join("report.json").exists());
}

#[test]
fn bgzf_embedded_eof_does_not_truncate_and_magic_detection_works() {
    let (d, data) = fixture();
    let split = data.len() / 2;
    let bytes = [
        bgzf(&data.as_bytes()[..split]),
        bgzf(b""),
        bgzf(&data.as_bytes()[split..]),
        bgzf(b""),
    ]
    .concat();
    fs::write(d.path().join("blocks.data"), bytes).unwrap();
    success(
        command(d.path())
            .args(["-i", "blocks.data", "-o", "out.fq"])
            .output()
            .unwrap(),
    );
    assert_eq!(fs::read_to_string(d.path().join("out.fq")).unwrap(), data);
    assert_eq!(
        report(d.path())["summary"]["before_filtering"]["total_reads"],
        30
    );
}

#[test]
fn compressed_output_is_finished_and_round_trips() {
    let (d, data) = fixture();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq.gz"])
            .output()
            .unwrap(),
    );
    success(
        command(d.path())
            .args(["-i", "out.fq.gz", "-o", "roundtrip.fq"])
            .output()
            .unwrap(),
    );
    assert_eq!(
        fs::read_to_string(d.path().join("roundtrip.fq")).unwrap(),
        data
    );
}

#[test]
fn overlap_correction_uses_encoded_phred_scores() {
    let (d, _) = fixture();
    let seq = "CGTAGCTACGATGCTAGCTAGGCTACGATCGTACGTGCATGCTAGCTAGTCAGTGCATCGATGCATCGTA";
    let rc: String = seq
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect();
    let mut wrong = rc.as_bytes().to_vec();
    let p = 20;
    wrong[p] = if wrong[p] == b'A' { b'C' } else { b'A' };
    let mut qualities = vec![b'I'; seq.len()];
    qualities[p] = b'!';
    fs::write(d.path().join("in.fq"), record("pair/1", seq, 'I')).unwrap();
    fs::write(
        d.path().join("in2.fq"),
        format!(
            "@pair/2\n{}\n+\n{}\n",
            String::from_utf8(wrong).unwrap(),
            String::from_utf8(qualities).unwrap()
        ),
    )
    .unwrap();
    success(
        command(d.path())
            .args([
                "-i", "in.fq", "-I", "in2.fq", "-o", "out1.fq", "-O", "out2.fq", "-c",
            ])
            .output()
            .unwrap(),
    );
    assert_eq!(
        fs::read_to_string(d.path().join("out2.fq")).unwrap(),
        record("pair/2", &rc, 'I')
    );
    assert_eq!(report(d.path())["filtering_result"]["corrected_bases"], 1);
    assert_eq!(report(d.path())["filtering_result"]["corrected_reads"], 1);
}

#[test]
fn stdin_stdout_interleaved_is_clean_fastq() {
    let (d, data) = fixture();
    let mut child = command(d.path())
        .args([
            "--stdin",
            "--stdout",
            "--interleaved_in",
            "-w",
            "3",
            "--batch_bytes",
            "1",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(data.as_bytes())
        .unwrap();
    let output = success(child.wait_with_output().unwrap());
    assert_eq!(output.stdout, data.as_bytes());
    assert_eq!(
        report(d.path())["summary"]["before_filtering"]["total_reads"],
        30
    );
}

#[test]
fn fragmented_stdin_matches_file_output_and_reports() {
    let (d, data) = fixture();
    let data = data.repeat(100);
    fs::write(d.path().join("in.fq"), &data).unwrap();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "expected.fq", "-f", "3"])
            .output()
            .unwrap(),
    );
    let expected = fs::read(d.path().join("expected.fq")).unwrap();
    let expected_report = report(d.path());
    for gzip in [false, true] {
        let input = if gzip {
            let split = data.len() / 2;
            [
                bgzf(&data.as_bytes()[..split]),
                bgzf(b""),
                bgzf(&data.as_bytes()[split..]),
            ]
            .concat()
        } else {
            data.as_bytes().to_vec()
        };
        let mut child = command(d.path())
            .args([
                "--stdin",
                "--stdout",
                "-f",
                "3",
                "-w",
                "4",
                "--batch_bytes",
                "4096",
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut stdin = child.stdin.take().unwrap();
        let producer = std::thread::spawn(move || {
            for part in input.chunks(37) {
                stdin.write_all(part).unwrap();
            }
        });
        let output = success(child.wait_with_output().unwrap());
        producer.join().unwrap();
        assert_eq!(output.stdout, expected);
        let actual = report(d.path());
        for section in [
            "summary",
            "filtering_result",
            "read1_before_filtering",
            "read1_after_filtering",
        ] {
            assert_eq!(
                actual[section], expected_report[section],
                "{section}, gzip={gzip}"
            );
        }
    }
}

#[test]
fn stdin_read_limit_does_not_wait_for_producer_eof() {
    let (d, data) = fixture();
    for interleaved in [false, true] {
        let mut cmd = command(d.path());
        cmd.args(["--stdin", "--stdout", "--reads_to_process", "3", "-w", "4"]);
        if interleaved {
            cmd.arg("--interleaved_in");
        }
        let mut child = cmd
            .stdin(Stdio::piped())
            .stdout(fs::File::create(d.path().join("out.fq")).unwrap())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let mut stdin = child.stdin.take().unwrap();
        stdin.write_all(data.as_bytes()).unwrap();
        // Leave the producer open: reaching the requested limit must be sufficient.
        assert!(bounded_exit(child).success());
        assert_eq!(
            fs::read_to_string(d.path().join("out.fq")).unwrap(),
            data.lines()
                .take(if interleaved { 24 } else { 12 })
                .map(|s| format!("{s}\n"))
                .collect::<String>()
        );
        drop(stdin);
    }
}

#[test]
fn stdin_stops_on_closed_output_without_waiting_for_producer() {
    for interleaved in [false, true] {
        let d = TempDir::new().unwrap();
        let mut cmd = command(d.path());
        cmd.args(["--stdin", "--stdout", "-w", "4", "--batch_bytes", "65536"]);
        if interleaved {
            cmd.arg("--interleaved_in");
        }
        let mut child = cmd
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        drop(child.stdout.take());
        let mut stdin = child.stdin.take().unwrap();
        let producer = std::thread::spawn(move || {
            let data = record("read", &"ACGT".repeat(38), 'I').repeat(4000);
            let _ = stdin.write_all(data.as_bytes());
            stdin // Keep the input open even if the child stops reading.
        });
        assert!(!bounded_exit(child).success());
        drop(producer.join().unwrap());
        assert!(!d.path().join("report.json").exists());
    }
}

#[test]
fn synchronous_paired_blocks_match_pipelined_output() {
    let (d, data) = fixture();
    let second = data.replace("@read", "@mate");
    fs::write(d.path().join("mate.fq"), second).unwrap();
    let mut expected = None;
    for backlog in [1, 4] {
        success(
            command(d.path())
                .args([
                    "-i",
                    "in.fq",
                    "-I",
                    "mate.fq",
                    "-o",
                    "out1.fq",
                    "-O",
                    "out2.fq",
                    "-w",
                    "4",
                    "--batch_bytes",
                    "1000",
                    "--max_backlog",
                    &backlog.to_string(),
                ])
                .output()
                .unwrap(),
        );
        let result = (
            fs::read(d.path().join("out1.fq")).unwrap(),
            fs::read(d.path().join("out2.fq")).unwrap(),
            report(d.path()),
        );
        if let Some(previous) = &expected {
            assert_eq!(&result, previous);
        } else {
            expected = Some(result);
        }
    }
}

#[test]
fn merge_stdout_omits_unmerged_mates_in_both_output_paths() {
    let d = TempDir::new().unwrap();
    // Both reads pass filtering but cannot overlap their reverse complements.
    let first = record("pair/1", &"A".repeat(80), 'I');
    let second = record("pair/2", &"A".repeat(80), 'I');
    fs::write(d.path().join("r1.fq"), &first).unwrap();
    fs::write(d.path().join("r2.fq"), &second).unwrap();
    for general_path in [false, true] {
        for include in [false, true] {
            let mut cmd = command(d.path());
            cmd.args(["-i", "r1.fq", "-I", "r2.fq", "--stdout", "-m", "-w", "4"]);
            if general_path {
                cmd.args(["--failed_out", "failed.fq"]);
            }
            if include {
                cmd.arg("--include_unmerged");
            }
            let output = success(cmd.output().unwrap());
            assert_eq!(
                output.stdout,
                if include {
                    format!("{first}{second}").into_bytes()
                } else {
                    Vec::new()
                }
            );
        }
    }
}

#[test]
fn quality_off_disables_n_and_average_filters() {
    let (d, _) = fixture();
    let data = record("read", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", '!');
    fs::write(d.path().join("in.fq"), &data).unwrap();
    success(
        command(d.path())
            .args([
                "--in1", "in.fq", "--out1", "out.fq", "-Q", "-e", "40", "-n", "0",
            ])
            .output()
            .unwrap(),
    );
    assert_eq!(fs::read_to_string(d.path().join("out.fq")).unwrap(), data);
}

#[test]
fn qc_only_requires_no_read_output() {
    let (d, _) = fixture();
    success(command(d.path()).args(["-i", "in.fq"]).output().unwrap());
    assert_eq!(
        report(d.path())["summary"]["after_filtering"]["total_reads"],
        30
    );
}

#[test]
fn phred64_is_converted_before_statistics_and_output() {
    let (d, _) = fixture();
    let seq = "ACGTGTCAGTACGATCGTAGCTAGCTACGATC";
    fs::write(d.path().join("in.fq"), record("read", seq, 'h')).unwrap();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq", "-6"])
            .output()
            .unwrap(),
    );
    assert_eq!(
        fs::read_to_string(d.path().join("out.fq")).unwrap(),
        record("read", seq, 'I')
    );
    assert_eq!(
        report(d.path())["summary"]["before_filtering"]["q40_bases"],
        seq.len()
    );
}

#[test]
fn explicit_zero_r2_overrides_inherited_trimming() {
    let (d, _) = fixture();
    fs::copy(d.path().join("in.fq"), d.path().join("in2.fq")).unwrap();
    success(
        command(d.path())
            .args([
                "-i", "in.fq", "-I", "in2.fq", "-o", "out.fq", "-O", "out2.fq", "-f", "5", "-F",
                "0", "-b", "20", "-B", "0",
            ])
            .output()
            .unwrap(),
    );
    assert_eq!(
        report(d.path())["read1_after_filtering"]["total_bases"],
        600
    );
    assert_eq!(
        report(d.path())["read2_after_filtering"]["total_bases"],
        960
    );
}

#[test]
fn read_limit_counts_pairs() {
    let (d, _) = fixture();
    fs::copy(d.path().join("in.fq"), d.path().join("in2.fq")).unwrap();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-I", "in2.fq", "--reads_to_process", "7"])
            .output()
            .unwrap(),
    );
    assert_eq!(
        report(d.path())["summary"]["before_filtering"]["total_reads"],
        14
    );
}

#[test]
fn umi_is_at_front_skip_follows_it_and_tag_precedes_comment() {
    let (d, _) = fixture();
    let seq = "ACGTGTCAGTACGATCGTAGCTAGCTACGATC";
    fs::write(
        d.path().join("in.fq"),
        record("read 1:N:0:ACGT+TGCA", seq, 'I'),
    )
    .unwrap();
    success(
        command(d.path())
            .args([
                "-i",
                "in.fq",
                "-o",
                "out.fq",
                "-U",
                "--umi_loc",
                "read1",
                "--umi_len",
                "4",
                "--umi_skip",
                "2",
                "--umi_prefix",
                "UMI",
                "--umi_delim",
                "_",
            ])
            .output()
            .unwrap(),
    );
    assert_eq!(
        fs::read_to_string(d.path().join("out.fq")).unwrap(),
        record("read_UMI_ACGT 1:N:0:ACGT+TGCA", &seq[6..], 'I')
    );
}

#[test]
fn index_umi_does_not_require_length_or_add_a_default_prefix() {
    let (d, _) = fixture();
    let seq = "ACGTGTCAGTACGATCGTAGCTAGCTACGATC";
    fs::write(
        d.path().join("in.fq"),
        record("read 1:N:0:ACGT+TGCA", seq, 'I'),
    )
    .unwrap();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq", "-U", "--umi_loc", "index1"])
            .output()
            .unwrap(),
    );
    assert_eq!(
        fs::read_to_string(d.path().join("out.fq")).unwrap(),
        record("read:ACGT 1:N:0:ACGT+TGCA", seq, 'I')
    );
}

#[test]
fn parallel_small_batches_preserve_order() {
    let (d, data) = fixture();
    success(
        command(d.path())
            .args([
                "-i",
                "in.fq",
                "-o",
                "out.fq",
                "-w",
                "3",
                "--batch_bytes",
                "1",
            ])
            .output()
            .unwrap(),
    );
    assert_eq!(fs::read_to_string(d.path().join("out.fq")).unwrap(), data);
}

#[test]
fn split_by_files_writes_every_record() {
    let (d, data) = fixture();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq", "-s", "3"])
            .output()
            .unwrap(),
    );
    let combined = (1..=3)
        .map(|i| fs::read_to_string(d.path().join(format!("{i:04}.out.fq"))).unwrap())
        .collect::<String>();
    assert_eq!(combined, data);
}

#[test]
fn split_by_lines_preserves_record_boundaries() {
    let (d, data) = fixture();
    let data = data.repeat(20);
    fs::write(d.path().join("in.fq"), &data).unwrap();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq", "-S", "1000"])
            .output()
            .unwrap(),
    );
    let combined = (1..=3)
        .map(|i| {
            let s = fs::read_to_string(d.path().join(format!("{i:04}.out.fq"))).unwrap();
            assert_eq!(s.lines().count() % 4, 0);
            s
        })
        .collect::<String>();
    assert_eq!(combined, data);
}

#[test]
fn dont_overwrite_checks_all_paths_before_opening_any() {
    let (d, data) = fixture();
    fs::write(d.path().join("existing.html"), "keep").unwrap();
    let output = command(d.path())
        .args([
            "-i",
            "in.fq",
            "-o",
            "out.fq",
            "--html",
            "existing.html",
            "--dont_overwrite",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(!d.path().join("out.fq").exists());
    assert_eq!(
        fs::read_to_string(d.path().join("existing.html")).unwrap(),
        "keep"
    );
    assert_eq!(fs::read_to_string(d.path().join("in.fq")).unwrap(), data);
}

#[test]
fn dont_overwrite_protects_numbered_split_files() {
    let (d, _) = fixture();
    fs::write(d.path().join("0001.out.fq"), "keep").unwrap();
    let result = command(d.path())
        .args(["-i", "in.fq", "-o", "out.fq", "-s", "3", "--dont_overwrite"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(
        fs::read_to_string(d.path().join("0001.out.fq")).unwrap(),
        "keep"
    );
    assert!(!d.path().join("0002.out.fq").exists());
}

#[test]
fn split_output_cannot_alias_a_report() {
    let (d, data) = fixture();
    fs::write(d.path().join("report.json"), "keep").unwrap();
    fs::hard_link(d.path().join("report.json"), d.path().join("0001.out.fq")).unwrap();
    let result = command(d.path())
        .args(["-i", "in.fq", "-o", "out.fq", "-s", "3"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(
        fs::read_to_string(d.path().join("report.json")).unwrap(),
        "keep"
    );
    assert!(!d.path().join("0002.out.fq").exists());
    assert_eq!(fs::read_to_string(d.path().join("in.fq")).unwrap(), data);
}

#[test]
fn split_outputs_cannot_alias_each_other() {
    for other in ["0002.out.fq", "0001.out2.fq"] {
        let (d, data) = fixture();
        fs::write(d.path().join("mate.fq"), &data).unwrap();
        fs::write(d.path().join("0001.out.fq"), "keep").unwrap();
        fs::hard_link(d.path().join("0001.out.fq"), d.path().join(other)).unwrap();
        let result = command(d.path())
            .args([
                "-i", "in.fq", "-I", "mate.fq", "-o", "out.fq", "-O", "out2.fq", "-s", "3",
            ])
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert_eq!(
            fs::read_to_string(d.path().join("0001.out.fq")).unwrap(),
            "keep"
        );
        assert!(!d.path().join("report.json").exists());
    }
}

#[cfg(unix)]
#[test]
fn dangling_split_symlink_cannot_alias_a_new_output() {
    for (link, target) in [
        ("0001.out.fq", "report.json"),
        ("0001.out.fq", "0002.out.fq"),
        ("0001.out.fq", "0001.out2.fq"),
        ("0001.out2.fq", "0001.out.fq"),
    ] {
        let (d, data) = fixture();
        fs::write(d.path().join("mate.fq"), &data).unwrap();
        std::os::unix::fs::symlink(target, d.path().join(link)).unwrap();
        let result = command(d.path())
            .args([
                "-i", "in.fq", "-I", "mate.fq", "-o", "out.fq", "-O", "out2.fq", "-S", "1000",
            ])
            .output()
            .unwrap();
        assert!(!result.status.success(), "{link} -> {target}");
        assert!(!d.path().join(target).exists());
        assert!(!d.path().join("report.json").exists());
    }
}

#[test]
fn input_cannot_be_an_output_or_report() {
    let (d, data) = fixture();
    for options in [["-o", "in.fq"], ["-j", "in.fq"], ["-h", "in.fq"]] {
        assert!(!command(d.path())
            .args(["-i", "in.fq"])
            .args(options)
            .output()
            .unwrap()
            .status
            .success());
        assert_eq!(fs::read_to_string(d.path().join("in.fq")).unwrap(), data);
    }
}

#[test]
fn malformed_and_truncated_reads_are_reported() {
    let (d, _) = fixture();
    for text in [
        "@read\nACGT\n+\nIII\n",
        "@read\nACGT\n+\nIIIII\n",
        "bad\nACGT\n+\nIIII\n",
    ] {
        fs::write(d.path().join("in.fq"), text).unwrap();
        assert!(!command(d.path())
            .args(["-i", "in.fq", "-w", "3"])
            .output()
            .unwrap()
            .status
            .success());
    }
}

#[test]
fn odd_interleaved_and_unequal_pair_counts_fail() {
    let (d, data) = fixture();
    fs::write(
        d.path().join("short.fq"),
        data.lines().take(4).collect::<Vec<_>>().join("\n") + "\n",
    )
    .unwrap();
    for options in [
        vec!["-i", "short.fq", "--interleaved_in"],
        vec!["-i", "in.fq", "-I", "short.fq"],
        vec!["-i", "short.fq", "-I", "in.fq"],
    ] {
        let output = command(d.path())
            .args(options)
            .args(["-w", "3"])
            .output()
            .unwrap();
        assert!(!output.status.success());
    }
}

#[test]
fn invalid_options_fail_before_creating_output() {
    let (d, _) = fixture();
    for options in [
        vec!["-w", "0"],
        vec!["--cut_window_size", "0"],
        vec!["-u", "101"],
        vec!["-S", "1001"],
        vec!["-m"],
    ] {
        assert!(!command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq"])
            .args(options)
            .output()
            .unwrap()
            .status
            .success());
        assert!(!d.path().join("out.fq").exists());
    }
}

#[test]
fn empty_input_has_zero_counts_and_finite_rates() {
    let (d, _) = fixture();
    fs::write(d.path().join("in.fq"), "").unwrap();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-o", "out.fq"])
            .output()
            .unwrap(),
    );
    assert_eq!(
        report(d.path())["summary"]["after_filtering"]["total_reads"],
        0
    );
    assert_eq!(
        report(d.path())["summary"]["after_filtering"]["q30_rate"],
        0.0
    );
}

#[test]
fn html_title_is_escaped() {
    let (d, _) = fixture();
    success(
        command(d.path())
            .args(["-i", "in.fq", "-R", "<script>alert('x')</script>"])
            .output()
            .unwrap(),
    );
    let html = fs::read_to_string(d.path().join("report.html")).unwrap();
    assert!(html.contains("&lt;script&gt;alert(&#39;x&#39;)&lt;/script&gt;"));
    assert!(!html.contains("<script>alert('x')</script>"));
}

#[test]
fn cli_retains_kebab_case_aliases_and_fastp_version_flag() {
    let (d, _) = fixture();
    success(
        command(d.path())
            .args([
                "--input",
                "in.fq",
                "--output",
                "out.fq",
                "--length-required",
                "20",
                "--threads",
                "1",
            ])
            .output()
            .unwrap(),
    );
    let output = success(
        Command::new(assert_cmd::cargo::cargo_bin!("fasterp"))
            .arg("-v")
            .output()
            .unwrap(),
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn reports_and_read_order_are_stable_across_workers_and_batch_boundaries() {
    let d = TempDir::new().unwrap();
    let mut input1 = String::new();
    let mut input2 = String::new();
    for i in 0..1400 {
        let mut seed = (i % 97 + 1) as u64;
        let prefix: String = (0..24)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(seed >> 32) as usize % 4] as char
            })
            .collect();
        let seq = prefix + "ACGTGTCAGTACGATCGTAGCTAGCTACGATCA";
        let mate: String = seq
            .chars()
            .rev()
            .map(|b| match b {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                _ => 'C',
            })
            .collect();
        input1.push_str(&record(
            &format!("pair{i}/1"),
            &seq,
            if i % 7 == 0 { '!' } else { 'I' },
        ));
        input2.push_str(&record(
            &format!("pair{i}/2"),
            &mate,
            if i % 11 == 0 { '!' } else { 'H' },
        ));
    }
    fs::write(d.path().join("in1.fq"), input1).unwrap();
    fs::write(d.path().join("in2.fq"), input2).unwrap();
    for options in [
        vec![],
        vec!["-D", "--dup_calc_accuracy", "1"],
        vec!["-m", "--merged_out", "merged.fq"],
        vec!["-m", "--merged_out", "merged.fq", "--include_unmerged"],
    ] {
        let mut expected = None;
        for (workers, batch, backlog) in [
            ("1", "1048576", "1"),
            ("4", "4096", "2"),
            ("8", "6144", "2"),
        ] {
            success(
                command(d.path())
                    .args([
                        "-i",
                        "in1.fq",
                        "-I",
                        "in2.fq",
                        "-o",
                        "out1.fq",
                        "-O",
                        "out2.fq",
                        "-p",
                        "-P",
                        "7",
                        "-w",
                        workers,
                        "--batch_bytes",
                        batch,
                        "--max_backlog",
                        backlog,
                    ])
                    .args(&options)
                    .output()
                    .unwrap(),
            );
            let json = report(d.path());
            assert!(!json["read1_before_filtering"]["overrepresented_sequences"]
                .as_object()
                .unwrap()
                .is_empty());
            let outputs: Vec<_> = ["out1.fq", "out2.fq", "merged.fq"]
                .iter()
                .map(|name| fs::read(d.path().join(name)).unwrap_or_default())
                .collect();
            if let Some((expected_json, expected_outputs)) = &expected {
                assert_eq!(
                    &json, expected_json,
                    "workers={workers}, options={options:?}"
                );
                assert_eq!(
                    &outputs, expected_outputs,
                    "workers={workers}, options={options:?}"
                );
            } else {
                expected = Some((json, outputs));
            }
        }
    }
}
