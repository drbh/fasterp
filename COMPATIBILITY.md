# fastp compatibility

The native CLI targets [fastp 1.4.0](https://github.com/OpenGene/fastp/releases/tag/v1.4.0). It accepts fastp's underscore option names and fasterp's hyphenated aliases.

## Supported behavior

- Adapter detection and trimming, quality and length filtering, poly-G/poly-X trimming, UMI extraction, overlap correction and merging.
- Single-end, paired-end and interleaved input; stdin/stdout; failed, unpaired, merged and split output.
- Phred+64 conversion, processing limits, duplication and overrepresented-sequence reports.
- Concatenated gzip input, including BGZF blocks and empty members.

FASTQ validation rejects truncated records, mismatched sequence/quality lengths and unequal mate counts. Output paths are checked before writing to prevent input or report files from being overwritten by another output.

## Differences to know

| Area | fasterp behavior |
|---|---|
| Reports | Defaults to `fasterp.json` and `fasterp.html`; styling, version and extra fields differ. |
| Insert sizes | Profiles all processed pairs; fastp samples one worker. Histograms can differ. |
| Split files | Preserves input order; records may be allocated to different files than in fastp. |
| Gzip output | Uses concatenated members for parallel local output. Compare decompressed reads; compression levels and file sizes differ between backends. |
| Memory budget | `--max_memory` controls processing batches, not total process memory. |
| Extensions | HTTP, S3, Zstandard and batch controls are fasterp features. |
| WebAssembly | Uses the existing buffer API; the full CLI compatibility path is native. |

Native Rust callers constructing `OutputWriter::ParGzip` use `fasterp::gzip::Gzip` as the format type. The `SplitWriter::new` constructor is unchanged.

## Testing

204 comparisons against fastp 1.4.0 passed on both macOS ARM64 and Ubuntu ARM64, at 1, 4 and 8 processing workers. They check ordered FASTQ records, counters and report curves. Regression tests also cover parsing boundaries, gzip members, cancellation and output-file collisions.

The WASM host tests and wasm32 release build pass. Browser execution, Linux x86 and release cross-compilation have not been verified. See [CONTRIBUTING.md](CONTRIBUTING.md) to run the checks.
