//! The engine's library mode: typed errors that the public
//! [`crate::Error`] maps from, input from readkit batches, and the in-memory
//! columnar destination that fills readkit batches instead of files.
use super::block::{ColumnSlice, Records};
use super::{Failure, Read};
use anyhow::Result;
use readkit_batch::{Batch as RkBatch, BatchBuilder, Pool, Read as RkRead, Segments};
use std::fmt;
use std::sync::Arc;

/// Malformed input, with the 1-based record number when known. Displays as
/// the message the command line has always printed.
#[derive(Debug)]
pub(crate) struct InputError {
    record: Option<usize>,
    message: String,
}

impl InputError {
    pub(crate) fn new(record: Option<usize>, message: impl Into<String>) -> Self {
        Self {
            record,
            message: message.into(),
        }
    }
    /// `"{what} at record {record}"`.
    pub(crate) fn at(record: usize, what: &str) -> Self {
        Self::new(Some(record), format!("{what} at record {record}"))
    }
    pub(crate) fn truncated(record: usize) -> Self {
        Self::new(Some(record), format!("Truncated FASTQ record {record}"))
    }
    pub(crate) fn odd_interleaved() -> Self {
        Self::new(None, "Odd read count in interleaved FASTQ")
    }
    pub(crate) fn mate_counts() -> Self {
        Self::new(None, "Paired input files have different read counts")
    }
    pub(crate) fn record(&self) -> Option<usize> {
        self.record
    }
}

impl fmt::Display for InputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for InputError {}

/// Invalid or inconsistent options.
#[derive(Debug)]
pub(crate) struct OptionsError(String);

impl OptionsError {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
    /// Marks an option check's error as an options error, keeping its text.
    pub(crate) fn wrap(error: anyhow::Error) -> anyhow::Error {
        if error.is::<Self>() {
            return error;
        }
        Self(format!("{error:#}")).into()
    }
}

impl fmt::Display for OptionsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for OptionsError {}

/// The error of a `from_source` source, kept whole as the cause.
pub(crate) struct UpstreamError(pub(crate) Box<dyn std::error::Error + Send + Sync + 'static>);

impl fmt::Debug for UpstreamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("UpstreamError").field(&self.0).finish()
    }
}

impl fmt::Display for UpstreamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "upstream source failed: {}", self.0)
    }
}

impl std::error::Error for UpstreamError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&*self.0)
    }
}

/// Batches handed to the engine without FASTQ text. Errors are
/// [`UpstreamError`]s or [`InputError`]s.
pub(crate) trait ColumnSource: Send {
    fn next(&mut self) -> Result<Option<RkBatch<Segments>>>;
}

/// The engine's view of a column source: checks that every batch has the
/// run's segment count and naming.
pub(crate) struct ColumnInput {
    pub(crate) source: Box<dyn ColumnSource>,
    pub(crate) segments: usize,
    pub(crate) names: bool,
    /// Records before the current batch, for error positions.
    pub(crate) records: usize,
}

impl ColumnInput {
    /// The next nonempty batch.
    pub(super) fn next(&mut self) -> Result<Option<Arc<RkBatch<Segments>>>> {
        loop {
            let Some(batch) = self.source.next()? else {
                return Ok(None);
            };
            if batch.segments() != self.segments {
                return Err(InputError::new(
                    Some(self.records + 1),
                    format!(
                        "Batch with {} segments in a run of {}",
                        batch.segments(),
                        self.segments
                    ),
                )
                .into());
            }
            if batch.is_empty() {
                continue;
            }
            if batch.has_names() != self.names {
                return Err(InputError::new(
                    Some(self.records + 1),
                    "Batches must all have read names or all have none",
                )
                .into());
            }
            self.records += batch.len();
            return Ok(Some(Arc::new(batch)));
        }
    }
}

/// Records `range` of every segment of `batch`; the second is empty for
/// single-end batches.
pub(super) fn records(batch: &Arc<RkBatch<Segments>>, start: usize, end: usize) -> [Records; 2] {
    std::array::from_fn(|segment| {
        if segment < batch.segments() {
            Records::from_columns(ColumnSlice {
                batch: Arc::clone(batch),
                segment,
                start,
                end,
            })
        } else {
            Records::default()
        }
    })
}

/// A record of a column batch as an engine read, borrowing every field. The
/// name has no FASTQ `@` (the engine runs with bare names) and is empty
/// when the batch has none; the separator is a bare `+`.
pub(super) fn read_at(batch: &RkBatch<Segments>, segment: usize, index: usize) -> Read<'_> {
    let column = batch.segment(segment);
    let read = column.read(index);
    Read {
        name: read.name().unwrap_or_default().into(),
        seq: read.seq().into(),
        plus: (&b"+"[..]).into(),
        qual: read.qual().into(),
        ..Read::default()
    }
}

/// Placeholder paths that stand for the in-memory destinations, so the
/// engine routes library output exactly as it routes files. They contain
/// `://` so no file system check touches them.
pub(crate) const IN1: &str = "fasterp-library://in1";
pub(crate) const IN2: &str = "fasterp-library://in2";
pub(crate) const OUT1: &str = "fasterp-library://out1";
pub(crate) const OUT2: &str = "fasterp-library://out2";
pub(crate) const MERGED: &str = "fasterp-library://merged";
pub(crate) const FAILED: &str = "fasterp-library://failed";
pub(crate) const UNPAIRED1: &str = "fasterp-library://unpaired1";
pub(crate) const UNPAIRED2: &str = "fasterp-library://unpaired2";

/// Library settings the engine takes besides its `Args`.
pub(crate) struct Library {
    /// Buffers of every output batch.
    pub(crate) pool: Pool,
    /// Segments of the main output (1 or 2).
    pub(crate) segments: usize,
    /// Whether output batches carry names (input from a batch without
    /// names has none).
    pub(crate) names: bool,
    /// `--adapter_fasta` sequences.
    pub(crate) adapters: Option<Vec<Vec<u8>>>,
    /// `--filter_by_index1` and `--filter_by_index2` sequences.
    pub(crate) index: [Option<Vec<Vec<u8>>>; 2],
}

/// What workers need to build column output.
#[derive(Clone)]
pub(super) struct Columns {
    pub(super) pool: Pool,
    pub(super) segments: usize,
    pub(super) names: bool,
    /// Names carry FASTQ's `@` (text input), which batches leave out.
    pub(super) strip_at: bool,
    pub(super) merge: bool,
}

impl Columns {
    fn builder(&self, segments: usize, records: usize) -> BatchBuilder<Segments> {
        BatchBuilder::with_segments(&self.pool, segments, records)
    }

    fn as_read<'a>(&self, r: &'a Read<'_>) -> RkRead<'a> {
        let (seq, qual) = (r.seq.as_ref(), r.qual.as_ref());
        if !self.names {
            return RkRead::new(seq, qual);
        }
        let name = r.name.as_ref();
        let name = if self.strip_at {
            name.strip_prefix(b"@").unwrap_or(name)
        } else {
            name
        };
        RkRead::named(name, seq, qual)
    }

    fn push(&self, builder: &mut BatchBuilder<Segments>, reads: &[&Read], tag: Option<u16>) {
        let record: [RkRead; 2] =
            std::array::from_fn(|i| self.as_read(reads[i.min(reads.len() - 1)]));
        // Reads reaching output passed validation, so their sequence and
        // quality lengths match; a step stays far below the 2 GiB limit.
        builder
            .push(&record[..reads.len()])
            .expect("validated reads fit a batch");
        if let Some(tag) = tag {
            builder.set_tag(tag);
        }
    }
}

/// Output of one worker chunk (the packed path): passing reads and merged
/// reads, appended to the step's batches in input order.
pub(super) struct ChunkColumns {
    columns: Columns,
    main: BatchBuilder<Segments>,
    merged: Option<BatchBuilder<Segments>>,
}

impl ChunkColumns {
    pub(super) fn new(columns: &Columns, records: usize) -> Self {
        Self {
            main: columns.builder(columns.segments, records),
            merged: columns.merge.then(|| columns.builder(1, records)),
            columns: columns.clone(),
        }
    }
    pub(super) fn main(&mut self, reads: &[Option<Read>; 2]) {
        let reads: Vec<&Read> = reads.iter().flatten().collect();
        self.columns.push(&mut self.main, &reads, None);
    }
    pub(super) fn merged(&mut self, read: &Read) {
        if let Some(merged) = &mut self.merged {
            self.columns.push(merged, &[read], None);
        }
    }
}

/// Batches of one step: the passing reads and fasterp's own sets.
pub(crate) struct Collected {
    pub(crate) main: RkBatch<Segments>,
    pub(crate) failed: RkBatch<Segments>,
    pub(crate) unpaired: RkBatch<Segments>,
    pub(crate) merged: RkBatch<Segments>,
}

/// The in-memory destination: fills one step's batches.
pub(super) struct Destination {
    columns: Columns,
    main: BatchBuilder<Segments>,
    /// Failed, unpaired and merged reads, when kept.
    sets: [Option<BatchBuilder<Segments>>; 3],
    keep: [bool; 3],
    /// Empty batches of sets that are not kept, without recycled buffers.
    empty: Pool,
    hint: usize,
}

const FAILED_SET: usize = 0;
const UNPAIRED_SET: usize = 1;
const MERGED_SET: usize = 2;

impl Destination {
    pub(super) fn new(columns: Columns, keep_failed: bool, keep_unpaired: bool) -> Self {
        let keep = [keep_failed, keep_unpaired, columns.merge];
        Self {
            main: columns.builder(columns.segments, 0),
            sets: keep.map(|keep| keep.then(|| columns.builder(1, 0))),
            keep,
            columns,
            empty: Pool::with_max_idle(0),
            hint: 0,
        }
    }

    /// A passing pair (or single read) for the main output.
    pub(super) fn main(&mut self, reads: &[Option<Read>; 2]) {
        let reads: Vec<&Read> = reads.iter().flatten().collect();
        self.columns.push(&mut self.main, &reads, None);
    }

    /// One read for the destination that `path` stands for. Failed reads are
    /// tagged with their reason and segment, unpaired reads with their
    /// segment (see [`crate::Failure`]).
    pub(super) fn write(&mut self, path: &str, read: &Read, failure: Option<(Failure, usize)>) {
        let (set, tag) = match path {
            MERGED => (MERGED_SET, None),
            UNPAIRED1 => (UNPAIRED_SET, Some(0)),
            UNPAIRED2 => (UNPAIRED_SET, Some(1 << 8)),
            FAILED => {
                let (reason, segment) = failure.unwrap_or((Failure::Pass, 0));
                (FAILED_SET, Some(reason as u16 | (segment as u16) << 8))
            }
            _ => return,
        };
        if let Some(builder) = &mut self.sets[set] {
            self.columns.push(builder, &[read], tag);
        }
    }

    /// Appends a worker chunk's output.
    pub(super) fn append(&mut self, chunk: ChunkColumns) {
        append(&mut self.main, &chunk.main.finish());
        if let (Some(dst), Some(src)) = (&mut self.sets[MERGED_SET], chunk.merged) {
            append(dst, &src.finish());
        }
    }

    /// The step's batches; starts the next step's.
    pub(super) fn take(&mut self) -> Collected {
        self.hint = self.main.len();
        let next = self.columns.builder(self.columns.segments, self.hint);
        let main = std::mem::replace(&mut self.main, next).finish();
        let [failed, unpaired, merged] = std::array::from_fn(|i| {
            let next = self.keep[i].then(|| self.columns.builder(1, 0));
            match std::mem::replace(&mut self.sets[i], next) {
                Some(builder) => builder.finish(),
                None => BatchBuilder::with_segments(&self.empty, 1, 0).finish(),
            }
        });
        Collected {
            main,
            failed,
            unpaired,
            merged,
        }
    }
}

/// Copies `src`'s records onto `dst` (same segment count and naming).
fn append(dst: &mut BatchBuilder<Segments>, src: &RkBatch<Segments>) {
    let k = src.segments();
    let columns: Vec<_> = src.columns().collect();
    for i in 0..src.len() {
        let record: [RkRead; 2] = std::array::from_fn(|s| columns[s.min(k - 1)].read(i));
        dst.push(&record[..k]).expect("chunk records fit a batch");
    }
}
