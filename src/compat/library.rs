//! The engine's library mode: typed errors that the public
//! [`crate::Error`] maps from, input from readkit batches, and the in-memory
//! columnar destination that fills readkit batches instead of files.
use super::block::{ColumnSlice, Records};
use super::{Failure, Read};
use anyhow::Result;
use rayon::prelude::*;
use readkit_batch::{
    Batch as RkBatch, BatchBuilder, BatchParts, Pool, PooledVec, Read as RkRead, SegmentParts,
    Segments,
};
use std::fmt;
use std::mem::MaybeUninit;
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
    /// Buffers of worker chunks' columns, which are joined into batches
    /// from `pool`.
    pub(super) scratch: Pool,
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
/// reads. The chunks of a step are joined into the step's batches on the
/// worker pool ([`StepColumns::join`]).
pub(super) struct ChunkColumns {
    columns: Columns,
    main: BatchBuilder<Segments>,
    merged: Option<BatchBuilder<Segments>>,
}

impl ChunkColumns {
    pub(super) fn new(columns: &Columns, records: usize) -> Self {
        let builder = |k| BatchBuilder::with_segments(&columns.scratch, k, records);
        Self {
            main: builder(columns.segments),
            merged: columns.merge.then(|| builder(1)),
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

/// One step's packed output: the passing and merged reads of all its worker
/// chunks, in input order, as finished batches.
pub(super) struct StepColumns {
    main: RkBatch<Segments>,
    merged: Option<RkBatch<Segments>>,
}

impl StepColumns {
    /// Joins a step's chunks, in order, into batches from the run's pool.
    /// Call on the worker pool: each chunk's columns are copied into their
    /// own region of the step's buffers in parallel, so the consumer thread
    /// copies nothing.
    pub(super) fn join(chunks: Vec<ChunkColumns>) -> Self {
        let pool = chunks.first().map(|c| c.columns.pool.clone());
        let (main, merged): (Vec<_>, Vec<_>) = chunks
            .into_iter()
            .map(|c| (c.main.finish(), c.merged.map(BatchBuilder::finish)))
            .unzip();
        let pool = pool.expect("a step has at least one chunk");
        let merged: Option<Vec<_>> = merged.into_iter().collect();
        Self {
            main: join(&pool, main),
            merged: merged.map(|m| join(&pool, m)),
        }
    }
}

/// One chunk's share of one column of a join: the column's values (or the
/// offsets after the chunk's first record) still to be written.
type Region<'a, T> = &'a mut [MaybeUninit<T>];

/// One chunk's share of one segment of a join.
struct Regions<'a> {
    offsets: Region<'a, i32>,
    seqs: Region<'a, u8>,
    quals: Region<'a, u8>,
    names: Option<(Region<'a, i32>, Region<'a, u8>)>,
}

/// One chunk of a join and, per segment, its regions and where its sequence
/// and name bytes start in the joined columns.
type Task<'a> = (&'a RkBatch<Segments>, Vec<(Regions<'a>, [i32; 2])>);

/// Splits `n` elements off the front of `slice`.
fn take_front<'a, T>(slice: &mut &'a mut [T], n: usize) -> &'a mut [T] {
    let (head, tail) = std::mem::take(slice).split_at_mut(n);
    *slice = tail;
    head
}

/// The first `len` spare elements of `v`, which must have room for them.
fn spare<T: readkit_batch::Element>(v: &mut PooledVec<T>, len: usize) -> Region<'_, T> {
    &mut v.vec_mut().spare_capacity_mut()[..len]
}

/// Bytes of the sequence (or name) column of segment `s` of `batch`.
fn bytes(batch: &RkBatch<Segments>, s: usize, names: bool) -> usize {
    let c = batch.segment(s);
    if names {
        c.names().map_or(0, |n| n.data().len())
    } else {
        c.seqs().data().len()
    }
}

/// Concatenates batches of the same segment count and naming (chunks
/// without records may lack names) into one batch from `pool`, copying
/// each batch into its own region in parallel.
fn join(pool: &Pool, batches: Vec<RkBatch<Segments>>) -> RkBatch<Segments> {
    let k = batches.first().map_or(1, RkBatch::segments);
    let batches: Vec<_> = batches.into_iter().filter(|b| !b.is_empty()).collect();
    if batches.is_empty() {
        return BatchBuilder::with_segments(pool, k, 0).finish();
    }
    let named = batches[0].has_names();
    let records: usize = batches.iter().map(|b| b.len()).sum();
    let total = |s, names| batches.iter().map(|b| bytes(b, s, names)).sum::<usize>();
    let sizes: Vec<[usize; 2]> = (0..k).map(|s| [total(s, false), total(s, true)]).collect();
    let mut out: Vec<SegmentParts> = sizes
        .iter()
        .map(|&[seqs, names]| SegmentParts {
            offsets: pool.take(records + 1),
            seqs: pool.take(seqs),
            quals: pool.take(seqs),
            names: named.then(|| (pool.take(records + 1), pool.take(names))),
        })
        .collect();
    let mut tasks: Vec<Task> = batches.iter().map(|b| (b, Vec::with_capacity(k))).collect();
    for (s, (seg, &[seq_bytes, name_bytes])) in out.iter_mut().zip(&sizes).enumerate() {
        // The offsets' first entry is the batch's start; chunks write the rest.
        seg.offsets.vec_mut().push(0);
        let mut offsets = spare(&mut seg.offsets, records);
        let mut seqs = spare(&mut seg.seqs, seq_bytes);
        let mut quals = spare(&mut seg.quals, seq_bytes);
        let mut names = seg.names.as_mut().map(|(o, v)| {
            o.vec_mut().push(0);
            (spare(o, records), spare(v, name_bytes))
        });
        let mut start = [0i32; 2];
        for (batch, regions) in &mut tasks {
            let (n, seq, name) = (batch.len(), bytes(batch, s, false), bytes(batch, s, true));
            regions.push((
                Regions {
                    offsets: take_front(&mut offsets, n),
                    seqs: take_front(&mut seqs, seq),
                    quals: take_front(&mut quals, seq),
                    names: names
                        .as_mut()
                        .map(|(o, v)| (take_front(o, n), take_front(v, name))),
                },
                start,
            ));
            // A step stays far below the 2 GiB offset limit.
            start = [start[0] + seq as i32, start[1] + name as i32];
        }
    }
    tasks.into_par_iter().for_each(|(batch, regions)| {
        for (s, (r, [seq_start, name_start])) in regions.into_iter().enumerate() {
            let c = batch.segment(s);
            r.seqs.write_copy_of_slice(c.seqs().data());
            r.quals.write_copy_of_slice(c.quals().data());
            rebase(r.offsets, c.seqs().offsets(), seq_start);
            match (r.names, c.names()) {
                (Some((offsets, bytes)), Some(names)) => {
                    bytes.write_copy_of_slice(names.data());
                    rebase(offsets, names.offsets(), name_start);
                }
                (None, None) => {}
                // Before any length is set, so nothing uninitialized is read.
                _ => panic!("worker chunks differ in naming"),
            }
        }
    });
    for (seg, &[seqs, names]) in out.iter_mut().zip(&sizes) {
        // SAFETY: the regions split each column's first `records` offsets
        // (after the pushed 0) and its `seqs`/`names` bytes into one disjoint
        // part per chunk, sized by that chunk's own columns, and every task
        // wrote its whole part: `write_copy_of_slice` writes all of it (and
        // panics on a length mismatch) and `rebase` writes `n` offsets for a
        // chunk of `n` records. So every element up to the new length is
        // initialized.
        unsafe {
            seg.offsets.vec_mut().set_len(records + 1);
            seg.seqs.vec_mut().set_len(seqs);
            seg.quals.vec_mut().set_len(seqs);
            if let Some((o, v)) = &mut seg.names {
                o.vec_mut().set_len(records + 1);
                v.vec_mut().set_len(names);
            }
        }
    }
    RkBatch::from_parts(BatchParts {
        segments: out,
        tags: None,
    })
    .expect("joined chunk columns form a valid batch")
}

/// Writes the ends of the `dst.len()` values `src` describes, moved to
/// start at `start`.
fn rebase(dst: Region<'_, i32>, src: &[i32], start: i32) {
    assert_eq!(dst.len() + 1, src.len(), "one offset per record");
    let shift = start - src[0];
    for (d, &o) in dst.iter_mut().zip(&src[1..]) {
        d.write(o + shift);
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
    /// The step's packed output, joined on the worker pool.
    joined: Option<StepColumns>,
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
            joined: None,
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

    /// The step's packed output, handed over without copying. Packed steps
    /// write nothing else to the main or merged output.
    pub(super) fn put(&mut self, step: StepColumns) {
        debug_assert!(self.main.is_empty() && self.joined.is_none());
        self.joined = Some(step);
    }

    /// The step's batches; starts the next step's.
    pub(super) fn take(&mut self) -> Collected {
        let (main, mut joined_merged) = match self.joined.take() {
            Some(StepColumns { main, merged }) => (main, merged),
            None => {
                let next = self.columns.builder(self.columns.segments, self.main.len());
                (std::mem::replace(&mut self.main, next).finish(), None)
            }
        };
        let [failed, unpaired, merged] = std::array::from_fn(|i| {
            if i == MERGED_SET {
                if let Some(merged) = joined_merged.take() {
                    return merged;
                }
            }
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A batch of `k`-segment records with `seqs[i]` as every segment's
    /// sequence of record `i` (reversed in segment 1), named or not.
    fn batch(pool: &Pool, k: usize, seqs: &[&[u8]], named: bool) -> RkBatch<Segments> {
        let mut b = BatchBuilder::with_segments(pool, k, seqs.len());
        for (i, seq) in seqs.iter().enumerate() {
            let rev: Vec<u8> = seq.iter().rev().copied().collect();
            let qual = vec![b'!' + i as u8; seq.len()];
            let name = format!("r{i}/{}", seq.len());
            let mk = |s: &[u8]| (s.to_vec(), qual[..s.len()].to_vec());
            let segs = [mk(seq), mk(&rev)];
            let record: Vec<RkRead> = segs[..k]
                .iter()
                .map(|(s, q)| {
                    if named {
                        RkRead::named(name.as_bytes(), s, q)
                    } else {
                        RkRead::new(s, q)
                    }
                })
                .collect();
            b.push(&record).unwrap();
        }
        b.finish()
    }

    type Record = Vec<(Option<Vec<u8>>, Vec<u8>, Vec<u8>)>;

    fn records(batch: &RkBatch<Segments>) -> Vec<Record> {
        batch
            .spots()
            .map(|spot| {
                spot.iter()
                    .map(|r| {
                        (
                            r.name().map(<[u8]>::to_vec),
                            r.seq().to_vec(),
                            r.qual().to_vec(),
                        )
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn join_concatenates_chunks_in_order() {
        let pool = Pool::new();
        let threads = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap();
        for k in [1, 2] {
            for named in [true, false] {
                let parts = [
                    &[&b"ACGT"[..], b"", b"GGGCC"][..],
                    &[],
                    &[b"T"],
                    &[b"", b""],
                    &[b"ACGTACGTAC", b"N"],
                ];
                let chunks: Vec<_> = parts.iter().map(|p| batch(&pool, k, p, named)).collect();
                let want: Vec<_> = chunks.iter().flat_map(records).collect();
                let joined = threads.install(|| join(&pool, chunks));
                assert_eq!(joined.segments(), k);
                assert_eq!(joined.has_names(), named);
                assert_eq!(records(&joined), want);
                for column in joined.columns() {
                    assert_eq!(column.seqs().offsets()[0], 0);
                }
            }
        }
    }

    #[test]
    fn join_of_empty_chunks_is_empty() {
        let pool = Pool::new();
        let joined = join(
            &pool,
            vec![batch(&pool, 2, &[], true), batch(&pool, 2, &[], true)],
        );
        assert_eq!((joined.len(), joined.segments()), (0, 2));
        let one = join(&pool, vec![batch(&pool, 1, &[b"AC"], true)]);
        assert_eq!(records(&one), records(&batch(&pool, 1, &[b"AC"], true)));
    }
}
