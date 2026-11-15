//! A library run: the command line's engine on readers or batch sources,
//! returning batches.
use super::{Batch, Error, ErrorKind, Options};
use crate::compat::library::{ColumnInput, ColumnSource, UpstreamError};
use crate::compat::{Engine, Input};
use readkit_batch::{Batch as Records, Layout, Paired, Pool, Read, Segments, Single};
use readkit_io::{BatchSource, ReadAhead};
use serde_json::Value;
use std::io::{BufRead, Write};
use std::marker::PhantomData;

/// Idle buffers the run's pool keeps per kind: enough for every buffer in
/// flight (worker chunks, queued batches, the batch the caller holds), so
/// the steady state reuses buffers instead of allocating.
const POOL_IDLE: usize = 1024;

/// One processing run over one input: fastp 1.4.0's trimming and filtering,
/// the same engine and the same results as the `fasterp` command line, with
/// batches in place of output files.
///
/// Pull batches with [`next_batch`](Self::next_batch) (or iterate); each
/// holds the passing reads of one input batch, in input order, plus
/// fasterp's own sets (failed, unpaired, merged reads). Then call
/// [`finish`](Self::finish) for the report. A run is itself a readkit
/// [`BatchSource`] of the passing reads, so it can feed the next stage
/// without copying.
///
/// Input is read ahead, and the next batch processed, on other threads
/// while the caller works on the current one (unless `--max_backlog 1`).
/// Dropping a run stops and joins those threads.
///
/// ```no_run
/// # fn main() -> Result<(), fasterp::Error> {
/// use std::{fs::File, io::BufReader};
/// let options = fasterp::Options::new().threads(8).cut_tail(20);
/// let r1 = BufReader::new(File::open("r1.fq").unwrap());
/// let r2 = BufReader::new(File::open("r2.fq").unwrap());
/// let mut run = fasterp::Run::paired(&options, r1, r2)?;
/// while let Some(batch) = run.next_batch()? {
///     // One contiguous buffer of read 1 bases plus offsets, no parsing.
///     let bases = batch.read1().seqs().data();
///     # let _ = bases;
/// }
/// let report = run.finish()?;
/// println!("{} of {} reads passed", report.passed_reads(), report.total_reads());
/// # Ok(())
/// # }
/// ```
pub struct Run<L> {
    engine: Option<Engine>,
    pool: Pool,
    batches: u64,
    _layout: PhantomData<fn() -> L>,
}

/// How far a run has come.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Progress {
    /// Input records processed (pairs for paired input).
    pub records: u64,
    /// Reads that passed (both mates count).
    pub passed_reads: u64,
    /// Batches returned.
    pub batches: u64,
}

/// The report of a finished run: fastp's JSON report, the one the command
/// line writes with `-j`.
#[derive(Clone, Debug)]
pub struct Report {
    json: Value,
}

impl Report {
    /// Reads before filtering (both mates count).
    pub fn total_reads(&self) -> u64 {
        self.json["summary"]["before_filtering"]["total_reads"]
            .as_u64()
            .unwrap_or(0)
    }

    /// Reads that passed the filters (both mates count).
    pub fn passed_reads(&self) -> u64 {
        self.json["filtering_result"]["passed_filter_reads"]
            .as_u64()
            .unwrap_or(0)
    }

    /// The report in fastp's JSON schema.
    pub fn json(&self) -> &Value {
        &self.json
    }

    /// The report in fastp's JSON schema.
    pub fn into_json(self) -> Value {
        self.json
    }
}

impl Run<Single> {
    /// A run over single-end FASTQ.
    pub fn single(options: &Options, input: impl BufRead + Send + 'static) -> Result<Self, Error> {
        Self::start(
            options,
            1,
            false,
            true,
            Input::Streams(Box::new(input), None),
        )
    }
}

impl Run<Paired> {
    /// A run over paired FASTQ, read 1 and read 2 in two inputs.
    pub fn paired(
        options: &Options,
        read1: impl BufRead + Send + 'static,
        read2: impl BufRead + Send + 'static,
    ) -> Result<Self, Error> {
        let input = Input::Streams(Box::new(read1), Some(Box::new(read2)));
        Self::start(options, 2, false, true, input)
    }

    /// A run over paired FASTQ interleaved in one input (read 1, then read 2
    /// of each pair).
    pub fn interleaved(
        options: &Options,
        input: impl BufRead + Send + 'static,
    ) -> Result<Self, Error> {
        Self::start(
            options,
            2,
            true,
            true,
            Input::Streams(Box::new(input), None),
        )
    }
}

impl<L: Layout> Run<L> {
    /// A run over batches from `source`, with no FASTQ text in between:
    /// sequences, qualities and names are processed where the source put
    /// them.
    ///
    /// The source is read ahead by at most one batch on a thread of its
    /// own (readkit's convention; `--max_backlog 1` turns that off). Its
    /// errors end the run as [`ErrorKind::Upstream`] with the source's
    /// error as `source()`. Records must have one or two segments (for
    /// [`Segments`], fixed by the first batch; more is
    /// [`ErrorKind::Unsupported`]) and names in every batch or in none.
    /// Batches without names give batches without names.
    pub fn from_source<S>(options: &Options, source: S) -> Result<Self, Error>
    where
        S: BatchSource<L> + Send + 'static,
    {
        let read_ahead = options.args.max_backlog != Some(1);
        let mut source = ReadAhead::with_read_ahead(source, read_ahead);
        // The segment count and naming come from the first batch.
        let first = loop {
            match source.next_batch() {
                Ok(Some(batch)) if batch.is_empty() => continue,
                Ok(batch) => break batch,
                Err(error) => {
                    return Err(Error::from_engine(UpstreamError(Box::new(error)).into()))
                }
            }
        };
        let segments = L::SEGMENTS
            .or(first.as_ref().map(Records::segments))
            .unwrap_or(1);
        if !(1..=2).contains(&segments) {
            return Err(Error::new(
                ErrorKind::Unsupported,
                format!("records with {segments} segments; fasterp processes 1 or 2"),
            ));
        }
        let names = first.as_ref().is_none_or(Records::has_names);
        let input = ColumnInput {
            source: Box::new(Adapter {
                first: first.map(Records::into_segments),
                source,
            }),
            segments,
            names,
            records: 0,
        };
        Self::start(options, segments, false, names, Input::Columns(input))
    }

    fn start(
        options: &Options,
        segments: usize,
        interleaved: bool,
        names: bool,
        input: Input,
    ) -> Result<Self, Error> {
        let pool = Pool::with_max_idle(POOL_IDLE);
        let (args, library) = options.engine(segments, interleaved, &pool, names);
        let engine = Engine::open(args, input, Some(library)).map_err(Error::from_engine)?;
        Ok(Self {
            engine: Some(engine),
            pool,
            batches: 0,
            _layout: PhantomData,
        })
    }

    /// Processes the next input batch. `None` at the end of input, and after
    /// an error has been returned.
    pub fn next_batch(&mut self) -> Result<Option<Batch<L>>, Error> {
        let Some(engine) = self.engine.as_mut() else {
            return Ok(None);
        };
        if !engine.step().map_err(Error::from_engine)? {
            return Ok(None);
        }
        let collected = engine.take_output().expect("a library engine collects");
        self.batches += 1;
        Ok(Some(Batch::from_collected(collected)))
    }

    /// How far the run has come.
    pub fn progress(&self) -> Progress {
        let (records, passed) = self.engine.as_ref().map_or((0, 0), Engine::progress);
        Progress {
            records: records as u64,
            passed_reads: passed as u64,
            batches: self.batches,
        }
    }

    /// The pool the run's batches draw their buffers from. Dropped batches
    /// return their buffers here; [`Pool::stats`] shows the reuse.
    pub fn pool(&self) -> &Pool {
        &self.pool
    }

    /// Processes the rest of the input, discarding its batches, and returns
    /// the report. Fails if the input was malformed, including when
    /// [`next_batch`](Self::next_batch) already returned that error.
    pub fn finish(mut self) -> Result<Report, Error> {
        let mut engine = self.engine.take().expect("finish runs once");
        while engine.step().map_err(Error::from_engine)? {
            engine.take_output();
        }
        let finished = engine.finish().map_err(Error::from_engine)?;
        Ok(Report {
            json: finished.report,
        })
    }

    fn write_all(mut self, writers: &mut [&mut dyn Write]) -> Result<Report, Error> {
        let mut text = Vec::new();
        while let Some(batch) = self.next_batch()? {
            let records = batch.records();
            for (column, writer) in records.columns().zip(writers.iter_mut()) {
                text.clear();
                for i in 0..column.len() {
                    write_fastq(&mut text, column.read(i));
                }
                writer.write_all(&text).map_err(io_error)?;
            }
        }
        for writer in writers {
            writer.flush().map_err(io_error)?;
        }
        self.finish()
    }
}

impl Run<Single> {
    /// Writes the passing reads to `output` as FASTQ (what the command line
    /// writes to `--out1`, with bare `+` lines) and returns the report.
    pub fn write_to(self, mut output: impl Write) -> Result<Report, Error> {
        self.write_all(&mut [&mut output])
    }
}

impl Run<Paired> {
    /// Writes the passing pairs to `read1` and `read2` as FASTQ (what the
    /// command line writes to `--out1` and `--out2`, with bare `+` lines)
    /// and returns the report.
    pub fn write_to(self, mut read1: impl Write, mut read2: impl Write) -> Result<Report, Error> {
        self.write_all(&mut [&mut read1, &mut read2])
    }
}

impl<L: Layout> Iterator for Run<L> {
    type Item = Result<Batch<L>, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_batch().transpose()
    }
}

/// The run as a stage: yields the passing reads of each batch and drops
/// fasterp's other sets.
impl<L: Layout> BatchSource<L> for Run<L> {
    type Error = Error;

    fn next_batch(&mut self) -> Result<Option<Records<L>>, Error> {
        Ok(Run::next_batch(self)?.map(Batch::into_records))
    }
}

impl<L: Layout> std::fmt::Debug for Run<L> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Run")
            .field("progress", &self.progress())
            .finish()
    }
}

/// A `BatchSource` as the engine's column source, its first batch already
/// taken to learn the layout.
struct Adapter<L: Layout, S: BatchSource<L>> {
    first: Option<Records<Segments>>,
    source: ReadAhead<L, S>,
}

impl<L, S> ColumnSource for Adapter<L, S>
where
    L: Layout,
    S: BatchSource<L> + Send,
{
    fn next(&mut self) -> anyhow::Result<Option<Records<Segments>>> {
        if let Some(first) = self.first.take() {
            return Ok(Some(first));
        }
        match self.source.next_batch() {
            Ok(batch) => Ok(batch.map(Records::into_segments)),
            Err(error) => Err(UpstreamError(Box::new(error)).into()),
        }
    }
}

/// One record as FASTQ text, with a bare `+` line.
pub(crate) fn write_fastq(out: &mut Vec<u8>, read: Read<'_>) {
    out.push(b'@');
    out.extend_from_slice(read.name().unwrap_or_default());
    out.push(b'\n');
    out.extend_from_slice(read.seq());
    out.extend_from_slice(b"\n+\n");
    out.extend_from_slice(read.qual());
    out.push(b'\n');
}

fn io_error(error: std::io::Error) -> Error {
    Error::from_engine(error.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use readkit_batch::BatchBuilder;
    use readkit_io::{iter_source, source_fn};
    use std::io::Cursor;

    fn assert_send<T: Send>() {}

    #[test]
    fn runs_can_move_between_threads() {
        assert_send::<Run<Single>>();
        assert_send::<Run<Paired>>();
        assert_send::<Run<Segments>>();
        assert_send::<Batch<Paired>>();
    }

    const SEQ: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTAC";
    const QUAL: &[u8] = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

    fn fastq(reads: usize) -> Cursor<Vec<u8>> {
        let mut fastq = Vec::new();
        for i in 0..reads {
            fastq.extend_from_slice(format!("@r{i}\n").as_bytes());
            fastq.extend_from_slice(SEQ);
            fastq.extend_from_slice(b"\n+\n");
            fastq.extend_from_slice(QUAL);
            fastq.push(b'\n');
        }
        Cursor::new(fastq)
    }

    fn options() -> Options {
        Options::from_fastp_options(["--disable_adapter_trimming", "--batch_bytes", "4096"])
            .unwrap()
    }

    #[test]
    fn single_end_batches_hold_every_passing_read() {
        let mut run = Run::single(&options(), fastq(1000)).unwrap();
        let mut names = Vec::new();
        while let Some(batch) = run.next_batch().unwrap() {
            assert!(batch.failed().is_empty() && batch.merged().is_empty());
            for read in batch.reads() {
                assert_eq!((read.seq(), read.qual()), (SEQ, QUAL));
                names.push(read.name().unwrap().to_vec());
            }
        }
        assert_eq!(names.len(), 1000);
        assert_eq!(names[999], b"r999");
        assert!(run.progress().batches > 1);
        let report = run.finish().unwrap();
        assert_eq!((report.total_reads(), report.passed_reads()), (1000, 1000));
    }

    #[test]
    fn malformed_input_is_an_input_error_with_its_record() {
        let mut text = fastq(10).into_inner();
        text.extend_from_slice(b"r\nACGT\n+\nIIII\n");
        let error = Run::single(&options(), Cursor::new(text))
            .and_then(|run| run.finish())
            .unwrap_err();
        assert_eq!(error.kind(), ErrorKind::InvalidInput);
        assert_eq!(error.record(), Some(11));
        assert_eq!(error.to_string(), "Invalid FASTQ header at record 11");
    }

    #[test]
    fn invalid_options_are_reported_when_a_run_starts() {
        let options = Options::new().adapters(["ACG"]);
        let error = Run::single(&options, fastq(10)).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::InvalidOptions);
        let error = Run::single(&Options::new().merge(), fastq(10)).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::InvalidOptions);
    }

    /// A batch of `n` reads of `k` segments, named or not.
    fn columns(pool: &Pool, k: usize, n: usize, names: bool) -> Records<Segments> {
        let mut b = BatchBuilder::<Segments>::with_segments(pool, k, n);
        for i in 0..n {
            let name = format!("c{i}");
            let read = if names {
                Read::named(name.as_bytes(), SEQ, QUAL)
            } else {
                Read::new(SEQ, QUAL)
            };
            b.push(&vec![read; k]).unwrap();
        }
        b.finish()
    }

    #[derive(Debug)]
    struct SourceFailed;
    impl std::fmt::Display for SourceFailed {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("decoder failed")
        }
    }
    impl std::error::Error for SourceFailed {}

    #[test]
    fn source_errors_surface_with_their_cause() {
        let pool = Pool::new();
        for read_ahead in [true, false] {
            let mut sent = 0;
            let source = {
                let pool = pool.clone();
                source_fn(move || {
                    sent += 1;
                    if sent > 3 {
                        return Err(SourceFailed);
                    }
                    Ok(Some(
                        columns(&pool, 2, 100, true)
                            .try_into_layout::<Paired>()
                            .ok()
                            .unwrap(),
                    ))
                })
            };
            let options = if read_ahead {
                options()
            } else {
                Options::from_fastp_options(["--max_backlog", "1"]).unwrap()
            };
            let mut run = Run::from_source(&options, source).unwrap();
            let mut pairs = 0;
            let error = loop {
                match run.next_batch() {
                    Ok(Some(batch)) => pairs += batch.len(),
                    Ok(None) => panic!("the source's error was lost"),
                    Err(error) => break error,
                }
            };
            assert_eq!(pairs, 300);
            assert_eq!(error.kind(), ErrorKind::Upstream);
            let cause = std::error::Error::source(&error).unwrap();
            assert!(cause.is::<SourceFailed>(), "{cause:?}");
            assert!(run.next_batch().unwrap().is_none());
            assert!(run.finish().is_err());
        }
        // A source that fails at once fails the constructor.
        let source = source_fn(|| Err::<Option<Records<Single>>, _>(SourceFailed));
        let error = Run::from_source(&options(), source).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Upstream);
    }

    #[test]
    fn column_input_gives_the_text_input_results() {
        let pool = Pool::new();
        for names in [true, false] {
            let batches: Vec<_> = (0..5)
                .map(|_| Ok::<_, SourceFailed>(columns(&pool, 1, 200, names)))
                .collect();
            let mut run = Run::<Segments>::from_source(&options(), iter_source(batches)).unwrap();
            let mut n = 0;
            while let Some(batch) = run.next_batch().unwrap() {
                assert_eq!(batch.records().segments(), 1);
                assert_eq!(batch.records().has_names(), names);
                for spot in batch.records().spots() {
                    assert_eq!(spot.get(0).seq(), SEQ);
                    n += 1;
                }
            }
            assert_eq!(n, 1000);
            assert_eq!(run.finish().unwrap().passed_reads(), 1000);
        }
    }

    #[test]
    fn more_than_two_segments_are_unsupported() {
        let pool = Pool::new();
        let source = iter_source([Ok::<_, SourceFailed>(columns(&pool, 3, 10, true))]);
        let error = Run::<Segments>::from_source(&options(), source).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Unsupported);
    }

    #[test]
    fn a_segment_count_change_is_an_input_error() {
        let pool = Pool::new();
        let source = iter_source([
            Ok::<_, SourceFailed>(columns(&pool, 1, 10, true)),
            Ok(columns(&pool, 2, 10, true)),
        ]);
        let error = Run::<Segments>::from_source(&options(), source)
            .and_then(|run| run.finish())
            .unwrap_err();
        assert_eq!(error.kind(), ErrorKind::InvalidInput);
    }

    /// Dropping a run partway stops its threads instead of waiting for the
    /// rest of the input.
    #[test]
    fn dropping_a_run_partway_returns() {
        let pool = Pool::new();
        let mut run = Run::single(&options(), fastq(50_000)).unwrap();
        assert!(run.next_batch().unwrap().is_some());
        drop(run);
        let source = {
            let pool = pool.clone();
            source_fn(move || Ok::<_, SourceFailed>(Some(columns(&pool, 1, 100, true))))
        };
        // An endless source: only dropping can end this run.
        let mut run = Run::<Segments>::from_source(&options(), source).unwrap();
        for _ in 0..3 {
            assert!(run.next_batch().unwrap().is_some());
        }
        drop(run);
    }

    #[test]
    fn a_run_is_a_batch_source_of_passing_reads() {
        let mut run = Run::single(&options(), fastq(500)).unwrap();
        let mut n = 0;
        while let Some(batch) = BatchSource::next_batch(&mut run).unwrap() {
            n += batch.len();
        }
        assert_eq!(n, 500);
    }

    #[test]
    fn index_filters_take_sequences_in_memory() {
        let mut text = Vec::new();
        for i in 0..100 {
            let index = if i % 4 == 0 { "AAAAAAAA" } else { "CCCCCCCC" };
            text.extend_from_slice(format!("@r{i} 1:N:0:{index}\n").as_bytes());
            text.extend_from_slice(SEQ);
            text.extend_from_slice(b"\n+\n");
            text.extend_from_slice(QUAL);
            text.push(b'\n');
        }
        let options = options().filter_by_index1(["AAAAAAAA", " "]);
        let report = Run::single(&options, Cursor::new(text))
            .unwrap()
            .finish()
            .unwrap();
        assert_eq!(report.passed_reads(), 75);
    }

    #[test]
    fn write_to_writes_the_passing_reads_as_fastq() {
        let mut out = Vec::new();
        let report = Run::single(&options(), fastq(20))
            .unwrap()
            .write_to(&mut out)
            .unwrap();
        assert_eq!(out, fastq(20).into_inner());
        assert_eq!(report.passed_reads(), 20);
    }
}
