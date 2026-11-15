//! Processed batches: the passing reads plus fasterp's own record sets.
use crate::compat::library::Collected;
use readkit_batch::{Batch as Records, Column, Layout, Paired, Read, Segments, Single};
use std::any::Any;

/// Why a read failed, as carried in the low byte of the tags of
/// [`Batch::failed`].
///
/// | tag (low byte) | reason | fastp `--failed_out` tag |
/// | --- | --- | --- |
/// | 0 | [`MateFailed`](Self::MateFailed) | `paired_read_is_failing` |
/// | 1 | [`TooManyN`](Self::TooManyN) | `failed_too_many_n_bases` |
/// | 2 | [`TooShort`](Self::TooShort) | `failed_too_short` |
/// | 3 | [`TooLong`](Self::TooLong) | `failed_too_long` |
/// | 4 | [`LowQuality`](Self::LowQuality) | `failed_quality_filter` |
/// | 5 | [`LowComplexity`](Self::LowComplexity) | `failed_low_complexity` |
/// | 6 | [`AdapterDimer`](Self::AdapterDimer) | `failed_adapter_dimer` |
///
/// The high byte is the segment the read came from: 0 for read 1 (and
/// single-end reads), 1 for read 2. So `0x0102` is a read 2 that was too
/// short. Tags of [`Batch::unpaired`] use the same layout with reason 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Failure {
    /// The read passed, but its mate failed and unpaired reads are not
    /// kept ([`Options::keep_unpaired`](crate::Options::keep_unpaired)).
    MateFailed,
    /// More N bases than the limit.
    TooManyN,
    /// Shorter than the required length (including reads trimmed away).
    TooShort,
    /// Longer than the length limit.
    TooLong,
    /// Too many low quality bases, or mean quality too low.
    LowQuality,
    /// Low sequence complexity.
    LowComplexity,
    /// An adapter dimer: adapters trimmed and little sequence left.
    AdapterDimer,
}

const FAILURES: [Failure; 7] = [
    Failure::MateFailed,
    Failure::TooManyN,
    Failure::TooShort,
    Failure::TooLong,
    Failure::LowQuality,
    Failure::LowComplexity,
    Failure::AdapterDimer,
];

impl Failure {
    /// The reason in a tag of [`Batch::failed`] (its low byte).
    pub fn from_tag(tag: u16) -> Option<Self> {
        FAILURES.get(usize::from(tag & 0xff)).copied()
    }

    /// The segment in a tag of [`Batch::failed`] or [`Batch::unpaired`]
    /// (its high byte): 0 for read 1, 1 for read 2.
    pub fn segment_of_tag(tag: u16) -> usize {
        usize::from(tag >> 8)
    }

    /// This reason's tag value for a read of segment 0.
    pub fn tag(self) -> u16 {
        FAILURES.iter().position(|&f| f == self).unwrap_or(0) as u16
    }

    /// The tag fastp appends to the read name in `--failed_out`.
    pub fn fastp_name(self) -> &'static str {
        match self {
            Self::MateFailed => "paired_read_is_failing",
            Self::TooManyN => "failed_too_many_n_bases",
            Self::TooShort => "failed_too_short",
            Self::TooLong => "failed_too_long",
            Self::LowQuality => "failed_quality_filter",
            Self::LowComplexity => "failed_low_complexity",
            Self::AdapterDimer => "failed_adapter_dimer",
        }
    }
}

/// One processed batch: the reads that passed, as a readkit batch, plus the
/// reads fasterp sets aside. Every set is a readkit batch too, so all of it
/// can be handed on without copying.
///
/// The passing reads are what the command line writes to `--out1` (and
/// `--out2`), in the same order. Names are the FASTQ header without `@`;
/// the `+` line is not kept.
pub struct Batch<L> {
    records: Records<L>,
    failed: Records<Single>,
    unpaired: Records<Single>,
    merged: Records<Single>,
}

impl<L: Layout> Batch<L> {
    pub(crate) fn from_collected(c: Collected) -> Self {
        Self {
            records: relabel(c.main),
            failed: relabel(c.failed),
            unpaired: relabel(c.unpaired),
            merged: relabel(c.merged),
        }
    }

    /// The passing reads.
    pub fn records(&self) -> &Records<L> {
        &self.records
    }

    /// The passing reads, dropping the other sets.
    pub fn into_records(self) -> Records<L> {
        self.records
    }

    /// All four sets: passing, failed, unpaired and merged reads.
    pub fn into_parts(
        self,
    ) -> (
        Records<L>,
        Records<Single>,
        Records<Single>,
        Records<Single>,
    ) {
        (self.records, self.failed, self.unpaired, self.merged)
    }

    /// Number of passing records (pairs for paired input).
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether no record passed.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Reads that failed a filter, one segment per record, in the order
    /// fastp writes `--failed_out`. Empty unless
    /// [`Options::keep_failed`](crate::Options::keep_failed). Each read's
    /// tag holds its [`Failure`] and segment; see [`Failure`] for the values.
    ///
    /// As in fastp, a pair is listed only when one mate passed; the mate
    /// that passed is listed with [`Failure::MateFailed`] unless unpaired
    /// reads are kept.
    pub fn failed(&self) -> &Records<Single> {
        &self.failed
    }

    /// The failed reads with their reason and segment.
    pub fn failures(&self) -> impl ExactSizeIterator<Item = (Read<'_>, Failure, usize)> + '_ {
        let tags = self.failed.tags().unwrap_or(&[]);
        self.failed.reads().zip(tags).map(|(read, &tag)| {
            let reason = Failure::from_tag(tag).unwrap_or(Failure::MateFailed);
            (read, reason, Failure::segment_of_tag(tag))
        })
    }

    /// Passing reads whose mate failed (`--unpaired1`/`--unpaired2`), each
    /// tagged with its segment in the high byte (`0x0000` read 1, `0x0100`
    /// read 2). Empty unless [`Options::keep_unpaired`](crate::Options::keep_unpaired).
    pub fn unpaired(&self) -> &Records<Single> {
        &self.unpaired
    }

    /// Merged pairs (`--merged_out`). Empty unless
    /// [`Options::merge`](crate::Options::merge).
    pub fn merged(&self) -> &Records<Single> {
        &self.merged
    }
}

impl Batch<Single> {
    /// The passing reads' columns.
    pub fn column(&self) -> Column<'_> {
        self.records.column()
    }

    /// Iterates over the passing reads.
    pub fn reads(&self) -> impl ExactSizeIterator<Item = Read<'_>> + '_ {
        self.records.reads()
    }
}

impl Batch<Paired> {
    /// Columns of read 1 of the passing pairs.
    pub fn read1(&self) -> Column<'_> {
        self.records.r1()
    }

    /// Columns of read 2 of the passing pairs.
    pub fn read2(&self) -> Column<'_> {
        self.records.r2()
    }

    /// Iterates over the passing pairs.
    pub fn pairs(&self) -> impl ExactSizeIterator<Item = (Read<'_>, Read<'_>)> + '_ {
        self.records.pairs()
    }
}

impl<L: Layout> std::fmt::Debug for Batch<L> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Batch")
            .field("records", &self.records.len())
            .field("failed", &self.failed.len())
            .field("unpaired", &self.unpaired.len())
            .field("merged", &self.merged.len())
            .finish()
    }
}

/// The engine builds runtime-`k` batches of the run's segment count; give
/// them the run's layout. readkit converts only to fixed layouts, so the
/// layout is matched through `Any`.
pub(crate) fn relabel<L: Layout>(batch: Records<Segments>) -> Records<L> {
    let any: Box<dyn Any> = match L::SEGMENTS {
        None => Box::new(batch),
        Some(1) => Box::new(batch.try_into_layout::<Single>().expect("one segment")),
        Some(2) => Box::new(batch.try_into_layout::<Paired>().expect("two segments")),
        Some(k) => unreachable!("no layout has {k} segments"),
    };
    *any.downcast::<Records<L>>()
        .expect("the layout's own batch type")
}
