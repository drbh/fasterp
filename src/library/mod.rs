//! fasterp as a library: the command line's fastp 1.4.0 compatible engine
//! on readers or readkit batch sources, returning readkit batches instead
//! of writing files.
//!
//! - [`Options`]: fastp's options, with in-memory replacements for the ones
//!   that name files.
//! - [`Run`]: one run over single-end ([`Run::single`]), paired
//!   ([`Run::paired`], [`Run::interleaved`]) or column input
//!   ([`Run::from_source`]). It yields [`Batch`]es and is itself a readkit
//!   [`BatchSource`](readkit_io::BatchSource).
//! - [`Batch`]: the passing reads as a [`readkit_batch::Batch`], plus the
//!   failed (with a [`Failure`] tag), unpaired and merged reads.
//! - [`Report`]: fastp's JSON report.
//! - [`Error`]: typed errors with an [`ErrorKind`].
//!
//! Results are the command line's for the same options and input: the
//! passing reads are what it writes to `--out1`/`--out2`, and the report is
//! what it writes with `-j`.
mod batch;
mod error;
mod options;
mod run;

pub use batch::{Batch, Failure};
pub use error::{Error, ErrorKind};
pub use options::Options;
pub use run::{Progress, Report, Run};
