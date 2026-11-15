//! The library's error type.
use crate::compat::{InputError, OptionsError, UpstreamError};
use std::fmt;

type BoxError = Box<dyn std::error::Error + Send + Sync + 'static>;

/// What went wrong. Match on this rather than on messages.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// The options are invalid or inconsistent, or name a file (library
    /// runs take readers and return batches).
    InvalidOptions,
    /// The input is malformed: not FASTQ, truncated, mates out of step, bad
    /// quality bytes. [`Error::record`] gives the position when known.
    InvalidInput,
    /// Reading the input failed; the cause is [`source()`](std::error::Error::source).
    Io,
    /// The [`BatchSource`](readkit_io::BatchSource) given to
    /// [`Run::from_source`](crate::Run::from_source) failed; its error is
    /// `source()`.
    Upstream,
    /// The input has a layout fasterp cannot process (more than two
    /// segments per record).
    Unsupported,
    /// A bug in fasterp.
    Internal,
}

impl ErrorKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::InvalidOptions => "invalid options",
            Self::InvalidInput => "invalid input",
            Self::Io => "input could not be read",
            Self::Upstream => "upstream source failed",
            Self::Unsupported => "unsupported input layout",
            Self::Internal => "internal error",
        }
    }
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// An error from a library run: a [kind](ErrorKind), a message, the record
/// it concerns for input errors, and the underlying cause as
/// [`source()`](std::error::Error::source).
pub struct Error {
    kind: ErrorKind,
    message: String,
    record: Option<u64>,
    cause: Option<BoxError>,
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
            record: None,
            cause: None,
        }
    }

    pub(crate) fn options(message: impl Into<String>) -> Self {
        Self::new(ErrorKind::InvalidOptions, message)
    }

    /// Maps an engine error, whose internals use anyhow, to its kind.
    pub(crate) fn from_engine(error: anyhow::Error) -> Self {
        let error = match error.downcast::<UpstreamError>() {
            Ok(upstream) => {
                return Self {
                    cause: Some(upstream.0),
                    ..Self::new(ErrorKind::Upstream, ErrorKind::Upstream.as_str())
                }
            }
            Err(error) => error,
        };
        let message = format!("{error:#}");
        if let Some(input) = error.downcast_ref::<InputError>() {
            return Self {
                record: input.record().map(|r| r as u64),
                ..Self::new(ErrorKind::InvalidInput, message)
            };
        }
        if error.is::<OptionsError>() {
            return Self::options(message);
        }
        if error.chain().any(|e| e.is::<std::io::Error>()) {
            let cause: BoxError = match error.downcast::<std::io::Error>() {
                Ok(io) => Box::new(io),
                Err(other) => other.into(),
            };
            return Self {
                cause: Some(cause),
                ..Self::new(ErrorKind::Io, message)
            };
        }
        Self::new(ErrorKind::Internal, message)
    }

    /// The error kind.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// For [`ErrorKind::InvalidInput`]: the 1-based record (read, or pair
    /// for paired files) the error was found at, when known.
    pub fn record(&self) -> Option<u64> {
        self.record
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Error")
            .field("kind", &self.kind)
            .field("message", &self.message)
            .field("record", &self.record)
            .field("cause", &self.cause)
            .finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.cause.as_deref().map(|e| e as _)
    }
}
