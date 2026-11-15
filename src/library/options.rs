//! Processing options for library runs.
use super::Error;
use crate::compat::library::{self as engine, Library};
use crate::Args;
use clap::parser::ValueSource;
use clap::{CommandFactory, FromArgMatches, Parser};
use std::ffi::OsString;

/// Options that name files, and what a library run offers instead. Library
/// runs read from readers or batch sources and return batches and the
/// report, so they write no file.
const PATH_OPTIONS: &[(&str, &str)] = &[
    (
        "input",
        "pass the input to Run::single, Run::paired, Run::interleaved or Run::from_source",
    ),
    ("input2", "pass both inputs to Run::paired"),
    ("stdin", "pass std::io::stdin().lock() to a Run constructor"),
    ("interleaved_in", "use Run::interleaved"),
    (
        "output",
        "passing reads are returned as batches: Run::next_batch, or Run::write_to",
    ),
    (
        "output2",
        "passing reads are returned as batches: Run::next_batch, or Run::write_to",
    ),
    (
        "stdout",
        "write the batches yourself, or use Run::write_to(std::io::stdout())",
    ),
    ("split", "split the returned batches yourself"),
    ("split_by_lines", "split the returned batches yourself"),
    ("split_prefix_digits", "split the returned batches yourself"),
    ("merged_out", "use Options::merge() and Batch::merged()"),
    (
        "unpaired1",
        "use Options::keep_unpaired() and Batch::unpaired()",
    ),
    (
        "unpaired2",
        "use Options::keep_unpaired() and Batch::unpaired()",
    ),
    (
        "failed_out",
        "use Options::keep_failed() and Batch::failed()",
    ),
    (
        "overlapped_out",
        "overlapped reads have no in-memory replacement yet",
    ),
    ("json", "Run::finish returns the report: Report::json()"),
    (
        "html",
        "Run::finish returns the report as JSON: Report::json()",
    ),
    (
        "report_title",
        "Run::finish returns the report as JSON: Report::json()",
    ),
    ("adapter_fasta", "use Options::adapters(sequences)"),
    (
        "filter_by_index1",
        "use Options::filter_by_index1(sequences)",
    ),
    (
        "filter_by_index2",
        "use Options::filter_by_index2(sequences)",
    ),
];

/// Processing options, starting from fastp's defaults.
///
/// The methods cover the common options; [`Options::from_fastp_options`]
/// takes any other fastp command line option. Options that read or write
/// files are replaced by in-memory equivalents: [`adapters`](Self::adapters)
/// for `--adapter_fasta`, [`filter_by_index1`](Self::filter_by_index1) and
/// [`filter_by_index2`](Self::filter_by_index2) for the index filter files,
/// and [`keep_failed`](Self::keep_failed), [`keep_unpaired`](Self::keep_unpaired)
/// and [`merge`](Self::merge) for the side outputs.
///
/// ```
/// let options = fasterp::Options::new()
///     .threads(4)
///     .cut_tail(20)
///     .length_required(30)
///     .keep_failed();
/// ```
#[derive(Clone, Debug)]
pub struct Options {
    pub(crate) args: Args,
    adapters: Option<Vec<Vec<u8>>>,
    index: [Option<Vec<Vec<u8>>>; 2],
    pub(crate) keep_failed: bool,
    pub(crate) keep_unpaired: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self::new()
    }
}

impl Options {
    /// fastp's defaults.
    pub fn new() -> Self {
        Self {
            args: Args::parse_from(["fasterp"]),
            adapters: None,
            index: [None, None],
            keep_failed: false,
            keep_unpaired: false,
        }
    }

    /// Options from fastp command line options, such as
    /// `["--trim_front1", "5", "--cut_right"]`.
    ///
    /// Options that name files (inputs, outputs, side outputs, reports,
    /// `--stdin`, `--stdout`, `--split`, adapter and index files) are
    /// rejected with [`ErrorKind::InvalidOptions`](crate::ErrorKind::InvalidOptions);
    /// the message names the in-memory replacement.
    pub fn from_fastp_options<I, T>(options: I) -> Result<Self, Error>
    where
        I: IntoIterator<Item = T>,
        T: Into<OsString> + Clone,
    {
        let argv =
            std::iter::once(OsString::from("fasterp")).chain(options.into_iter().map(Into::into));
        let matches = Args::command()
            .try_get_matches_from(argv)
            .map_err(|e| Error::options(e.to_string().trim_end().to_string()))?;
        for (id, instead) in PATH_OPTIONS {
            if matches.value_source(id) == Some(ValueSource::CommandLine) {
                return Err(Error::options(format!(
                    "{} is not available in library mode: {instead}",
                    flag_name(id)
                )));
            }
        }
        let args = Args::from_arg_matches(&matches).map_err(|e| Error::options(e.to_string()))?;
        Ok(Self {
            args,
            ..Self::new()
        })
    }

    /// Worker threads, 1 to 16 (default: every CPU, up to 16). Output does not
    /// depend on the count.
    pub fn threads(mut self, threads: usize) -> Self {
        self.args.threads = Some(threads);
        self
    }

    /// Trim the 5' end while the mean quality of the window is below
    /// `quality` (`--cut_front --cut_front_mean_quality`).
    pub fn cut_front(mut self, quality: u8) -> Self {
        self.args.cut_front = true;
        self.args.cut_front_mean_quality = Some(quality);
        self
    }

    /// Trim the 3' end while the mean quality of the window is below
    /// `quality` (`--cut_tail --cut_tail_mean_quality`).
    pub fn cut_tail(mut self, quality: u8) -> Self {
        self.args.cut_tail = true;
        self.args.cut_tail_mean_quality = Some(quality);
        self
    }

    /// Window size for `cut_front` and `cut_tail` (default: 4).
    pub fn cut_window_size(mut self, bases: usize) -> Self {
        self.args.cut_window_size = bases;
        self
    }

    /// Bases to trim from the front of every read (`--trim_front1`; read 2
    /// follows read 1 unless `--trim_front2` is set).
    pub fn trim_front(mut self, bases: usize) -> Self {
        self.args.trim_front1 = bases;
        self
    }

    /// Bases to trim from the tail of every read (`--trim_tail1`; read 2
    /// follows read 1 unless `--trim_tail2` is set).
    pub fn trim_tail(mut self, bases: usize) -> Self {
        self.args.trim_tail1 = bases;
        self
    }

    /// Discard reads shorter than this after trimming (default: 15).
    pub fn length_required(mut self, bases: usize) -> Self {
        self.args.length_required = bases;
        self
    }

    /// Discard reads longer than this (default: 0, no limit).
    pub fn length_limit(mut self, bases: usize) -> Self {
        self.args.length_limit = bases;
        self
    }

    /// Bases below this Phred quality are unqualified (default: 15).
    pub fn qualified_quality(mut self, phred: u8) -> Self {
        self.args.qualified_quality_phred = phred;
        self
    }

    /// Discard reads with more than this percentage of unqualified bases
    /// (default: 40).
    pub fn unqualified_percent_limit(mut self, percent: usize) -> Self {
        self.args.unqualified_percent_limit = percent;
        self
    }

    /// Discard reads whose mean quality is below this (default: 0, off).
    pub fn average_quality(mut self, phred: u8) -> Self {
        self.args.average_qual = phred;
        self
    }

    /// Discard reads with more N bases than this (default: 5).
    pub fn n_base_limit(mut self, bases: usize) -> Self {
        self.args.n_base_limit = bases;
        self
    }

    /// Adapter for read 1 instead of auto-detection.
    pub fn adapter(mut self, sequence: impl Into<String>) -> Self {
        self.args.adapter_sequence = Some(sequence.into());
        self
    }

    /// Adapter for read 2.
    pub fn adapter_r2(mut self, sequence: impl Into<String>) -> Self {
        self.args.adapter_sequence_r2 = Some(sequence.into());
        self
    }

    /// More adapters to trim from every read, in place of `--adapter_fasta`:
    /// one sequence per entry, at least six A/C/G/T bases each (checked when
    /// a run starts). Lower case is accepted.
    pub fn adapters<I, S>(mut self, sequences: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.adapters = Some(
            sequences
                .into_iter()
                .map(|s| s.as_ref().trim().to_ascii_uppercase().into_bytes())
                .collect(),
        );
        self
    }

    /// Detect adapters for paired input too, not only from read overlap.
    pub fn detect_adapter_for_pe(mut self) -> Self {
        self.args.detect_adapter_for_pe = true;
        self
    }

    /// Turn off adapter trimming.
    pub fn disable_adapter_trimming(mut self) -> Self {
        self.args.disable_adapter_trimming = true;
        self
    }

    /// Turn off the quality, mean quality and N base filters.
    pub fn disable_quality_filtering(mut self) -> Self {
        self.args.disable_quality_filtering = true;
        self
    }

    /// Turn off the length filters.
    pub fn disable_length_filtering(mut self) -> Self {
        self.args.disable_length_filtering = true;
        self
    }

    /// Trim poly-G tails for every read, not only for two-color platforms.
    pub fn trim_poly_g(mut self) -> Self {
        self.args.trim_poly_g = true;
        self
    }

    /// Never trim poly-G tails.
    pub fn disable_trim_poly_g(mut self) -> Self {
        self.args.disable_trim_poly_g = true;
        self
    }

    /// Trim homopolymer tails of any base.
    pub fn trim_poly_x(mut self) -> Self {
        self.args.trim_poly_x = true;
        self
    }

    /// Discard reads whose complexity is below `percent` (fastp's default
    /// threshold is 30).
    pub fn low_complexity_filter(mut self, percent: usize) -> Self {
        self.args.low_complexity_filter = true;
        self.args.complexity_threshold = percent;
        self
    }

    /// Correct mismatched bases in the overlap of paired reads.
    pub fn correction(mut self) -> Self {
        self.args.correction = true;
        self
    }

    /// Merge overlapping pairs (`--merge`, which also turns on correction).
    /// Merged reads are in [`Batch::merged`](crate::Batch::merged); pairs
    /// that do not merge stay in the main batch.
    pub fn merge(mut self) -> Self {
        self.args.merge = true;
        self
    }

    /// Remove duplicated reads.
    pub fn dedup(mut self) -> Self {
        self.args.dedup = true;
        self
    }

    /// Skip the duplication rate in the report, which saves up to about
    /// 1 GiB of memory.
    pub fn dont_eval_duplication(mut self) -> Self {
        self.args.dont_eval_duplication = true;
        self
    }

    /// Input qualities are Phred+64; output is Phred+33.
    pub fn phred64(mut self) -> Self {
        self.args.phred64 = true;
        self
    }

    /// Stop after this many reads or pairs (default: 0, all).
    pub fn reads_to_process(mut self, reads: usize) -> Self {
        self.args.reads_to_process = reads;
        self
    }

    /// Drop reads whose read 1 index matches one of `sequences`, in place of
    /// `--filter_by_index1` (mismatches allowed: `--filter_by_index_threshold`).
    pub fn filter_by_index1<I, S>(mut self, sequences: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.index[0] = Some(index_list(sequences));
        self
    }

    /// Drop pairs whose read 2 index matches one of `sequences`, in place of
    /// `--filter_by_index2`.
    pub fn filter_by_index2<I, S>(mut self, sequences: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.index[1] = Some(index_list(sequences));
        self
    }

    /// Keep reads that fail a filter, in [`Batch::failed`](crate::Batch::failed)
    /// with their reason as the tag (`--failed_out`). Processing takes a
    /// slower path that keeps every record until it is placed.
    pub fn keep_failed(mut self) -> Self {
        self.keep_failed = true;
        self
    }

    /// Keep the passing mate of a pair whose other mate failed, in
    /// [`Batch::unpaired`](crate::Batch::unpaired) (`--unpaired1` and
    /// `--unpaired2`). Like [`keep_failed`](Self::keep_failed), this takes
    /// the slower path.
    pub fn keep_unpaired(mut self) -> Self {
        self.keep_unpaired = true;
        self
    }

    /// The engine's arguments and settings for a run with `segments` (1 or
    /// 2) segments per record.
    pub(crate) fn engine(
        &self,
        segments: usize,
        interleaved: bool,
        pool: &readkit_batch::Pool,
        names: bool,
    ) -> (Args, Library) {
        let mut a = self.args.clone();
        let pe = segments == 2;
        a.input = engine::IN1.into();
        a.input2 = (pe && !interleaved).then(|| engine::IN2.into());
        a.interleaved_in = interleaved;
        a.output = engine::OUT1.into();
        a.output2 = pe.then(|| engine::OUT2.into());
        a.merged_out = a.merge.then(|| engine::MERGED.into());
        a.failed_out = self.keep_failed.then(|| engine::FAILED.into());
        a.unpaired1 = self.keep_unpaired.then(|| engine::UNPAIRED1.into());
        a.unpaired2 = self.keep_unpaired.then(|| engine::UNPAIRED2.into());
        a.json.clear();
        a.html.clear();
        let library = Library {
            pool: pool.clone(),
            segments,
            names,
            adapters: self.adapters.clone(),
            index: self.index.clone(),
        };
        (a, library)
    }
}

fn index_list<I: IntoIterator<Item = S>, S: AsRef<str>>(sequences: I) -> Vec<Vec<u8>> {
    sequences
        .into_iter()
        .map(|s| s.as_ref().trim().as_bytes().to_vec())
        .filter(|s| !s.is_empty())
        .collect()
}

/// `--long` for an option id.
fn flag_name(id: &str) -> String {
    Args::command()
        .get_arguments()
        .find(|a| a.get_id() == id)
        .and_then(|a| a.get_long())
        .map_or_else(|| id.to_string(), |long| format!("--{long}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each method sets what its fastp option sets.
    #[test]
    fn methods_match_their_fastp_options() {
        let cases: Vec<(Options, &[&str])> = vec![
            (Options::new().threads(3), &["-w", "3"]),
            (
                Options::new().cut_front(18),
                &["--cut_front", "--cut_front_mean_quality", "18"],
            ),
            (
                Options::new().cut_tail(25),
                &["--cut_tail", "--cut_tail_mean_quality", "25"],
            ),
            (
                Options::new().cut_window_size(6),
                &["--cut_window_size", "6"],
            ),
            (Options::new().trim_front(5), &["--trim_front1", "5"]),
            (Options::new().trim_tail(2), &["--trim_tail1", "2"]),
            (
                Options::new().length_required(30),
                &["--length_required", "30"],
            ),
            (Options::new().length_limit(200), &["--length_limit", "200"]),
            (
                Options::new().qualified_quality(20),
                &["--qualified_quality_phred", "20"],
            ),
            (
                Options::new().unqualified_percent_limit(30),
                &["--unqualified_percent_limit", "30"],
            ),
            (
                Options::new().average_quality(22),
                &["--average_qual", "22"],
            ),
            (Options::new().n_base_limit(2), &["--n_base_limit", "2"]),
            (
                Options::new().adapter("AGATCGGAAGAGC"),
                &["--adapter_sequence", "AGATCGGAAGAGC"],
            ),
            (
                Options::new().adapter_r2("AGATCGGAAGAGC"),
                &["--adapter_sequence_r2", "AGATCGGAAGAGC"],
            ),
            (
                Options::new().detect_adapter_for_pe(),
                &["--detect_adapter_for_pe"],
            ),
            (
                Options::new().disable_adapter_trimming(),
                &["--disable_adapter_trimming"],
            ),
            (
                Options::new().disable_quality_filtering(),
                &["--disable_quality_filtering"],
            ),
            (
                Options::new().disable_length_filtering(),
                &["--disable_length_filtering"],
            ),
            (Options::new().trim_poly_g(), &["--trim_poly_g"]),
            (
                Options::new().disable_trim_poly_g(),
                &["--disable_trim_poly_g"],
            ),
            (Options::new().trim_poly_x(), &["--trim_poly_x"]),
            (
                Options::new().low_complexity_filter(40),
                &["--low_complexity_filter", "--complexity_threshold", "40"],
            ),
            (Options::new().correction(), &["--correction"]),
            (Options::new().merge(), &["--merge"]),
            (Options::new().dedup(), &["--dedup"]),
            (
                Options::new().dont_eval_duplication(),
                &["--dont_eval_duplication"],
            ),
            (Options::new().phred64(), &["--phred64"]),
            (
                Options::new().reads_to_process(100),
                &["--reads_to_process", "100"],
            ),
        ];
        let defaults = format!("{:?}", Options::new().args);
        for (built, flags) in cases {
            let parsed = Options::from_fastp_options(flags).unwrap();
            let built = format!("{:?}", built.args);
            assert_eq!(built, format!("{:?}", parsed.args), "{flags:?}");
            assert_ne!(built, defaults, "{flags:?} changes nothing");
        }
    }

    /// Every option that names a file is refused, naming what replaces it.
    #[test]
    fn path_options_are_rejected_with_their_replacement() {
        let cases: &[(&[&str], &str)] = &[
            (&["-i", "in.fq"], "Run::single"),
            (&["--in2", "in2.fq"], "Run::paired"),
            (&["--stdin"], "stdin"),
            (&["--interleaved_in"], "Run::interleaved"),
            (&["-o", "out.fq"], "Run::next_batch"),
            (&["-O", "out2.fq"], "Run::next_batch"),
            (&["--stdout"], "Run::write_to"),
            (&["--split", "2"], "split"),
            (&["--split_by_lines", "4000"], "split"),
            (&["--split_prefix_digits", "3"], "split"),
            (&["--merged_out", "m.fq"], "Options::merge()"),
            (&["--unpaired1", "u.fq"], "Options::keep_unpaired()"),
            (&["--unpaired2", "u.fq"], "Options::keep_unpaired()"),
            (&["--failed_out", "f.fq"], "Options::keep_failed()"),
            (&["--overlapped_out", "o.fq"], "no in-memory"),
            (&["-j", "r.json"], "Report::json()"),
            (&["-h", "r.html"], "Report::json()"),
            (&["--report_title", "x"], "Report::json()"),
            (&["--adapter_fasta", "a.fa"], "Options::adapters"),
            (
                &["--filter_by_index1", "i.txt"],
                "Options::filter_by_index1",
            ),
            (
                &["--filter_by_index2", "i.txt"],
                "Options::filter_by_index2",
            ),
        ];
        for (flags, replacement) in cases {
            let error = Options::from_fastp_options(*flags).unwrap_err();
            assert_eq!(error.kind(), crate::ErrorKind::InvalidOptions, "{flags:?}");
            let message = error.to_string();
            assert!(
                message.contains("not available in library mode") && message.contains(replacement),
                "{flags:?}: {message}"
            );
        }
        let error = Options::from_fastp_options(["--no_such_option"]).unwrap_err();
        assert_eq!(error.kind(), crate::ErrorKind::InvalidOptions);
        // Defaults that name files (`-j fasterp.json`) are not rejected.
        assert!(Options::from_fastp_options(["-w", "2", "--max_backlog", "1"]).is_ok());
    }
}
