//! In-memory stand-in for `io` when building without the `native` feature
//! (WebAssembly). It exposes the same functions the engine uses, backed by a
//! map from path to bytes, so `compat::run` reads and writes "files" that the
//! caller provides and collects. Gzip is detected from magic bytes on input and
//! from the `.gz` extension on output, as natively.
use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Cursor, Write};
use std::sync::{Arc, Mutex, MutexGuard};

static FILES: Mutex<BTreeMap<String, Arc<[u8]>>> = Mutex::new(BTreeMap::new());
static RUN: Mutex<()> = Mutex::new(());

fn files() -> MutexGuard<'static, BTreeMap<String, Arc<[u8]>>> {
    FILES.lock().unwrap_or_else(|e| e.into_inner())
}

/// Run `f` with exactly `inputs` as the file system and return every file it
/// wrote, inputs excluded. Runs are serialized because the store is global.
pub fn with_files<T>(
    inputs: Vec<(String, Vec<u8>)>,
    f: impl FnOnce() -> Result<T>,
) -> (Result<T>, BTreeMap<String, Vec<u8>>) {
    let _run = RUN.lock().unwrap_or_else(|e| e.into_inner());
    let names: Vec<String> = inputs.iter().map(|(name, _)| name.clone()).collect();
    *files() = inputs
        .into_iter()
        .map(|(name, data)| (name, Arc::from(data)))
        .collect();
    let result = f();
    let outputs = std::mem::take(&mut *files())
        .into_iter()
        .filter(|(name, _)| !names.contains(name))
        .map(|(name, data)| (name, data.to_vec()))
        .collect();
    (result, outputs)
}

fn get(path: &str) -> Result<Arc<[u8]>> {
    files()
        .get(path)
        .cloned()
        .with_context(|| format!("No such input: {path}"))
}

pub(crate) fn is_s3_path(_path: &str) -> bool {
    false
}

pub(crate) fn set_split_decoding(_enabled: bool) {}

pub(crate) fn binned_gzip_input(_path: &str) -> bool {
    false
}

pub(crate) fn threads_per_output(_outputs: usize) -> usize {
    1
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CompressionFormat {
    None,
    Gzip,
}

impl CompressionFormat {
    pub(crate) fn from_path(path: &str) -> Self {
        let lower = path.to_lowercase();
        if path != "-"
            && [".gz", ".gzip", ".bgz", ".bgzf"]
                .iter()
                .any(|e| lower.ends_with(e))
        {
            CompressionFormat::Gzip
        } else {
            CompressionFormat::None
        }
    }
}

pub(crate) fn open_input(path: &str) -> Result<Box<dyn BufRead + Send>> {
    let data = get(path)?;
    if data.starts_with(&[0x1f, 0x8b]) {
        let decoder = flate2::read::MultiGzDecoder::new(Cursor::new(data));
        Ok(Box::new(BufReader::with_capacity(1 << 20, decoder)))
    } else if data.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
        bail!("zstd input is not supported in this build: {path}")
    } else {
        Ok(Box::new(Cursor::new(data)))
    }
}

/// A file that is stored under its path when dropped.
pub struct MemFile {
    path: String,
    data: Vec<u8>,
}

impl MemFile {
    fn new(path: &str) -> Self {
        files().insert(path.into(), Arc::from(Vec::new()));
        Self {
            path: path.into(),
            data: Vec::new(),
        }
    }
}

impl Write for MemFile {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.data.extend_from_slice(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Drop for MemFile {
    fn drop(&mut self) {
        files().insert(
            std::mem::take(&mut self.path),
            Arc::from(std::mem::take(&mut self.data)),
        );
    }
}

pub enum OutputWriter {
    Plain(MemFile),
    Gzip(GzEncoder<MemFile>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Gzip(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Gzip(w) => w.flush(),
        }
    }
}

impl OutputWriter {
    pub fn finish(self) -> Result<()> {
        if let OutputWriter::Gzip(w) = self {
            w.finish().context("Failed to finish gzip encoding")?;
        }
        Ok(())
    }
}

pub(crate) fn open_output_with_threads(
    path: &str,
    compression_level: Option<u32>,
    _parallel: bool,
    _compression_threads: usize,
) -> Result<OutputWriter> {
    let file = MemFile::new(path);
    Ok(match CompressionFormat::from_path(path) {
        CompressionFormat::Gzip => OutputWriter::Gzip(GzEncoder::new(
            file,
            Compression::new(compression_level.unwrap_or(4).min(9)),
        )),
        CompressionFormat::None => OutputWriter::Plain(file),
    })
}

pub(crate) fn create_file(path: &str) -> Result<Box<dyn Write>> {
    Ok(Box::new(MemFile::new(path)))
}

pub(crate) fn read_to_string(path: &str) -> Result<String> {
    String::from_utf8(get(path)?.to_vec()).with_context(|| format!("{path} is not UTF-8"))
}

pub(crate) fn write_file(path: &str, data: &[u8]) -> Result<()> {
    files().insert(path.into(), Arc::from(data.to_vec()));
    Ok(())
}
