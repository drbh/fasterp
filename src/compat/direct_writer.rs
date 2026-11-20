//! Write an uncompressed local output on its own thread. The engine thread
//! hands over finished buffers instead of copying them into a `BufWriter`
//! and waiting for `write` to return, so a slow or throttled disk no longer
//! stops the engine from dispatching more work. Buffers reach the file in the
//! order they were sent.
use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Sender};
use std::fs::File;
use std::io::{self, Write};
use std::thread::{self, JoinHandle};

/// Small writes are gathered into buffers of about this size before they are
/// sent to the writer thread.
const BATCH: usize = 1 << 20;
/// Buffers waiting for the writer thread. Packed chunks are usually well
/// under a megabyte, so this holds a few batches of chunks.
const QUEUE: usize = 64;

pub(super) struct DirectWriter {
    pending: Vec<u8>,
    queue: Option<Sender<Vec<u8>>>,
    writer: Option<JoinHandle<io::Result<()>>>,
}

impl DirectWriter {
    /// `recycle` receives each buffer after it is written.
    pub(super) fn create(path: &str, recycle: fn(Vec<u8>)) -> Result<Self> {
        let mut file =
            File::create(path).with_context(|| format!("Failed to create output file: {path}"))?;
        let (queue, buffers) = bounded::<Vec<u8>>(QUEUE);
        let writer = thread::Builder::new().spawn(move || {
            for buffer in buffers {
                file.write_all(&buffer)?;
                recycle(buffer);
            }
            file.flush()
        })?;
        Ok(Self {
            pending: Vec::new(),
            queue: Some(queue),
            writer: Some(writer),
        })
    }

    /// The writer thread only stops early after an error, which `join` returns.
    fn join(&mut self) -> Result<()> {
        self.queue.take();
        match self.writer.take().map(JoinHandle::join) {
            Some(Ok(result)) => Ok(result?),
            Some(Err(_)) => Err(anyhow!("output writer thread panicked")),
            None => Ok(()),
        }
    }

    fn send(&mut self, buffer: Vec<u8>) -> Result<()> {
        if buffer.is_empty() {
            return Ok(());
        }
        let Some(queue) = &self.queue else {
            return Err(anyhow!("output writer already stopped"));
        };
        if queue.send(buffer).is_err() {
            self.join()?;
            return Err(anyhow!("output writer stopped"));
        }
        Ok(())
    }

    /// Queue a whole buffer without copying it, after any gathered writes.
    pub(super) fn write_owned(&mut self, buffer: Vec<u8>) -> Result<()> {
        let pending = std::mem::take(&mut self.pending);
        self.send(pending)?;
        self.send(buffer)
    }

    pub(super) fn write(&mut self, bytes: &[u8]) -> Result<()> {
        if self.pending.capacity() == 0 {
            self.pending.reserve(BATCH);
        }
        self.pending.extend_from_slice(bytes);
        if self.pending.len() >= BATCH {
            let pending = std::mem::take(&mut self.pending);
            self.send(pending)?;
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<()> {
        let pending = std::mem::take(&mut self.pending);
        self.send(pending)?;
        self.join()
    }
}

impl Drop for DirectWriter {
    fn drop(&mut self) {
        let pending = std::mem::take(&mut self.pending);
        let _ = self.send(pending);
        let _ = self.join();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_write_order_across_owned_and_gathered_writes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.fq");
        let mut writer = DirectWriter::create(path.to_str().unwrap(), drop).unwrap();
        let mut expected = Vec::new();
        for i in 0..2000usize {
            let bytes: Vec<u8> = (0..(i * 37) % 3000).map(|j| (i + j) as u8).collect();
            expected.extend_from_slice(&bytes);
            if i % 3 == 0 {
                writer.write_owned(bytes).unwrap();
            } else {
                writer.write(&bytes).unwrap();
            }
        }
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), expected);
    }

    #[test]
    fn writes_pending_bytes_on_drop() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.fq");
        let mut writer = DirectWriter::create(path.to_str().unwrap(), drop).unwrap();
        writer.write(b"@r\nA\n+\nI\n").unwrap();
        drop(writer);
        assert_eq!(std::fs::read(&path).unwrap(), b"@r\nA\n+\nI\n");
    }

    #[test]
    fn reports_write_errors() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.fq");
        std::fs::write(&path, b"").unwrap();
        let mut writer = DirectWriter::create(path.to_str().unwrap(), drop).unwrap();
        // Swap in a read-only handle so every write fails.
        let (queue, buffers) = bounded::<Vec<u8>>(QUEUE);
        writer.queue.take();
        writer.writer.take().unwrap().join().unwrap().unwrap();
        let mut file = File::open(&path).unwrap();
        writer.queue = Some(queue);
        writer.writer = Some(thread::spawn(move || {
            for buffer in buffers {
                file.write_all(&buffer)?;
            }
            Ok(())
        }));
        let mut result = Ok(());
        for _ in 0..QUEUE * 4 {
            result = writer.write_owned(vec![b'A'; 16]);
            if result.is_err() {
                break;
            }
        }
        assert!(result.is_err() || writer.finish().is_err());
    }
}
