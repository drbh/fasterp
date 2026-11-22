//! Share the existing compression workers between two local gzip outputs.
//! Submission order is also write order. A bounded queue holds the result
//! for each independent member, including its destination mate.
use crate::gzip::Gzip;
use anyhow::{Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use gzp::{Compression, FormatSpec};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, Write};
use std::thread::{self, JoinHandle};

struct Job {
    data: Vec<u8>,
    reply: Sender<io::Result<Vec<u8>>>,
}

pub(super) struct PairedGzip {
    files: [File; 2],
    pending: [Vec<u8>; 2],
    jobs: Option<Sender<Job>>,
    workers: Vec<JoinHandle<()>>,
    results: VecDeque<(usize, Receiver<io::Result<Vec<u8>>>)>,
    capacity: usize,
}

impl PairedGzip {
    pub(super) fn new(paths: [&str; 2], level: u32, threads: usize) -> Result<Self> {
        let open = |path| {
            File::create(path).with_context(|| format!("Failed to create output file: {path}"))
        };
        let files = [open(paths[0])?, open(paths[1])?];
        let capacity = 2 * threads;
        let (jobs, work) = bounded::<Job>(capacity);
        let mut workers = Vec::with_capacity(threads);
        let level = Compression::new(level);
        for _ in 0..threads {
            let work = work.clone();
            let mut encoder = Gzip.create_compressor(level)?;
            workers.push(thread::Builder::new().spawn(move || {
                for job in work {
                    let encoded = Gzip
                        .encode(&job.data, &mut encoder, level, None, false)
                        .map_err(io::Error::other);
                    let _ = job.reply.send(encoded);
                }
            })?);
        }
        Ok(Self {
            files,
            pending: Default::default(),
            jobs: Some(jobs),
            workers,
            results: VecDeque::with_capacity(capacity),
            capacity,
        })
    }

    fn write_next(&mut self) -> Result<()> {
        let (mate, ready) = self.results.pop_front().unwrap();
        let data = ready
            .recv()
            .context("gzip worker exited before finishing a member")??;
        self.files[mate].write_all(&data)?;
        Ok(())
    }

    fn submit(&mut self, mate: usize) -> Result<()> {
        if self.results.len() == self.capacity {
            self.write_next()?;
        }
        let (reply, ready) = bounded(1);
        self.jobs
            .as_ref()
            .unwrap()
            .send(Job {
                data: std::mem::take(&mut self.pending[mate]),
                reply,
            })
            .map_err(|_| io::Error::other("gzip workers exited"))?;
        self.results.push_back((mate, ready));
        Ok(())
    }

    pub(super) fn write(&mut self, mate: usize, mut data: &[u8]) -> Result<()> {
        while !data.is_empty() {
            if self.pending[mate].len() == Gzip::DEFAULT_BUFSIZE {
                self.submit(mate)?;
            }
            let n = (Gzip::DEFAULT_BUFSIZE - self.pending[mate].len()).min(data.len());
            self.pending[mate].extend_from_slice(&data[..n]);
            data = &data[n..];
        }
        Ok(())
    }

    fn finish_members(&mut self, explicit_flush: bool) -> Result<()> {
        for mate in 0..2 {
            self.submit(mate)?;
            // The existing writer explicitly flushes its pending member before
            // finish on success. On early drop it only finishes that member.
            if explicit_flush {
                self.submit(mate)?;
            }
        }
        while !self.results.is_empty() {
            self.write_next()?;
        }
        for file in &mut self.files {
            file.flush()?;
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<()> {
        let result = self.finish_members(true);
        self.jobs.take();
        result
    }
}

impl Drop for PairedGzip {
    fn drop(&mut self) {
        if self.jobs.is_some() {
            let _ = self.finish_members(false);
            self.jobs.take();
        }
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::open_output_with_threads;

    #[test]
    fn matches_independent_writers_at_member_boundaries_and_on_drop() {
        let dir = tempfile::tempdir().unwrap();
        let paths = ["one.gz", "two.gz"].map(|name| dir.path().join(name));
        let names = paths.each_ref().map(|p| p.to_str().unwrap());
        let limit = Gzip::DEFAULT_BUFSIZE;
        let data: Vec<_> = (0..limit * 5 + 3).map(|i| (i ^ (i >> 7)) as u8).collect();
        for close in [false, true] {
            for (size1, size2) in [
                (0, 0),
                (0, 1),
                (limit - 1, limit + 1),
                (limit * 2, data.len()),
            ] {
                let mut pair = PairedGzip::new(names, 4, 2).unwrap();
                for (i, size) in [size1, size2].into_iter().enumerate() {
                    for chunk in data[..size].chunks(8193) {
                        pair.write(i, chunk).unwrap();
                    }
                }
                if close {
                    pair.finish().unwrap();
                } else {
                    drop(pair);
                }
                for (i, size) in [size1, size2].into_iter().enumerate() {
                    let path = dir.path().join("expected.gz");
                    let mut expected =
                        open_output_with_threads(path.to_str().unwrap(), Some(4), true, 1).unwrap();
                    expected.write_all(&data[..size]).unwrap();
                    if close {
                        expected.finish().unwrap();
                    } else {
                        drop(expected);
                    }
                    assert_eq!(
                        std::fs::read(&path).unwrap(),
                        std::fs::read(&paths[i]).unwrap()
                    );
                }
            }
        }
    }

    #[test]
    fn reports_output_errors_and_shuts_down_workers() {
        let dir = tempfile::tempdir().unwrap();
        let paths = ["one.gz", "two.gz"].map(|name| dir.path().join(name));
        let mut pair =
            PairedGzip::new(paths.each_ref().map(|p| p.to_str().unwrap()), 4, 2).unwrap();
        pair.files[1] = File::open(&paths[1]).unwrap();
        pair.write(1, b"reads").unwrap();
        assert!(pair.finish().is_err());
    }
}
