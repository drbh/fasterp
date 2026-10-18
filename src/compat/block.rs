//! Bounded FASTQ batches. Fields borrow their batch until worker output is packed.
//! The slow boundary path accepts the same strict four-line format as Reader.
use super::{Read, Reader};
use anyhow::{ensure, Result};
use std::borrow::Cow;

#[derive(Default)]
pub(super) struct Records {
    data: Vec<u8>,
    ends: Vec<[usize; 4]>,
}

impl Records {
    pub(super) fn len(&self) -> usize {
        self.ends.len()
    }

    fn field(&self, index: usize, field: usize) -> &[u8] {
        let start = if field > 0 {
            self.ends[index][field - 1]
        } else if index > 0 {
            self.ends[index - 1][3]
        } else {
            0
        };
        let line = &self.data[start..self.ends[index][field]];
        let line = line.strip_suffix(b"\n").unwrap_or(line);
        line.strip_suffix(b"\r").unwrap_or(line)
    }

    pub(super) fn sequence(&self, index: usize) -> &[u8] {
        self.field(index, 1)
    }
    pub(super) fn name(&self, index: usize) -> &[u8] {
        self.field(index, 0)
    }

    pub(super) fn truncate(&mut self, count: usize) {
        if count < self.len() {
            self.data.truncate(if count == 0 {
                0
            } else {
                self.ends[count - 1][3]
            });
            self.ends.truncate(count);
        }
    }

    pub(super) fn read(&self, index: usize, ordinal: usize, phred64: bool) -> Result<Read<'_>> {
        let mut start = if index == 0 {
            0
        } else {
            self.ends[index - 1][3]
        };
        let fields = self.ends[index].map(|end| {
            let line = &self.data[start..end];
            start = end;
            let line = line.strip_suffix(b"\n").unwrap_or(line);
            Cow::Borrowed(line.strip_suffix(b"\r").unwrap_or(line))
        });
        let [name, seq, plus, qual] = fields;
        let mut read = Read {
            name,
            seq,
            plus,
            qual,
            ..Read::default()
        };
        read.validate(ordinal + 1, phred64)?;
        Ok(read)
    }

    fn read_one(&mut self, reader: &mut Reader) -> Result<bool> {
        let mut ends = [0; 4];
        for (line, end) in ends.iter_mut().enumerate() {
            let n = reader.inner.read_until(b'\n', &mut self.data)?;
            if n == 0 && line == 0 {
                return Ok(false);
            }
            ensure!(n > 0, "Truncated FASTQ record {}", reader.record + 1);
            *end = self.data.len();
        }
        self.ends.push(ends);
        reader.record += 1;
        Ok(true)
    }
}

#[inline]
fn next_ends(newlines: &mut impl Iterator<Item = usize>) -> Option<[usize; 4]> {
    Some([
        newlines.next()? + 1,
        newlines.next()? + 1,
        newlines.next()? + 1,
        newlines.next()? + 1,
    ])
}

pub(super) fn read_batch(
    first: &mut Reader,
    second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
) -> Result<Option<[Records; 2]>> {
    read_batch_inner::<false, false>(first, second, budget, limit, usize::MAX)
}

pub(super) fn read_interleaved_batch(
    first: &mut Reader,
    budget: usize,
    pairs: usize,
) -> Result<Option<[Records; 2]>> {
    read_batch_inner::<false, true>(first, None, budget, pairs, usize::MAX)
}

pub(super) fn read_sample_batch(
    first: &mut Reader,
    second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
    base_limit: usize,
) -> Result<Option<[Records; 2]>> {
    read_batch_inner::<true, false>(first, second, budget, limit, base_limit)
}

fn read_batch_inner<const SAMPLE: bool, const INTERLEAVED: bool>(
    first: &mut Reader,
    mut second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
    base_limit: usize,
) -> Result<Option<[Records; 2]>> {
    let mut blocks: [Records; 2] = std::array::from_fn(|_| Records::default());
    let mut bases = 0;
    let stride = if INTERLEAVED { 2 } else { 1 };
    while blocks[0].len() / stride < limit
        && blocks.iter().map(|b| b.data.len()).sum::<usize>() < budget
        && (!SAMPLE || bases < base_limit)
    {
        // Scan each existing input buffer once, then copy a contiguous prefix.
        // Select paired prefixes together so asymmetric read lengths still obey
        // the byte budget (with the same one-pair overshoot as the owned reader).
        let remaining = limit - blocks[0].len() / stride;
        let b1 = first.inner.fill_buf()?;
        let b2 = match second.as_mut() {
            Some(r) => Some(r.inner.fill_buf()?),
            None => None,
        };
        let mut newlines1 = memchr::memchr_iter(b'\n', b1);
        let mut newlines2 = b2.map(|b| memchr::memchr_iter(b'\n', b));
        let base = [blocks[0].data.len(), blocks[1].data.len()];
        let used = base[0] + base[1];
        let mut count = 0;
        let mut consumed = [0; 2];
        let mut full = false;
        while count < remaining {
            let Some(ends1) = next_ends(&mut newlines1) else {
                break;
            };
            let ends2 = if INTERLEAVED {
                let Some(ends) = next_ends(&mut newlines1) else {
                    break;
                };
                Some(ends)
            } else if let Some(lines) = &mut newlines2 {
                let Some(ends) = next_ends(lines) else {
                    break;
                };
                Some(ends)
            } else {
                None
            };
            let next_bytes = if INTERLEAVED {
                ends2.unwrap()[3]
            } else {
                ends1[3] + ends2.map_or(0, |e| e[3])
            };
            if used + next_bytes > budget && blocks[0].len() > 0 {
                full = true;
                break;
            }
            consumed[0] = ends1[3];
            blocks[0].ends.push(ends1.map(|n| n + base[0]));
            if let Some(ends) = ends2 {
                let mate = usize::from(!INTERLEAVED);
                consumed[mate] = ends[3];
                blocks[mate].ends.push(ends.map(|n| n + base[mate]));
            }
            count += 1;
            if SAMPLE {
                let seq = &b1[ends1[0]..ends1[1] - 1];
                bases += seq.strip_suffix(b"\r").unwrap_or(seq).len();
                if bases >= base_limit {
                    full = true;
                    break;
                }
            }
        }
        if count > 0 {
            blocks[0].data.extend_from_slice(&b1[..consumed[0]]);
            if let Some(b) = b2 {
                blocks[1].data.extend_from_slice(&b[..consumed[1]]);
            }
            first.inner.consume(consumed[0]);
            first.record += count * stride;
            if let Some(r) = second.as_mut() {
                r.inner.consume(consumed[1]);
                r.record += count;
            }
            if full {
                break;
            }
        } else if full {
            break;
        } else {
            // A record crosses a buffer boundary, or EOF needs checking. This
            // path also preserves a nonempty final quality line without LF.
            let has_first = blocks[0].read_one(first)?;
            if INTERLEAVED && has_first {
                ensure!(
                    blocks[0].read_one(first)?,
                    "Odd read count in interleaved FASTQ"
                );
            }
            if let Some(r) = second.as_mut() {
                let has_second = blocks[1].read_one(r)?;
                ensure!(
                    has_first == has_second,
                    "Paired input files have different read counts"
                );
            }
            if !has_first {
                break;
            }
            if SAMPLE {
                bases += blocks[0].sequence(blocks[0].len() - stride).len();
            }
        }
    }
    Ok((blocks[0].len() > 0).then_some(blocks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn reader(data: &[u8], capacity: usize, phred64: bool) -> Reader {
        Reader {
            inner: Box::new(BufReader::with_capacity(
                capacity,
                Cursor::new(data.to_vec()),
            )),
            record: 0,
            phred64,
        }
    }

    fn normalized(data: &[u8], phred64: bool) -> Vec<u8> {
        let mut r = reader(data, 17, phred64);
        let mut out = Vec::new();
        while let Some(read) = r.next().unwrap() {
            read.write(&mut out, None).unwrap();
        }
        out
    }

    #[test]
    fn interleaved_blocks_match_owned_records_at_every_buffer_boundary() {
        let data = format!(
            "@first/1\r\nACG\r\n+first\r\nIII\r\n@first/2\n{}\n+\n{}\n@empty/1\n\n+\n\n@empty/2\nT\n+\nH",
            "G".repeat(97), "I".repeat(97)
        );
        for phred64 in [false, true] {
            let expected = normalized(data.as_bytes(), phred64);
            for capacity in 1..=data.len() + 1 {
                for budget in [1, 31, 256, 20000] {
                    let mut source = reader(data.as_bytes(), capacity, phred64);
                    let mut actual = Vec::new();
                    let mut ordinal = 0;
                    while let Some(batch) =
                        read_interleaved_batch(&mut source, budget, usize::MAX).unwrap()
                    {
                        assert_eq!(batch[0].len() % 2, 0);
                        assert_eq!(batch[1].len(), 0);
                        for index in 0..batch[0].len() {
                            batch[0]
                                .read(index, ordinal, phred64)
                                .unwrap()
                                .write(&mut actual, None)
                                .unwrap();
                            ordinal += 1;
                        }
                    }
                    assert_eq!(
                        actual, expected,
                        "capacity={capacity} budget={budget} phred64={phred64}"
                    );
                }
            }
        }
    }

    #[test]
    fn interleaved_blocks_count_limits_in_pairs_and_reject_dangling_mates() {
        let data = b"@a/1\nAC\n+\nII\n@a/2\nG\n+\nI\n@b/1\nT\n+\nH\n";
        for capacity in 1..=data.len() + 1 {
            let mut source = reader(data, capacity, false);
            let first = read_interleaved_batch(&mut source, 1024, 1)
                .unwrap()
                .unwrap();
            assert_eq!(first[0].len(), 2);
            assert_eq!(first[0].name(0), b"@a/1");
            assert_eq!(first[0].name(1), b"@a/2");
            assert!(read_interleaved_batch(&mut source, 1024, 1).is_err());
        }
    }

    #[test]
    fn sample_base_limit_matches_owned_prefix_without_consuming_extra_records() {
        let data = b"@a\r\nAC\r\n+\r\nII\r\n@b\nACGTA\n+\nIIIII\n@c\nT\n+\nI";
        for capacity in 1..=data.len() + 1 {
            for base_limit in 1..=10 {
                let mut expected = reader(data, capacity, false);
                let mut names = Vec::new();
                let mut bases = 0;
                while bases < base_limit {
                    let Some(r) = expected.next().unwrap() else {
                        break;
                    };
                    bases += r.seq.len();
                    names.push(r.name.into_owned());
                }
                let mut source = reader(data, capacity, false);
                let batch = read_sample_batch(&mut source, None, 1024, 99, base_limit)
                    .unwrap()
                    .unwrap();
                assert_eq!(
                    (0..batch[0].len())
                        .map(|i| batch[0].name(i).to_vec())
                        .collect::<Vec<_>>(),
                    names
                );
                assert_eq!(
                    source.next().unwrap().map(|r| r.name),
                    expected.next().unwrap().map(|r| r.name)
                );
            }
        }
    }

    #[test]
    fn blocks_match_owned_reader_at_boundaries_and_asymmetric_lengths() {
        let first = b"@one note\r\nACGTN\r\n+one\r\n@ABCD\r\n@empty\n\n+\n\n@last\nAC\n+\nHH";
        let second = format!(
            "@one\n{}\n+\n{}\n@empty\r\nA\r\n+note\r\nH\r\n@last\nT\n+\nH",
            "G".repeat(4097),
            "I".repeat(4097)
        );
        for phred64 in [false, true] {
            let expected = [
                normalized(first, phred64),
                normalized(second.as_bytes(), phred64),
            ];
            for capacity in 1..=first.len() + 1 {
                for budget in [1, 31, 256, 20000] {
                    let mut r1 = reader(first, capacity, phred64);
                    let mut r2 = reader(second.as_bytes(), capacity * 2 + 1, phred64);
                    let mut actual = [Vec::new(), Vec::new()];
                    let mut offset = 0;
                    while let Some(batch) =
                        read_batch(&mut r1, Some(&mut r2), budget, usize::MAX).unwrap()
                    {
                        assert_eq!(batch[0].len(), batch[1].len());
                        for mate in 0..2 {
                            for i in 0..batch[mate].len() {
                                batch[mate]
                                    .read(i, offset + i, phred64)
                                    .unwrap()
                                    .write(&mut actual[mate], None)
                                    .unwrap();
                            }
                        }
                        offset += batch[0].len();
                    }
                    assert_eq!(
                        actual, expected,
                        "capacity={capacity}, budget={budget}, phred64={phred64}"
                    );
                }
            }
        }
    }

    #[test]
    fn blocks_stop_at_read_limit_and_reject_unbalanced_mates() {
        let data = b"@one\nAC\n+\nII\n@two\nT\n+\nI\n";
        for capacity in 1..=data.len() + 1 {
            let mut first = reader(data, capacity, false);
            let mut second = reader(&data[..13], capacity, false);
            let batch = read_batch(&mut first, Some(&mut second), 1024, 1)
                .unwrap()
                .unwrap();
            assert_eq!(batch[0].len(), 1);
            assert_eq!(batch[1].len(), 1);
            assert!(read_batch(&mut first, Some(&mut second), 1024, usize::MAX).is_err());
        }
    }

    #[test]
    fn blocks_reject_the_same_malformed_records_as_owned_reader() {
        for data in [
            &b"wrong\nA\n+\nI\n"[..],
            &b"@r\nA\nwrong\nI\n"[..],
            &b"@r\nAC\n+\nI\n"[..],
            &b"@r\nA\n+\n \n"[..],
            &b"@r\nA\n+\n\xff\n"[..],
            &b"@r\nA\n+\n"[..],
            &b"@r\n"[..],
        ] {
            let expected = reader(data, 7, false).next().unwrap_err().to_string();
            for capacity in 1..=data.len() + 1 {
                let mut source = reader(data, capacity, false);
                let result = read_batch(&mut source, None, 1024, usize::MAX)
                    .and_then(|batch| batch.unwrap()[0].read(0, 0, false).map(|_| ()));
                assert_eq!(
                    result.unwrap_err().to_string(),
                    expected,
                    "capacity={capacity}, data={data:?}"
                );
            }
        }
    }
}
