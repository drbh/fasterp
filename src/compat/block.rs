//! Bounded FASTQ batches. Fields borrow their batch until worker output is packed.
//! The slow boundary path accepts the same strict four-line format as Reader,
//! including skipping blank lines between records.
mod newlines;

use super::{Read, Reader};
use anyhow::{ensure, Result};
use std::borrow::Cow;
use std::sync::Mutex;

#[derive(Default)]
pub(super) struct Records {
    data: Vec<u8>,
    ends: Vec<[usize; 4]>,
}

/// Buffers of dropped batches, reused by later batches so their pages stay
/// mapped. Freed large blocks go back to the kernel and every new batch
/// faulted them in again. A few batches are in flight at once, so a few
/// spares are enough.
static SPARE: Mutex<Vec<Records>> = Mutex::new(Vec::new());
const SPARES: usize = 8;

impl Drop for Records {
    fn drop(&mut self) {
        if self.data.capacity() == 0 {
            return;
        }
        let mut spare = SPARE.lock().unwrap_or_else(|e| e.into_inner());
        if spare.len() < SPARES {
            let mut data = std::mem::take(&mut self.data);
            let mut ends = std::mem::take(&mut self.ends);
            data.clear();
            ends.clear();
            spare.push(Records { data, ends });
        }
    }
}

impl Records {
    fn reused() -> Self {
        let spare = SPARE.lock().unwrap_or_else(|e| e.into_inner()).pop();
        spare.unwrap_or_default()
    }

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
        let record_start = start;
        let fields = self.ends[index].map(|end| {
            let line = &self.data[start..end];
            start = end;
            let line = line.strip_suffix(b"\n").unwrap_or(line);
            Cow::Borrowed(line.strip_suffix(b"\r").unwrap_or(line))
        });
        let [name, seq, plus, qual] = fields;
        let record = &self.data[record_start..start];
        let verbatim = (record.ends_with(b"\n")
            && record.len() == name.len() + seq.len() + plus.len() + qual.len() + 4)
            .then_some(record);
        let mut read = Read {
            verbatim,
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
        reader.skip_blank_lines()?;
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

pub(super) fn read_sample_batch<const INTERLEAVED: bool>(
    first: &mut Reader,
    second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
    base_limit: usize,
) -> Result<Option<[Records; 2]>> {
    read_batch_inner::<true, INTERLEAVED>(first, second, budget, limit, base_limit)
}

fn read_batch_inner<const SAMPLE: bool, const INTERLEAVED: bool>(
    first: &mut Reader,
    mut second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
    base_limit: usize,
) -> Result<Option<[Records; 2]>> {
    let mut blocks: [Records; 2] = std::array::from_fn(|_| Records::reused());
    // A new buffer grown by doubling copies itself several times before it
    // holds a batch, so give it the byte budget it is expected to fill.
    // Reused buffers already have about that capacity.
    let share = if second.is_some() { budget / 2 } else { budget };
    for block in &mut blocks[..1 + usize::from(second.is_some())] {
        if block.data.capacity() == 0 {
            block.data.reserve_exact(share.min(16 << 20));
        }
    }
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
        let mut newlines1 = newlines::offsets(b1);
        let mut newlines2 = b2.map(newlines::offsets);
        let base = [blocks[0].data.len(), blocks[1].data.len()];
        let used = base[0] + base[1];
        let mut count = 0;
        let mut consumed = [0; 2];
        let mut full = false;
        while count < remaining {
            // Leave any line that is not a header, such as a blank line
            // between records, to the slow path.
            if b1.get(consumed[0]) != Some(&b'@')
                || b2.is_some_and(|b| b.get(consumed[1]) != Some(&b'@'))
            {
                break;
            }
            let Some(ends1) = next_ends(&mut newlines1) else {
                break;
            };
            let ends2 = if INTERLEAVED {
                if b1.get(ends1[3]) != Some(&b'@') {
                    break;
                }
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
            let has_first = read_slow::<INTERLEAVED>(&mut blocks, first, second.as_deref_mut())
                .map_err(|error| {
                    let records = [first.record, second.as_ref().map_or(0, |r| r.record)];
                    earliest_error(&blocks, records, first.phred64).unwrap_or(error)
                })?;
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

/// Read one record from each input (a pair from interleaved input).
fn read_slow<const INTERLEAVED: bool>(
    blocks: &mut [Records; 2],
    first: &mut Reader,
    second: Option<&mut Reader>,
) -> Result<bool> {
    let has_first = blocks[0].read_one(first)?;
    if INTERLEAVED && has_first {
        ensure!(
            blocks[0].read_one(first)?,
            "Odd read count in interleaved FASTQ"
        );
    }
    if let Some(r) = second {
        let has_second = blocks[1].read_one(r)?;
        ensure!(
            has_first == has_second,
            "Paired input files have different read counts"
        );
    }
    Ok(has_first)
}

/// Report records in order, as Reader does. A blank line inside a record
/// shifts every later line, which otherwise surfaces as truncation at EOF
/// instead of an error at the record that holds it.
fn earliest_error(
    blocks: &[Records; 2],
    records: [usize; 2],
    phred64: bool,
) -> Option<anyhow::Error> {
    let first = records[0] - blocks[0].len();
    let second = records[1] - blocks[1].len();
    (0..blocks[0].len().max(blocks[1].len())).find_map(|i| {
        let check = |mate: usize, start: usize| {
            (i < blocks[mate].len())
                .then(|| blocks[mate].read(i, start + i, phred64).err())
                .flatten()
        };
        check(0, first).or_else(|| check(1, second))
    })
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
    fn verbatim_records_require_matching_fields_and_all_four_newlines() {
        for (data, eligible) in [
            (&b"@a\nACGT\n+a\nIIII\n"[..], true),
            (b"@a\r\nACGT\n+a\nIIII\n", false),
            (b"@a\r\nACGT\n+a\nIIII", false),
        ] {
            let mut source = reader(data, 64, false);
            let records = read_batch(&mut source, None, 1024, 1).unwrap().unwrap();
            let r = records[0].read(0, 0, false).unwrap();
            assert_eq!(r.verbatim.is_some(), eligible);
            for (mut r, tag, expected) in [
                (r.clone(), None, &b"@a\nACG\n+a\nIII\n"[..]),
                (r.clone(), Some("tag"), b"@a tag\nACG\n+a\nIII\n"),
            ] {
                r.resize(3);
                let mut out = Vec::new();
                r.write(&mut out, tag).unwrap();
                assert_eq!(out, expected);
            }
            let mut renamed = r.clone();
            renamed.name_mut().push(b'b');
            let mut out = Vec::new();
            renamed.write(&mut out, None).unwrap();
            assert_eq!(out, b"@ab\nACGT\n+a\nIIII\n");
            let mut out = Vec::new();
            r.write(&mut out, Some("tag")).unwrap();
            assert_eq!(out, b"@a tag\nACGT\n+a\nIIII\n");
            assert!(records[0].read(0, 0, true).unwrap().verbatim.is_none());
        }
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
    fn interleaved_sample_budget_counts_first_mates_and_keeps_whole_pairs() {
        let data = b"@a/1\nAC\n+\nII\n@a/2\nAAAAA\n+\nIIIII\n@b/1\nGT\n+\nII\n@b/2\nC\n+\nI\n@c/1\nA\n+\nI\n@c/2\nC\n+\nI\n";
        for capacity in 1..=data.len() {
            let mut source = reader(data, capacity, false);
            let block = read_sample_batch::<true>(&mut source, None, 1024, 99, 4)
                .unwrap()
                .unwrap();
            assert_eq!(block[0].len(), 4);
            assert_eq!(block[0].name(3), b"@b/2");
            assert_eq!(source.next().unwrap().unwrap().name.as_ref(), b"@c/1");
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
                let batch = read_sample_batch::<false>(&mut source, None, 1024, 99, base_limit)
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

    fn batch_names(batch: &[Records; 2]) -> [Vec<Vec<u8>>; 2] {
        batch
            .each_ref()
            .map(|b| (0..b.len()).map(|i| b.name(i).to_vec()).collect())
    }

    #[test]
    fn blocks_skip_blank_lines_between_records_in_every_mode() {
        let blank = format!("\n \t\r\n{}\n", " ".repeat(40));
        let records = [
            "@a\nAC\n+\nII\n",
            "@b\nT\n+\nI\n",
            "@c\nG\n+\nH\n",
            "@d\nA\n+\nI",
        ];
        let clean = records.concat();
        let messy = format!("{blank}{}", records.join(&blank)) + "\n\n  ";
        let names = |n: &[&str]| n.iter().map(|n| n.as_bytes().to_vec()).collect::<Vec<_>>();
        let all = names(&["@a", "@b", "@c", "@d"]);
        for capacity in 1..=messy.len() + 1 {
            for budget in [1, 31, 20000] {
                let mut collected: [Vec<Vec<u8>>; 2] = Default::default();
                let mut first = reader(messy.as_bytes(), capacity, false);
                let mut second = reader(clean.as_bytes(), capacity, false);
                while let Some(batch) =
                    read_batch(&mut first, Some(&mut second), budget, usize::MAX).unwrap()
                {
                    let [a, b] = batch_names(&batch);
                    collected[0].extend(a);
                    collected[1].extend(b);
                }
                assert_eq!(collected, [all.clone(), all.clone()], "capacity={capacity}");
                assert_eq!((first.record, second.record), (4, 4));

                let mut interleaved = Vec::new();
                let mut source = reader(messy.as_bytes(), capacity, false);
                while let Some(batch) =
                    read_interleaved_batch(&mut source, budget, usize::MAX).unwrap()
                {
                    interleaved.extend(batch_names(&batch)[0].clone());
                }
                assert_eq!(interleaved, all, "capacity={capacity}");

                let mut source = reader(messy.as_bytes(), capacity, false);
                let batch = read_sample_batch::<false>(&mut source, None, 20000, 99, 3)
                    .unwrap()
                    .unwrap();
                assert_eq!(batch_names(&batch)[0], names(&["@a", "@b"]));
                assert_eq!(source.next().unwrap().unwrap().name.as_ref(), b"@c");
            }
        }
        for data in [&b""[..], b"\n", b"  \n\t\r\n  "] {
            assert!(read_batch(&mut reader(data, 1, false), None, 1024, 9)
                .unwrap()
                .is_none());
        }
    }

    #[test]
    fn blocks_report_blank_lines_inside_records_like_owned_reader() {
        for data in [
            &b"\n@a\nA\n+\nI\n\n@b\n\nA\n+\nI\n@c\nA\n+\nI\n"[..],
            b"@a\nA\n+\nI\n@b\nA\n\n+\nI\n",
            b"@a\nA\n+\nI\n@b\nA\n+\n\nI\n@c\nA\n+\nI\n",
            b"@a\nA\n+\nI\n\n  @b\nA\n+\nI\n",
        ] {
            let mut owned = reader(data, 7, false);
            let expected = loop {
                match owned.next() {
                    Ok(Some(_)) => {}
                    Ok(None) => panic!("{data:?} parsed"),
                    Err(error) => break error.to_string(),
                }
            };
            for capacity in 1..=data.len() + 1 {
                let mut source = reader(data, capacity, false);
                let mut ordinal = 0;
                let actual = loop {
                    match read_batch(&mut source, None, 1024, usize::MAX) {
                        Ok(Some(batch)) => {
                            for i in 0..batch[0].len() {
                                if let Err(error) = batch[0].read(i, ordinal, false) {
                                    panic!("{data:?} capacity={capacity}: late {error}");
                                }
                                ordinal += 1;
                            }
                        }
                        Ok(None) => panic!("{data:?} parsed"),
                        Err(error) => break error.to_string(),
                    }
                };
                assert_eq!(actual, expected, "capacity={capacity}, data={data:?}");
            }
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
