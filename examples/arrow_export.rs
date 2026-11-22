//! fasterp's output as Arrow, without copying.
//!
//! Each processed batch becomes an arrow-rs `RecordBatch` whose buffers are
//! fasterp's own; when Arrow drops them (here, or in Python after a C Data
//! Interface export) they go back to the run's pool.
//!
//! ```text
//! cargo run --release --features arrow --example arrow_export
//! ```
use fasterp::readkit_batch::arrow::arrow_array::{Array, BinaryArray, StructArray};
use fasterp::{Options, Run};
use std::io::Cursor;

fn main() -> Result<(), fasterp::Error> {
    // A small single-end FASTQ: the same read 1000 times, with a poor tail.
    let mut fastq = Vec::new();
    for i in 0..1000 {
        fastq.extend_from_slice(format!("@read{i}\n").as_bytes());
        fastq.extend_from_slice(b"ACGTTGCAAGGCTTACCGATACGTTGCAAGGCTTACCGATNNNNN\n+\n");
        fastq.extend_from_slice(b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#####\n");
    }
    let options = Options::new().cut_tail(20);
    let mut run = Run::single(&options, Cursor::new(fastq))?;

    while let Some(batch) = run.next_batch()? {
        // Zero-copy: the RecordBatch's buffers are the batch's columns.
        let table = batch.into_records().into_arrow();
        println!("schema: {}", table.schema());
        let seg0 = table.column(0).as_any().downcast_ref::<StructArray>();
        let seqs = seg0.and_then(|s| s.column_by_name("seq"));
        let seqs = seqs.and_then(|c| c.as_any().downcast_ref::<BinaryArray>());
        if let Some(seqs) = seqs {
            let first = String::from_utf8_lossy(seqs.value(0));
            println!("{} rows; first trimmed read: {first}", seqs.len());
        }
    }
    println!(
        "pool after Arrow released the buffers: {:?}",
        run.pool().stats()
    );
    let report = run.finish()?;
    println!("{} reads passed", report.passed_reads());
    Ok(())
}
