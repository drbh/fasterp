#!/usr/bin/env python3
"""
Basic example of using fasterp's record-based API.

This demonstrates how to process FASTQ files and access individual
records in memory without writing output files.
"""

import fasterp
import os

# Path to test data
# input_file = "../SRR27128759_1_part.fastq"
# input_file = "./SRR27128759_1_part.fastq"
input_file = "SRR5808766/SRR5808766_1.fastq.gz"

if not os.path.exists(input_file):
    print(f"Error: Test file '{input_file}' not found")
    print("Please run this script from the python_example directory")
    exit(1)

print("=" * 70)
print("fasterp Basic Record Processing Example")
print("=" * 70)

# Process FASTQ file and get records in memory
print(f"\nProcessing: {input_file}")
print("Filters: min_length=30, n_base_limit=5")

result = fasterp.process_records(
    input=input_file,
    min_length=30,
    n_base_limit=5,
)

print(f"\n{'Results':^70}")
print("-" * 70)
print(f"Total reads processed: {result.total_reads}")
print(f"Reads passed:          {result.passed_reads}")
print(f"Reads failed:          {result.failed_reads}")
print(f"Records in memory:     {len(result.records)}")

# Display quality metrics
print(f"\n{'Quality Metrics':^70}")
print("-" * 70)
print(f"Total bases before:    {result.total_bases_before:,}")
print(f"Total bases after:     {result.total_bases_after:,}")
print(f"Q20 rate:              {result.q20_rate:.2%}")
print(f"Q30 rate:              {result.q30_rate:.2%}")
print(f"GC content:            {result.gc_content:.2%}")
print(f"Duplication rate:      {result.duplication_rate:.2%}")

# Show first 3 records
print(f"\n{'First 3 Records':^70}")
print("-" * 70)
for i, record in enumerate(result.records[:3], 1):
    print(f"\nRecord {i}:")
    print(f"  Header:   {record.header}")
    print(f"  Sequence: {record.sequence[:60]}...")
    print(f"  Quality:  {record.quality[:60]}...")
    print(f"  Length:   {record.length} bp")
    print(f"  Status:   {'PASS' if record.passed else 'FAIL'}")

# Calculate sequence statistics
print(f"\n{'Sequence Statistics':^70}")
print("-" * 70)
lengths = [record.length for record in result.records]
if lengths:
    print(f"Min length:   {min(lengths)} bp")
    print(f"Max length:   {max(lengths)} bp")
    print(f"Mean length:  {sum(lengths) / len(lengths):.1f} bp")

# Calculate GC content per record
print(f"\n{'GC Content Distribution (first 10 reads)':^70}")
print("-" * 70)
for i, record in enumerate(result.records[:10], 1):
    gc_count = record.sequence.count('G') + record.sequence.count('C')
    gc_pct = (gc_count / record.length * 100) if record.length > 0 else 0
    print(f"  Read {i:2d}: {gc_pct:5.1f}% GC")

print("\n" + "=" * 70)
print("Done! All records are available in memory as Python objects.")
print("=" * 70)
