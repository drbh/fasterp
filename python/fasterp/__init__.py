"""
fasterp - High-performance FASTQ preprocessing

A Python library for ultra-fast FASTQ file processing, built on Rust for maximum performance.
fasterp is 2-10x faster than fastp while maintaining the same functionality.

Examples:
    File-based processing:
    >>> import fasterp
    >>> result = fasterp.process(
    ...     input="reads.fq.gz",
    ...     output="filtered.fq.gz",
    ...     min_length=15,
    ...     average_qual=20
    ... )
    >>> print(f"Passed: {result.passed_reads}, Failed: {result.failed_reads}")

    Record-based processing (in-memory):
    >>> result = fasterp.process_records(
    ...     input="reads.fq.gz",
    ...     min_length=15,
    ...     n_base_limit=5
    ... )
    >>> for record in result.records:
    ...     print(f">{record.header}")
    ...     print(f"Sequence: {record.sequence}")
    ...     print(f"Quality: {record.quality}")

    Processing from bytes:
    >>> with open("reads.fq.gz", "rb") as f:
    ...     data = f.read()
    >>> result = fasterp.process_records(
    ...     input_bytes=data,
    ...     min_length=20
    ... )
    >>> sequences = [r.sequence for r in result.records]
"""

from .fasterp import (
    ProcessResult,
    FastqRecord,
    ProcessRecordsResult,
    process,
    process_records,
)

__all__ = [
    "ProcessResult",
    "FastqRecord",
    "ProcessRecordsResult",
    "process",
    "process_records",
]
__version__ = "0.2.1"
