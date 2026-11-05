//! Native command-line interface. Compatibility target: OpenGene/fastp v1.4.0.
use anyhow::Result;
use clap::Parser;

// The system allocator returns freed blocks of a few hundred KiB to the
// kernel, and the gzip writer frees one per block, so their pages were
// faulted in again. mimalloc keeps them and is also faster for small
// allocations: 6% faster on paired-end NovaSeq data.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> Result<()> {
    fasterp::compat::run(fasterp::Args::parse())
}
