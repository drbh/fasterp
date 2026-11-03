//! Native command-line interface. Compatibility target: OpenGene/fastp v1.4.0.
use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    fasterp::compat::run(fasterp::Args::parse())
}
