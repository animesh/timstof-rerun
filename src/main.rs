// timstof-rerun/src/main.rs
//
// Pure-Rust timsTOF viewer.  No Bruker SDK.  No Python.
//
// Pipeline:
//   .d/analysis.tdf       (SQLite)  -> frame metadata, precursor info
//   .d/analysis.tdf_bin   (binary)  -> raw TOF indices + intensities
//   .d/chromatography-data.sqlite   -> BPC, TIC, nanoElute traces
//
// Reading via rustdf (from the rustims monorepo).
// Visualisation via the Rerun Rust SDK.
//
// Subcommands:
//   rerun   - stream 5D data (RT, m/z, 1/K0, intensity, MS-type) to Rerun viewer
//   sql     - dump all analysis.tdf tables as TSV
//   chrom   - dump all chromatography traces as TSV
//   mgf     - write MS/MS spectra in MGF format

use anyhow::Result;
use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

mod chrom;
mod heatmap;
mod mgf;
mod rerun_view;
mod sql_dump;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser)]
#[command(
    name    = "timstof-rerun",
    version = env!("CARGO_PKG_VERSION"),
    about   = "timsTOF 5D raw data viewer - no Bruker SDK required",
    long_about = "\
Reads Bruker timsTOF .d directories using rustdf (rustims framework, \
J. Proteome Res. 2025) and streams data directly to the Rerun viewer.\n\
No Bruker SDK, no Python, no intermediate files required."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Stream 5D ion map + chromatograms to Rerun viewer (live or .rrd file)
    Rerun {
        /// Path to the .d directory
        #[arg(value_name = "FILE.d")]
        path: PathBuf,

        /// Max MS1 frames to visualise (0 = all, default 2000)
        #[arg(long, default_value = "2000")]
        max_frames: usize,

        /// Top N peaks per frame for heatmap (default 2000)
        #[arg(long, default_value = "2000")]
        top_n: usize,

        /// Save to .rrd file instead of launching live viewer
        #[arg(long, value_name = "OUT.rrd")]
        save: Option<PathBuf>,

        /// Connect to running Rerun viewer at this address (default: spawn new)
        #[arg(long, value_name = "HOST:PORT")]
        connect: Option<String>,
    },

    /// Dump every table in analysis.tdf as TSV files
    Sql {
        /// Path to the .d directory
        #[arg(value_name = "FILE.d")]
        path: PathBuf,

        /// Output directory (default: <run>_tdf/)
        #[arg(long)]
        out_dir: Option<PathBuf>,
    },

    /// Dump all chromatography traces (BPC, TIC, nanoElute) as TSV files
    Chrom {
        /// Path to the .d directory
        #[arg(value_name = "FILE.d")]
        path: PathBuf,

        /// Output directory (default: <run>_chrom/)
        #[arg(long)]
        out_dir: Option<PathBuf>,
    },

    /// Write MS/MS spectra in MGF format
    Mgf {
        /// Path to the .d directory
        #[arg(value_name = "FILE.d")]
        path: PathBuf,

        /// Output MGF file (default: <run>.mgf)
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or("info"),
    )
    .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Rerun { path, max_frames, top_n, save, connect } => {
            info!("Streaming to Rerun: {}", path.display());
            rerun_view::run(&path, max_frames, top_n, save.as_deref(), connect.as_deref())?;
        }
        Commands::Sql { path, out_dir } => {
            info!("Dumping SQL tables: {}", path.display());
            sql_dump::run(&path, out_dir.as_deref())?;
        }
        Commands::Chrom { path, out_dir } => {
            info!("Dumping chromatography traces: {}", path.display());
            chrom::run(&path, out_dir.as_deref())?;
        }
        Commands::Mgf { path, out } => {
            info!("Writing MGF: {}", path.display());
            mgf::run(&path, out.as_deref())?;
        }
    }

    Ok(())
}
