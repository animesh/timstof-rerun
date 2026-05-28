// src/rerun_view.rs
//
// Stream timsTOF 5D data directly to Rerun.
//
// 5D mapping:
//   RT         -> Rerun timeline  "rt_seconds" (duration)
//   m/z        -> X axis of frame heatmap image (px column)
//   1/K0       -> Y axis of frame heatmap image (px row, flipped)
//   Intensity  -> pixel colour  (plasma LUT, log-scaled)
//   MS type    -> separate BPC/TIC time series entities
//
// Data flow (rustdf, no SDK):
//   TimsDataset::new(path)          open .d directory
//   dataset.get_frame(i)            -> TimsFrame (all IMS scans)
//   frame.mz_values                 Vec<f64>  (calibrated from SQLite polynomial)
//   frame.scan_numbers              Vec<u32>  (IMS scan index)
//   frame.intensities               Vec<u32>
//   frame.to_one_over_k0_vec()      Vec<f64>  (1/K0 from scan via SQLite poly)
//   frame.retention_time            f64 (seconds)
//   frame.ms_type                   MsType::{MS1, MS2}

use anyhow::Result;
use log::info;
use rayon::prelude::*;
use rustdf::data::dataset::TimsDataset;
use rustdf::data::meta::MsType;
use std::path::Path;

use crate::heatmap::{make_heatmap_rgb, MZ_BINS, MZ_MAX, MZ_MIN, MOB_BINS, MOB_MAX, MOB_MIN};
use crate::chrom::read_chrom_trace;

// Rerun Rust SDK (0.22.x)
use rerun::{
    RecordingStreamBuilder,
    Points2D, Image, Scalar, SeriesLine, TextDocument,
};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const TIMELINE: &str = "rt_seconds";

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

pub fn run(
    d_path: &Path,
    max_frames: usize,
    top_n: usize,
    save: Option<&Path>,
    connect: Option<&str>,
) -> Result<()> {
    // Build Rerun recording stream
    let mut builder = RecordingStreamBuilder::new("timstof-rerun");

    let rec = match (save, connect) {
        (Some(p), _) => {
            info!("Saving to {}", p.display());
            builder.save(p)?
        }
        (_, Some(addr)) => {
            info!("Connecting to {}", addr);
            builder.connect_grpc_opts(addr, Default::default())?
        }
        _ => {
            info!("Spawning Rerun viewer");
            builder.spawn()?
        }
    };

    // Open dataset via rustdf (reads SQLite + mmap's tdf_bin)
    let dataset = TimsDataset::new(d_path.to_str().unwrap())?;
    let n_frames = dataset.get_frame_count();
    info!("{} total frames in {}", n_frames, d_path.display());

    // Log chromatograms from chromatography-data.sqlite
    log_chromatograms(&rec, d_path)?;

    // Log layout blueprint
    log_blueprint(&rec)?;

    // Collect MS1 frame indices
    let ms1_indices: Vec<usize> = (0..n_frames)
        .filter(|&i| {
            dataset
                .get_frame(i)
                .map(|f| matches!(f.ms_type, MsType::MS1))
                .unwrap_or(false)
        })
        .collect();

    let take = if max_frames == 0 {
        ms1_indices.len()
    } else {
        ms1_indices.len().min(max_frames)
    };

    info!("Logging {} / {} MS1 frames (top_n={})", take, ms1_indices.len(), top_n);

    for (n, &frame_idx) in ms1_indices[..take].iter().enumerate() {
        let frame = match dataset.get_frame(frame_idx) {
            Ok(f)  => f,
            Err(e) => { log::warn!("skip frame {}: {}", frame_idx, e); continue; }
        };

        let rt = frame.retention_time;
        rec.set_time_seconds(TIMELINE, rt);

        // Get calibrated 1/K0 for all peaks in this frame
        let k0_vals: Vec<f64> = frame.to_one_over_k0_vec();

        // Build arrays: mz, k0, intensity
        let mz_arr:    &[f64] = &frame.mz_values;
        let inten_arr: &[u32] = &frame.intensities;

        // Apply top-N filter by intensity
        let (mz_f, mob_f, int_f) = if mz_arr.len() > top_n {
            let mut indexed: Vec<usize> = (0..mz_arr.len()).collect();
            indexed.sort_unstable_by(|&a, &b| inten_arr[b].cmp(&inten_arr[a]));
            indexed.truncate(top_n);
            let mz_f:  Vec<f32> = indexed.iter().map(|&i| mz_arr[i]  as f32).collect();
            let mob_f: Vec<f32> = indexed.iter().map(|&i| k0_vals[i] as f32).collect();
            let int_f: Vec<f32> = indexed.iter().map(|&i| inten_arr[i] as f32).collect();
            (mz_f, mob_f, int_f)
        } else {
            let mz_f:  Vec<f32> = mz_arr.iter().map(|&v| v as f32).collect();
            let mob_f: Vec<f32> = k0_vals.iter().map(|&v| v as f32).collect();
            let int_f: Vec<f32> = inten_arr.iter().map(|&v| v as f32).collect();
            (mz_f, mob_f, int_f)
        };

        // Frame heatmap image
        let rgb = make_heatmap_rgb(&mz_f, &mob_f, &int_f);
        rec.log(
            "frame_heatmap",
            &Image::from_rgb24(rgb, [MZ_BINS as u32, MOB_BINS as u32]),
        )?;

        // Mobilogram: (x=intensity, y=1/K0) - filter zeros
        let mob_pts: Vec<[f32; 2]> = mob_f
            .iter()
            .zip(int_f.iter())
            .filter(|(_, &i)| i > 0.0)
            .map(|(&mob, &i)| [i.ln_1p(), mob])
            .collect();
        if !mob_pts.is_empty() {
            rec.log("mobilogram", &Points2D::new(&mob_pts))?;
        }

        // Summed spectrum: (x=m/z, y=intensity)
        let spec_pts: Vec<[f32; 2]> = mz_f
            .iter()
            .zip(int_f.iter())
            .filter(|(_, &i)| i > 0.0)
            .map(|(&mz, &i)| [mz, i.ln_1p()])
            .collect();
        if !spec_pts.is_empty() {
            rec.log("spectrum", &Points2D::new(&spec_pts))?;
        }

        if (n + 1) % 200 == 0 {
            info!("  {} frames logged...", n + 1);
        }
    }

    info!("Done. {} MS1 frames logged.", take);

    if save.is_none() && connect.is_none() {
        info!("Viewer running - press Ctrl+C to exit");
        loop {
            std::thread::sleep(std::time::Duration::from_secs(1));
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Chromatogram logging
// ---------------------------------------------------------------------------

fn log_chromatograms(rec: &rerun::RecordingStream, d_path: &Path) -> Result<()> {
    let chrom_db = d_path.join("chromatography-data.sqlite");
    if !chrom_db.exists() {
        log::warn!("chromatography-data.sqlite not found - skipping chromatograms");
        return Ok(());
    }

    // Style hints (static)
    rec.log_static("chromatogram/BPC",
        &SeriesLine::new().with_color([0, 180, 255]).with_width(1.0))?;
    rec.log_static("chromatogram/TIC",
        &SeriesLine::new().with_color([255, 140, 0]).with_width(1.0))?;

    // Trace 18 = BPC +-MS, Trace 20 = TIC +-MS
    for (trace_id, entity) in &[(18u32, "chromatogram/BPC"), (20u32, "chromatogram/TIC")] {
        match read_chrom_trace(&chrom_db, *trace_id) {
            Ok((times, intensities)) => {
                let rt0 = times.first().copied().unwrap_or(0.0);
                for (t, v) in times.iter().zip(intensities.iter()) {
                    rec.set_time_seconds(TIMELINE, t - rt0);
                    rec.log(*entity, &Scalar(*v as f64))?;
                }
                info!("  {} : {} points", entity, times.len());
            }
            Err(e) => log::warn!("  trace {}: {}", trace_id, e),
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Blueprint as static TextDocument explaining the layout
// (The Rerun Rust SDK 0.22 does not yet have a typed blueprint API;
//  the viewer auto-lays out panels. We log a help text instead.)
// ---------------------------------------------------------------------------

fn log_blueprint(rec: &rerun::RecordingStream) -> Result<()> {
    rec.log_static(
        "info/layout",
        &TextDocument::new(
            "## timsTOF 5D Viewer  (timstof-rerun)\n\n\
            **Top row**: BPC (blue) and TIC (orange) chromatograms.\n\n\
            **Frame Heatmap**: one MS1 frame at the current RT cursor.\n\
            - X = m/z  (200-1700 Da, left->right)\n\
            - Y = 1/K0 (0.60-1.80 Vs/cm², bottom->top)\n\
            - Colour = log(intensity): black=0, purple=low, red=mid, yellow=high\n\
            → **Drag the timeline scrubber** to the chromatographic peak region.\n\n\
            **Mobilogram**: ln(I) vs 1/K0 for the selected frame.\n\
            **Spectrum**: m/z vs ln(I) for the selected frame.\n\n\
            Data read by rustdf (rustims, J. Proteome Res. 2025). No Bruker SDK.",
        )
        .with_media_type(rerun::MediaType::markdown()),
    )?;
    Ok(())
}
