// src/rerun_view.rs
// API targets rerun 0.27.x
// Key changes from 0.22:
//   connect_tcp_opts -> connect_grpc (0.27 uses gRPC by default)
//   Image::from_rgb24 still works in 0.27
//   Scalar::new still works in 0.27

use anyhow::Result;
use log::info;
use rustdf::data::dataset::TimsDataset;
use rustdf::data::handle::TimsData;
use std::path::Path;

use crate::heatmap::{make_heatmap_rgb, MZ_BINS, MOB_BINS};
use crate::chrom::read_chrom_trace_by_desc;

use rerun::{
    RecordingStreamBuilder,
    Points2D, Image, Scalars, SeriesLines, TextDocument,
};

const TIMELINE: &str = "rt_seconds";
const MS1_NUMERIC: i32 = 0;

pub fn run(
    d_path: &Path,
    max_frames: usize,
    top_n: usize,
    save: Option<&Path>,
    connect: Option<&str>,
) -> Result<()> {
    let builder = RecordingStreamBuilder::new("timstof-rerun");

    let rec = match (save, connect) {
        (Some(p), _) => {
            info!("Saving to {}", p.display());
            builder.save(p)?
        }
        (_, Some(addr)) => {
            info!("Connecting to {}", addr);
            // 0.27: connect_grpc takes a URL string directly
            builder.connect_grpc_opts(addr)?
        }
        _ => {
            info!("Spawning Rerun viewer");
            builder.spawn()?
        }
    };

    let d_str = d_path.to_str().unwrap();
    let dataset = TimsDataset::new("", d_str, false, false);
    let n_frames = dataset.loader.get_frame_count() as usize;
    info!("{} total frames in {}", n_frames, d_path.display());

    log_chromatograms(&rec, d_path)?;
    log_layout_info(&rec)?;

    // Collect MS1 frame IDs
    let ms1_ids: Vec<u32> = (1..=n_frames as u32)
        .filter(|&fid| {
            dataset.loader.get_frame(fid).ms_type.ms_type_numeric() == MS1_NUMERIC
        })
        .collect();

    let take = if max_frames == 0 { ms1_ids.len() } else { ms1_ids.len().min(max_frames) };
    info!("Logging {} / {} MS1 frames  top_n={}", take, ms1_ids.len(), top_n);

    // Parallel: build images first, then send sequentially
    // rayon prep + sequential gRPC is faster than per-frame serial
    use rayon::prelude::*;

    let frame_data: Vec<(f64, Vec<u8>, Vec<[f32;2]>, Vec<[f32;2]>)> =
        ms1_ids[..take].par_iter().map(|&fid| {
            let frame  = dataset.loader.get_frame(fid);
            let rt     = frame.ims_frame.retention_time;
            let mz     = &*frame.ims_frame.mz;
            let mob    = &*frame.ims_frame.mobility;
            let inten  = &*frame.ims_frame.intensity;

            let (mz_f, mob_f, int_f) = top_n_filter(mz, mob, inten, top_n);

            let rgb = make_heatmap_rgb(&mz_f, &mob_f, &int_f);

            let mob_pts: Vec<[f32;2]> = int_f.iter().zip(mob_f.iter())
                .filter(|(&i, _)| i > 0.0)
                .map(|(&i, &k)| [i.ln_1p(), k])
                .collect();

            let spec_pts: Vec<[f32;2]> = mz_f.iter().zip(int_f.iter())
                .filter(|(_, &i)| i > 0.0)
                .map(|(&m, &i)| [m, i.ln_1p()])
                .collect();

            (rt, rgb, mob_pts, spec_pts)
        }).collect();

    info!("Frame prep done, sending to viewer...");

    for (n, (rt, rgb, mob_pts, spec_pts)) in frame_data.into_iter().enumerate() {
        rec.set_duration_secs(TIMELINE, rt);

        rec.log("frame_heatmap",
            &Image::from_rgb24(rgb, [MZ_BINS as u32, MOB_BINS as u32]))?;

        if !mob_pts.is_empty() {
            rec.log("mobilogram", &Points2D::new(&mob_pts))?;
        }
        if !spec_pts.is_empty() {
            rec.log("spectrum", &Points2D::new(&spec_pts))?;
        }

        if (n + 1) % 200 == 0 {
            info!("  {} frames sent...", n + 1);
        }
    }

    info!("Done. {} MS1 frames logged.", take);

    if save.is_none() && connect.is_none() {
        info!("Viewer running - Ctrl+C to exit");
        loop { std::thread::sleep(std::time::Duration::from_secs(1)); }
    }
    Ok(())
}

fn top_n_filter(
    mz: &[f64], mob: &[f64], inten: &[f64], top_n: usize,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let n = mz.len().min(mob.len()).min(inten.len());
    if n == 0 { return (vec![], vec![], vec![]); }
    if n <= top_n {
        return (
            mz[..n].iter().map(|&v| v as f32).collect(),
            mob[..n].iter().map(|&v| v as f32).collect(),
            inten[..n].iter().map(|&v| v as f32).collect(),
        );
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_unstable_by(|&a, &b| inten[b].partial_cmp(&inten[a]).unwrap());
    idx.truncate(top_n);
    (
        idx.iter().map(|&i| mz[i]    as f32).collect(),
        idx.iter().map(|&i| mob[i]   as f32).collect(),
        idx.iter().map(|&i| inten[i] as f32).collect(),
    )
}

fn log_chromatograms(rec: &rerun::RecordingStream, d_path: &Path) -> Result<()> {
    let chrom_db = d_path.join("chromatography-data.sqlite");
    if !chrom_db.exists() {
        log::warn!("chromatography-data.sqlite not found - skipping");
        return Ok(());
    }

    rec.log_static("chromatogram/BPC",
        &SeriesLines::new().with_colors([(0_u8, 180_u8, 255_u8)]).with_widths([1.0]))?;
    rec.log_static("chromatogram/TIC",
        &SeriesLines::new().with_colors([(255_u8, 140_u8, 0_u8)]).with_widths([1.0]))?;

    // Dynamic lookup: find BPC and TIC traces by description substring
    // (trace IDs vary between acquisition methods and HyStar versions)
    for (desc_fragment, entity) in &[("BPC", "chromatogram/BPC"), ("TIC", "chromatogram/TIC")] {
        match read_chrom_trace_by_desc(&chrom_db, desc_fragment) {
            Ok(Some((times, intens))) => {
                let rt0 = times.first().copied().unwrap_or(0.0);
                for (t, v) in times.iter().zip(intens.iter()) {
                    rec.set_duration_secs(TIMELINE, t - rt0);
                    rec.log(*entity, &Scalars::new([*v as f64]))?;
                }
                info!("  {} : {} points", entity, times.len());
            }
            Ok(None) => log::warn!("  no trace matching '{}' found", desc_fragment),
            Err(e)   => log::warn!("  {}: {}", desc_fragment, e),
        }
    }
    Ok(())
}

fn log_layout_info(rec: &rerun::RecordingStream) -> Result<()> {
    rec.log_static("info/layout",
        &TextDocument::new(
            "## timsTOF 5D Viewer (timstof-rerun)\n\n\
            **BPC/TIC** top row: full-run chromatograms.\n\
            **Frame Heatmap**: x=m/z  y=1/K0  colour=log(I)\n\
            -> Drag the timeline scrubber to the peak region (30-60 min).\n\
            **Mobilogram**: ln(I) vs 1/K0.  **Spectrum**: m/z vs ln(I).\n\
            No Bruker SDK - calibration via rustdf polynomial fit."
        ).with_media_type(rerun::MediaType::markdown()),
    )?;
    Ok(())
}
