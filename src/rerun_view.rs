// src/rerun_view.rs  (rerun 0.33)
//
// Architecture (mzPeak-inspired):
//   Per frame: log Points3D(mz, mobility, log_intensity) with intensity-mapped colour.
//   This gives a true rotatable 3D scatter.  The Rerun timeline = RT axis.
//   Dragging the scrubber to a BPC peak shows the corresponding ion-mobility scatter.
//   MS2 precursor m/z values are logged as Scalar markers on a separate entity.
//
// No pre-built Images.  No 20 GB TSV.  No RAM pre-load of all frames.
// Each frame is built and sent immediately, streaming to the viewer.

use anyhow::Result;
use log::info;
use rustdf::data::dataset::TimsDataset;
use rustdf::data::handle::TimsData;
use rustdf::data::meta::read_dda_precursor_meta;
use std::collections::HashMap;
use std::path::Path;

use crate::heatmap::{make_heatmap_rgb, MZ_BINS, MOB_BINS};
use crate::chrom::read_chrom_trace_by_desc;

use rerun::{
    RecordingStreamBuilder,
    Points3D, Points2D, Image, Scalars, SeriesLines, TextDocument,
};

const TIMELINE:     &str = "rt_seconds";
const MS1_NUMERIC:  i32  = 0;
const MS2_DDA:      i32  = 8;

pub fn run(
    d_path: &Path,
    max_frames: usize,
    top_n: usize,
    save: Option<&Path>,
    connect: Option<&str>,
) -> Result<()> {
    let builder = RecordingStreamBuilder::new("timstof-rerun");

    let rec = match (save, connect) {
        (Some(p), _) => { info!("Saving to {}", p.display()); builder.save(p)? }
        (_, Some(addr)) => { info!("Connecting to {}", addr); builder.connect_grpc_opts(addr)? }
        _ => { info!("Spawning Rerun viewer"); builder.spawn()? }
    };

    let d_str   = d_path.to_str().unwrap();
    let dataset = TimsDataset::new("", d_str, false, false);
    let n       = dataset.loader.get_frame_count() as usize;
    info!("{} frames in {}", n, d_path.display());

    // Load precursor metadata for MS2 markers
    let prec_by_frame = load_precursors(d_str);

    log_chromatograms(&rec, d_path)?;
    log_ms2_markers(&rec, d_str, &prec_by_frame)?;
    log_layout_info(&rec)?;

    // MS1 frames only
    let ms1_ids: Vec<u32> = (1..=n as u32)
        .filter(|&fid| dataset.loader.get_frame(fid).ms_type.ms_type_numeric() == MS1_NUMERIC)
        .collect();

    let take = if max_frames == 0 { ms1_ids.len() } else { ms1_ids.len().min(max_frames) };
    info!("Streaming {} / {} MS1 frames  top_n={}", take, ms1_ids.len(), top_n);

    for (n_done, &fid) in ms1_ids[..take].iter().enumerate() {
        let frame  = dataset.loader.get_frame(fid);
        let rt     = frame.ims_frame.retention_time;
        let mz     = &*frame.ims_frame.mz;
        let mob    = &*frame.ims_frame.mobility;
        let inten  = &*frame.ims_frame.intensity;

        rec.set_duration_secs(TIMELINE, rt);

        let (mz_f, mob_f, int_f) = top_n_filter(mz, mob, inten, top_n);
        if mz_f.is_empty() { continue; }

        // --- 3D ion map: (mz, mobility, log_intensity) with plasma colour ---
        let log_max = int_f.iter().cloned().fold(0.0f32, f32::max).ln_1p().max(1.0);
        let positions: Vec<[f32; 3]> = mz_f.iter().zip(mob_f.iter()).zip(int_f.iter())
            .map(|((&m, &k), &i)| [m, k, i.ln_1p()])
            .collect();
        let colours: Vec<rerun::Color> = int_f.iter()
            .map(|&i| plasma_colour(i.ln_1p() / log_max))
            .collect();
        let radii = vec![0.002f32; positions.len()];

        rec.log("ion_map",
            &Points3D::new(&positions)
                .with_colors(colours)
                .with_radii(radii))?;

        // --- Frame heatmap image (m/z vs 1/K0, plasma LUT) ---
        let rgb = make_heatmap_rgb(&mz_f, &mob_f, &int_f);
        rec.log("frame_heatmap",
            &Image::from_rgb24(rgb, [MZ_BINS as u32, MOB_BINS as u32]))?;

        // --- Mobilogram: (ln(I), 1/K0) ---
        let mob_pts: Vec<[f32; 2]> = int_f.iter().zip(mob_f.iter())
            .filter(|(&i, _)| i > 0.0)
            .map(|(&i, &k)| [i.ln_1p(), k])
            .collect();
        if !mob_pts.is_empty() {
            rec.log("mobilogram", &Points2D::new(&mob_pts))?;
        }

        // --- Summed spectrum: (mz, ln(I)) ---
        let spec_pts: Vec<[f32; 2]> = mz_f.iter().zip(int_f.iter())
            .filter(|(_, &i)| i > 0.0)
            .map(|(&m, &i)| [m, i.ln_1p()])
            .collect();
        if !spec_pts.is_empty() {
            rec.log("spectrum", &Points2D::new(&spec_pts))?;
        }

        // --- MS2 precursors for this frame (if any) ---
        if let Some(precs) = prec_by_frame.get(&(fid as i64)) {
            let prec_pts: Vec<[f32; 2]> = precs.iter()
                .filter_map(|p| {
                    let mz = p.precursor_mz_monoisotopic.unwrap_or(p.precursor_mz_average);
                    if mz > 0.0 { Some([mz as f32, 0.0f32]) } else { None }
                })
                .collect();
            if !prec_pts.is_empty() {
                rec.log("ms2_precursors",
                    &Points2D::new(&prec_pts)
                        .with_colors([rerun::Color::from_rgb(255, 80, 0)]))?;
            }
        }

        if (n_done + 1) % 200 == 0 {
            info!("  {} frames...", n_done + 1);
        }
    }

    info!("Done. {} frames.", take);
    if save.is_none() && connect.is_none() {
        info!("Viewer running - drag the RT scrubber to a BPC peak (30-60 min)");
        info!("The ion_map, frame_heatmap, mobilogram, spectrum all update automatically.");
        loop { std::thread::sleep(std::time::Duration::from_secs(1)); }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Plasma colour for a normalised value t in [0,1]
// ---------------------------------------------------------------------------
fn plasma_colour(t: f32) -> rerun::Color {
    let t = t.clamp(0.0, 1.0);
    let r = (t * 3.0 - 0.5).clamp(0.0, 1.0);
    let g = (t * 2.5 - 1.2).clamp(0.0, 1.0);
    let b = (0.8 - t * 1.5).clamp(0.0, 1.0);
    // Force black for empty bins
    if r == 0.0 && g == 0.0 && b == 0.0 {
        return rerun::Color::from_rgb(0, 0, 0);
    }
    rerun::Color::from_rgb(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    )
}

// ---------------------------------------------------------------------------
// Top-N filter
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Precursor metadata
// ---------------------------------------------------------------------------
fn load_precursors(d_str: &str) -> HashMap<i64, Vec<rustdf::data::meta::DDAPrecursorMeta>> {
    match read_dda_precursor_meta(d_str).map_err(|e| anyhow::anyhow!("{e}")) {
        Ok(precs) => {
            let mut map: HashMap<i64, Vec<_>> = HashMap::new();
            for p in precs { map.entry(p.precursor_frame_id).or_default().push(p); }
            map
        }
        Err(e) => { log::warn!("Could not load precursors: {}", e); HashMap::new() }
    }
}

/// Log MS2 precursor m/z as Scalar on a separate entity
/// so you can see which RT positions have MS/MS fragmentation.
fn log_ms2_markers(
    rec: &rerun::RecordingStream,
    d_str: &str,
    prec_by_frame: &HashMap<i64, Vec<rustdf::data::meta::DDAPrecursorMeta>>,
) -> Result<()> {
    use rusqlite::{Connection, OpenFlags};
    let tdf = format!("{}/analysis.tdf", d_str);
    let conn = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)?;

    rec.log_static("chromatogram/MS2_count",
        &SeriesLines::new()
            .with_colors([rerun::Color::from_rgb(200, 0, 200)])
            .with_widths([0.5_f32]))?;

    let mut stmt = conn.prepare(
        "SELECT Time, MsMsType FROM Frames WHERE MsMsType != 0 ORDER BY Id")?;
    let rows: Vec<(f64, i32)> = stmt.query_map([], |r| Ok((r.get(0)?, r.get(1)?)))?
        .filter_map(|r| r.ok()).collect();

    let rt0 = rows.first().map(|(t,_)| *t).unwrap_or(0.0);
    for (rt, _mstype) in &rows {
        rec.set_duration_secs(TIMELINE, rt - rt0);
        // Log count of precursors at this RT as a scalar (so MS2 density shows in chromatogram area)
        let n_prec = prec_by_frame.get(&0).map(|v| v.len()).unwrap_or(1) as f64;
        rec.log("chromatogram/MS2_count", &Scalars::new([n_prec]))?;
    }
    info!("  MS2 markers: {} frames", rows.len());
    Ok(())
}

fn log_chromatograms(rec: &rerun::RecordingStream, d_path: &Path) -> Result<()> {
    let chrom_db = d_path.join("chromatography-data.sqlite");
    if !chrom_db.exists() {
        log::warn!("chromatography-data.sqlite not found");
        return Ok(());
    }
    rec.log_static("chromatogram/BPC",
        &SeriesLines::new()
            .with_colors([rerun::Color::from_rgb(0, 180, 255)])
            .with_widths([1.0_f32]))?;
    rec.log_static("chromatogram/TIC",
        &SeriesLines::new()
            .with_colors([rerun::Color::from_rgb(255, 140, 0)])
            .with_widths([1.0_f32]))?;

    for (frag, entity) in &[("BPC","chromatogram/BPC"),("TIC","chromatogram/TIC")] {
        if let Ok(Some((times, intens))) = crate::chrom::read_chrom_trace_by_desc(&chrom_db, frag) {
            let rt0 = times.first().copied().unwrap_or(0.0);
            for (t, v) in times.iter().zip(intens.iter()) {
                rec.set_duration_secs(TIMELINE, t - rt0);
                rec.log(*entity, &Scalars::new([*v as f64]))?;
            }
            info!("  {} : {} points", entity, times.len());
        }
    }
    Ok(())
}

fn log_layout_info(rec: &rerun::RecordingStream) -> Result<()> {
    rec.log_static("info/layout", &TextDocument::new(
        "## timsTOF 5D Viewer\n\n\
        **5 dimensions, all linked to the timeline:**\n\
        - RT (X axis of BPC/TIC) = timeline scrubber\n\
        - m/z = X axis of ion_map and spectrum\n\
        - 1/K0 (mobility) = Y axis of ion_map and mobilogram\n\
        - Intensity = colour + Z of ion_map\n\
        - MS type = separate BPC (MS1) vs MS2_count entities\n\n\
        **How to interact:**\n\
        1. Drag the RT scrubber to a BPC peak (try 30-60 min region)\n\
        2. ion_map, frame_heatmap, mobilogram and spectrum all update\n\
        3. In the ion_map 3D view: left-drag to rotate, scroll to zoom\n\
        4. MS2 precursor m/z shown as orange dots in the spectrum panel\n\n\
        Calibration: rustdf polynomial fit (no Bruker SDK)."
    ).with_media_type(rerun::MediaType::markdown()))?;
    Ok(())
}
