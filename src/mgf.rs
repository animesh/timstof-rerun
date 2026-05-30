// src/mgf.rs
use anyhow::Result;
use log::info;
use rustdf::data::dataset::TimsDataset;
use rustdf::data::handle::TimsData;
use rustdf::data::meta::read_dda_precursor_meta;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

// ms_type_numeric(): 0=Precursor(MS1), 8=FragmentDda(MS2), 9=FragmentDia
const MS2_DDA: i32 = 8;

pub fn run(d_path: &Path, out: Option<&Path>) -> Result<()> {
    let base     = d_path.file_name().unwrap_or_default().to_string_lossy().to_string();
    let default_out = PathBuf::from(format!("{}.mgf", base));
    let out_path = out.unwrap_or(&default_out);

    let d_str = d_path.to_str().unwrap();
    let dataset = TimsDataset::new("", d_str, false, false);
    let n = dataset.loader.get_frame_count() as usize;
    info!("{} frames", n);

    // read_dda_precursor_meta takes the .d folder path (opens its own connection)
    let precursors = read_dda_precursor_meta(d_str)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    info!("{} DDA precursors", precursors.len());

    // Index precursors by frame_id
    let mut prec_by_frame: HashMap<i64, Vec<_>> = HashMap::new();
    for p in &precursors {
        prec_by_frame.entry(p.precursor_frame_id).or_default().push(p);
    }

    let file = std::fs::File::create(out_path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "# MGF from {} (timstof-rerun, rustdf - no Bruker SDK)", d_str)?;
    writeln!(w)?;

    let mut ms2_written = 0usize;
    let mut skipped     = 0usize;

    for fid in 1..=n as u32 {
        let frame = dataset.loader.get_frame(fid);
        // Use integer comparison to avoid mscore version conflicts
        if frame.ms_type.ms_type_numeric() != MS2_DDA {
            continue;
        }

        let rt    = frame.ims_frame.retention_time;
        let mz    = &*frame.ims_frame.mz;
        let inten = &*frame.ims_frame.intensity;
        let mob   = &*frame.ims_frame.mobility;

        let precs = match prec_by_frame.get(&(fid as i64)) {
            Some(v) => v,
            None    => { skipped += 1; continue; }
        };

        for prec in precs {
            let prec_mz = match prec.precursor_mz_monoisotopic {
                Some(m) if m > 0.0 => m,
                _ => {
                    let avg = prec.precursor_mz_average;
                    if avg <= 0.0 { skipped += 1; continue; }
                    avg
                }
            };

            let avg_mob = if !mob.is_empty() {
                mob.iter().sum::<f64>() / mob.len() as f64
            } else { 0.0 };

            writeln!(w, "BEGIN IONS")?;
            writeln!(w, "TITLE=Frame_{}_Precursor_{}", fid, prec.precursor_id)?;
            writeln!(w, "RTINSECONDS={:.6}", rt)?;
            writeln!(w, "MOBILITY={:.6}", avg_mob)?;
            if let Some(z) = prec.precursor_charge {
                if z > 0 { writeln!(w, "CHARGE={}+", z)?; }
            }
            write!(w, "PEPMASS={:.6}", prec_mz)?;
            if prec.precursor_total_intensity > 0.0 {
                write!(w, " {:.0}", prec.precursor_total_intensity)?;
            }
            writeln!(w)?;

            for (m, i) in mz.iter().zip(inten.iter()) {
                if *i <= 0.0 { continue; }
                writeln!(w, "{:.6} {:.0}", m, i)?;
            }
            writeln!(w, "END IONS")?;
            writeln!(w)?;
            ms2_written += 1;
        }

        if ms2_written % 10_000 == 0 && ms2_written > 0 {
            info!("  {} spectra...", ms2_written);
        }
    }

    info!("Written: {} spectra, {} skipped", ms2_written, skipped);
    info!("Output: {}", out_path.display());
    Ok(())
}
