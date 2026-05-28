// src/mgf.rs
//
// Write MS/MS spectra in MGF format using rustdf.
// rustdf provides calibrated m/z and 1/K0 directly from the SQLite polynomials.

use anyhow::Result;
use log::info;
use rustdf::data::dataset::TimsDataset;
use rustdf::data::meta::MsType;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn run(d_path: &Path, out: Option<&Path>) -> Result<()> {
    let base = d_path.file_name().unwrap_or_default().to_string_lossy().to_string();
    let default_out = PathBuf::from(format!("{}.mgf", d_path.display()));
    let out_path = out.unwrap_or(&default_out);

    let dataset = TimsDataset::new(d_path.to_str().unwrap())?;
    let n = dataset.get_frame_count();
    info!("{} frames", n);

    let file = std::fs::File::create(out_path)?;
    let mut w = BufWriter::new(file);

    writeln!(w, "# MGF from {} via timstof-rerun (rustdf)", d_path.display())?;
    writeln!(w, "# ZERO-FILTER ONLY")?;
    writeln!(w)?;

    let mut ms2_count = 0usize;
    let mut skipped   = 0usize;

    // rustdf provides get_all_ms2_spectra() which returns fully assembled
    // spectra with precursor info already merged from Precursors + PasefFrameMsMsInfo
    let spectra = dataset.get_all_ms2_spectra()?;
    info!("{} MS/MS spectra", spectra.len());

    for spec in &spectra {
        let prec_mz = spec.precursor_mz;
        if prec_mz <= 0.0 { skipped += 1; continue; }

        writeln!(w, "BEGIN IONS")?;
        writeln!(w, "TITLE=Frame_{}_Precursor_{}", spec.frame_id, spec.precursor_id)?;
        writeln!(w, "RTINSECONDS={:.6}", spec.retention_time)?;
        writeln!(w, "MOBILITY={:.6}", spec.one_over_k0)?;
        if spec.precursor_charge > 0 {
            writeln!(w, "CHARGE={}+", spec.precursor_charge)?;
        }
        write!(w, "PEPMASS={:.6}", prec_mz)?;
        if spec.precursor_intensity > 0.0 {
            write!(w, " {:.0}", spec.precursor_intensity)?;
        }
        if spec.collision_energy > 0.0 {
            writeln!(w)?;
            write!(w, "COLLISION_ENERGY={:.1}", spec.collision_energy)?;
        }
        writeln!(w)?;

        for (mz, i) in spec.mz_values.iter().zip(spec.intensities.iter()) {
            if *i == 0 { continue; }
            writeln!(w, "{:.6} {}", mz, i)?;
        }

        writeln!(w, "END IONS")?;
        writeln!(w)?;
        ms2_count += 1;

        if ms2_count % 10_000 == 0 {
            info!("  {} spectra written...", ms2_count);
        }
    }

    info!("Written: {} spectra, {} skipped (no precursor m/z)", ms2_count, skipped);
    info!("Output: {}", out_path.display());
    Ok(())
}
