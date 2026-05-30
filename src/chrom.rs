// src/chrom.rs
//
// Read chromatography traces from chromatography-data.sqlite.
// Blob encoding: Times = little-endian f64, Intensities = little-endian f32.
// (Same logic as timsread.cpp -chrom, rewritten in Rust.)
//
// Trace IDs:
//   18 = BPC +-MS          19 = TIC +-AllMS/MS       20 = TIC +-MS
//    1 = Flow A             2 = Pressure A             7 = Solvent B %
//   14 = Column Temp       15 = Valve I               16 = Valve T

use anyhow::{Context, Result};
use rusqlite::{Connection, OpenFlags};
use std::path::{Path, PathBuf};

pub struct TraceInfo {
    pub id:          i64,
    pub description: String,
    pub instrument:  String,
    pub unit:        i64,
}

/// Read a single trace's (times, intensities) from the sqlite file.
pub fn read_chrom_trace(db_path: &Path, trace_id: u32) -> Result<(Vec<f64>, Vec<f32>)> {
    let conn = Connection::open_with_flags(db_path, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .with_context(|| format!("opening {}", db_path.display()))?;

    let mut stmt = conn.prepare(
        "SELECT Times, Intensities FROM TraceChunks WHERE Trace=? ORDER BY rowid",
    )?;

    let mut all_times:   Vec<f64> = Vec::new();
    let mut all_intens:  Vec<f32> = Vec::new();

    let rows = stmt.query_map([trace_id], |row| {
        let t_blob: Vec<u8> = row.get(0)?;
        let i_blob: Vec<u8> = row.get(1)?;
        Ok((t_blob, i_blob))
    })?;

    for row in rows {
        let (t_blob, i_blob) = row?;
        // Times: f64 LE
        for chunk in t_blob.chunks_exact(8) {
            all_times.push(f64::from_le_bytes(chunk.try_into().unwrap()));
        }
        // Intensities: f32 LE
        for chunk in i_blob.chunks_exact(4) {
            all_intens.push(f32::from_le_bytes(chunk.try_into().unwrap()));
        }
    }

    Ok((all_times, all_intens))
}

fn unit_name(unit: i64) -> &'static str {
    match unit {
        2 => "mL_per_min",
        3 => "bar",
        4 => "uL",
        5 => "degC",
        6 => "counts",
        7 => "deg",
        _ => "unknown",
    }
}

fn sanitise(desc: &str) -> String {
    desc.chars()
        .map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' })
        .collect::<String>()
        .split('_')
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join("_")
}

/// Dump all traces from a chromatography sqlite file to a directory.
pub fn dump_to_dir(db_path: &Path, out_dir: &Path, prefix: &str) -> Result<()> {
    let conn = Connection::open_with_flags(db_path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;

    // Load trace metadata
    let mut sources: Vec<TraceInfo> = Vec::new();
    {
        let mut stmt = conn.prepare(
            "SELECT Id, Description, Instrument, Unit FROM TraceSources ORDER BY Id",
        )?;
        let rows = stmt.query_map([], |row| {
            Ok(TraceInfo {
                id:          row.get(0)?,
                description: row.get(1)?,
                instrument:  row.get(2)?,
                unit:        row.get(3)?,
            })
        })?;
        for r in rows { sources.push(r?); }
    }

    // Which traces have data?
    let mut has_data: std::collections::HashSet<i64> = Default::default();
    {
        let mut stmt = conn.prepare("SELECT DISTINCT Trace FROM TraceChunks")?;
        let rows = stmt.query_map([], |row| row.get::<_, i64>(0))?;
        for r in rows { has_data.insert(r?); }
    }

    std::fs::create_dir_all(out_dir)?;

    for src in &sources {
        if !has_data.contains(&src.id) { continue; }

        let (times, intens) = read_chrom_trace(db_path, src.id as u32)?;
        if times.is_empty() { continue; }

        let fname = out_dir.join(format!(
            "{}_{}_trace{}.txt",
            prefix,
            sanitise(&src.description),
            src.id
        ));
        let mut f = std::fs::File::create(&fname)?;
        use std::io::Write;
        writeln!(f, "# Trace_ID: {}", src.id)?;
        writeln!(f, "# Description: {}", src.description)?;
        writeln!(f, "# Instrument: {}", src.instrument)?;
        writeln!(f, "# Unit: {}", unit_name(src.unit))?;
        writeln!(f, "Time_seconds\tTime_minutes\tValue")?;
        let n = times.len().min(intens.len());
        for i in 0..n {
            writeln!(f, "{:.6}\t{:.6}\t{:.6}", times[i], times[i] / 60.0, intens[i])?;
        }
        println!("  [{:2}] {:30}  {} pts -> {}", src.id, src.description, n, fname.display());
    }

    Ok(())
}

/// CLI entry point for `chrom` subcommand.
pub fn run(d_path: &Path, out_dir: Option<&Path>) -> Result<()> {
    let base = d_path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let default_out = PathBuf::from(format!("{}_chrom", d_path.display()));
    let out = out_dir.unwrap_or(&default_out);

    for (suffix, label) in &[("", "run"), ("-pre", "pre")] {
        let fname = format!("chromatography-data{}.sqlite", suffix);
        let db = d_path.join(&fname);
        if db.exists() {
            println!("\n--- {} ({}) ---", fname, label);
            dump_to_dir(&db, out, &format!("{}_{}", base, label))?;
        } else {
            println!("Not found: {} (skipping)", db.display());
        }
    }
    println!("\nOutput: {}", out.display());
    Ok(())
}

/// Find a trace by description substring (case-insensitive) and read its data.
/// Returns Ok(None) if no matching trace found or trace has no data.
pub fn read_chrom_trace_by_desc(
    db_path: &Path,
    desc_fragment: &str,
) -> anyhow::Result<Option<(Vec<f64>, Vec<f32>)>> {
    use rusqlite::{Connection, OpenFlags};
    let conn = Connection::open_with_flags(db_path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;

    // Find the trace ID whose description contains desc_fragment
    // Prefer MS trace (Instrument like "timsTOF") over HPLC traces
    let mut stmt = conn.prepare(
        "SELECT ts.Id FROM TraceSources ts \
         WHERE UPPER(ts.Description) LIKE UPPER(?1) \
         AND EXISTS (SELECT 1 FROM TraceChunks tc WHERE tc.Trace = ts.Id) \
         ORDER BY \
           CASE WHEN UPPER(ts.Instrument) LIKE '%TIMS%' THEN 0 ELSE 1 END, \
           ts.Id \
         LIMIT 1"
    )?;
    let pattern = format!("%{}%", desc_fragment);
    let trace_id: Option<u32> = stmt
        .query_map([&pattern], |row| row.get::<_, u32>(0))?
        .next()
        .transpose()?;

    match trace_id {
        None => Ok(None),
        Some(tid) => {
            let result = read_chrom_trace(db_path, tid)?;
            if result.0.is_empty() { Ok(None) } else { Ok(Some(result)) }
        }
    }
}
