// src/sql_dump.rs
//
// Dump every table and view in analysis.tdf as a TSV file.
// Dynamic: works for any TDF version (ddaPASEF, diaPASEF, DIA-PASEF).

use anyhow::Result;
use rusqlite::{Connection, OpenFlags, types::ValueRef};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Format a float in "shortest useful" form:
/// integers print without decimal, small/large values use scientific notation.
fn fmt_float(f: f64) -> String {
    if f == 0.0 {
        return "0".to_string();
    }
    let abs = f.abs();
    // Use scientific notation for very small or very large values
    if abs < 1e-4 || abs >= 1e10 {
        format!("{:.6e}", f)
    } else {
        // Use enough decimal places to be useful
        let s = format!("{:.10}", f);
        // Strip trailing zeros after decimal point
        let s = s.trim_end_matches('0');
        let s = s.trim_end_matches('.');
        s.to_string()
    }
}

fn cell_to_string(v: ValueRef) -> String {
    match v {
        ValueRef::Null        => String::new(),
        ValueRef::Integer(i)  => i.to_string(),
        ValueRef::Real(f)     => fmt_float(f),
        ValueRef::Text(s)     => String::from_utf8_lossy(s).to_string(),
        ValueRef::Blob(b)     => format!("[BLOB {} bytes]", b.len()),
    }
}

pub fn run(d_path: &Path, out_dir: Option<&Path>) -> Result<()> {
    let tdf = d_path.join("analysis.tdf");
    let base = d_path.file_name().unwrap_or_default().to_string_lossy().to_string();
    let default_out = PathBuf::from(format!("{}_tdf", d_path.display()));
    let out = out_dir.unwrap_or(&default_out);
    std::fs::create_dir_all(out)?;

    let conn = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)?;

    let mut stmt = conn.prepare(
        "SELECT name, type FROM sqlite_master \
         WHERE type IN ('table','view') ORDER BY type DESC, name",
    )?;
    let objects: Vec<(String, String)> = stmt
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .filter_map(|r| r.ok())
        .collect();

    println!("Found {} tables/views in {}", objects.len(), tdf.display());

    for (name, kind) in &objects {
        let path = out.join(format!("{}_{}.txt", base, name));
        let mut f = std::fs::File::create(&path)?;

        let mut stmt = conn.prepare(&format!("SELECT * FROM [{}]", name))?;
        let ncols = stmt.column_count();
        if ncols == 0 { continue; }

        let header: Vec<String> = (0..ncols)
            .map(|i| stmt.column_name(i).unwrap_or("?").to_string())
            .collect();
        writeln!(f, "{}", header.join("\t"))?;

        let mut nrows = 0i64;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let cells: Vec<String> = (0..ncols)
                .map(|i| cell_to_string(row.get_ref(i).unwrap_or(ValueRef::Null)))
                .collect();
            writeln!(f, "{}", cells.join("\t"))?;
            nrows += 1;
        }

        println!("  {:6} {:<30} {:>10} rows -> {}", kind, name, nrows, path.display());
    }

    println!("\nOutput: {}", out.display());
    Ok(())
}
