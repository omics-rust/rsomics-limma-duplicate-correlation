use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

fn open(path: &Path) -> Result<BufReader<File>> {
    let f = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    Ok(BufReader::new(f))
}

/// Parse a numeric cell, rejecting non-finite values loudly. Rust's f64 parser
/// accepts `inf`/`nan` literals; limma drops such observations, but a silently
/// accepted NaN would poison the per-gene REML objective at every rho and leak a
/// garbage correlation into the consensus, so we fail fast instead.
fn parse_f64(s: &str) -> Result<f64> {
    let t = s.trim();
    let v = t
        .parse::<f64>()
        .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric value '{t}'")))?;
    if !v.is_finite() {
        return Err(RsomicsError::InvalidInput(format!(
            "non-finite value '{t}'"
        )));
    }
    Ok(v)
}

pub struct Expr {
    pub samples: Vec<String>,
    pub genes: Vec<String>,
    /// row-major [gene][sample]
    pub y: Vec<Vec<f64>>,
}

pub fn read_expr(path: &Path) -> Result<Expr> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty expression matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let samples: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if samples.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "expression matrix needs at least one sample column".into(),
        ));
    }
    let mut genes = Vec::new();
    let mut y = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing gene id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != samples.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "gene '{gene}' has {} values, header declares {} samples",
                row.len(),
                samples.len()
            )));
        }
        genes.push(gene.to_string());
        y.push(row);
    }
    if genes.is_empty() {
        return Err(RsomicsError::InvalidInput("no genes in matrix".into()));
    }
    Ok(Expr { samples, genes, y })
}

pub struct Design {
    pub coef_names: Vec<String>,
    /// row-major [sample][coef]
    pub x: Vec<Vec<f64>>,
}

/// Design TSV: first column = sample id, header first cell may be empty or a
/// label; remaining columns are coefficient names with numeric model-matrix
/// entries, one row per sample in sample order.
pub fn read_design(path: &Path) -> Result<Design> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty design matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let coef_names: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if coef_names.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "design matrix needs at least one coefficient column".into(),
        ));
    }
    let mut x = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let id = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing design row id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != coef_names.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "design row '{id}' has {} values, header declares {} coefficients",
                row.len(),
                coef_names.len()
            )));
        }
        x.push(row);
    }
    Ok(Design { coef_names, x })
}

/// Block TSV: first column = sample id, second column = block label. The label
/// may be any string; samples sharing a label form one block. A header line is
/// consumed unconditionally (first row is treated as a header).
pub fn read_block(path: &Path) -> Result<Vec<usize>> {
    let mut lines = open(path)?.lines();
    lines.next();
    let mut labels: Vec<String> = Vec::new();
    let mut ids: Vec<usize> = Vec::new();
    let mut seen: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let label = line
            .split('\t')
            .nth(1)
            .ok_or_else(|| RsomicsError::InvalidInput("block row missing label column".into()))?
            .trim()
            .to_string();
        let next = seen.len();
        let id = *seen.entry(label.clone()).or_insert(next);
        if id == labels.len() {
            labels.push(label);
        }
        ids.push(id);
    }
    if ids.is_empty() {
        return Err(RsomicsError::InvalidInput("no block assignments".into()));
    }
    Ok(ids)
}
