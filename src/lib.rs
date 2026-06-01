//! limma duplicateCorrelation: the consensus intra-block correlation.
//!
//! Clean-room reimplementation following Smyth, Michaud & Scott (2005),
//! Bioinformatics 21:2067, doi:10.1093/bioinformatics/bti270. No limma (GPL)
//! source was consulted; the method follows the published paper and is
//! validated black-box against the limma binary.

mod block;
mod matrix;
mod reml;

use std::io::{BufWriter, Write};
use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

pub use matrix::{read_block, read_design, read_expr};

pub const UPPER: f64 = 0.99;
pub const TRIM: f64 = 0.15;

pub struct Options<'a> {
    pub expr: &'a Path,
    pub design: &'a Path,
    pub block: &'a Path,
    pub threads: usize,
}

pub struct Results {
    pub consensus: f64,
    /// per-gene atanh(rho); non-estimable genes are None and excluded from the consensus
    pub atanh: Vec<Option<f64>>,
    pub genes: Vec<String>,
    pub n_used: usize,
}

pub fn run(opts: &Options) -> Result<Results> {
    let expr = read_expr(opts.expr)?;
    let design = read_design(opts.design)?;
    let block_ids = read_block(opts.block)?;

    let nsamp = expr.samples.len();
    if design.x.len() != nsamp {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows, expression has {nsamp} samples",
            design.x.len()
        )));
    }
    if block_ids.len() != nsamp {
        return Err(RsomicsError::InvalidInput(format!(
            "block factor has {} entries, expression has {nsamp} samples",
            block_ids.len()
        )));
    }

    let d = block::Design::new(&design.x, &block_ids);
    if d.n <= d.p {
        return Err(RsomicsError::InvalidInput(format!(
            "{} samples <= {} coefficients (no residual df)",
            d.n, d.p
        )));
    }
    let lower = -1.0 / (d.max_size as f64 - 1.0) + 0.01;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build()
        .map_err(|e| RsomicsError::InvalidInput(e.to_string()))?;

    let atanh: Vec<Option<f64>> = pool.install(|| {
        expr.y
            .par_iter()
            .map(|row| {
                let g = d.gene_stats(&design.x, row);
                reml::gene_rho(&d, &g, lower, UPPER).map(|rho| rho.atanh())
            })
            .collect()
    });

    let used: Vec<f64> = atanh
        .iter()
        .filter_map(|&v| v)
        .filter(|v| v.is_finite())
        .collect();
    let consensus = (trimmed_mean(used.clone(), TRIM)).tanh();

    Ok(Results {
        consensus,
        atanh,
        genes: expr.genes,
        n_used: used.len(),
    })
}

/// R's `mean(x, trim)`: drop floor(n*trim) from each tail, mean the rest.
fn trimmed_mean(mut v: Vec<f64>, trim: f64) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lo = (n as f64 * trim).floor() as usize;
    let hi = n - lo;
    let slice = &v[lo..hi];
    slice.iter().sum::<f64>() / slice.len() as f64
}

pub fn write_results(res: &Results, out: &mut dyn Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    let mut fmt = ryu::Buffer::new();
    writeln!(w, "consensus.correlation\t{}", fmt.format(res.consensus))
        .map_err(RsomicsError::Io)?;
    writeln!(w, "gene\tatanh.correlation").map_err(RsomicsError::Io)?;
    let mut line = String::with_capacity(64);
    for (gene, a) in res.genes.iter().zip(&res.atanh) {
        line.clear();
        line.push_str(gene);
        line.push('\t');
        match a {
            Some(v) => line.push_str(fmt.format(*v)),
            None => line.push_str("NA"),
        }
        line.push('\n');
        w.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trimmed_mean_matches_r() {
        // R: mean(c(1,2,3,4,5,6,7,8,9,10), trim=0.15) drops floor(10*0.15)=1 each side
        let v = vec![1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];
        let m = trimmed_mean(v, 0.15);
        assert!((m - 5.5).abs() < 1e-12);
    }

    #[test]
    fn trimmed_mean_no_trim() {
        let v = vec![1., 2., 3., 4., 5., 6.];
        // floor(6*0.15)=0 -> plain mean
        assert!((trimmed_mean(v, 0.15) - 3.5).abs() < 1e-12);
    }
}
