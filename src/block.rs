//! Per-block sufficient statistics for the compound-symmetry REML.
//!
//! The design X is shared across genes, so each block's X'X and column sums are
//! computed once (`BlockDesign`); the y-dependent terms (`GeneStats`) are filled
//! per gene, keeping the per-gene work cheap and embarrassingly parallel.

pub struct BlockDesign {
    pub size: usize,
    rows: Vec<usize>,
    /// p×p row-major X'X within the block
    pub xx: Vec<f64>,
    /// length-p X column sums within the block
    pub sx: Vec<f64>,
}

pub struct Design {
    pub blocks: Vec<BlockDesign>,
    pub n: usize,
    pub p: usize,
    pub max_size: usize,
}

/// y-dependent block statistics for one gene, aligned 1:1 with `Design.blocks`.
pub struct GeneStats {
    /// length-p X'y per block, concatenated
    pub xy: Vec<f64>,
    /// sum of y per block
    pub sy: Vec<f64>,
    /// y'y per block
    pub yy: Vec<f64>,
}

impl Design {
    pub fn new(x: &[Vec<f64>], block_ids: &[usize]) -> Design {
        let p = x[0].len();
        let nblk = block_ids.iter().copied().max().map_or(0, |m| m + 1);
        let mut rows: Vec<Vec<usize>> = vec![Vec::new(); nblk];
        for (i, &b) in block_ids.iter().enumerate() {
            rows[b].push(i);
        }
        let mut blocks = Vec::with_capacity(nblk);
        let mut max_size = 0;
        for r in rows {
            if r.is_empty() {
                continue;
            }
            max_size = max_size.max(r.len());
            let mut xx = vec![0.0f64; p * p];
            let mut sx = vec![0.0f64; p];
            for &i in &r {
                let xi = &x[i];
                for j in 0..p {
                    sx[j] += xi[j];
                    for k in 0..p {
                        xx[j * p + k] += xi[j] * xi[k];
                    }
                }
            }
            blocks.push(BlockDesign {
                size: r.len(),
                rows: r,
                xx,
                sx,
            });
        }
        Design {
            blocks,
            n: block_ids.len(),
            p,
            max_size,
        }
    }

    pub fn gene_stats(&self, x: &[Vec<f64>], y: &[f64]) -> GeneStats {
        let p = self.p;
        let nb = self.blocks.len();
        let mut xy = vec![0.0f64; p * nb];
        let mut sy = vec![0.0f64; nb];
        let mut yy = vec![0.0f64; nb];
        for (bi, blk) in self.blocks.iter().enumerate() {
            for &i in &blk.rows {
                let yi = y[i];
                sy[bi] += yi;
                yy[bi] += yi * yi;
                let xi = &x[i];
                for j in 0..p {
                    xy[bi * p + j] += xi[j] * yi;
                }
            }
        }
        GeneStats { xy, sy, yy }
    }
}
