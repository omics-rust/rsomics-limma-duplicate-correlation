//! Per-gene REML of the compound-symmetry intra-block correlation.
//!
//! Smyth, Michaud & Scott (2005), Bioinformatics 21:2067, doi:10.1093/bioinformatics/bti270.
//! Each block of size m has within-block covariance proportional to
//! (1-rho) I + rho J. The residual REML log-likelihood profiled over the scale
//! variance is maximised over rho on a single shared interval; per-gene rho is
//! mapped through atanh and the consensus is the (trimmed) atanh-mean.

use crate::block::{Design, GeneStats};

/// Profiled REML terms at a single rho: the residual sum of squares, residual
/// df, and the two log-determinants entering the objective.
///
/// With a = 1/(1-rho) and g = rho/(1+(m-1)rho), the compound-symmetry inverse
/// is Vinv = a(I - g J), so each block contraction is dot products plus a
/// rank-one correction from the within-block X and y sums.
struct RemlTerms {
    rss: f64,
    df: f64,
    logdet_v: f64,
    logdet_xtvix: f64,
}

fn reml_terms(rho: f64, d: &Design, g: &GeneStats) -> RemlTerms {
    let p = d.p;
    let mut xtvix = vec![0.0f64; p * p];
    let mut xtviy = vec![0.0f64; p];
    let mut ytviy = 0.0f64;
    let mut logdet_v = 0.0f64;

    for (bi, blk) in d.blocks.iter().enumerate() {
        let m = blk.size as f64;
        let a = 1.0 / (1.0 - rho);
        let gg = rho / (1.0 + (m - 1.0) * rho);
        let sy = g.sy[bi];

        ytviy += a * (g.yy[bi] - gg * sy * sy);
        for j in 0..p {
            let sxj = blk.sx[j];
            xtviy[j] += a * (g.xy[bi * p + j] - gg * sxj * sy);
            for k in 0..p {
                xtvix[j * p + k] += a * (blk.xx[j * p + k] - gg * sxj * blk.sx[k]);
            }
        }
        logdet_v += (m - 1.0) * (1.0 - rho).ln() + (1.0 + (m - 1.0) * rho).ln();
    }

    let (chol, logdet_xtvix) = cholesky_logdet(&xtvix, p);
    let beta = chol_solve(&chol, &xtviy, p);
    let mut rss = ytviy;
    for j in 0..p {
        rss -= xtviy[j] * beta[j];
    }
    RemlTerms {
        rss,
        df: (d.n - p) as f64,
        logdet_v,
        logdet_xtvix,
    }
}

/// Negative REML profile log-likelihood at a single rho, up to a constant.
fn neg_reml(rho: f64, d: &Design, g: &GeneStats) -> f64 {
    let t = reml_terms(rho, d, g);
    let s2 = t.rss / t.df;
    0.5 * (t.df * s2.ln() + t.logdet_v + t.logdet_xtvix)
}

/// Lower-triangular Cholesky of a row-major p×p SPD matrix, returning the
/// factor and 2·sum(log diag) (= log determinant of the input).
fn cholesky_logdet(m: &[f64], p: usize) -> (Vec<f64>, f64) {
    let mut l = vec![0.0f64; p * p];
    let mut logdet = 0.0f64;
    for i in 0..p {
        for j in 0..=i {
            let mut s = m[i * p + j];
            for k in 0..j {
                s -= l[i * p + k] * l[j * p + k];
            }
            if i == j {
                let dd = s.max(1e-300).sqrt();
                l[i * p + j] = dd;
                logdet += 2.0 * dd.ln();
            } else {
                l[i * p + j] = s / l[j * p + j];
            }
        }
    }
    (l, logdet)
}

fn chol_solve(l: &[f64], rhs: &[f64], p: usize) -> Vec<f64> {
    let mut y = vec![0.0f64; p];
    for i in 0..p {
        let mut s = rhs[i];
        for k in 0..i {
            s -= l[i * p + k] * y[k];
        }
        y[i] = s / l[i * p + i];
    }
    let mut x = vec![0.0f64; p];
    for i in (0..p).rev() {
        let mut s = y[i];
        for k in (i + 1)..p {
            s -= l[k * p + i] * x[k];
        }
        x[i] = s / l[i * p + i];
    }
    x
}

/// REML estimate of rho for one gene, or None if the gene is non-estimable.
///
/// A residual sum of squares that has collapsed to zero (a constant / all-zero
/// gene the design fits exactly) drives the profile objective to -inf at every
/// rho, and Brent then drifts to a boundary and reports a spurious extreme
/// correlation. limma's mixedModel2Fit reports such genes as NA; we mirror that
/// by returning None, which excludes the gene from the trimmed-mean consensus.
pub fn gene_rho(d: &Design, g: &GeneStats, lower: f64, upper: f64) -> Option<f64> {
    let rho = brent_min(lower, upper, |r| neg_reml(r, d, g));
    if !rho.is_finite() {
        return None;
    }
    let t = reml_terms(rho, d, g);
    let obj = 0.5 * (t.df * (t.rss / t.df).ln() + t.logdet_v + t.logdet_xtvix);
    if t.rss <= 0.0 || !obj.is_finite() {
        return None;
    }
    Some(rho)
}

/// Brent's combined golden-section / parabolic-interpolation minimiser.
/// limma drives `optimize` to the true REML optimum (a tight tolerance), so we
/// match it rather than R's loose default of eps^0.25.
fn brent_min(ax: f64, bx: f64, mut f: impl FnMut(f64) -> f64) -> f64 {
    let tol = 1e-10;
    let c = 0.5 * (3.0 - 5.0f64.sqrt());
    let eps = f64::EPSILON.sqrt();

    let mut a = ax;
    let mut b = bx;
    let mut v = a + c * (b - a);
    let mut w = v;
    let mut x = v;
    let mut fx = f(x);
    let mut fv = fx;
    let mut fw = fx;
    let mut d = 0.0f64;
    let mut e = 0.0f64;

    loop {
        let xm = 0.5 * (a + b);
        let tol1 = eps * x.abs() + tol / 3.0;
        let tol2 = 2.0 * tol1;
        if (x - xm).abs() <= tol2 - 0.5 * (b - a) {
            break;
        }

        let mut use_golden = true;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let etemp = e;
            e = d;
            if p.abs() < (0.5 * q * etemp).abs() && p > q * (a - x) && p < q * (b - x) {
                d = p / q;
                let u = x + d;
                if (u - a) < tol2 || (b - u) < tol2 {
                    d = if xm - x >= 0.0 { tol1 } else { -tol1 };
                }
                use_golden = false;
            }
        }
        if use_golden {
            e = if x >= xm { a - x } else { b - x };
            d = c * e;
        }

        let u = if d.abs() >= tol1 {
            x + d
        } else if d > 0.0 {
            x + tol1
        } else {
            x - tol1
        };
        let fu = f(u);

        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    x
}
