//! Differential compat against limma duplicateCorrelation.
//!
//! - `golden_*` always runs: ours vs a committed R-captured result.
//! - `live_r_*` runs only when an Rscript with limma is found; it regenerates
//!   the oracle and diffs against ours (loud-skip otherwise).

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

// The deliverable — the consensus correlation — is value-exact.
const CONSENSUS_EPS: f64 = 1e-6;
// Per-gene atanh: limma's REML `optimize` (GPL) does not always converge to the
// true optimum, so a handful of genes drift from our tight optimum at the ~1e-3
// level; the trimmed-mean consensus absorbs that. Median deviation is ~1e-7.
const PERGENE_EPS: f64 = 2e-3;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-limma-duplicate-correlation"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn manifest(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

/// Parsed result: the consensus plus the per-gene atanh map.
struct Parsed {
    consensus: f64,
    atanh: BTreeMap<String, Option<f64>>,
}

fn parse(text: &str) -> Parsed {
    let mut lines = text.lines();
    let head = lines.next().unwrap();
    let consensus: f64 = head.split('\t').nth(1).unwrap().trim().parse().unwrap();
    lines.next(); // "gene\tatanh.correlation"
    let mut atanh = BTreeMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f.next().unwrap().to_string();
        let raw = f.next().unwrap().trim();
        let v = if raw == "NA" {
            None
        } else {
            Some(raw.parse::<f64>().unwrap())
        };
        atanh.insert(gene, v);
    }
    Parsed { consensus, atanh }
}

fn assert_close(a: &Parsed, b: &Parsed, label: &str) {
    let crel = (a.consensus - b.consensus).abs() / b.consensus.abs().max(1e-9);
    assert!(
        crel < CONSENSUS_EPS,
        "{label}: consensus ours={} ref={} rel={crel:e}",
        a.consensus,
        b.consensus
    );
    assert_eq!(a.atanh.len(), b.atanh.len(), "{label}: gene count mismatch");
    let mut max_rel = 0.0f64;
    for (gene, x) in &a.atanh {
        let y = b
            .atanh
            .get(gene)
            .unwrap_or_else(|| panic!("{label}: missing gene {gene}"));
        match (x, y) {
            (Some(vx), Some(vy)) => {
                let rel = (vx - vy).abs() / vy.abs().max(1e-9);
                max_rel = max_rel.max(rel);
                assert!(
                    rel < PERGENE_EPS,
                    "{label}: {gene} ours={vx} ref={vy} rel={rel:e}"
                );
            }
            (None, None) => {}
            _ => panic!("{label}: {gene} NA mismatch ours={x:?} ref={y:?}"),
        }
    }
    eprintln!("{label}: consensus rel={crel:e}, max per-gene rel={max_rel:e}");
}

fn run_ours_on(expr: &str, design: &str, block: &str) -> std::process::Output {
    Command::new(ours())
        .arg(golden(expr))
        .args(["--design", golden(design).to_str().unwrap()])
        .args(["--block", golden(block).to_str().unwrap()])
        .output()
        .unwrap()
}

fn run_ours() -> String {
    let out = run_ours_on("expr.tsv", "design.tsv", "blocks.tsv");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn golden_consensus() {
    let ours_out = run_ours();
    let expected = std::fs::read_to_string(golden("result.expected.tsv")).unwrap();
    assert_close(
        &parse(&ours_out),
        &parse(&expected),
        "duplicateCorrelation (golden)",
    );
}

/// A non-estimable gene (here an exactly all-zero row, which the design fits
/// with residual sum of squares 0) must be reported NA and dropped from the
/// consensus, matching limma's mixedModel2Fit. Without the guard the profile
/// REML objective is -inf at every rho, Brent drifts to the upper clamp, and the
/// gene leaks a spurious rho=0.99 into the trimmed-mean consensus. Expected
/// values captured from limma 3.62.1 duplicateCorrelation.
#[test]
fn golden_degenerate_gene_is_na() {
    let out = run_ours_on(
        "degenerate_expr.tsv",
        "degenerate_design.tsv",
        "degenerate_blocks.tsv",
    );
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let ours = parse(&String::from_utf8(out.stdout).unwrap());
    assert_eq!(
        ours.atanh.get("g_zero"),
        Some(&None),
        "all-zero gene must be NA (excluded), got {:?}",
        ours.atanh.get("g_zero")
    );
    let expected = std::fs::read_to_string(golden("degenerate_result.expected.tsv")).unwrap();
    assert_close(
        &ours,
        &parse(&expected),
        "duplicateCorrelation (degenerate)",
    );
}

/// Non-finite input cells (`inf`/`nan`, which Rust's f64 parser silently
/// accepts) must be rejected loudly rather than poisoning the REML objective.
#[test]
fn nonfinite_cell_is_rejected() {
    let dir = std::env::temp_dir().join(format!("dupcor_inf_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let expr = dir.join("expr.tsv");
    std::fs::write(
        &expr,
        "gene\ts001\ts002\ts003\ts004\ts005\ts006\n\
         g1\t1\t2\tinf\t4\t5\t6\n\
         g2\t2\t1\t3\t2\t4\t3\n",
    )
    .unwrap();
    let design = dir.join("design.tsv");
    std::fs::write(
        &design,
        "sample\tIntercept\tGroup\n\
         s001\t1\t0\ns002\t1\t1\ns003\t1\t0\ns004\t1\t1\ns005\t1\t0\ns006\t1\t1\n",
    )
    .unwrap();
    let block = dir.join("blocks.tsv");
    std::fs::write(
        &block,
        "sample\tblock\n\
         s001\tb1\ns002\tb1\ns003\tb1\ns004\tb2\ns005\tb2\ns006\tb2\n",
    )
    .unwrap();
    let out = Command::new(ours())
        .arg(&expr)
        .args(["--design", design.to_str().unwrap()])
        .args(["--block", block.to_str().unwrap()])
        .output()
        .unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    assert!(!out.status.success(), "must fail on a non-finite cell");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("non-finite"),
        "error should name the non-finite cell, got: {err}"
    );
}

fn rscript() -> Option<String> {
    let conda = format!(
        "{}/miniconda3/envs/r-bioc/bin/Rscript",
        std::env::var("HOME").unwrap_or_default()
    );
    for cand in [conda.as_str(), "Rscript"] {
        let ok = Command::new(cand)
            .args([
                "-e",
                "if(!requireNamespace('limma',quietly=TRUE)) quit(status=1)",
            ])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(cand.to_string());
        }
    }
    None
}

#[test]
fn live_r_consensus() {
    let Some(rs) = rscript() else {
        eprintln!("SKIP live_r_consensus: no Rscript with limma found");
        return;
    };
    let scratch = std::env::temp_dir();
    let r_out = scratch.join(format!("dupcor_r_{}.tsv", std::process::id()));
    let oracle = Command::new(&rs)
        .arg(manifest("tests/dupcor_oracle.R"))
        .arg(golden("expr.tsv"))
        .arg(golden("design.tsv"))
        .arg(golden("blocks.tsv"))
        .arg(&r_out)
        .output()
        .unwrap();
    assert!(
        oracle.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&oracle.stderr)
    );
    let ours_out = run_ours();
    let r_ref = std::fs::read_to_string(&r_out).unwrap();
    let _ = std::fs::remove_file(&r_out);
    assert_close(
        &parse(&ours_out),
        &parse(&r_ref),
        "duplicateCorrelation (live R)",
    );
}
