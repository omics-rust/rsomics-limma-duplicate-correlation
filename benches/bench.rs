use std::hint::black_box;
use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_limma_duplicate_correlation::{Options, run};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bench_dupcor(c: &mut Criterion) {
    let expr = fixture("expr.tsv");
    let design = fixture("design.tsv");
    let block = fixture("blocks.tsv");
    if !expr.exists() {
        return;
    }
    c.bench_function("duplicate_correlation", |b| {
        b.iter(|| {
            let opts = Options {
                expr: &expr,
                design: &design,
                block: &block,
                threads: 1,
            };
            black_box(run(&opts).unwrap());
        })
    });
}

criterion_group!(benches, bench_dupcor);
criterion_main!(benches);
