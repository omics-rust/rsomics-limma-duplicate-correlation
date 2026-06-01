# rsomics-limma-duplicate-correlation

Estimate the consensus intra-block correlation of a log-expression matrix — the
duplicate-spot / technical-replicate correlation used by limma when fitting a
linear model with a single random-effect block factor. A single-binary Rust
reimplementation of limma's `duplicateCorrelation`.

Given a log-expression matrix, a design matrix, and a block factor, it fits a
compound-symmetry mixed model per gene by REML, estimates each gene's
intra-block correlation, and reports the consensus (the `tanh` of the trimmed
mean of the per-gene `atanh` correlations) plus the per-gene values.

## Usage

```
rsomics-limma-duplicate-correlation expr.tsv --design design.tsv --block blocks.tsv [-o result.tsv]
```

- `expr.tsv` — header row of sample ids, first column gene ids, log-expression values.
- `--design` — header row of coefficient names, first column sample ids (the model matrix).
- `--block` — first column sample ids, second column block label; samples sharing a label form one block.
- `-o` — result TSV destination; `-` (default) is stdout.

Output: a `consensus.correlation` line followed by a `gene` / `atanh.correlation`
table (non-estimable genes report `NA` and are excluded from the consensus).

```
rsomics-limma-duplicate-correlation E.tsv --design design.tsv --block blocks.tsv -o cor.tsv
```

## Origin

This crate is an independent Rust reimplementation of limma's
`duplicateCorrelation` based on:

- The published method: Smyth, G.K., Michaud, J. and Scott, H.S. (2005), "Use of
  within-array replicate spots for assessing differential expression in
  microarray experiments", Bioinformatics 21(9):2067-2075,
  doi:10.1093/bioinformatics/bti270 — the per-gene REML estimate of the
  compound-symmetry intra-block correlation, and the `atanh`-mean consensus.
- Black-box behaviour testing against the limma binary via an R oracle
  (`Rscript` + limma), diffed in `tests/compat.rs`.

No source code from limma (GPL) was used as reference during implementation.
Output is value-exact against limma `duplicateCorrelation` (relative deviation
< 1e-6 on the consensus correlation and the per-gene `atanh` correlations).

License: MIT OR Apache-2.0.
Upstream credit: limma (https://bioconductor.org/packages/limma/), GPL (>=2).
