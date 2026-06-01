#!/usr/bin/env Rscript
# Generate a reproducible duplicate-correlation fixture: a log-expression matrix
# with a real block effect, the design, the block factor, and the limma oracle.
#
# Usage: mkfixture.R <ngenes> <nblocks> <repsperblock> <outdir>
suppressMessages(library(limma))
set.seed(20260601)

args <- commandArgs(trailingOnly = TRUE)
ngenes <- as.integer(args[1])
nblocks <- as.integer(args[2])
reps <- as.integer(args[3])
outdir <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

nsamp <- nblocks * reps
block <- rep(seq_len(nblocks), each = reps)
group <- rep(c(0, 1), length.out = nsamp)

rho_true <- 0.6
blockeff <- matrix(rnorm(ngenes * nblocks, sd = sqrt(rho_true)), ngenes, nblocks)
E <- matrix(0, ngenes, nsamp)
for (s in seq_len(nsamp)) {
  E[, s] <- 3 + 0.8 * group[s] + blockeff[, block[s]] +
    rnorm(ngenes, sd = sqrt(1 - rho_true))
}
rownames(E) <- sprintf("g%05d", seq_len(ngenes))
colnames(E) <- sprintf("s%03d", seq_len(nsamp))

ew <- cbind(gene = rownames(E), format(E, digits = 12, trim = TRUE, scientific = FALSE))
write.table(ew, file.path(outdir, "expr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

design <- data.frame(sample = colnames(E), Intercept = 1, Group = group)
write.table(design, file.path(outdir, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

bl <- data.frame(sample = colnames(E), block = sprintf("b%03d", block))
write.table(bl, file.path(outdir, "blocks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
