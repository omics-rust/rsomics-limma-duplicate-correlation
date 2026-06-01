#!/usr/bin/env Rscript
# duplicateCorrelation oracle: read a log-expression matrix TSV (header = sample
# ids, col 1 = gene ids), a design matrix TSV (header = coefficient names, col 1
# = sample ids), and a block factor TSV (col 1 = sample id, col 2 = block
# label), then write the consensus correlation and per-gene atanh correlations.
#
# Usage: dupcor_oracle.R <expr.tsv> <design.tsv> <block.tsv> <out.tsv>
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
E <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
design <- as.matrix(read.delim(args[2], row.names = 1, check.names = FALSE))
blk <- read.delim(args[3], check.names = FALSE)
out_path <- args[4]

dc <- duplicateCorrelation(E, design, block = blk[[2]])

con <- file(out_path, "w")
g <- function(x) formatC(x, digits = 15, format = "g", flag = "")
writeLines(paste0("consensus.correlation\t", g(dc$consensus.correlation)), con)
writeLines("gene\tatanh.correlation", con)
ac <- dc$atanh.correlations
for (i in seq_len(nrow(E))) {
  v <- if (is.finite(ac[i])) g(ac[i]) else "NA"
  writeLines(paste(rownames(E)[i], v, sep = "\t"), con)
}
close(con)
