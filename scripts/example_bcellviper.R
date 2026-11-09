# End-to-end example on public data (Bioconductor package bcellViper): no private data needed.
#   Rscript scripts/example_bcellviper.R [outdir]
#
# Steps: ARACNe-AP-style network file -> aracne2regulon -> pruneRegulon -> viper -> msviper (+ shadow)
#        -> a small target-shuffled null as a sanity control.
#
# NOTE: a real ARACNe-AP run needs hundreds of samples and an HPC (see docs/aracne_ap_hpc.md).
# To keep this example runnable anywhere, the "network" below is SYNTHETIC: it is written from the
# regulon shipped with bcellViper, in the same file layout ARACNe-AP produces (network.txt with a
# header: Regulator, Target, MI, pvalue). Replace `net_file` and `expr` with your own to reuse the code.

suppressPackageStartupMessages({
  library(viper)
  library(bcellViper)
})

outdir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "example_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data(bcellViper)                      # gives `dset` (ExpressionSet, 6249 genes x 211 samples) and `regulon`
expr <- exprs(dset)                   # gene x sample matrix, gene symbols as rownames
group <- pData(dset)$description      # sample labels (CB, CC, N, ...)

# ---- 1. synthetic ARACNe-AP-style network.txt ---------------------------------------------------
net <- do.call(rbind, lapply(names(regulon), function(r) {
  data.frame(Regulator = r, Target = names(regulon[[r]]$tfmode),
             MI = regulon[[r]]$likelihood, pvalue = 1e-9, stringsAsFactors = FALSE)
}))
net_file <- file.path(outdir, "network.txt")
write.table(net, net_file, sep = "\t", quote = FALSE, row.names = FALSE)   # header + 4 columns, as ARACNe-AP

# ---- 2. regulon: aracne2regulon -> pruneRegulon ------------------------------------------------
# format = "3col" reads EVERY line and takes the first three columns, so the header and the p-value
# column must be removed first (a file with a header fails).
net_in <- read.delim(net_file, stringsAsFactors = FALSE)
net3 <- file.path(outdir, "network_3col.txt")
write.table(net_in[, c("Regulator", "Target", "MI")], net3, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# `expr` must be the matrix the network was built from: mode of regulation is the sign of the
# regulator-target correlation in it, the weight comes from MI.
reg <- aracne2regulon(net3, expr, format = "3col", verbose = FALSE)
# Lab choice (NOT the default): exactly the 50 strongest targets by MI. Default is adaptive = TRUE.
reg <- pruneRegulon(reg, cutoff = 50, adaptive = FALSE)
saveRDS(reg, file.path(outdir, "regulon.rds"))
cat("regulons:", length(reg), "| targets per regulon (median):",
    median(sapply(reg, function(z) length(z$tfmode))), "\n")

# ---- 3. per-sample activity: viper --------------------------------------------------------------
# method = "none" expects a signature that is already standardised. Here: each gene z-scored over
# all samples (in a real project use a reference cohort). For a raw expression matrix use "scale".
sig <- t(scale(t(expr)))
# Lab choices: minsize 25 (= viper default), eset.filter = FALSE (NOT the default TRUE).
act <- viper(sig, reg, method = "none", minsize = 25, eset.filter = FALSE, verbose = FALSE)
saveRDS(act, file.path(outdir, "activity.rds"))
cat("activity matrix:", paste(dim(act), collapse = " x "), "(regulators x samples)\n")

# ---- 4. master regulators between two groups: msviper -------------------------------------------
A <- expr[, group %in% c("CB", "CC")]     # germinal-centre B cells
B <- expr[, group == "N"]                 # naive B cells
tt <- rowTtest(A, B)
# msviper needs a NAMED NUMERIC VECTOR (not the list that rowTtest returns): signed z from the p-values
ges <- (qnorm(tt$p.value / 2, lower.tail = FALSE) * sign(tt$statistic))[, 1]
null <- ttestNull(A, B, per = 1000, repos = TRUE, verbose = FALSE)   # sample-permutation null
mrs <- msviper(ges, reg, null, minsize = 25, verbose = FALSE)
# all tested regulators (NES, p, FDR, regulon size), unfiltered
write.table(summary(mrs, mrs = length(mrs$regulon)), file.path(outdir, "mr_CBCC_vs_N.tsv"),
            sep = "\t", quote = FALSE)
# shadow analysis: flags pairs of MRs whose enrichment overlaps (see ?shadow for how to read a pair).
# With shadow() applied, summary() returns a list: $msviper.results and $Shadow.pairs.
shd <- summary(shadow(mrs, regulators = 25, verbose = FALSE), mrs = 10)
writeLines(shd$Shadow.pairs, file.path(outdir, "mr_CBCC_vs_N_shadow_pairs.txt"))
print(shd$msviper.results)

# ---- 5. control: target-shuffled null -----------------------------------------------------------
# One permutation of the gene labels applied to all targets; a regulator landing in its own regulon
# is dropped. Real activities should be far more structured than shuffled ones.
set.seed(1)
perm <- setNames(sample(rownames(expr)), rownames(expr))
shuf <- lapply(names(reg), function(r) {
  z <- reg[[r]]; nt <- unname(perm[names(z$tfmode)]); keep <- nt != r & !duplicated(nt)
  list(tfmode = setNames(z$tfmode[keep], nt[keep]), likelihood = z$likelihood[keep])
})
names(shuf) <- names(reg); class(shuf) <- "regulon"
act0 <- viper(sig, shuf, method = "none", minsize = 25, eset.filter = FALSE, verbose = FALSE)
cat(sprintf("median |NES|, real regulons: %.2f | target-shuffled: %.2f\n",
            median(abs(act)), median(abs(act0))))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("outputs in", outdir, "\n")
