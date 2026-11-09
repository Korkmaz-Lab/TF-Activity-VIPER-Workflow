# TF-Activity-VIPER-Workflow

A short lab quick-reference for inferring transcription factor (TF) activity from bulk RNA-seq with **VIPER**, using regulons from either **ARACNe-AP** (data-driven) or **DoRothEA** (curated). It is not a full tutorial: for details see the original papers and the package documentation (references below). What this page tries to get right is the choice of options and the places where the defaults are not what we use.

---

## Overview

Expression of a TF does not tell whether the TF is active. VIPER treats a TF's regulon (the genes it controls) as a readout of its activity: if the targets are coordinately up- or down-regulated in a sample (or a contrast), the TF is inferred to be active or suppressed.

```
RNA-seq expression matrix
        │
        ├──► ARACNe-AP ──► network.txt ──► aracne2regulon() ──► pruneRegulon() ─┐
        │                                                                        ├──► regulon
        └──► DoRothEA / CollecTRI (curated) ──► df2regulon() ───────────────────┘
                                                                                 │
                                                          viper() / msviper()
                                                                 │
                                         ┌───────────────────────┴──────────────────┐
                                  per-sample activity                   master regulators (MR)
                                  matrix (NES)                          between two groups (NES)
```

**Which regulon?**

| | ARACNe-AP (Option A) | DoRothEA / CollecTRI (Option B) |
|---|---|---|
| Source | Inferred from your own cohort | Curated (literature, ChIP-seq, ...) |
| Needs | A homogeneous cohort, typically ≥ 100 samples, and an HPC for the bootstraps | Nothing but the R package |
| Strength | Context-specific; covers many regulators (TFs, co-TFs, ...) | Quick; well-characterised TFs; context-agnostic |
| Weakness | Depends on the training cohort; costly to run | Few targets per TF, so small `minsize`; not tissue-specific |

Running both is a useful cross-check; where they disagree, look at the tissue match of the network first.

## Quick start (runs without private data)

[`scripts/example_bcellviper.R`](scripts/example_bcellviper.R) runs the whole chain on the public `bcellViper` data: network file → `aracne2regulon` → `pruneRegulon` → `viper` → `msviper` + `shadow` → a target-shuffled null. The "network" in it is synthetic (written from the regulon shipped with `bcellViper`) because a real ARACNe-AP run needs hundreds of samples; swap in your own `network.txt` and expression matrix.

```bash
Rscript scripts/example_bcellviper.R example_out
```

(`viper` prints some `number of iterations=` lines while it fits the regulon; this is harmless.)

---

## Stage 1: regulon

### Option A: ARACNe-AP

ARACNe-AP infers a network from a large expression matrix with mutual information (MI) and the data processing inequality (DPI). Run it as a command-line Java tool: three steps (MI threshold, 100 bootstraps as a job array, consolidation). Notes on Java 8 / Singularity and SLURM templates are in [`docs/aracne_ap_hpc.md`](docs/aracne_ap_hpc.md).

> *Example (our earlier TCGA-BRCA analysis; VST-normalised, 1,113 tumour and 113 normal samples, separate networks):* the tumour network had 267,016 edges over 1,471 regulators and the normal network 223,545 edges over 1,460 regulators.

The output `network.txt` has a **header and four columns** (`Regulator`, `Target`, `MI`, `pvalue`). `aracne2regulon(format = "3col")` reads every line and takes the first three columns, so **the header and the p-value column must be removed first**; passing `network.txt` directly fails.

```r
library(viper)

net <- read.delim("network.txt")                                   # Regulator, Target, MI, pvalue
write.table(net[, c("Regulator", "Target", "MI")], "network_3col.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# expr: gene x sample matrix (gene symbols as rownames) = the matrix the network was built from
regulon <- aracne2regulon("network_3col.txt", expr, format = "3col", verbose = FALSE)
regulon <- pruneRegulon(regulon, cutoff = 50, adaptive = FALSE)    # see the parameter table
saveRDS(regulon, "regulon.rds")
```

`expr` is used to set the **mode of regulation** (sign of the regulator-target correlation in that matrix) and must contain all regulators and targets; use the network's own input matrix. The MI gives the weight.

### Option B: DoRothEA / CollecTRI (curated)

DoRothEA is a curated regulon with confidence levels A–E. In `dorothea` 1.22.0 it holds 1,333 TFs over all levels; levels A+B (recommended for most analyses) have 118 TFs and 6,408 interactions. These are small regulons (median 20.5 targets at A+B), so the default `minsize = 25` would drop 66 of the 118 TFs (52 pass); use a small `minsize` (5 here, all 118 pass). CollecTRI (via `decoupleR`) is a larger curated collection; the lab project uses `minsize = 5` for it.

```r
library(dorothea)
library(dplyr)          # for filter(); without it stats::filter is called and fails
library(viper)

regulon <- df2regulon(filter(dorothea_hs, confidence %in% c("A", "B")))
# later: viper(sig, regulon, ..., minsize = 5)
```

---

## Stage 2: VIPER

### Per-sample activity: `viper()`

Output: regulator × sample matrix of NES (positive = inferred active, negative = suppressed). It can be used for clustering, differential activity, or as features in a classifier.

```r
# sig: gene x sample signature; regulon: from Stage 1
activity <- viper(sig, regulon, method = "none", minsize = 25,
                  eset.filter = FALSE, verbose = FALSE)
```

- `method` sets how the input is standardised. `"none"` (default) assumes the signature is already standardised, e.g. a z-score of each sample against a reference cohort (what we do). For a raw expression matrix use `"scale"`.
- There is no `nes.method` argument (the earlier README used one and the call fails with `unused argument`). `nes` (logical, default `TRUE`) is the real one.
- `minsize`: default 25 (use this for ARACNe regulons); small for curated regulons (Option B).
- `eset.filter = FALSE` is **not** the default (`TRUE`): each sample is ranked over all measured genes instead of only the regulon genes.

### Master regulators between groups: `msviper()`

```r
# groupA, groupB: gene x sample matrices (e.g. columns of expr)
tt   <- rowTtest(groupA, groupB)
ges  <- (qnorm(tt$p.value / 2, lower.tail = FALSE) * sign(tt$statistic))[, 1]   # NAMED numeric vector
null <- ttestNull(groupA, groupB, per = 1000, repos = TRUE)                     # sample permutations

mrs <- msviper(ges, regulon, null, minsize = 25)
summary(mrs, mrs = 20)                       # top 20 regulators
mrs <- shadow(mrs, regulators = 25)          # flags MRs whose enrichment overlaps another MR's
```

`msviper()` needs a named numeric vector for `ges`, not the list returned by `rowTtest()`. Without the null model the p-values are not sample-permutation based and NES are strongly inflated when there are few samples. After `shadow()`, `summary()` returns a list (`$msviper.results`, `$Shadow.pairs`).

## Parameters

Package defaults come from viper 1.44.0. **Non-default** values are deliberate choices.

| Parameter | Value used here | Default? | Why |
|---|---|---|---|
| ARACNe-AP MI p-value | `1E-8` | yes | |
| ARACNe-AP DPI tolerance | `0` | yes | |
| ARACNe-AP consolidation | Bonferroni | yes | |
| ARACNe-AP bootstraps | 100 (seeds 1-100) | no (common practice) | stable consolidated network |
| `pruneRegulon(cutoff)` | 50 | yes | |
| `pruneRegulon(adaptive)` | `FALSE` | **no** (default `TRUE`) | keeps exactly the 50 strongest targets by MI; `TRUE` keeps more. Check sensitivity with `cutoff = 100`. |
| `viper(minsize)` | 25 (ARACNe); 5 (DoRothEA / CollecTRI) | 25 | small curated regulons |
| `viper(eset.filter)` | `FALSE` | **no** (default `TRUE`) | rank each sample over all measured genes |
| `viper(method)` | `"none"` with a pre-standardised signature | yes | use `"scale"` for raw expression |
| `msviper` null | 1000 sample permutations | no (null optional) | permutation p-values |

## Controls

Activity scores are indirect. Before interpreting a master-regulator list:

1. **Positive controls**, written down before looking at results (e.g. in breast: ESR1, FOXA1, GATA3 higher in luminal than basal). If they do not separate in the expected direction, check network, regulon and signature and do not interpret the list.
2. **Target-shuffled null:** permute the gene labels of the regulon targets (one permutation for all regulons, fixed seed, repeated e.g. 100 times) and score again; real activities should be clearly more structured than shuffled ones. A short version is in the example script.
3. **Sensitivity:** `pruneRegulon(cutoff = 50)` versus `100`; and `eset.filter` `FALSE` versus `TRUE` (correlation of activities per regulator).
4. **Baseline:** compare with a pre-built or curated regulon scored on the same signature, and with the regulator's own mRNA level.

## Files

Suggested names: `network.txt` (consolidated ARACNe-AP output), `regulon.rds` (pruned regulon), `activity.rds` (regulator × sample NES), `mr_<comparison>.tsv` (all tested regulators: NES, p, FDR, regulon size). Keep the parameters you used (and the `sessionInfo()`) next to the results.

## Dependencies

Tested with R 4.5.2 and Bioconductor 3.22: `viper` 1.44.0, `bcellViper` 1.46.0, `dorothea` 1.22.0, `decoupleR` 2.16.0, `dplyr` 1.2.1. (Bioconductor 3.23, R 4.6, ships `viper` 1.46.0 and `bcellViper` 1.48.0; not tested here.)

```r
install.packages("BiocManager")
BiocManager::install(c("viper", "bcellViper", "dorothea", "dplyr"))   # decoupleR only for CollecTRI
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")           # record versions with the results
```

ARACNe-AP is a separate Java tool: <https://github.com/califano-lab/ARACNe-AP> (Java 8; see [`docs/aracne_ap_hpc.md`](docs/aracne_ap_hpc.md)).

## References

- Alvarez MJ et al. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. *Nat Genet* 48(8):838-847. doi:10.1038/ng.3593 [VIPER]
- Lachmann A et al. (2016). ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. *Bioinformatics* 32(14):2233-2235. doi:10.1093/bioinformatics/btw216 [ARACNe-AP]
- Margolin AA et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics* 7(Suppl 1):S7. doi:10.1186/1471-2105-7-S1-S7 [ARACNE]
- Garcia-Alonso L et al. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Res* 29(8):1363-1375. doi:10.1101/gr.240663.118 [DoRothEA]
- Badia-I-Mompel P et al. (2022). decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinform Adv* 2(1):vbac016. doi:10.1093/bioadv/vbac016 [decoupleR]
- Müller-Dott S et al. (2023). Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. *Nucleic Acids Res* 51(20):10934-10949. doi:10.1093/nar/gkad841 [CollecTRI]
- Ding H et al. (2018). Quantitative assessment of protein activity in orphan tissues and single cells using the metaVIPER algorithm. *Nat Commun* 9:1471. doi:10.1038/s41467-018-03843-3 [metaVIPER]
