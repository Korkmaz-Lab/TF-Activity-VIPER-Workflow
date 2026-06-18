# TF-Activity-VIPER-Workflow

# ARACNe + VIPER: Transcription Factor Activity Inference

A methods reference for inferring transcription factor (TF) activity from bulk RNA-seq using network-based approaches. This covers two complementary strategies: **ARACNe-AP** (data-driven regulon inference) and **DoRothEA** (curated regulon), both evaluated through **VIPER**.

---

## Overview

Gene expression levels alone do not tell which transcription factors are driving a phenotype. A TF can be active without changing its own mRNA level, and conversely a highly expressed TF may be inhibited post-translationally. VIPER addresses this by treating each TF's regulon (the set of genes it controls) as a molecular readout of TF activity — if a TF's downstream targets are coordinately up- or downregulated, the TF is inferred to be active or suppressed.

This workflow has two stages:

1. **Build the regulon** — either from a co-expression network (ARACNe) or from a curated database (DoRothEA)
2. **Run VIPER** — score TF activity per sample against a gene expression signature

---

## Workflow

```
RNA-seq expression matrix
        │
        ├──► ARACNe-AP ──► data-driven regulon (.txt)
        │                         │
        └──► DoRothEA  ──► curated regulon (R package)
                                  │
                             aracne2regulon()
                                  │
                              viper() / msviper()
                                  │
                    ┌─────────────┴─────────────┐
              per-sample TF               master regulator
              activity matrix (NES)       analysis (msviper NES)
```

---

## Stage 1 — Regulon Construction

### Option A: ARACNe-AP (data-driven)

ARACNe-AP infers a gene regulatory network from a large expression dataset using mutual information and data processing inequality (DPI) to retain only direct interactions.

**Requirements:** A large, homogeneous expression cohort (typically n ≥ 100 samples). TCGA disease-specific matrices are a common choice.

> *Example:* We used TCGA BRCA VST-normalized expression (1,226 samples) to infer separate tumor and normal networks. The tumor network contained ~267,000 interactions across 1,471 TFs; the normal network ~224,000 interactions across 1,460 TFs.

ARACNe-AP is run as a command-line Java tool (see [ARACNe-AP GitHub](https://github.com/califano-lab/ARACNe-AP)). The output is a tab-delimited network file.

Loading the network in R:

```r
library(viper)

regulon <- aracne2regulon(
  afile   = "tumor_network.txt",
  eset    = tcga_expr_matrix,   # used to estimate mode of regulation
  format  = "3col"
)
```

The expression matrix passed to `aracne2regulon()` is used solely to determine the sign (activation vs. repression) of each edge — it does not need to be the same dataset you will later run VIPER on.

### Option B: DoRothEA (curated)

DoRothEA provides a literature- and ChIP-seq-curated regulon for ~500 TFs, stratified by confidence (A–E). Confidence levels A and B are recommended for most analyses.

```r
library(dorothea)
library(viper)

regulon <- df2regulon(
  filter(dorothea_hs, confidence %in% c("A", "B"))
)
```

DoRothEA is smaller (A+B: ~118 TFs) but higher precision. It works well when the biological question involves well-characterised TFs and you do not have a tissue-matched expression dataset for ARACNe.

---

## Stage 2 — VIPER

### Per-sample TF activity

`viper()` scores TF activity for each sample in an expression matrix. The output is a TF × sample matrix of normalized enrichment scores (NES), where positive NES = TF is active, negative NES = TF is suppressed.

```r
tf_activity <- viper(
  eset    = query_expr_matrix,   # genes × samples
  regulon = regulon,
  nes.method = "none",
  method  = "scale",
  minsize = 4,
  verbose = FALSE
)
```

The resulting matrix can be used as features for classification, clustering, or differential activity testing, in place of or alongside raw gene expression.

### Master regulator analysis

`msviper()` tests which TFs best explain a gene expression signature — for example, the differential expression results from a 2D vs. 3D comparison, or a disease vs. normal contrast.

```r
# Build a gene expression signature (e.g. from DESeq2 results)
signature <- rowTtest(
  x = expr_condition_A,
  y = expr_condition_B
)

# Run msviper
mrs <- msviper(
  ges     = signature,
  regulon = regulon,
  minsize = 4,
  verbose = FALSE
)

summary(mrs, mrs = 20)   # top 20 master regulators
```

The output gives each TF a NES and p-value representing how coherently its regulon is enriched at the top or bottom of the expression signature.

---

## DoRothEA as Complementary Evidence

A practical use of running both is cross-validation: findings that appear in both the data-driven ARACNe network (built from an independent cohort) and the curated DoRothEA regulon are more likely to reflect real regulatory biology rather than network inference artefacts.

Differences between the two are also informative. ARACNe networks are tissue- and context-specific, so TFs that are strong master regulators in one tissue context may drop out when the network is inferred from a different one. DoRothEA, being context-agnostic, can recover these when the tissue match is poor.

---

## Output Interpretation

| Output | Interpretation |
|---|---|
| `viper()` NES > 0 | TF is inferred active in that sample |
| `viper()` NES < 0 | TF is inferred suppressed |
| `msviper()` NES > 0 | TF drives the upregulated end of the signature |
| `msviper()` NES < 0 | TF drives the downregulated end |
| Large \|NES\|, small p-value | Strong, coherent regulon enrichment |

NES from `viper()` can be used directly as feature matrices for machine learning classifiers (replacing or supplementing gene-level features). NES from `msviper()` identifies which TFs are most mechanistically relevant to the phenotype contrast being studied.

---

## Key Parameters

| Parameter | Function | Notes |
|---|---|---|
| `minsize` | Minimum regulon size to test | 4–15 typical; smaller = more TFs tested, more noise |
| `confidence` (DoRothEA) | A–E confidence tiers | A+B recommended; C adds coverage at lower precision |
| DPI tolerance (ARACNe) | Edge pruning stringency | Default 0 is conservative; relax for sparse networks |

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("viper", "dorothea"))

# ARACNe-AP: Java tool, separate installation
# https://github.com/califano-lab/ARACNe-AP
```

---

## References

- Alvarez MJ et al. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. *Nature Genetics*, 48(8), 838–847. [VIPER]
- Badia-i-Mompel P et al. (2023). Gene regulatory network inference from scRNA-seq with DoRothEA. *Bioinformatics*. [DoRothEA]
- Margolin AA et al. (2006). ARACNE: An algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7. [ARACNe]
