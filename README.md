# Linear-Mixed-Models-Recover-Immune-Pathway

R pipeline for linear mixed model–based recovery of immune pathway signals in γδ T cell RNA-seq.

**Manuscript Title**  
Linear Mixed Models Recover Immune Pathway Signatures Masked by Inter-individual Variance in Sorted γδ T Cell RNA-seq: IFN-α and B Cell Activation as Correlates of HIV Neutralization Breadth

**Author**  
Dohoon Kim  
PromptGenix LLC  
Corresponding author: dkim@promptgenix.org  

---

## Overview

This repository provides a unified R pipeline for transcriptomic and pathway-level analysis using linear mixed-effects models.  
The workflow integrates differential expression, gene set enrichment, and validation across independent cohorts.

---

## How to Run

Place the required input files in `~/Downloads`, then run:

    source("RV217_UNIFIED_PIPELINE_GitHub.R")

---

## Input Data

- final_metadata.csv  
  Sample metadata (reconstructed from GSE271442 supplementary data)

- GSE271442_Merged_with_Symbols.csv  
  TPM expression matrix (31 MB) with gene symbols  
  (unzip csv.zip before use)

---

## Public Datasets

- Primary dataset: GSE271442  
- Validation datasets: GSE124731, GSE24081  

All datasets are publicly available through NCBI GEO.

---

## Repository Contents

- RV217_UNIFIED_PIPELINE_GitHub.R  
  Compact reproducible version of the unified analysis pipeline and figure-generation code

---

## Analysis Workflow

- metadata harmonization  
- expression matrix processing and filtering  
- differential expression analysis (limma)  
- linear mixed-effects modeling (LMM)  
- pathway analysis (fGSEA, GSVA)  
- cross-cohort validation  
- manuscript figure generation  

---

## Output

The pipeline generates:

- RV217_UNIFIED_WORKSPACE.RData (analysis objects)  
- publication-ready figures (PDF format)  

---

## Data Availability

The RNA-seq dataset used in this study (GSE271442) is publicly available in NCBI GEO:  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271442  

Validation datasets (GSE124731 and GSE24081) are also publicly available in GEO.  

Original data processing details are available in:  
https://pmc.ncbi.nlm.nih.gov/articles/PMC11805418/  

---

## Citation

If you use this repository, please cite the associated manuscript and the original GEO datasets.
