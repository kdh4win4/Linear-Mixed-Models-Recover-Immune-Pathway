# Linear-Mixed-Models-Recover-Immune-Pathway

Repository for Bioinformatics(Oxford) journal data availability.
"Linear Mixed Models Recover Immune Pathway Signatures Masked by Inter-individual Variance in Sorted γδ T Cell RNA-seq: IFN-α and B Cell Activation as Correlates of HIV Neutralization Breadth"

Dohoon Kim  |
PromptGenix LLC  |
Corresponding author: dkim@promptgenix.org


## Overview

This repository contains the R code used for the unified analysis pipeline and figure generation scripts associated with the manuscript on recovery of immune pathway signals using linear mixed models.

## Public datasets

- **Primary dataset:** GSE271442  
- **Validation datasets:** GSE124731, GSE24081  

All datasets are publicly available through NCBI GEO.

## Repository contents

- `RV217_UNIFIED_PIPELINE_GitHub_min.R`  
  Compact repository version of the final unified R pipeline and figure-generation code.

## Expected local input files

Place the following in `~/Downloads` before running the script:

- `final_metadata.csv`
- `GSE271442_Merged_with_Symbols.csv`

The script downloads or accesses GEO resources for the validation datasets as needed.

## Main analysis components

- metadata harmonization
- expression matrix loading and filtering
- limma differential expression
- linear mixed-effects modeling
- fGSEA and GSVA pathway analysis
- cross-cohort and validation analyses
- manuscript figure generation

## Output

The script saves analysis objects to:

- `RV217_UNIFIED_WORKSPACE.RData`

and writes manuscript figure files as PDF outputs.

## Data availability statement

The GSE271442 RNA-seq dataset analyzed in this study is publicly available in NCBI Gene Expression Omnibus:  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271442

Validation datasets GSE124731 and GSE24081 are also publicly available in GEO.

Original data processing methods and supplementary materials for the primary dataset publication are available here:  
https://pmc.ncbi.nlm.nih.gov/articles/PMC11805418/

## Citation

If you use this repository, please cite the associated manuscript and the original public datasets.
