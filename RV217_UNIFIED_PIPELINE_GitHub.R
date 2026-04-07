################################################################################
# RV217 γδ T Cell LMM-GSEA Unified Pipeline — FINAL
# Author:  Dohoon Kim | PromptGenix LLC
# Dataset: GSE271442 (Griffith et al. 2025, PLoS Pathogens)
# original file name : RV217_UNIFIED_PIPELINE_FINAL_v7_EN_v3.R
#
# Single-file end-to-end analysis pipeline.
#
# Required input files (~/Downloads):
#   1. final_metadata.csv                    — metadata (82 samples)
#   2. GSE271442_Merged_with_Symbols.csv     — expression data (with gene symbols)
#   3. GSE124731/                            — GEO supplementary (auto-download)
#   4. GSE24081                              — auto-download from GEO
#
# Datasets:
#   GSE271442  — primary: RV217 sorted γδ T cells, neutralization breadth
#   GSE124731  — methodological validation: innate T cells (Vδ1/NK), low-input
#   GSE24081   — biological validation: HIV Controller vs Progressor
#
# Object naming:
#   rv217_*     → GSE271442
#   gse124731_* → GSE124731
#   gse24081_*  → GSE24081
#   shared_*    → common (gene sets, targets, etc.)
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 0: Setup
# ══════════════════════════════════════════════════════════════════════════════
# NOTE: Do NOT use library(tidyverse) — causes conflicts with fgsea/patchwork/cowplot.
# Load packages individually; use dplyr:: namespace explicitly.

cat("══════════════════════════════════════════\n")
cat("SECTION 0: Environment Setup\n")
cat("══════════════════════════════════════════\n")

setwd("~/Downloads/")
set.seed(2025)

# CRAN
library(dplyr); library(tidyr); library(stringr); library(tibble)
library(readr); library(ggplot2); library(ggrepel); library(patchwork)
library(scales); library(RColorBrewer); library(lme4); library(lmerTest)
library(broom.mixed); library(jsonlite)

# Bioconductor
library(Biobase); library(GEOquery); library(limma); library(edgeR)
library(fgsea); library(msigdbr); library(GSVA); library(BiocParallel)
library(ComplexHeatmap); library(circlize); library(org.Hs.eg.db)
library(AnnotationDbi)

nCores <- max(1L, parallel::detectCores() - 2L)
cat(sprintf("[Setup] Packages loaded | %d cores\n\n", nCores))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Metadata + Expression (GSE271442)
# ══════════════════════════════════════════════════════════════════════════════
# GEO series matrix returns 0 expression rows for GSE271442.
# Load from local CSV files instead.

cat("══════════════════════════════════════════\n")
cat("SECTION 1: GSE271442 Metadata + TPM\n")
cat("══════════════════════════════════════════\n")

# 1-1. Metadata
rv217_meta <- read_csv("final_metadata.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    TimeClean = dplyr::case_when(
      Time %in% c("Pre-Infection", "Pre-Peak") ~ "Pre", TRUE ~ Time),
    TimeF    = factor(TimeClean, levels = c("Pre", "Peak", "SetPoint", "Chronic")),
    BroadF   = factor(Broad_vs_Nonbroad, levels = c("Non-Broad", "Broad")),
    SubjectF = factor(sampleID),
    # Region: first digit of sampleID — 1/2/3=Africa, 4=Thailand
    Region   = dplyr::case_when(
      grepl("^[123]", as.character(sampleID)) ~ "Africa",
      grepl("^4",     as.character(sampleID)) ~ "Thailand",
      TRUE ~ NA_character_),
    RegionF      = factor(Region, levels = c("Africa", "Thailand")),
    Sex          = ifelse(Region == "Africa", "Female", "TransgenderMale"),
    SexF         = factor(Sex, levels = c("Female", "TransgenderMale")),
    HIVsubtype_g = ifelse(Region == "Africa", "DiverseAfrican", "CRF01_AE"),
    HIVsubtypeF  = factor(HIVsubtype_g, levels = c("DiverseAfrican", "CRF01_AE"))
  )

cat("  Cohort structure:\n  Region × Broad (unique subjects):\n")
rv217_meta %>%
  dplyr::distinct(sampleID, Region, Broad_vs_Nonbroad) %>%
  dplyr::count(Region, Broad_vs_Nonbroad) %>%
  tidyr::pivot_wider(names_from = Broad_vs_Nonbroad, values_from = n) %>%
  print()

# 1-2. TPM expression matrix
# GSE271442_Merged_with_Symbols.csv: RowName, Symbol, + 82 sample columns
cat("\n[TPM] Loading expression data...\n")
rv217_merged    <- read.csv("GSE271442_Merged_with_Symbols.csv")
rv217_expr_cols <- rv217_meta$BNB   # 82 sample column names

rv217_tpm_mat <- rv217_merged %>%
  dplyr::select(RowName, Symbol, all_of(rv217_expr_cols)) %>%
  dplyr::filter(!is.na(Symbol), Symbol != "", Symbol != "NA") %>%
  dplyr::group_by(Symbol) %>%
  dplyr::slice_max(rowMeans(across(all_of(rv217_expr_cols))), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("Symbol") %>%
  dplyr::select(all_of(rv217_expr_cols)) %>%
  as.matrix()

# Low-expression filter: TPM > 1 in >= 30% of samples
rv217_keep  <- rowSums(rv217_tpm_mat > 1) >= (ncol(rv217_tpm_mat) * 0.30)
rv217_tpm_filt <- rv217_tpm_mat[rv217_keep, ]
rv217_ltpm  <- log2(rv217_tpm_filt + 1)

stopifnot("Column mismatch!" = all(colnames(rv217_ltpm) == rv217_meta$BNB))
cat(sprintf("  TPM matrix: %d genes × %d samples\n", nrow(rv217_ltpm), ncol(rv217_ltpm)))
cat(sprintf("  Value range: [%.2f, %.2f]\n", min(rv217_ltpm), max(rv217_ltpm)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: MSigDB Gene Sets
# ══════════════════════════════════════════════════════════════════════════════
# msigdbr v10+: use collection/subcollection (not deprecated category/subcategory)
# fGSEA: Hallmark + C2 Reactome/KEGG + C5 GO:BP/MF + C7 ImmuneSigDB
# GSVA:  C5 GO:BP (immune keyword filter) + C5 GO:MF (full)

cat("\n══════════════════════════════════════════\n")
cat("SECTION 2: MSigDB Gene Sets\n")
cat("══════════════════════════════════════════\n")

# fGSEA: all collections
shared_gs_hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>% split(.$gs_name) %>% lapply(`[[`, "gene_symbol")
shared_gs_c2_reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% split(.$gs_name) %>% lapply(`[[`, "gene_symbol")
shared_gs_c2_kegg <- msigdbr(species = "Homo sapiens",collection = "C2",  subcollection = "CP:KEGG_LEGACY") %>% split(.$gs_name) %>% lapply(`[[`, "gene_symbol")

shared_c5bp_all <- msigdbr(species = "Homo sapiens",  collection = "C5", subcollection = "GO:BP")
shared_c5mf_all <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")

shared_gs_c5bp_fgsea <- shared_c5bp_all %>% split(.$gs_name) %>% lapply(`[[`, "gene_symbol")
shared_gs_c5mf_fgsea <- shared_c5mf_all %>% split(.$gs_name) %>% lapply(`[[`, "gene_symbol")

shared_gs_c7 <- msigdbr(species = "Homo sapiens", collection = "C7",  subcollection = "IMMUNESIGDB") %>% split(.$gs_name) %>% lapply(`[[`, "gene_symbol")

shared_all_gene_sets <- c(shared_gs_hallmark, shared_gs_c2_reactome, shared_gs_c2_kegg, shared_gs_c5bp_fgsea,shared_gs_c5mf_fgsea, shared_gs_c7)
cat(sprintf("  fGSEA total: %d gene sets\n", length(shared_all_gene_sets)))

# GSVA: C5 GO:BP immune-filtered + GO:MF full
shared_immune_keywords <- c(
  "t_cell", "b_cell", "nk_cell", "lymphocyte", "leukocyte",
  "gamma_delta", "natural_killer", "tcr", "cytotoxic", "cytotoxicity",
  "killing", "apoptosis", "cytokine", "chemokine", "interleukin", "interferon",
  "immune_activation", "immune_response", "immune_effector",
  "lymphocyte_activation", "cell_activation", "adaptive_immune", "innate_immune",
  "antibody", "immunoglobulin", "humoral", "b_cell_activation",
  "b_cell_differentiation", "germinal_center", "somatic_hypermutation",
  "class_switch", "viral", "antiviral", "virus", "defense_response_to_virus",
  "response_to_virus", "response_to_interferon", "inflammatory", "inflammation",
  "nf_kb", "fc_receptor", "fcgr", "phagocytosis", "opsonization",
  "cell_cycle", "g2_m", "mitotic", "glycolysis", "oxidative_phosphorylation",
  "hiv", "retrovirus", "syncytium", "translation", "ribosom"
)
shared_immune_pattern <- paste(shared_immune_keywords, collapse = "|")

ms2list <- function(df) split(df$gene_symbol, df$gs_name)

shared_gs_c5bp_immune <- shared_c5bp_all %>%
  dplyr::mutate(gs_name_lower = tolower(gs_name)) %>%
  dplyr::filter(grepl(shared_immune_pattern, gs_name_lower)) %>%
  ms2list()

shared_gs_c5mf_gsva   <- ms2list(shared_c5mf_all)
shared_gs_c5_for_gsva <- c(shared_gs_c5bp_immune, shared_gs_c5mf_gsva)
cat(sprintf("  GSVA: %d sets (BP immune:%d + MF:%d)\n",
            length(shared_gs_c5_for_gsva),
            length(shared_gs_c5bp_immune), length(shared_gs_c5mf_gsva)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: limma — Naïve vs Region-corrected
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 3: limma DE (Naive vs Corrected)\n")
cat("══════════════════════════════════════════\n")

# 3-A: Naive (BroadF only, TimeF covariate)
rv217_design_naive <- model.matrix(~ TimeF + BroadF, data = rv217_meta)
rv217_fit_naive    <- lmFit(rv217_ltpm, rv217_design_naive) %>% eBayes()
rv217_naive_ranks  <- topTable(rv217_fit_naive, coef = "BroadFBroad",
                                number = Inf, sort.by = "t") %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::select(Symbol, t_naive = t, logFC_naive = logFC,
                P_naive = P.Value, adjP_naive = adj.P.Val) %>%
  tibble::as_tibble()

# 3-B: Region-corrected limma
rv217_design_corr <- model.matrix(~ TimeF + BroadF + RegionF, data = rv217_meta)
rv217_fit_corr    <- lmFit(rv217_ltpm, rv217_design_corr) %>% eBayes()
rv217_corrected_ranks <- topTable(rv217_fit_corr, coef = "BroadFBroad",
                                   number = Inf, sort.by = "t") %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::select(Symbol, t_corrected = t, logFC_corrected = logFC,
                P_corrected = P.Value, adjP_corrected = adj.P.Val) %>%
  tibble::as_tibble()

rv217_std_ranks <- rv217_corrected_ranks %>% dplyr::rename(t_std = t_corrected)

cat(sprintf("  Naive DE (adjP<0.05):     %d\n", sum(rv217_naive_ranks$adjP_naive < 0.05)))
cat(sprintf("  Corrected DE (adjP<0.05): %d\n", sum(rv217_corrected_ranks$adjP_corrected < 0.05)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: LMM with Region covariate
# ══════════════════════════════════════════════════════════════════════════════
# Model: log2(TPM+1) ~ TimeF + BroadF + RegionF + TimeF:BroadF + (1|SubjectF)

cat("\n══════════════════════════════════════════\n")
cat("SECTION 4: LMM gene-level DE\n")
cat("══════════════════════════════════════════\n")

rv217_fit_lmm_gene <- function(expr_vec, df) {
  df$y <- expr_vec
  tryCatch({
    m <- lmer(y ~ TimeF + BroadF + RegionF + TimeF:BroadF + (1 | SubjectF),
              data = df, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    tbl <- coef(summary(m))
    c(t_broad  = tbl["BroadFBroad", "t value"],
      t_region = if ("RegionFThailand" %in% rownames(tbl))
        tbl["RegionFThailand", "t value"] else NA_real_)
  }, error = function(e) c(t_broad = NA_real_, t_region = NA_real_))
}

cat(sprintf("  LMM on %d genes (%d cores)...\n", nrow(rv217_ltpm), nCores))
cat("  Estimated time: 15-25 min\n")

# rv217_t_list <- bplapply(
#   seq_len(nrow(rv217_ltpm)),
#   function(i) rv217_fit_lmm_gene(rv217_ltpm[i, ], rv217_meta),
#   BPPARAM = MulticoreParam(workers = nCores, progressbar = TRUE)
# )
# rv217_t_mat <- do.call(rbind, rv217_t_list)
# if failed above

#Fallback: stable single-thread loop (if bplapply freezes on macOS)
n_genes <- nrow(rv217_ltpm); rv217_t_list <- vector("list", n_genes); t0 <- Sys.time()
for (i in seq_len(n_genes)) {
  rv217_t_list[[i]] <- rv217_fit_lmm_gene(rv217_ltpm[i, ], rv217_meta)
  if (i %% 500 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  %d/%d (%.0f%%) | %.1f min | ~%.1f min left\n",
                i, n_genes, i/n_genes*100, elapsed, elapsed/i*(n_genes-i)))
  }
}
rv217_t_mat <- do.call(rbind, rv217_t_list)

rv217_lmm_ranks <- tibble::tibble(
  Symbol   = rownames(rv217_ltpm),
  t_lmm    = rv217_t_mat[, "t_broad"],
  t_region = rv217_t_mat[, "t_region"]
) %>%
  dplyr::filter(!is.na(t_lmm)) %>%
  dplyr::arrange(desc(t_lmm))

rv217_n_de_lmm <- sum(abs(rv217_lmm_ranks$t_lmm) > 2, na.rm = TRUE)
cat(sprintf("\n  LMM DE genes (|t|>2): %d\n", rv217_n_de_lmm))
cat(sprintf("  vs Naive: %.1fx | vs Corrected: %.1fx\n",
            rv217_n_de_lmm / max(sum(rv217_naive_ranks$adjP_naive < 0.05), 1),
            rv217_n_de_lmm / max(sum(rv217_corrected_ranks$adjP_corrected < 0.05), 1)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Timepoint-specific + Stratified + Interaction
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 5: Timepoint + Stratified analyses\n")
cat("══════════════════════════════════════════\n")

rv217_timepoints   <- c("Peak", "SetPoint", "Chronic")
rv217_lmm_tp_ranks <- list()

for (tp in rv217_timepoints) {
  tp_meta <- rv217_meta %>% dplyr::filter(TimeClean == tp)
  tp_ltpm <- rv217_ltpm[, tp_meta$BNB, drop = FALSE]
  d <- model.matrix(~ BroadF + RegionF, data = tp_meta)
  f <- lmFit(tp_ltpm, d) %>% eBayes()
  rv217_lmm_tp_ranks[[tp]] <- topTable(f, coef = "BroadFBroad",
                                         number = Inf, sort.by = "t") %>%
    tibble::rownames_to_column("Symbol") %>%
    dplyr::select(Symbol, t_lmm = t, logFC, P.Value, adj.P.Val) %>%
    tibble::as_tibble()
  cat(sprintf("  %s: adjP<0.05=%d\n", tp, sum(rv217_lmm_tp_ranks[[tp]]$adj.P.Val < 0.05)))
}

# Stratified (Africa-only, Thailand-only)
rv217_strat_ranks <- list()
for (region in c("Africa", "Thailand")) {
  reg_meta <- rv217_meta %>% dplyr::filter(Region == region)
  reg_ltpm <- rv217_ltpm[, reg_meta$BNB, drop = FALSE]
  d <- model.matrix(~ TimeF + BroadF, data = reg_meta)
  f <- lmFit(reg_ltpm, d) %>% eBayes()
  rv217_strat_ranks[[region]] <- topTable(f, coef = "BroadFBroad",
                                            number = Inf, sort.by = "t") %>%
    tibble::rownames_to_column("Symbol") %>%
    dplyr::select(Symbol, t_strat = t, logFC, P.Value, adj.P.Val) %>%
    tibble::as_tibble()
  cat(sprintf("  Stratified %s: adjP<0.05=%d\n", region,
              sum(rv217_strat_ranks[[region]]$adj.P.Val < 0.05)))
}

# Interaction (Region × Broad)
rv217_design_int <- model.matrix(~ TimeF + BroadF + RegionF + BroadF:RegionF, data = rv217_meta)
rv217_fit_int    <- lmFit(rv217_ltpm, rv217_design_int) %>% eBayes()
rv217_int_coef   <- grep("BroadFBroad:RegionFThailand", colnames(rv217_design_int), value = TRUE)
rv217_interaction_ranks <- if (length(rv217_int_coef) > 0) {
  topTable(rv217_fit_int, coef = rv217_int_coef, number = Inf, sort.by = "t") %>%
    tibble::rownames_to_column("Symbol") %>%
    dplyr::select(Symbol, t_int = t, logFC_int = logFC,
                  P_int = P.Value, adjP_int = adj.P.Val) %>%
    tibble::as_tibble()
} else NULL


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: fGSEA — All comparisons
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 6: fGSEA\n")
cat("══════════════════════════════════════════\n")

rv217_run_fgsea <- function(ranked_tbl, stat_col = "t_lmm", label = "", n_perm = 10000) {
  vec <- tibble::deframe(ranked_tbl %>% dplyr::select(Symbol, all_of(stat_col)))
  vec <- vec[!duplicated(names(vec))]
  if (length(vec) < 100) { cat(sprintf("  Skip %s — too few genes\n", label)); return(NULL) }
  res <- fgsea(pathways = shared_all_gene_sets, stats = vec,
               minSize = 15, maxSize = 500, nPermSimple = n_perm, nproc = nCores) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::mutate(
      analysis   = label,
      collection = dplyr::case_when(
        grepl("^HALLMARK_", pathway) ~ "Hallmark",
        grepl("^REACTOME_", pathway) ~ "Reactome",
        grepl("^KEGG_",     pathway) ~ "KEGG",
        grepl("^GOBP_",     pathway) ~ "GO:BP",
        grepl("^GOMF_",     pathway) ~ "GO:MF",
        TRUE                         ~ "ImmuneSigDB"))
  cat(sprintf("  %s: padj<0.05=%d | padj<0.25=%d\n", label,
              sum(res$padj < 0.05, na.rm = TRUE), sum(res$padj < 0.25, na.rm = TRUE)))
  res
}

cat("  6-A: Naive fGSEA\n")
rv217_gsea_naive <- rv217_run_fgsea(
  rv217_naive_ranks %>% dplyr::rename(t_lmm = t_naive), label = "naive")

cat("  6-B: Corrected fGSEA\n")
rv217_gsea_corrected <- rv217_run_fgsea(
  rv217_std_ranks, stat_col = "t_std", label = "corrected_limma")

cat("  6-C: LMM fGSEA\n")
rv217_gsea_lmm <- rv217_run_fgsea(rv217_lmm_ranks, stat_col = "t_lmm", label = "lmm_region")

# Timepoint-specific fGSEA
# (Peak, SetPoint, Chronic)
rv217_gsea_tp <- list()
for (tp in rv217_timepoints) {
  cat(sprintf("  6-D: %s fGSEA\n", tp))
  rv217_gsea_tp[[tp]] <- rv217_run_fgsea(
    rv217_lmm_tp_ranks[[tp]], stat_col = "t_lmm", label = paste0("tp_", tp))
}

# Stratified fGSEA
rv217_gsea_strat <- list()
for (region in c("Africa", "Thailand")) {
  cat(sprintf("  6-E: Stratified %s\n", region))
  rv217_gsea_strat[[region]] <- rv217_run_fgsea(
    rv217_strat_ranks[[region]] %>% dplyr::rename(t_lmm = t_strat),
    label = paste0("strat_", region))
}

# Cross-cohort concordance
rv217_gsea_cross <- NULL
if (!is.null(rv217_gsea_strat[["Africa"]]) && !is.null(rv217_gsea_strat[["Thailand"]])) {
  rv217_gsea_cross <- rv217_gsea_strat[["Africa"]] %>%
    dplyr::select(pathway, collection, NES_africa = NES, padj_africa = padj) %>%
    dplyr::inner_join(
      rv217_gsea_strat[["Thailand"]] %>% dplyr::select(pathway, NES_thailand = NES, padj_thailand = padj),
      by = "pathway") %>%
    dplyr::mutate(
      concordant    = sign(NES_africa) == sign(NES_thailand),
      both_sig      = padj_africa < 0.25 & padj_thailand < 0.25,
      africa_only   = padj_africa < 0.25 & padj_thailand >= 0.25,
      thailand_only = padj_thailand < 0.25 & padj_africa >= 0.25)
}

# NES 3-way comparison (SetPoint)
rv217_nes_comparison <- NULL
if (!is.null(rv217_gsea_tp[["SetPoint"]]) && !is.null(rv217_gsea_naive) && !is.null(rv217_gsea_corrected)) {
  rv217_nes_comparison <- rv217_gsea_tp[["SetPoint"]] %>%
    dplyr::select(pathway, collection, NES_lmm = NES, padj_lmm = padj) %>%
    dplyr::inner_join(rv217_gsea_naive %>% dplyr::select(pathway, NES_naive = NES, padj_naive = padj),
                      by = "pathway") %>%
    dplyr::inner_join(rv217_gsea_corrected %>% dplyr::select(pathway, NES_corrected = NES, padj_corrected = padj),
                      by = "pathway") %>%
    dplyr::mutate(
      NES_delta = NES_lmm - NES_naive,
      status = dplyr::case_when(
        abs(NES_naive) < 1.3 & abs(NES_lmm) >= 1.5 & padj_lmm < 0.1 & NES_lmm > 0 ~ "rescued_up",
        abs(NES_naive) < 1.3 & abs(NES_lmm) >= 1.5 & padj_lmm < 0.1 & NES_lmm < 0 ~ "rescued_down",
        padj_naive < 0.05 & padj_lmm < 0.05 ~ "confirmed",
        padj_naive >= 0.05 & padj_lmm < 0.05 ~ "new",
        TRUE ~ "unchanged")) %>%
    dplyr::arrange(desc(abs(NES_delta)))
  cat("\n  NES status:\n"); print(table(rv217_nes_comparison$status))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: GSVA — C5 GO:BP (immune-filtered) + GO:MF
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 7: GSVA\n")
cat("══════════════════════════════════════════\n")

rv217_gsva_param  <- gsvaParam(exprData = rv217_ltpm, geneSets = shared_gs_c5_for_gsva,
                                kcdf = "Gaussian", minSize = 10, maxSize = 500)
rv217_gsva_scores <- gsva(rv217_gsva_param, BPPARAM = MulticoreParam(workers = nCores, progressbar = FALSE))
cat(sprintf("  GSVA: %d pathways × %d samples\n", nrow(rv217_gsva_scores), ncol(rv217_gsva_scores)))

rv217_go_type_map <- tibble::tibble(
  pathway = rownames(rv217_gsva_scores),
  go_type = dplyr::case_when(
    grepl("^GOBP_", pathway) ~ "GO:BP",
    grepl("^GOMF_", pathway) ~ "GO:MF",
    TRUE ~ "Other"))

rv217_gsva_long <- rv217_gsva_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  tidyr::pivot_longer(-pathway, names_to = "BNB", values_to = "gsva_score") %>%
  dplyr::left_join(rv217_meta %>% dplyr::select(BNB, sampleID, TimeClean, Broad_vs_Nonbroad,
                                                 TimeF, BroadF, SubjectF, RegionF, Region),
                   by = "BNB") %>%
  dplyr::left_join(rv217_go_type_map, by = "pathway")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: Pathway-level LMM on GSVA scores
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 8: Pathway-level LMM on GSVA\n")
cat("══════════════════════════════════════════\n")

rv217_pathway_lmm <- rv217_gsva_long %>%
  dplyr::group_by(pathway, go_type) %>%
  dplyr::group_modify(function(df, key) {
    tryCatch({
      m <- lmer(gsva_score ~ TimeF * BroadF + RegionF + (1 | sampleID),
                data = df, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      tbl <- coef(summary(m))
      tibble::tibble(term = rownames(tbl), estimate = tbl[, "Estimate"],
                     se = tbl[, "Std. Error"], t_value = tbl[, "t value"],
                     p_value = tbl[, "Pr(>|t|)"])
    }, error = function(e) {
      tibble::tibble(term = NA_character_, estimate = NA_real_,
                     se = NA_real_, t_value = NA_real_, p_value = NA_real_)
    })
  }) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(term))

rv217_broad_effects <- rv217_pathway_lmm %>%
  dplyr::filter(grepl("BroadFBroad", term)) %>%
  dplyr::mutate(
    timepoint = dplyr::case_when(
      term == "BroadFBroad"              ~ "Overall",
      term == "TimeFPeak:BroadFBroad"     ~ "Peak",
      term == "TimeFSetPoint:BroadFBroad" ~ "SetPoint",
      term == "TimeFChronic:BroadFBroad"  ~ "Chronic",
      TRUE ~ term),
    pp_broad = pnorm(t_value),
    sig      = pp_broad >= 0.85 | pp_broad <= 0.15) %>%
  dplyr::arrange(desc(abs(t_value)))

rv217_region_effects <- rv217_pathway_lmm %>%
  dplyr::filter(grepl("RegionFThailand", term)) %>%
  dplyr::mutate(pp_thailand = pnorm(t_value),
                sig_region  = pp_thailand >= 0.85 | pp_thailand <= 0.15) %>%
  dplyr::select(pathway, go_type, estimate_region = estimate,
                t_region = t_value, pp_thailand, sig_region)

cat("  BroadF sig pathways (PP): ", sum(rv217_broad_effects$sig), "\n")
cat("  Region sig pathways (PP): ", sum(rv217_region_effects$sig_region), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9: GSE124731 — Methodological validation (Vδ1 vs NK)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 9: GSE124731 Methodological Validation\n")
cat("══════════════════════════════════════════\n")

# 9-1. Load GEO + supplementary files
gse124731_gse     <- getGEO("GSE124731", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse124731_supp_dir <- file.path(".", "GSE124731")
if (!dir.exists(gse124731_supp_dir)) getGEOSuppFiles("GSE124731", baseDir = ".")

gse124731_expr_raw <- read.table(
  gzfile(file.path(gse124731_supp_dir, "GSE124731_low_input_rnaseq_gene_normalized.txt.gz")),
  header = TRUE, sep = "\t", row.names = 1)
gse124731_meta_raw <- read.table(
  gzfile(file.path(gse124731_supp_dir, "GSE124731_low_input_rnaseq_meta_data.txt.gz")),
  header = TRUE, sep = "\t")
cat(sprintf("  Raw: %d genes × %d samples\n", nrow(gse124731_expr_raw), ncol(gse124731_expr_raw)))

# 9-2. Meta matching + filter
gse124731_expr_mat  <- as.matrix(gse124731_expr_raw)
gse124731_cn_clean  <- gsub("\\.", "_", colnames(gse124731_expr_mat))
gse124731_midx      <- match(gse124731_cn_clean, gse124731_meta_raw$sampleID)
gse124731_meta      <- gse124731_meta_raw[gse124731_midx, ]
gse124731_keep      <- gse124731_meta$cell_type %in% c("Vd1", "Vd2", "NK", "CD4", "CD8", "iNKT", "MAIT")
gse124731_expr_clean <- gse124731_expr_mat[, gse124731_keep]
gse124731_meta_clean <- gse124731_meta[gse124731_keep, ]
cat(sprintf("  Filtered: %d samples\n", ncol(gse124731_expr_clean)))
# All 79 samples are already Vd1/Vd2/NK/CD4/CD8/iNKT/MAIT — 79→79 expected.

# 9-3. Ensembl → Gene Symbol
cat("  Ensembl → Symbol...\n")
gse124731_symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = rownames(gse124731_expr_clean),
  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gse124731_has_sym  <- !is.na(gse124731_symbol_map)
gse124731_expr_sym <- gse124731_expr_clean[gse124731_has_sym, ]
rownames(gse124731_expr_sym) <- gse124731_symbol_map[gse124731_has_sym]
gse124731_expr_sym <- limma::avereps(gse124731_expr_sym, ID = rownames(gse124731_expr_sym))
cat(sprintf("  Mapped: %d → %d genes (%.0f%%)\n",
            sum(gse124731_has_sym), nrow(gse124731_expr_sym), mean(gse124731_has_sym)*100))

# 9-4. Vδ1 vs NK: limma + LMM
gse124731_vd1nk        <- gse124731_meta_clean$cell_type %in% c("Vd1", "NK")
gse124731_expr_vd1nk   <- gse124731_expr_sym[, gse124731_vd1nk]
gse124731_meta_vd1nk   <- gse124731_meta_clean[gse124731_vd1nk, ]
gse124731_meta_vd1nk$donor_id <- gse124731_meta_vd1nk$individual_id  # use existing column
cat(sprintf("  Vδ1: %d | NK: %d\n",
            sum(gse124731_meta_vd1nk$cell_type == "Vd1"),
            sum(gse124731_meta_vd1nk$cell_type == "NK")))

# limma
gse124731_grp    <- factor(gse124731_meta_vd1nk$cell_type, levels = c("NK", "Vd1"))
gse124731_design <- model.matrix(~ gse124731_grp)
gse124731_fit    <- eBayes(lmFit(gse124731_expr_vd1nk, gse124731_design))
gse124731_tt     <- topTable(gse124731_fit, coef = 2, number = Inf, sort.by = "none")
gse124731_n_de_limma  <- sum(gse124731_tt$adj.P.Val < 0.05, na.rm = TRUE)
cat(sprintf("  limma DE: %d\n", gse124731_n_de_limma))

gse124731_ranks_limma <- sort(
  setNames(gse124731_tt$t, rownames(gse124731_tt)) %>% .[!is.na(.)], decreasing = TRUE)

# LMM (donor random intercept)
cat("  Running LMM...\n")
gse124731_lmm_res <- lapply(rownames(gse124731_expr_vd1nk), function(g) {
  df <- data.frame(y = as.numeric(gse124731_expr_vd1nk[g, ]),
                   group = gse124731_meta_vd1nk$cell_type,
                   donor = gse124731_meta_vd1nk$donor_id)
  tryCatch({
    m <- lmer(y ~ group + (1|donor), data = df, REML = FALSE)
    ct <- summary(m)$coefficients
    if (nrow(ct) < 2) return(NULL)
    data.frame(gene = g, t_value = ct[2, "t value"], estimate = ct[2, "Estimate"])
  }, error = function(e) NULL)
})

gse124731_lmm_df    <- dplyr::bind_rows(gse124731_lmm_res)
gse124731_n_de_lmm  <- sum(abs(gse124731_lmm_df$t_value) > 2, na.rm = TRUE)
cat(sprintf("  LMM DE: %d (%.1fx vs limma)\n",
            gse124731_n_de_lmm, gse124731_n_de_lmm / max(gse124731_n_de_limma, 1)))

gse124731_lmm_clean <- gse124731_lmm_df %>%
  dplyr::filter(!is.na(gene)) %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(abs(t_value), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()
gse124731_ranks_lmm <- sort(
  setNames(gse124731_lmm_clean$t_value, gse124731_lmm_clean$gene) %>% .[!is.na(.)],
  decreasing = TRUE)

# 9-5. fGSEA on target pathways
shared_target_pathways <- c(
  "GOBP_B_CELL_ACTIVATION", "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
  "GOBP_CYTOPLASMIC_TRANSLATION", "GOBP_RIBOSOME_BIOGENESIS",
  "GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME", "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
  "GOBP_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "GOBP_DEFENSE_RESPONSE_TO_VIRUS",
  "GOBP_INNATE_IMMUNE_RESPONSE", "GOBP_HUMORAL_IMMUNE_RESPONSE",
  "GOBP_IMMUNOGLOBULIN_PRODUCTION", "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY", "GOMF_CYTOKINE_RECEPTOR_BINDING",
  "GOBP_CELL_KILLING", "GOBP_RRNA_PROCESSING",
  "GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS",
  "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA", "GOBP_REGULATION_OF_VIRAL_PROCESS",
  "GSE21063_WT_VS_NFATC1_KO_16H_ANTI_IGM_STIM_BCELL_DN"
)
shared_target_list <- shared_all_gene_sets[intersect(shared_target_pathways, names(shared_all_gene_sets))]

gse124731_fg_limma <- fgsea(shared_target_list, gse124731_ranks_limma, minSize = 5, maxSize = 500, nPermSimple = 2000)
gse124731_fg_lmm   <- fgsea(shared_target_list, gse124731_ranks_lmm,   minSize = 5, maxSize = 500, nPermSimple = 2000)

gse124731_comp <- gse124731_fg_limma %>%
  dplyr::select(pathway, NES_limma = NES, padj_limma = padj) %>%
  dplyr::left_join(gse124731_fg_lmm %>% dplyr::select(pathway, NES_lmm = NES, padj_lmm = padj),
                   by = "pathway") %>%
  dplyr::mutate(
    category = dplyr::case_when(
      padj_limma < 0.05 & padj_lmm < 0.05                          ~ "Both",
      padj_limma < 0.05 & (is.na(padj_lmm) | padj_lmm >= 0.05)    ~ "limma only",
      (is.na(padj_limma) | padj_limma >= 0.05) & padj_lmm < 0.05  ~ "LMM only",
      TRUE ~ "Neither"))

cat("\n  Pathway classification:\n"); print(table(gse124731_comp$category))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10: GSE24081 — Biological validation (Controller vs Progressor)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 10: GSE24081 Biological Validation\n")
cat("══════════════════════════════════════════\n")

gse24081_gse   <- getGEO("GSE24081", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse24081_eset  <- gse24081_gse[[1]]
gse24081_pdata <- pData(gse24081_eset)
gse24081_expr  <- exprs(gse24081_eset)
cat(sprintf("  Probes: %d | Samples: %d\n", nrow(gse24081_expr), ncol(gse24081_expr)))

# Group assignment
gse24081_pdata$hiv_group <- dplyr::case_when(
  grepl("controller|elite|LTNP|non.prog", gse24081_pdata$title, ignore.case = TRUE) ~ "Controller",
  grepl("progressor|viremic|chronic",     gse24081_pdata$title, ignore.case = TRUE) ~ "Progressor",
  grepl("uninfect|seroneg|healthy|negat", gse24081_pdata$title, ignore.case = TRUE) ~ "Control",
  TRUE ~ "Unknown")

# Fallback to disease state column if majority are Unknown
gse24081_ds_col <- grep("disease|state", colnames(gse24081_pdata), value = TRUE, ignore.case = TRUE)
if (sum(gse24081_pdata$hiv_group == "Unknown") > nrow(gse24081_pdata)/2 && length(gse24081_ds_col) > 0) {
  ds <- gse24081_pdata[[gse24081_ds_col[1]]]
  gse24081_pdata$hiv_group <- dplyr::case_when(
    grepl("controller|elite",    ds, ignore.case = TRUE) ~ "Controller",
    grepl("progressor|viremic",  ds, ignore.case = TRUE) ~ "Progressor",
    TRUE ~ "Unknown")
}
cat("  Groups:\n"); print(table(gse24081_pdata$hiv_group))

# Probe → gene collapse
gse24081_fdata  <- fData(gse24081_eset)
gse24081_sym_col <- grep("Gene.Symbol|gene_symbol|Symbol", colnames(gse24081_fdata),
                          value = TRUE, ignore.case = TRUE)[1]
gse24081_syms   <- gse24081_fdata[[gse24081_sym_col]]
gse24081_syms[gse24081_syms == "" | is.na(gse24081_syms)] <- NA
gse24081_has    <- !is.na(gse24081_syms)
gse24081_expr_gene <- limma::avereps(gse24081_expr[gse24081_has, ], ID = gse24081_syms[gse24081_has])
cat(sprintf("  Probes→Genes: %d → %d\n", nrow(gse24081_expr), nrow(gse24081_expr_gene)))

# DE
gse24081_grp1 <- "Controller"
gse24081_grp2 <- if ("Progressor" %in% gse24081_pdata$hiv_group) "Progressor" else "Control"
gse24081_keep  <- gse24081_pdata$hiv_group %in% c(gse24081_grp1, gse24081_grp2)
gse24081_group <- factor(gse24081_pdata$hiv_group[gse24081_keep], levels = c(gse24081_grp2, gse24081_grp1))
gse24081_design <- model.matrix(~ gse24081_group)
gse24081_fit    <- eBayes(lmFit(gse24081_expr_gene[, gse24081_keep], gse24081_design))
gse24081_tt     <- topTable(gse24081_fit, coef = 2, number = Inf, sort.by = "none")
cat(sprintf("  DE genes (adjP<0.05/0.1/0.2/0.3): %d / %d / %d / %d\n",
            sum(gse24081_tt$adj.P.Val < 0.05, na.rm=TRUE),
            sum(gse24081_tt$adj.P.Val < 0.10, na.rm=TRUE),
            sum(gse24081_tt$adj.P.Val < 0.20, na.rm=TRUE),
            sum(gse24081_tt$adj.P.Val < 0.30, na.rm=TRUE)))

gse24081_ranks <- sort(setNames(gse24081_tt$t, rownames(gse24081_tt)) %>% .[!is.na(.)], decreasing = TRUE)
gse24081_fg    <- fgsea(shared_target_list, gse24081_ranks, minSize = 5, maxSize = 500, nPermSimple = 2000)
cat(sprintf("  Target sig: %d / %d\n", sum(gse24081_fg$padj < 0.05, na.rm = TRUE), nrow(gse24081_fg)))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11: Manuscript numbers summary
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("MANUSCRIPT KEY NUMBERS\n")
cat("══════════════════════════════════════════\n")

cat("\n── GSE271442 Gene-level DE ──\n")
cat(sprintf("  Naive (adjP<0.05):     %d\n", sum(rv217_naive_ranks$adjP_naive < 0.05)))
cat(sprintf("  Corrected (adjP<0.05): %d\n", sum(rv217_corrected_ranks$adjP_corrected < 0.05)))
cat(sprintf("  LMM (|t|>2):           %d\n", rv217_n_de_lmm))

cat("\n── GSE271442 Pathway-level (padj<0.05) ──\n")
cat(sprintf("  Naive:     %d\n", sum(rv217_gsea_naive$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  Corrected: %d\n", sum(rv217_gsea_corrected$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  LMM:       %d\n", sum(rv217_gsea_lmm$padj < 0.05, na.rm = TRUE)))

cat("\n── GSE124731 (Vδ1 vs NK) ──\n")
cat(sprintf("  limma DE: %d | LMM DE: %d (%.1fx)\n",
            gse124731_n_de_limma, gse124731_n_de_lmm,
            gse124731_n_de_lmm / max(gse124731_n_de_limma, 1)))

cat("\n── Gene set counts ──\n")
cat("Hallmark:", length(shared_gs_hallmark), "\n")
cat("C2 Reactome:", length(shared_gs_c2_reactome), "\n")
cat("C2 KEGG:", length(shared_gs_c2_kegg), "\n")
cat("C5 GO:BP:", length(shared_gs_c5bp_fgsea), "\n")
cat("C5 GO:MF:", length(shared_gs_c5mf_fgsea), "\n")
cat("C7 ImmuneSigDB:", length(shared_gs_c7), "\n")
cat("Total fGSEA:", length(shared_all_gene_sets), "\n")
cat("GSVA GO:BP immune:", length(shared_gs_c5bp_immune), "\n")
cat("GSVA GO:MF:", length(shared_gs_c5mf_gsva), "\n")
cat("GSVA total:", length(shared_gs_c5_for_gsva), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12: Save workspace
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("SECTION 12: Save\n")
cat("══════════════════════════════════════════\n")

save(
  rv217_meta, rv217_ltpm,
  rv217_naive_ranks, rv217_corrected_ranks, rv217_std_ranks,
  rv217_lmm_ranks, rv217_lmm_tp_ranks, rv217_strat_ranks,
  rv217_interaction_ranks,
  rv217_gsea_naive, rv217_gsea_corrected, rv217_gsea_lmm,
  rv217_gsea_tp, rv217_gsea_strat, rv217_gsea_cross,
  rv217_nes_comparison,
  shared_gs_c5_for_gsva, shared_gs_c5bp_immune, shared_gs_c5mf_gsva,
  rv217_gsva_scores, rv217_gsva_long,
  rv217_broad_effects, rv217_region_effects, rv217_pathway_lmm,
  rv217_go_type_map,
  gse124731_expr_sym, gse124731_meta_clean,
  gse124731_tt, gse124731_lmm_df,
  gse124731_fg_limma, gse124731_fg_lmm, gse124731_comp,
  gse124731_n_de_limma, gse124731_n_de_lmm,
  gse24081_tt, gse24081_fg,
  shared_all_gene_sets, shared_target_list,
  file = "RV217_UNIFIED_WORKSPACE.RData"
)
cat("  Saved: RV217_UNIFIED_WORKSPACE.RData\n")
cat("\n══════════════════════════════════════════\n")
cat("UNIFIED PIPELINE COMPLETED!\n")
cat("══════════════════════════════════════════\n")


################################################################################
# FIGURES — uses rv217_* objects (run pipeline above, or load workspace)
################################################################################

setwd("~/Downloads/")

# Color palettes + helper
status_colors  <- c(rescued_up = "#7F77DD", rescued_down = "#D4537E",
                    confirmed  = "#1D9E75", new          = "#BA7517",
                    unchanged  = "#C8C7C0")
go_type_colors <- c("GO:BP" = "#1D9E75", "GO:MF" = "#BA7517")
region_colors  <- c("Africa" = "#D85A30", "Thailand" = "#7F77DD")
broad_colors   <- c("Non-Broad" = "#B4B2A9", "Broad" = "#534AB7")

clean_go <- function(x, width = 45) {
  x %>% gsub("^GOBP_|^GOMF_", "", .) %>% tolower() %>%
    gsub("_", " ", .) %>% str_wrap(width = width)
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Cross-cohort GSEA — Africa vs Thailand NES scatter
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(rv217_gsea_cross)) {
  p1_dat <- rv217_gsea_cross %>%
    dplyr::filter(padj_africa < 0.25 | padj_thailand < 0.25) %>%
    dplyr::filter(collection %in% c("GO:BP","GO:MF","Hallmark","Reactome")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        both_sig & concordant ~
          gsub("^HALLMARK_|^REACTOME_|^GOBP_|^GOMF_", "", pathway) %>%
          tolower() %>% gsub("_"," ",.) %>% stringr::str_wrap(width = 25),
        TRUE ~ ""),
      both_sig = as.logical(both_sig))

  p1 <- ggplot(p1_dat, aes(NES_africa, NES_thailand, color = collection, size = both_sig)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="gray60", linewidth=0.6) +
    geom_hline(yintercept=0, linewidth=0.4, color="gray70") +
    geom_vline(xintercept=0, linewidth=0.4, color="gray70") +
    geom_point(aes(alpha = both_sig)) +
    ggrepel::geom_text_repel(aes(label=label), size=3.8, lineheight=0.85,
                              max.overlaps=Inf, segment.size=0.3,
                              segment.color="gray60", box.padding=0.6) +
    scale_size_manual(values  = c("FALSE"=1.5, "TRUE"=3.5), guide="none") +
    scale_alpha_manual(values = c("FALSE"=0.35,"TRUE"=0.9), guide="none") +
    scale_color_brewer(palette="Set2") +
    facet_wrap(~ collection, ncol=2) +
    labs(title    = "Fig 1: Cross-cohort GSEA replication",
         subtitle = paste0("Africa (Female, A1/C/D)  vs.  Thailand (Transgender, CRF01_AE)\n",
                           "Diagonal = perfect agreement   |   Large points = sig in BOTH cohorts"),
         x = "NES \u2014 Africa only", y = "NES \u2014 Thailand only",
         color   = "Collection",
         caption = "Broad vs. Non-Broad | GSE271442 | Dohoon Kim \u00b7 PromptGenix LLC") +
    theme_minimal(base_size=15) +
    theme(plot.title=element_text(size=19,face="bold"),
          plot.subtitle=element_text(size=13,color="gray40"),
          plot.caption=element_text(size=11,color="gray50"),
          strip.text=element_text(size=13,face="bold"),
          axis.title=element_text(size=14), axis.text=element_text(size=12),
          legend.title=element_text(size=13), legend.text=element_text(size=12),
          legend.position="bottom")

  ggsave("Fig1_CrossCohort_GSEA_v1.pdf", p1, width=16, height=13)
  message("Saved: Fig1_CrossCohort_GSEA_v1.pdf")
}

toupper("atp synthesis")
p1_dat$pathway[grepl("ATP_SYNTHESIS", p1_dat$pathway)]
vv <- as.data.frame(p1_dat[grepl("ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT", p1_dat$pathway),])
tolower(vv$pathway)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Three-way NES comparison — Naive → Corrected → LMM (SetPoint)
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(rv217_nes_comparison)) {

  # Abbreviation fix for pathway labels
  fix_pathway_label <- function(x) {
    x %>% stringr::str_to_title() %>%
      stringr::str_replace_all("\\bRna\\b","RNA")   %>% stringr::str_replace_all("\\bDna\\b","DNA")   %>%
      stringr::str_replace_all("\\bMrna\\b","mRNA") %>% stringr::str_replace_all("\\bTrna\\b","tRNA") %>%
      stringr::str_replace_all("\\bRrna\\b","rRNA") %>% stringr::str_replace_all("\\bSnrna\\b","snRNA") %>%
      stringr::str_replace_all("\\bSnoRna\\b","snoRNA") %>% stringr::str_replace_all("\\bNcrna\\b","ncRNA") %>%
      stringr::str_replace_all("\\bFc Gamma\\b","FcγR")  %>% stringr::str_replace_all("\\bFc-Gamma\\b","FcγR") %>%
      stringr::str_replace_all("\\bFcgr\\b","FcγR")      %>% stringr::str_replace_all("\\bFc Gamma Receptor\\b","FcγR") %>%
      stringr::str_replace_all("\\bMhc\\b","MHC")   %>% stringr::str_replace_all("\\bNk\\b","NK")     %>%
      stringr::str_replace_all("\\bNf-?Kb\\b","NF-κB") %>% stringr::str_replace_all("\\bNfkb\\b","NF-κB") %>%
      stringr::str_replace_all("\\bTcr\\b","TCR")   %>% stringr::str_replace_all("\\bBcr\\b","BCR")   %>%
      stringr::str_replace_all("\\bIl\\b","IL")     %>% stringr::str_replace_all("\\bIfn\\b","IFN")   %>%
      stringr::str_replace_all("\\bTnf\\b","TNF")   %>% stringr::str_replace_all("\\bTgf\\b","TGF")   %>%
      stringr::str_replace_all("\\bGtp\\b","GTP")   %>% stringr::str_replace_all("\\bAtp\\b","ATP")   %>%
      stringr::str_replace_all("\\bAdp\\b","ADP")   %>% stringr::str_replace_all("\\bNad\\b","NAD")   %>%
      stringr::str_replace_all("\\bCoa\\b","CoA")   %>% stringr::str_replace_all("\\bAcyl-Coa\\b","Acyl-CoA") %>%
      stringr::str_replace_all("\\bMtoc\\b","MTOC") %>% stringr::str_replace_all("\\bAdcc\\b","ADCC") %>%
      stringr::str_replace_all("\\bHiv\\b","HIV")   %>% stringr::str_replace_all("\\bAids\\b","AIDS") %>%
      stringr::str_replace_all("\\bEr\\b","ER")
  }

  p2_dat <- rv217_nes_comparison %>%
    dplyr::filter(collection %in% c("GO:BP","GO:MF","Hallmark","Reactome"),
                  status != "unchanged") %>%
    dplyr::slice_max(abs(NES_lmm), n=25) %>%
    dplyr::mutate(
      label = clean_go(pathway, width=40),
      label = fix_pathway_label(label),
      label = reorder(label, NES_lmm))

  p2_long <- p2_dat %>%
    tidyr::pivot_longer(c(NES_naive, NES_corrected, NES_lmm),
                        names_to="method", values_to="NES") %>%
    dplyr::mutate(
      method = factor(method,
                      levels = c("NES_naive","NES_corrected","NES_lmm"),
                      labels = c("① Naïve\n(no correction)",
                                 "② Region-corrected\n(limma + RegionF)",
                                 "③ LMM + Region\n(best estimate)")))

  x_lim <- ceiling(max(abs(p2_long$NES), na.rm=TRUE) * 1.15 * 2) / 2

  collection_label <- p2_dat %>%
    dplyr::select(label, collection) %>%
    dplyr::mutate(method = factor("③ LMM + Region\n(best estimate)"), NES = x_lim * 0.97)

  p2 <- ggplot(p2_long, aes(x=NES, y=label, color=status)) +
    annotate("rect", xmin=-1.5, xmax=1.5, ymin=-Inf, ymax=Inf, fill="#EEEEEE", alpha=0.35) +
    geom_segment(aes(x=0, xend=NES, y=label, yend=label), linewidth=0.45, color="gray75") +
    geom_point(size=3.5, alpha=0.9) +
    geom_vline(xintercept=0, linewidth=0.6, color="gray30") +
    geom_vline(xintercept=c(-1.5,1.5), linetype="dashed", color="gray50", linewidth=0.5) +
    geom_text(data=collection_label, aes(x=NES, y=label, label=collection),
              hjust=1, size=3.0, color="gray50", fontface="italic", inherit.aes=FALSE) +
    scale_color_manual(values=status_colors, name="Pathway Status",
                       labels=c(confirmed="Confirmed (all methods)", new="New (LMM-only)",
                                rescued_up="Rescued ↑ (by LMM)", rescued_down="Rescued ↓ (by LMM)")) +
    scale_x_continuous(limits=c(-x_lim, x_lim), breaks=seq(-floor(x_lim), floor(x_lim), by=1)) +
    facet_wrap(~ method, nrow=1, scales="fixed") +
    labs(title    = "Fig 2: Pathway NES — Naïve  →  Region-corrected  →  LMM+Region",
         subtitle = "SetPoint VL (day 42) | Broad vs. Non-Broad | Top 25 by |NES_LMM| | C5 GO + Hallmark + Reactome",
         caption  = "GSE271442 | Dohoon Kim · PromptGenix LLC",
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_minimal(base_size=13) +
    theme(plot.margin=margin(8,6,6,6,"mm"), panel.spacing=unit(0.7,"cm"),
          strip.text=element_text(size=13,face="bold",color="gray20",margin=margin(t=3,b=3)),
          strip.background=element_rect(fill="gray95",color=NA),
          panel.grid.major.y=element_line(color="gray92",linewidth=0.3),
          panel.grid.major.x=element_line(color="gray88",linewidth=0.3),
          axis.text.y=element_text(size=10.5), axis.text.x=element_text(size=11),
          axis.title.x=element_text(size=12),
          plot.title=element_text(size=17,face="bold",margin=margin(b=4)),
          plot.subtitle=element_text(size=12,color="gray40",margin=margin(b=6)),
          plot.caption=element_text(size=10,color="gray55"),
          legend.position="bottom", legend.title=element_text(size=12,face="bold"),
          legend.text=element_text(size=11), legend.key.size=unit(0.5,"cm"),
          legend.margin=margin(t=2,unit="mm"))

  ggsave("Fig2_NES_threeway_final.pdf", p2, width=16, height=11)
  message("Saved: Fig2_NES_threeway_final.pdf")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: C5 GO GSVA trajectories — 4 groups
# ══════════════════════════════════════════════════════════════════════════════

rv217_top_go_setpoint <- rv217_broad_effects %>%
  dplyr::filter(timepoint == "SetPoint", sig) %>%
  dplyr::arrange(desc(abs(t_value))) %>%
  head(8) %>% dplyr::pull(pathway)

if (length(rv217_top_go_setpoint) > 0) {
  p3_dat <- rv217_gsva_long %>%
    dplyr::filter(pathway %in% rv217_top_go_setpoint, !is.na(TimeF)) %>%
    dplyr::mutate(pw_label = clean_go(pathway, width=26),
                  Group    = paste0(Broad_vs_Nonbroad, "\n", Region)) %>%
    dplyr::group_by(pw_label, TimeF, Broad_vs_Nonbroad, Region, Group) %>%
    dplyr::summarise(mean_score = mean(gsva_score, na.rm=TRUE),
                     se_score   = sd(gsva_score, na.rm=TRUE) / sqrt(dplyr::n()),
                     .groups="drop")

  p3 <- ggplot(p3_dat, aes(TimeF, mean_score, color=Broad_vs_Nonbroad, linetype=Region,
                             group=interaction(Broad_vs_Nonbroad, Region))) +
    geom_line(linewidth=1.2) + geom_point(size=3.2) +
    geom_errorbar(aes(ymin=mean_score-se_score, ymax=mean_score+se_score),
                  width=0.25, alpha=0.55, linewidth=0.7) +
    facet_wrap(~ pw_label, scales="free_y", ncol=2) +
    scale_color_manual(values=broad_colors) +
    scale_linetype_manual(values=c("Africa"="solid","Thailand"="dashed")) +
    scale_x_discrete(labels=c("Pre"="Pre\n(-185d)","Peak"="Peak\n(18d)",
                               "SetPoint"="Set\n(42d)","Chronic"="CHI\n(980d)")) +
    labs(title    = "Fig 3: C5 GO GSVA trajectories \u2014 4-group comparison",
         subtitle = paste0("Top GO pathways at SetPoint VL (region-corrected PP)\n",
                           "Solid = Africa (Female / A1-D)   |   Dashed = Thailand (TG / CRF01_AE)"),
         x=NULL, y="Mean GSVA score (\u00b1SE)",
         color="Neutralization", linetype="Region",
         caption="GSE271442 | Dohoon Kim \u00b7 PromptGenix LLC") +
    theme_minimal(base_size=14) +
    theme(strip.text=element_text(size=12,face="bold"), legend.position="top",
          legend.box="horizontal", legend.title=element_text(size=13),
          legend.text=element_text(size=12), legend.key.size=unit(0.8,"cm"),
          axis.title.y=element_text(size=13), axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=11), plot.title=element_text(size=18,face="bold"),
          plot.subtitle=element_text(size=12,color="gray40"),
          plot.caption=element_text(size=11,color="gray50"))

  ggsave("Fig3_GSVA_GO_trajectories.pdf", p3, width=14, height=14)
  message("Saved: Fig3_GSVA_GO_trajectories.pdf")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: C5 GO Posterior Probability Heatmap
# ══════════════════════════════════════════════════════════════════════════════

if (nrow(rv217_pp_top) > 0) {
  pp_mat <- rv217_pp_top %>%
    tibble::column_to_rownames("pathway") %>%
    dplyr::select(any_of(c("Peak","SetPoint","Chronic"))) %>%
    as.matrix()
  pp_mat[is.na(pp_mat)] <- 0.5
  rownames(pp_mat) <- rownames(pp_mat) %>%
    gsub("^GOBP_|^GOMF_","",.) %>% gsub("_"," ",.) %>%
    stringr::str_to_title() %>% stringr::str_trunc(100)

  go_type_vec <- rv217_pp_top$go_type
  row_ann <- rowAnnotation(
    GO_type = go_type_vec,
    col     = list(GO_type = c("GO:BP"="#1D9E75","GO:MF"="#BA7517")),
    annotation_name_gp = gpar(fontsize=13), annotation_label="GO type",
    simple_anno_size   = unit(0.6,"cm"))

  col_fun <- colorRamp2(c(0,0.15,0.5,0.85,1),
                         c("#D4537E","#F7C1C1","#F1EFE8","#CECBF6","#534AB7"))

  pdf("Fig4_PP_GO_heatmap.pdf", width=15, height=16)
  pushViewport(viewport(y=0.92, height=0.88, just="top"))
  draw(
    Heatmap(pp_mat, name="PP\nBroad>Non-Broad", col=col_fun,
            cluster_rows=TRUE, cluster_columns=FALSE,
            show_row_dend=TRUE, row_dend_width=unit(1.5,"cm"), row_dend_gp=gpar(lwd=1.0),
            row_split=go_type_vec,
            row_names_gp=gpar(fontsize=13), column_names_gp=gpar(fontsize=15,fontface="bold"),
            row_title_gp=gpar(fontsize=15,fontface="bold"), column_title=NULL,
            right_annotation=row_ann,
            cell_fun = function(j,i,x,y,width,height,fill) {
              v <- pp_mat[i,j]
              if (!is.na(v)) grid.text(sprintf("%.2f",v), x, y,
                gp=gpar(fontsize=11, col=ifelse(abs(v-0.5)>0.3,"white","gray35")))
            },
            heatmap_legend_param = list(
              title_gp=gpar(fontsize=13,fontface="bold"), labels_gp=gpar(fontsize=12),
              at=c(0,0.15,0.5,0.85,1),
              labels=c("0.00 (Non-Broad)","0.15","0.50 (uncertain)","0.85","1.00 (Broad)")),
            row_names_max_width=unit(11,"cm"), column_names_rot=0),
    heatmap_legend_side="right", annotation_legend_side="bottom", newpage=FALSE)
  upViewport()
  grid.text("Fig 4: Region-corrected Posterior Probability of Broad Neutralizer Effect",
            x=0.5, y=0.985, just="top", gp=gpar(fontsize=16,fontface="bold",col="gray15"))
  grid.text("C5 GO:BP and GO:MF  |  SetPoint-ranked  |  RV217 Cohort (GSE271442)",
            x=0.5, y=0.958, just="top", gp=gpar(fontsize=13,col="gray40"))
  grid.text("GSE271442 | Dohoon Kim · PromptGenix LLC",
            x=0.98, y=0.012, just="right", gp=gpar(fontsize=10,col="gray55"))
  dev.off()
  message("Saved: Fig4_PP_GO_heatmap_v2.pdf")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Region-specific GO pathways
# ══════════════════════════════════════════════════════════════════════════════

rv217_top_region_go <- rv217_region_effects %>%
  dplyr::filter(sig_region) %>% dplyr::arrange(desc(abs(t_region))) %>%
  head(8) %>% dplyr::pull(pathway)

if (length(rv217_top_region_go) > 0) {
  p5_dat <- rv217_gsva_long %>%
    dplyr::filter(pathway %in% rv217_top_region_go, !is.na(TimeF)) %>%
    dplyr::mutate(pw_label = clean_go(pathway, width=28)) %>%
    dplyr::group_by(pw_label, TimeF, Region) %>%
    dplyr::summarise(mean_score = mean(gsva_score, na.rm=TRUE),
                     se_score   = sd(gsva_score, na.rm=TRUE) / sqrt(dplyr::n()),
                     .groups="drop")

  p5 <- ggplot(p5_dat, aes(TimeF, mean_score, color=Region, group=Region)) +
    geom_line(linewidth=1.2) + geom_point(size=3.2) +
    geom_errorbar(aes(ymin=mean_score-se_score, ymax=mean_score+se_score),
                  width=0.22, alpha=0.55, linewidth=0.7) +
    facet_wrap(~ pw_label, scales="free_y", ncol=2) +
    scale_color_manual(values=region_colors) +
    scale_x_discrete(labels=c("Pre"="Pre","Peak"="Peak\n(18d)",
                               "SetPoint"="Set\n(42d)","Chronic"="CHI\n(980d)")) +
    labs(title    = "Fig 5: Pure Region-specific GO pathways",
         subtitle = paste0("Africa (Female / A1-D)  vs.  Thailand (Transgender / CRF01_AE)\n",
                           "After controlling for Broad/Non-Broad status"),
         x=NULL, y="Mean GSVA score (\u00b1SE)", color=NULL,
         caption="GSE271442 | Dohoon Kim \u00b7 PromptGenix LLC") +
    theme_minimal(base_size=14) +
    theme(strip.text=element_text(size=12,face="bold"), legend.position="top",
          legend.text=element_text(size=13), legend.key.size=unit(0.8,"cm"),
          axis.title.y=element_text(size=13), axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=11), plot.title=element_text(size=18,face="bold"),
          plot.subtitle=element_text(size=12,color="gray40"),
          plot.caption=element_text(size=11,color="gray50"))

  ggsave("Fig5_Region_GO_pathways.pdf", p5, width=14, height=12)
  message("Saved: Fig5_Region_GO_pathways.pdf")
}

message("\nAll 5 figures saved.")
message("  Fig1: 16x13  CrossCohort_GSEA")
message("  Fig2: 16x11  NES_threeway")
message("  Fig3: 14x14  GSVA_GO_trajectories")
message("  Fig4: 15x16  PP_GO_heatmap (PDF)")
message("  Fig5: 14x12  Region_GO_pathways")


################################################################################
# GSE124731 Validation Figure (Supplementary) — 3 panels
# Panel A: DE gene count bar chart (limma vs LMM)
# Panel B: NES scatter (limma vs LMM, rescued pathways labeled)
# Panel C: Significant pathway bar chart by biological theme
################################################################################

setwd("~/Downloads/")

# Theme mapping (target pathways → 4 themes)
theme_map <- dplyr::case_when(
  grepl("NFATC1|B_CELL|IMMUNOGLOBULIN|HUMORAL",       gse124731_comp$pathway) ~ "B cell / NFATC1",
  grepl("INTERFERON|VIRUS|ANTIVIRAL|INNATE_IMMUNE|VIRAL", gse124731_comp$pathway) ~ "IFN / Antiviral",
  grepl("RIBOSOM|TRANSLATION|RRNA",                    gse124731_comp$pathway) ~ "Ribosome / Translation",
  grepl("KILLER|KILLING|LYMPHOCYTE_MEDIATED|CYTOKINE_RECEPTOR", gse124731_comp$pathway) ~ "NK-like",
  TRUE ~ "Other")
gse124731_comp$theme <- theme_map

theme_colors <- c("B cell / NFATC1"="#E41A1C", "IFN / Antiviral"="#377EB8",
                   "Ribosome / Translation"="#4DAF4A", "NK-like"="#FF7F00", "Other"="#999999")

# Panel A: DE gene count bar chart
fig_a_df <- tibble::tibble(
  Method   = factor(c("limma","LMM"), levels=c("limma","LMM")),
  DE_genes = c(gse124731_n_de_limma, gse124731_n_de_lmm))

fig_a <- ggplot(fig_a_df, aes(x=Method, y=DE_genes, fill=Method)) +
  geom_col(width=0.6) +
  geom_text(aes(label=scales::comma(DE_genes)), vjust=-0.5, size=6) +
  scale_fill_manual(values=c("limma"="#B4B2A9","LMM"="#E41A1C")) +
  labs(title    = "A) DE Genes (V\u03b41 vs NK)",
       subtitle = sprintf("LMM: %.1fx more DE genes", gse124731_n_de_lmm/max(gse124731_n_de_limma,1)),
       y="DE genes", x=NULL) +
  theme_minimal(base_size=14) +
  theme(legend.position="none", plot.title=element_text(size=16,face="bold"),
        plot.subtitle=element_text(size=12,color="gray40"),
        axis.text=element_text(size=13), axis.title=element_text(size=13)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15)))

# Panel B: NES scatter (limma vs LMM)
fig_b_df <- gse124731_comp %>%
  dplyr::mutate(
    pathway_short = pathway %>%
      stringr::str_remove("GOBP_|GOMF_|GSE[0-9]+_") %>%
      stringr::str_replace_all("_"," ") %>% stringr::str_to_title() %>% stringr::str_trunc(32))

fig_b <- ggplot(fig_b_df, aes(x=NES_limma, y=NES_lmm, color=category)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
  geom_hline(yintercept=0, linewidth=0.3, color="grey70") +
  geom_vline(xintercept=0, linewidth=0.3, color="grey70") +
  geom_point(size=4, alpha=0.85) +
  ggrepel::geom_text_repel(
    data=fig_b_df %>% dplyr::filter(category %in% c("LMM only","limma only")),
    aes(label=pathway_short), size=4, max.overlaps=30, box.padding=0.5,
    segment.size=0.3, segment.color="grey50", fontface="bold", color="black") +
  ggrepel::geom_text_repel(
    data=fig_b_df %>% dplyr::filter(category=="Both", abs(NES_lmm)>2.5),
    aes(label=pathway_short), size=3.5, max.overlaps=20, box.padding=0.3, segment.size=0.2) +
  scale_color_manual(values=c("Both"="#377EB8","LMM only"="#E41A1C",
                               "limma only"="#FF7F00","Neither"="#CCCCCC"), name="Category") +
  labs(title="B) NES: limma vs LMM", subtitle="Diagonal = perfect agreement",
       x="NES (limma)", y="NES (LMM)") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(size=16,face="bold"), plot.subtitle=element_text(size=12,color="gray40"),
        axis.text=element_text(size=12), axis.title=element_text(size=13),
        legend.position="bottom", legend.title=element_text(size=12), legend.text=element_text(size=11))

# Panel C: Significant pathway bar chart by theme
fig_c_df <- gse124731_comp %>%
  dplyr::filter(padj_lmm < 0.05 | padj_limma < 0.05) %>%
  dplyr::mutate(
    pathway_short = pathway %>%
      stringr::str_remove("GOBP_|GOMF_|GSE[0-9]+_") %>%
      stringr::str_replace_all("_"," ") %>% stringr::str_to_title() %>% stringr::str_trunc(38),
    is_rescued = category == "LMM only") %>%
  dplyr::arrange(NES_lmm)
fig_c_df$pathway_short <- factor(fig_c_df$pathway_short, levels=fig_c_df$pathway_short)

fig_c <- ggplot(fig_c_df, aes(x=NES_lmm, y=pathway_short, fill=theme)) +
  geom_col(width=0.7) +
  geom_point(data=fig_c_df %>% dplyr::filter(is_rescued),
             aes(x=NES_lmm, y=pathway_short),
             shape=8, size=4, color="black", stroke=1.2, show.legend=FALSE) +
  geom_vline(xintercept=0, linewidth=0.4) +
  scale_fill_manual(values=theme_colors, name="Biological Theme") +
  labs(title="C) Significant Pathways", subtitle="\u2605 = LMM-rescued (limma missed)",
       x="NES (LMM)", y=NULL) +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(size=16,face="bold"), plot.subtitle=element_text(size=12,color="gray40"),
        axis.text.y=element_text(size=11), axis.text.x=element_text(size=12),
        axis.title=element_text(size=13), legend.position="bottom",
        legend.title=element_text(size=12), legend.text=element_text(size=11))

# Combine 3 panels
fig_validation <- fig_a + fig_b + fig_c +
  plot_layout(widths=c(0.8,1.3,1.5)) +
  plot_annotation(
    title    = "GSE124731: Methodological Validation (V\u03b41 vs NK, Low-Input RNA-seq)",
    subtitle = sprintf("limma DE: %s  |  LMM DE: %s (%.1fx)  |  IFN pathway rescued by LMM",
                       scales::comma(gse124731_n_de_limma), scales::comma(gse124731_n_de_lmm),
                       gse124731_n_de_lmm/max(gse124731_n_de_limma,1)),
    caption  = "GSE124731 | Dohoon Kim \u00b7 PromptGenix LLC",
    theme    = theme(plot.title=element_text(size=20,face="bold"),
                     plot.subtitle=element_text(size=14,color="gray30"),
                     plot.caption=element_text(size=11,color="gray50")))

ggsave("Fig_GSE124731_validation.pdf", fig_validation, width=22, height=10)
ggsave("Fig_GSE124731_panelA.pdf", fig_a, width=6,  height=7)
ggsave("Fig_GSE124731_panelB.pdf", fig_b, width=9,  height=8)
ggsave("Fig_GSE124731_panelC.pdf", fig_c, width=10, height=8)
message("Saved: Fig_GSE124731_validation.pdf + individual panels A, B, C")

# Summary
cat("\n\u2550\u2550 GSE124731 Validation Summary \u2550\u2550\n")
cat(sprintf("  limma DE genes: %d\n", gse124731_n_de_limma))
cat(sprintf("  LMM DE genes:   %d (%.1fx)\n", gse124731_n_de_lmm, gse124731_n_de_lmm/max(gse124731_n_de_limma,1)))
cat("\n  Pathway categories:\n"); print(table(gse124731_comp$category))
cat("\n  LMM-rescued pathway(s):\n")
gse124731_comp %>%
  dplyr::filter(category == "LMM only") %>%
  dplyr::select(pathway, NES_limma, padj_limma, NES_lmm, padj_lmm) %>%
  print()

