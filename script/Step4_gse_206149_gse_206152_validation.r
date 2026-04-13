#!/usr/bin/env Rscript

# =========================================================
# Validation script for SLIT bulk RNA-seq project
# Primary aim:
#   - validate the GSE206149 SLIT vs Placebo signature in GSE206152
# Secondary aims:
#   - test DEG overlap and direction concordance
#   - compare pathway-level consistency
#   - score the GSE206149 signature in the validation cohort
#
# Assumptions:
#   - You have already run the main GSE206149 analysis and saved:
#       results/GSE206149_SLIT_vs_Placebo.csv
#       results/GSE206149_SLIT_vs_Placebo_sig.csv
#       results/GSE206149_vst_matrix.csv (optional, not required here)
#   - GSE206152 contains a comparable bulk RNA-seq cohort with SLIT / placebo
#     or pre/post labels that can be mapped into group labels.
#
# Output:
#   - Validation DEG table, overlap summary, correlation plot, heatmaps,
#     and a simple signature score comparison.
# =========================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Packages
# -----------------------------
pkgs_cran <- c("dplyr", "tibble", "stringr", "readr", "ggplot2", "tidyr", "purrr")
pkgs_bioc <- c(
  "GEOquery", "DESeq2", "limma", "EnhancedVolcano", "pheatmap",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(GEOquery)
  library(DESeq2)
  library(limma)
  library(EnhancedVolcano)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
})

# -----------------------------
# 1) Paths
# -----------------------------
base_dir <- getwd()
train_accession <- "GSE206149"
valid_accession <- "GSE206152"

train_dir <- file.path(base_dir, train_accession)
valid_dir <- file.path(base_dir, valid_accession)

train_results <- file.path(train_dir, "results")
train_plots <- file.path(train_dir, "plots")

valid_project <- file.path(base_dir, paste0(valid_accession, "_validation"))
valid_meta_dir <- file.path(valid_project, "metadata")
valid_results <- file.path(valid_project, "results")
valid_plots <- file.path(valid_project, "plots")
for (d in c(valid_project, valid_meta_dir, valid_results, valid_plots)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# 2) Helpers
# -----------------------------
clean_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_png <- function(path, expr, width = 1400, height = 1100, res = 160) {
  png(path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

parse_characteristics <- function(df_row) {
  cols <- grep("^characteristics_ch1", names(df_row), value = TRUE)
  out <- list()
  if (length(cols) == 0) return(out)
  vals <- as.character(df_row[cols])
  vals <- vals[!is.na(vals) & vals != ""]
  if (length(vals) == 0) return(out)
  for (v in vals) {
    m <- regexec("^\\s*([^:]+?)\\s*:\\s*(.*?)\\s*$", v)
    parts <- regmatches(v, m)[[1]]
    if (length(parts) == 3) {
      key <- clean_name(parts[2])
      val <- trimws(parts[3])
      if (nzchar(key) && nzchar(val)) {
        out[[key]] <- unique(c(out[[key]], val))
      }
    }
  }
  out
}

collapse_characteristics <- function(meta) {
  parsed <- lapply(seq_len(nrow(meta)), function(i) parse_characteristics(meta[i, , drop = FALSE]))
  keys <- sort(unique(unlist(lapply(parsed, names))))
  if (length(keys) == 0) return(meta)
  for (k in keys) {
    meta[[k]] <- vapply(parsed, function(x) {
      if (k %in% names(x)) paste(unique(na.omit(x[[k]])), collapse = "; ") else NA_character_
    }, character(1))
  }
  meta
}

normalize_group <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("healthy", "healthy_control", "healthy control", "normal control", "control", "non-allergic")) return("Healthy")
  if (x %in% c("placebo")) return("Placebo")
  if (x %in% c("scit", "subcutaneous immunotherapy", "subcutaneous")) return("SCIT")
  if (x %in% c("slit", "sublingual immunotherapy", "sublingual")) return("SLIT")
  if (x %in% c("pre", "baseline", "before")) return("Pre")
  if (x %in% c("post", "after", "followup", "follow-up")) return("Post")
  NA_character_
}

infer_group_from_row <- function(x) {
  blob <- paste(x, collapse = " | ")
  blob <- tolower(blob)
  if (str_detect(blob, "healthy|non-allergic|normal control")) return("Healthy")
  if (str_detect(blob, "placebo")) return("Placebo")
  if (str_detect(blob, "scit|subcutaneous")) return("SCIT")
  if (str_detect(blob, "slit|sublingual")) return("SLIT")
  if (str_detect(blob, "pre|baseline|before")) return("Pre")
  if (str_detect(blob, "post|after|follow-up|followup")) return("Post")
  NA_character_
}

read_count_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    ext <- tolower(sub(".*\\.([A-Za-z0-9]+)\\.gz$", "\\1", path))
  }
  df <- switch(
    ext,
    csv = readr::read_csv(path, show_col_types = FALSE),
    tsv = readr::read_tsv(path, show_col_types = FALSE),
    txt = readr::read_tsv(path, show_col_types = FALSE),
    readr::read_delim(path, delim = "\t", show_col_types = FALSE)
  )
  as.data.frame(df)
}

pick_count_file <- function(root_dir, accession) {
  files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(csv|tsv|txt)(\\.gz)?$", files, ignore.case = TRUE)]
  if (length(files) == 0) stop("No candidate count files found in: ", root_dir)

  priority <- files[grepl("raw_counts|counts|expression|matrix|", basename(files), ignore.case = TRUE)]
  if (length(priority) == 0) priority <- files

  for (f in priority) {
    ok <- tryCatch({
      df <- read_count_table(f)
      if (ncol(df) < 3) return(FALSE)
      num_frac <- sapply(df[-1], function(z) mean(!is.na(suppressWarnings(as.numeric(as.character(z))))))
      any(num_frac > 0.7)
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(f)
  }
  stop("Could not identify a valid count table among candidates for ", accession)
}

make_sample_key <- function(title, geo_accession, sample_lib = NA_character_, subject_id = NA_character_) {
  vals <- c(title, geo_accession, sample_lib, subject_id)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return(NA_character_)
  clean_name(vals[1])
}

match_metadata_to_counts <- function(counts, meta) {
  cn <- colnames(counts)
  rn <- rownames(meta)

  common <- intersect(cn, rn)
  if (length(common) >= 4) {
    counts <- counts[, common, drop = FALSE]
    meta <- meta[common, , drop = FALSE]
    return(list(counts = counts, meta = meta, method = "direct"))
  }

  meta$sample_key <- mapply(
    make_sample_key,
    meta$title,
    meta$geo_accession,
    if ("sample_lib" %in% names(meta)) meta$sample_lib else NA_character_,
    if ("subject_id" %in% names(meta)) meta$subject_id else NA_character_,
    USE.NAMES = FALSE
  )

  common2 <- intersect(clean_name(cn), meta$sample_key)
  if (length(common2) >= 4) {
    map <- setNames(meta$geo_accession[match(common2, meta$sample_key)], common2)
    keep <- !is.na(match(clean_name(cn), common2))
    counts <- counts[, keep, drop = FALSE]
    colnames(counts) <- map[clean_name(colnames(counts))[keep]]
    meta <- meta[colnames(counts), , drop = FALSE]
    rownames(meta) <- meta$geo_accession
    return(list(counts = counts, meta = meta, method = "sample_key"))
  }

  mapped <- rep(NA_character_, length(cn))
  for (i in seq_along(cn)) {
    cclean <- clean_name(cn[i])
    hit <- meta$geo_accession[meta$sample_key == cclean]
    if (length(hit) >= 1 && !is.na(hit[1]) && nzchar(hit[1])) {
      mapped[i] <- hit[1]
      next
    }
    hit <- meta$geo_accession[str_detect(clean_name(meta$title), fixed(cclean))]
    if (length(hit) >= 1) mapped[i] <- hit[1]
  }

  keep <- !is.na(mapped)
  if (sum(keep) >= 4) {
    counts <- counts[, keep, drop = FALSE]
    colnames(counts) <- mapped[keep]
    meta <- meta[colnames(counts), , drop = FALSE]
    rownames(meta) <- meta$geo_accession
    return(list(counts = counts, meta = meta, method = "fallback"))
  }

  stop("Could not match count columns to metadata for validation dataset.")
}

orient_if_needed <- function(counts, meta) {
  cn <- colnames(counts)
  rn <- rownames(counts)
  sample_hits_cols <- sum(clean_name(cn) %in% clean_name(meta$geo_accession), na.rm = TRUE)
  sample_hits_rows <- sum(clean_name(rn) %in% clean_name(meta$geo_accession), na.rm = TRUE)
  if (sample_hits_rows > sample_hits_cols) {
    counts <- t(counts)
  }
  counts <- as.matrix(counts)
  mode(counts) <- "numeric"
  counts[is.na(counts)] <- 0
  counts
}

run_deseq <- function(counts, meta, design_formula) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = meta,
    design = design_formula
  )
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  dds
}

get_results_df <- function(dds, contrast_vec, out_csv = NULL) {
  res <- results(dds, contrast = contrast_vec)
  res_df <- as.data.frame(res) |> rownames_to_column("gene_id")
  if (!is.null(out_csv)) write.csv(res_df, out_csv, row.names = FALSE)
  res_df
}

zscore_rows <- function(mat) {
  t(scale(t(mat)))
}

signature_score <- function(vst_mat, genes) {
  genes <- intersect(genes, rownames(vst_mat))
  if (length(genes) < 5) return(NULL)
  colMeans(vst_mat[genes, , drop = FALSE], na.rm = TRUE)
}

load_gene_sets <- function() {
  list(
    immune_core = c("LEUKOCYTE_MIGRATION", "CHEMOTAXIS", "LYMPHOCYTE_ACTIVATION", "HUMORAL_IMMUNE_RESPONSE"),
    t_cell = c("T_CELL_ACTIVATION", "T_CELL_DIFFERENTIATION", "REGULATION_OF_T_CELL_ACTIVATION"),
    monocyte = c("INFLAMMATORY_RESPONSE", "MONOCYTE_ACTIVATION", "CYTOKINE_PRODUCTION"),
    antigen = c("ANTIGEN_PROCESSING_AND_PRESENTATION", "MHC_CLASS_II_PROTEIN_COMPLEX_ASSEMBLY"),
    remodeling = c("EXTRACELLULAR_MATRIX_ORGANIZATION", "CYTOSKELETON_ORGANIZATION", "CELL_MIGRATION")
  )
}

# -----------------------------
# 3) Load training results
# -----------------------------
train_sig_path <- file.path(train_results, paste0(train_accession, "_SLIT_vs_Placebo_sig.csv"))
train_all_path <- file.path(train_results, paste0(train_accession, "_SLIT_vs_Placebo.csv"))
train_vst_path <- file.path(train_results, paste0(train_accession, "_vst_matrix.csv"))

if (!file.exists(train_sig_path) || !file.exists(train_all_path)) {
  stop("Could not find training result files in ", train_results, ". Run the main GSE206149 analysis first.")
}

train_sig <- read.csv(train_sig_path)
train_all <- read.csv(train_all_path)
train_sig <- train_sig %>% mutate(direction = ifelse(log2FoldChange > 0, "up", "down"))
train_all <- train_all %>% mutate(direction = ifelse(log2FoldChange > 0, "up", "down"))
train_up <- train_sig$gene_id[train_sig$log2FoldChange > 0]
train_down <- train_sig$gene_id[train_sig$log2FoldChange < 0]
train_universe <- train_all$gene_id

# -----------------------------
# 4) Load validation GEO metadata
# -----------------------------
message("Loading validation GEO metadata for ", valid_accession, "...")
geo_obj <- getGEO(valid_accession, GSEMatrix = TRUE, getGPL = FALSE)
if (length(geo_obj) > 1) {
  message("Validation GEO has multiple ExpressionSets; using the first one.")
}
eset <- geo_obj[[1]]
meta <- pData(eset) |> as.data.frame() |> rownames_to_column("geo_row")
meta <- collapse_characteristics(meta)

meta$group <- NA_character_
for (i in seq_len(nrow(meta))) {
  vals <- unlist(meta[i, ], use.names = FALSE)
  meta$group[i] <- infer_group_from_row(vals)
}

# fallback to common metadata fields
for (nm in c("treatment", "group", "status", "condition", "timepoint")) {
  if (nm %in% names(meta)) {
    meta$group <- ifelse(is.na(meta$group), vapply(meta[[nm]], normalize_group, character(1)), meta$group)
  }
}

meta$tissue <- case_when(
  str_detect(tolower(paste(meta$title, meta$source_name_ch1, sep = " | ")), "nasal") ~ "Nasal",
  str_detect(tolower(paste(meta$title, meta$source_name_ch1, sep = " | ")), "pbmc|blood") ~ "PBMC",
  TRUE ~ "Unknown"
)

meta$group <- factor(meta$group, levels = c("Healthy", "Placebo", "SCIT", "SLIT", "Pre", "Post"))
meta$tissue <- factor(meta$tissue, levels = c("Nasal", "PBMC", "Unknown"))
rownames(meta) <- meta$geo_accession
write.csv(meta, file.path(valid_meta_dir, paste0(valid_accession, "_metadata_parsed.csv")), row.names = FALSE)

message("Validation group counts:")
print(table(meta$group, useNA = "ifany"))
print(table(meta$tissue, useNA = "ifany"))

# -----------------------------
# 5) Download validation files and import counts
# -----------------------------
message("Downloading validation supplementary files...")
getGEOSuppFiles(valid_accession, makeDirectory = FALSE, baseDir = valid_project)
count_path <- pick_count_file(valid_project, valid_accession)
message("Selected validation count file: ", count_path)

counts_df <- read_count_table(count_path)
message("Validation raw count table: ", nrow(counts_df), " x ", ncol(counts_df))

gene_col <- names(counts_df)[1]
counts_df[[gene_col]] <- gsub("\\..*$", "", as.character(counts_df[[gene_col]]))
counts <- counts_df %>% column_to_rownames(gene_col) %>% as.data.frame()
counts[] <- lapply(counts, function(x) suppressWarnings(as.numeric(as.character(x))))
counts[is.na(counts)] <- 0
counts <- orient_if_needed(counts, meta)
matched <- match_metadata_to_counts(counts, meta)
counts <- matched$counts
meta <- matched$meta
message("Validation sample matching method: ", matched$method)
message("Validation matched samples: ", ncol(counts))

counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
libsize <- colSums(counts, na.rm = TRUE)
zero_samples <- names(libsize)[libsize == 0]
if (length(zero_samples) > 0) {
  counts <- counts[, libsize > 0, drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
}
stopifnot(identical(colnames(counts), rownames(meta)))

# -----------------------------
# 6) Choose validation contrast + enforce SLIT vs Placebo and raw-count checks
# -----------------------------
# Keep only SLIT and Placebo samples (as requested)
meta <- meta[meta$group %in% c("Placebo", "SLIT"), , drop = FALSE]
counts <- counts[, rownames(meta), drop = FALSE]

# Ensure factor levels
meta$group <- factor(meta$group, levels = c("Placebo", "SLIT"))

# Use limma for validation because the matrix contains normalized/log-scale values.
# Keep only SLIT and Placebo samples (as requested)
meta <- meta[meta$group %in% c("Placebo", "SLIT"), , drop = FALSE]
counts <- counts[, rownames(meta), drop = FALSE]
meta$group <- factor(meta$group, levels = c("Placebo", "SLIT"))

contrast_name <- "SLIT_vs_Placebo"
message("Using validation contrast: ", contrast_name)
message("Validation matrix is not raw counts, so using limma.")

# Expression matrix for limma
expr_mat <- as.matrix(counts)
mode(expr_mat) <- "numeric"

# Optional PCA on normalized data
# prcomp(scale.=TRUE) fails if any row/column has zero variance, so remove them first.
expr_pca <- expr_mat
row_vars <- apply(expr_pca, 1, var, na.rm = TRUE)
expr_pca <- expr_pca[row_vars > 0, , drop = FALSE]

if (nrow(expr_pca) >= 2 && ncol(expr_pca) >= 2) {
  pca <- prcomp(t(expr_pca), scale. = TRUE)
  pca_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = meta[rownames(pca$x), "group"],
    tissue = meta[rownames(pca$x), "tissue"]
  )
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = tissue)) +
    geom_point(size = 3) +
    theme_bw(base_size = 12) +
    labs(title = paste0(valid_accession, " PCA (limma input)"), color = "Group", shape = "Tissue")
  ggsave(file.path(valid_plots, paste0(valid_accession, "_PCA.png")), p_pca, width = 8, height = 6, dpi = 300)
} else {
  message("Skipping PCA because too few variable genes remain after filtering.")
}

# limma differential expression
# Use group order Placebo -> SLIT so the coefficient is SLIT vs Placebo
if (sum(meta$group == "Placebo") < 2 || sum(meta$group == "SLIT") < 2) {
  stop("Need at least 2 samples per group for limma validation.")
}
design <- model.matrix(~ 0 + group, data = meta)
# Rename columns to the clean group labels expected by makeContrasts()
colnames(design) <- sub("^group", "", colnames(design))
fit <- lmFit(expr_mat, design)
contrast <- makeContrasts(SLIT - Placebo, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
res_df <- topTable(fit2, number = Inf, sort.by = "P") %>% rownames_to_column("gene_id")
write.csv(res_df, file.path(valid_results, paste0(valid_accession, "_", contrast_name, ".csv")), row.names = FALSE)

sig_df <- res_df %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 1)
write.csv(sig_df, file.path(valid_results, paste0(valid_accession, "_", contrast_name, "_sig.csv")), row.names = FALSE)

safe_png(file.path(valid_plots, paste0(valid_accession, "_", contrast_name, "_volcano.png")), {
  limma_volc <- res_df %>% mutate(pvalue = P.Value, log2FoldChange = logFC, padj = adj.P.Val)
  print(
    EnhancedVolcano(
      limma_volc,
      lab = ifelse(is.na(limma_volc$padj), "", limma_volc$gene_id),
      x = "log2FoldChange",
      y = "pvalue",
      title = paste0(valid_accession, " ", contrast_name, " (limma)"),
      pCutoff = 0.05,
      FCcutoff = 1
    )
  )
})

# Normalize column names for downstream validation steps
res_df$log2FoldChange <- res_df$logFC
res_df$padj <- res_df$adj.P.Val
res_df$pvalue <- res_df$P.Value
sig_df$log2FoldChange <- sig_df$logFC
sig_df$padj <- sig_df$adj.P.Val

# -----------------------------
# 9) DEG overlap and direction concordance
# -----------------------------
train_sig$direction <- ifelse(train_sig$log2FoldChange > 0, "up", "down")
sig_df$direction <- ifelse(sig_df$log2FoldChange > 0, "up", "down")

common_sig <- intersect(train_sig$gene_id, sig_df$gene_id)
common_all <- intersect(train_all$gene_id, res_df$gene_id)

# Direction consistency among overlapping significant genes
common_dir <- if (length(common_sig) > 0) {
  train_map <- train_sig %>% select(gene_id, train_lfc = log2FoldChange, train_padj = padj, train_dir = direction)
  valid_map <- sig_df %>% select(gene_id, valid_lfc = log2FoldChange, valid_padj = padj, valid_dir = direction)
  overlap_tbl <- inner_join(train_map, valid_map, by = "gene_id")
  overlap_tbl$concordant <- overlap_tbl$train_dir == overlap_tbl$valid_dir
  overlap_tbl
} else {
  tibble()
}

# Correlation of logFC across all common genes with numeric values
corr_tbl <- if (length(common_all) > 10) {
  tmp <- inner_join(
    train_all %>% select(gene_id, train_lfc = log2FoldChange),
    res_df %>% select(gene_id, valid_lfc = log2FoldChange),
    by = "gene_id"
  ) %>% filter(is.finite(train_lfc), is.finite(valid_lfc))
  tmp
} else {
  tibble()
}

lfc_cor <- NA_real_
p_cor <- NA_real_
if (nrow(corr_tbl) >= 10) {
  cor_out <- cor.test(corr_tbl$train_lfc, corr_tbl$valid_lfc, method = "spearman")
  lfc_cor <- unname(cor_out$estimate)
  p_cor <- cor_out$p.value
  p_cor_plot <- ggplot(corr_tbl, aes(train_lfc, valid_lfc)) +
    geom_point(size = 1.6, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Log2FC correlation: ", valid_accession, " vs ", train_accession),
      subtitle = paste0("Spearman rho = ", signif(lfc_cor, 3), ", p = ", signif(p_cor, 3)),
      x = paste0(train_accession, " log2FC"),
      y = paste0(valid_accession, " log2FC")
    )
  ggsave(file.path(valid_plots, paste0(valid_accession, "_logFC_correlation.png")), p_cor_plot, width = 7, height = 6, dpi = 300)
}

# Simple overlap summary and Fisher-style enrichment against the train universe
n_train_sig <- nrow(train_sig)
n_valid_sig <- nrow(sig_df)
n_overlap <- length(common_sig)

# universe is the shared gene set across both analyses
shared_universe <- intersect(train_universe, res_df$gene_id)
if (length(shared_universe) == 0) shared_universe <- unique(c(train_all$gene_id, res_df$gene_id))

# Fisher test on significant overlap in shared universe
train_sig_set <- shared_universe %in% train_sig$gene_id
valid_sig_set <- shared_universe %in% sig_df$gene_id
ft <- fisher.test(matrix(c(
  sum(train_sig_set & valid_sig_set),
  sum(train_sig_set & !valid_sig_set),
  sum(!train_sig_set & valid_sig_set),
  sum(!train_sig_set & !valid_sig_set)
), nrow = 2, byrow = TRUE))

# -----------------------------
# 10) Pathway-level validation
# -----------------------------
# Use GO-BP enrichment on each dataset and compare top terms.
map_ids_to_entrez <- function(ids) {
  ids <- unique(na.omit(as.character(ids)))
  ids <- gsub("\\..*$", "", ids)

  candidate_keytypes <- c("ENSEMBL", "SYMBOL", "ALIAS", "ENTREZID")
  best_keytype <- NULL
  best_hits <- 0

  for (kt in candidate_keytypes) {
    valid_keys <- tryCatch(AnnotationDbi::keys(org.Hs.eg.db, keytype = kt), error = function(e) character(0))
    hits <- sum(ids %in% valid_keys)
    if (hits > best_hits) {
      best_hits <- hits
      best_keytype <- kt
    }
  }

  if (is.null(best_keytype) || best_hits == 0) {
    message("No valid gene ID keytype found for mapping; skipping Entrez conversion.")
    return(character(0))
  }

  message("Using keytype for mapping: ", best_keytype, " (", best_hits, " matching IDs)")
  suppressWarnings({
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ids[ids %in% AnnotationDbi::keys(org.Hs.eg.db, keytype = best_keytype)],
      keytype = best_keytype,
      columns = c("ENTREZID")
    )
  })
  unique(na.omit(m$ENTREZID))
}

train_entrez <- map_ids_to_entrez(train_all$gene_id)
valid_entrez <- map_ids_to_entrez(res_df$gene_id)

train_ego <- NULL
valid_ego <- NULL
if (length(train_entrez) >= 10) {
  train_ego <- enrichGO(gene = train_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
}
if (length(valid_entrez) >= 10) {
  valid_ego <- enrichGO(gene = valid_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
}

if (!is.null(train_ego)) {
  write.csv(as.data.frame(train_ego), file.path(valid_results, paste0(train_accession, "_GO_for_validation.csv")), row.names = FALSE)
}
if (!is.null(valid_ego)) {
  write.csv(as.data.frame(valid_ego), file.path(valid_results, paste0(valid_accession, "_GO_BP.csv")), row.names = FALSE)
}

# Compare top pathway labels by overlap of descriptions
train_terms <- if (!is.null(train_ego) && nrow(as.data.frame(train_ego)) > 0) {
  as.data.frame(train_ego) %>% arrange(p.adjust) %>% slice_head(n = 20) %>% pull(Description)
} else character(0)
valid_terms <- if (!is.null(valid_ego) && nrow(as.data.frame(valid_ego)) > 0) {
  as.data.frame(valid_ego) %>% arrange(p.adjust) %>% slice_head(n = 20) %>% pull(Description)
} else character(0)
term_overlap <- intersect(train_terms, valid_terms)

# -----------------------------
# 11) Signature score validation
# -----------------------------
# Compute a simple sample-level score from the training signature genes.
# Positive score means the sample looks more like the SLIT-up signature.
train_up_score <- signature_score(expr_mat, train_up)
train_down_score <- signature_score(expr_mat, train_down)
if (!is.null(train_up_score) && !is.null(train_down_score)) {
  sig_score <- train_up_score - train_down_score
  score_df <- data.frame(
    sample = names(sig_score),
    score = as.numeric(sig_score),
    group = meta[names(sig_score), "group"]
  )
  p_score <- ggplot(score_df, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Training signature score in ", valid_accession),
      subtitle = "Higher score means the validation sample resembles the SLIT-up / SLIT-down pattern from GSE206149",
      x = "Validation group",
      y = "Signature score"
    )
  ggsave(file.path(valid_plots, paste0(valid_accession, "_signature_score.png")), p_score, width = 8, height = 6, dpi = 300)
}

# Optional heatmap for top overlapping genes
if (nrow(common_dir) >= 5) {
  top_overlap <- common_dir %>%
    mutate(score = abs(train_lfc) + abs(valid_lfc)) %>%
    arrange(desc(score)) %>%
    slice_head(n = min(30, n()))
  top_genes <- top_overlap$gene_id
  mat <- expr_mat[intersect(top_genes, rownames(expr_mat)), , drop = FALSE]
  if (nrow(mat) >= 2) {
    safe_png(file.path(valid_plots, paste0(valid_accession, "_top_overlap_heatmap.png")), {
      pheatmap(zscore_rows(mat), scale = "none", show_colnames = FALSE,
               main = paste0("Top overlap genes: ", train_accession, " vs ", valid_accession))
    })
  }
}

# -----------------------------
# 12) Write validation summary
# -----------------------------
summary_df <- tibble(
  metric = c(
    "validation_accession",
    "contrast_used",
    "validation_sig_genes",
    "train_sig_genes",
    "overlapping_sig_genes",
    "overlap_direction_concordance",
    "spearman_logFC_rho",
    "spearman_logFC_p",
    "fisher_overlap_p",
    "shared_GO_term_overlap"
  ),
  value = c(
    valid_accession,
    contrast_name,
    n_valid_sig,
    n_train_sig,
    n_overlap,
    if (nrow(common_dir) > 0) mean(common_dir$concordant) else NA_real_,
    lfc_cor,
    p_cor,
    ft$p.value,
    length(term_overlap)
  )
)
write.csv(summary_df, file.path(valid_results, paste0(valid_accession, "_validation_summary.csv")), row.names = FALSE)

# Also write a more detailed overlap table for figure building.
if (nrow(common_dir) > 0) {
  write.csv(common_dir, file.path(valid_results, paste0(valid_accession, "_overlapping_sig_genes.csv")), row.names = FALSE)
}

message("Validation complete.")
message("Outputs written to: ", valid_project)

source("fig2_4_combined.R")
fig5 <- plot_grid(
  make_panel(file.path(fig_dir, "GSE206152_PCA.png")),
  make_panel(file.path(fig_dir, "GSE206152_signature_score.png")),
  labels = c("A", "B"),          # ⭐ add labels
  label_size = 18,
  label_fontface = "bold",       # ⭐ match STS style
  label_x = 0.02,                # ⭐ move slightly inward
  label_y = 0.98,                # ⭐ top-left corner
  ncol = 2,
  rel_widths = c(1, 1)
)

ggsave(
  file.path(fig_dir, "Fig5_validation_combined.png"),
  fig5,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"   # ⭐ important
)
