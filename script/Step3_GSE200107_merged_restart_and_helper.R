#!/usr/bin/env Rscript

# =========================================================
# GSE200107 merged restart script
#
# What this script does:
#   1) load the saved Seurat object
#   2) join layers for Seurat v5 compatibility
#   3) create cluster markers (optional, but enabled here)
#   4) annotate clusters into broad cell types
#   5) save the annotated Seurat object
#   6) run paired pseudobulk Post vs Pre by cell type
#   7) provide a helper to convert ENSEMBL IDs to gene symbols
#
# Notes:
#   - Edit cluster_map after inspecting the marker CSV.
#   - The bulk ENSEMBL -> SYMBOL helper is optional and can be used
#     separately by calling convert_ensembl_to_symbol(...).
# =========================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Packages
# -----------------------------
pkgs_cran <- c("dplyr", "purrr", "tibble", "stringr", "ggplot2", "readr")
pkgs_bioc <- c("Seurat", "DESeq2", "Matrix", "AnnotationDbi", "org.Hs.eg.db")

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
  library(purrr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(Seurat)
  library(DESeq2)
  library(Matrix)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# -----------------------------
# 1) Paths
# -----------------------------
accession <- "GSE200107"
project_dir <- file.path(getwd(), accession)
results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "plots")

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Small helpers
# -----------------------------
as_chr1 <- function(x) {
  if (is.list(x)) {
    purrr::map_chr(x, ~ as.character(.x)[1])
  } else {
    as.character(x)
  }
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# -----------------------------
# 3) Load saved Seurat object
# -----------------------------
obj_path <- file.path(results_dir, paste0(accession, "_combined_seurat_annotated.rds"))
if (!file.exists(obj_path)) {
  stop("Cannot find Seurat object: ", obj_path)
}

combined <- readRDS(obj_path)
message("Loaded object: ", obj_path)
message("Cells: ", ncol(combined), " | Genes: ", nrow(combined))

# Seurat v5 compatibility
combined <- JoinLayers(combined)
DefaultAssay(combined) <- "RNA"

# -----------------------------
# 4) Create cluster markers
# -----------------------------
# Uses seurat_clusters as the identity class.
# This is useful for reviewing marker genes before relabeling clusters.
if ("seurat_clusters" %in% colnames(combined@meta.data)) {
  Idents(combined) <- "seurat_clusters"
  combined <- droplevels(combined)

  markers <- FindAllMarkers(
    combined,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  marker_file <- file.path(results_dir, paste0(accession, "_cluster_markers.csv"))
  write.csv(markers, marker_file, row.names = FALSE)
  message("Saved cluster markers: ", marker_file)

  if (nrow(markers) > 0) {
    top_markers <- markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)

    top_marker_file <- file.path(results_dir, paste0(accession, "_top10_markers_per_cluster.csv"))
    write.csv(top_markers, top_marker_file, row.names = FALSE)
    message("Saved top markers: ", top_marker_file)
  }
} else {
  stop("Could not find seurat_clusters in metadata.")
}

# -----------------------------
# 5) Manual cluster -> broad cell type map
# -----------------------------
# Edit this mapping after reviewing the marker file.
# Any cluster not listed will become "Other".
cluster_map <- c(
  "0" = "T_cell",
  "1" = "T_cell",
  "2" = "CD4_T",
  "3" = "CD8_T",
  "4" = "B_cell",
  "5" = "Monocyte",
  "6" = "NK_cell",
  "7" = "Dendritic",
  "8" = "Monocyte",
  "9" = "T_cell",
  "10" = "NK_cell",
  "11" = "B_cell",
  "12" = "Plasma",
  "13" = "Dendritic",
  "14" = "Rare",
  "15" = "Rare",
  "16" = "Rare",
  "17" = "Rare"
)

# If an auto annotation file exists, use it to fill in missing cluster labels.
auto_anno_file <- file.path(results_dir, "GSE200107_cluster_annotation_auto.csv")
if (file.exists(auto_anno_file)) {
  auto_anno <- read.csv(auto_anno_file, stringsAsFactors = FALSE)
  anno_col <- intersect(c("celltype", "annotation", "broad_celltype"), colnames(auto_anno))
  if (length(anno_col) >= 1 && "cluster" %in% colnames(auto_anno)) {
    auto_anno$cluster <- as.character(auto_anno$cluster)
    auto_anno[[anno_col[1]]] <- as.character(auto_anno[[anno_col[1]]])

    for (cl in unique(auto_anno$cluster)) {
      current <- if (cl %in% names(cluster_map)) unname(cluster_map[cl]) else NA_character_
      if (is.na(current) || current == "") {
        cluster_map[[cl]] <- auto_anno[[anno_col[1]]][match(cl, auto_anno$cluster)]
      }
    }
    message("Loaded auto annotation file and merged with manual map.")
  }
}

# Apply annotation to Seurat object.
cluster_ids <- as.character(Idents(combined))
combined$celltype <- unname(cluster_map[cluster_ids])
combined$celltype[is.na(combined$celltype) | combined$celltype == ""] <- "Other"
combined$broad_celltype <- combined$celltype

cat("\nCell type counts after annotation:\n")
print(table(combined$broad_celltype))

write.csv(
  as.data.frame(table(combined$broad_celltype)),
  file.path(results_dir, paste0(accession, "_broad_celltype_counts.csv")),
  row.names = FALSE
)

annotated_path <- file.path(results_dir, paste0(accession, "_combined_seurat_annotated_relabelled.rds"))
saveRDS(combined, annotated_path)
message("Saved annotated object: ", annotated_path)

# -----------------------------
# 6) Pseudobulk helper
# -----------------------------
safe_pseudobulk_deseq <- function(seu, ct, meta_pb, min_cells_per_sample = 20) {
  expr_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"

  meta_ct <- meta_pb %>% filter(broad_celltype == ct)
  if (nrow(meta_ct) == 0) return(NULL)

  # Keep only patients with both Pre and Post.
  keep_patients <- meta_ct %>%
    dplyr::count(patient_id, timepoint, name = "n_cells") %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(n_timepoints = dplyr::n_distinct(timepoint), .groups = "drop") %>%
    dplyr::filter(n_timepoints == 2) %>%
    dplyr::pull(patient_id)

  meta_ct <- meta_ct %>% filter(patient_id %in% keep_patients)
  if (nrow(meta_ct) == 0) return(NULL)

  common_cells <- intersect(colnames(expr_mat), rownames(meta_ct))
  if (length(common_cells) == 0) return(NULL)

  expr_mat <- expr_mat[, common_cells, drop = FALSE]
  meta_ct <- meta_ct[common_cells, , drop = FALSE]

  meta_ct$group_id <- paste(meta_ct$patient_id, meta_ct$timepoint, sep = "__")

  # Require enough cells per pseudobulk sample.
  group_sizes <- table(meta_ct$group_id)
  keep_groups <- names(group_sizes)[group_sizes >= min_cells_per_sample]
  meta_ct <- meta_ct %>% filter(group_id %in% keep_groups)
  if (nrow(meta_ct) == 0) return(NULL)

  common_cells <- intersect(colnames(expr_mat), rownames(meta_ct))
  expr_mat <- expr_mat[, common_cells, drop = FALSE]
  meta_ct <- meta_ct[common_cells, , drop = FALSE]

  group_levels <- unique(meta_ct$group_id)
  if (length(group_levels) < 4) return(NULL)

  pb_list <- lapply(group_levels, function(g) {
    cell_ids <- rownames(meta_ct)[meta_ct$group_id == g]
    Matrix::rowSums(expr_mat[, cell_ids, drop = FALSE])
  })

  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- group_levels
  rownames(pb_mat) <- rownames(expr_mat)

  pb_meta <- data.frame(
    sample = group_levels,
    patient_id = sub("__.*$", "", group_levels),
    timepoint = sub("^.*__", "", group_levels),
    row.names = group_levels,
    stringsAsFactors = FALSE
  )
  pb_meta$timepoint <- factor(pb_meta$timepoint, levels = c("Pre", "Post"))

  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_mat),
    colData = pb_meta,
    design = ~ patient_id + timepoint
  )

  keep_genes <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_genes, ]
  if (nrow(dds) == 0) return(NULL)

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("timepoint", "Post", "Pre"))
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

  safe_ct <- safe_name(ct)
  write.csv(
    res_df,
    file.path(results_dir, paste0(accession, "_", safe_ct, "_Post_vs_Pre.csv")),
    row.names = FALSE
  )
  write.csv(
    sig,
    file.path(results_dir, paste0(accession, "_", safe_ct, "_Post_vs_Pre_sig.csv")),
    row.names = FALSE
  )

  list(res = res_df, sig = sig, dds = dds)
}

# -----------------------------
# 7) Run pseudobulk Post vs Pre by cell type
# -----------------------------
meta_pb <- combined@meta.data %>%
  as.data.frame() %>%
  mutate(
    broad_celltype = as_chr1(broad_celltype),
    patient_id = as_chr1(patient_id),
    timepoint = as_chr1(timepoint)
  ) %>%
  filter(
    !is.na(broad_celltype), broad_celltype != "",
    !is.na(patient_id), patient_id != "",
    !is.na(timepoint), timepoint %in% c("Pre", "Post")
  )

celltypes_to_test <- meta_pb %>%
  dplyr::count(broad_celltype, name = "n_cells") %>%
  filter(n_cells >= 100) %>%
  pull(broad_celltype) %>%
  as.character()

message("Cell types entering pseudobulk:")
print(celltypes_to_test)

pb_results <- list()
summary_rows <- list()

for (ct in celltypes_to_test) {
  message("Pseudobulk testing: ", ct)
  out <- safe_pseudobulk_deseq(combined, ct, meta_pb, min_cells_per_sample = 20)
  if (!is.null(out)) {
    pb_results[[ct]] <- out
    summary_rows[[ct]] <- data.frame(
      celltype = ct,
      n_genes_tested = nrow(out$res),
      n_sig = nrow(out$sig),
      stringsAsFactors = FALSE
    )
  }
}

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  summary_path <- file.path(results_dir, paste0(accession, "_pseudobulk_summary.csv"))
  write.csv(summary_df, summary_path, row.names = FALSE)
  saveRDS(pb_results, file.path(results_dir, paste0(accession, "_pseudobulk_results.rds")))
  message("Saved pseudobulk summary: ", summary_path)
  message("Saved pseudobulk results RDS.")
} else {
  message("No pseudobulk results generated.")
}

# -----------------------------
# 8) Optional helper: ENSEMBL -> SYMBOL for bulk results
# -----------------------------
convert_ensembl_to_symbol <- function(input_csv, output_csv = NULL) {
  df <- read.csv(input_csv, stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(df)) {
    stop("Expected a gene_id column in: ", input_csv)
  }

  ens <- gsub("\\..*$", "", df$gene_id)
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ens),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL", "ENTREZID")
  )

  df2 <- df %>%
    mutate(ENSEMBL = ens) %>%
    left_join(mapping, by = "ENSEMBL")

  if (is.null(output_csv)) {
    output_csv <- sub("\\.csv$", "_symbol.csv", input_csv)
  }
  write.csv(df2, output_csv, row.names = FALSE)
  message("Saved symbol-converted file: ", output_csv)
  invisible(df2)
}

# Example usage:
# convert_ensembl_to_symbol(file.path("GSE206149", "results", "GSE206149_SLIT_vs_Placebo_sig.csv"))

cat("\n✅ Merge/restart script complete\n")
