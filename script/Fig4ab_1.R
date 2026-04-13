# =========================================================
# Fig 4A/4B: SLIT-responsive cell-type mechanism figures
#
# Fig 4A: Key immune genes heatmap across cell types
# Fig 4B: Pathway summary heatmap across cell types
#
# Input:
#   ./GSE200107/results/GSE200107_<CellType>_Post_vs_Pre.csv
#
# Output:
#   ./figures/Fig4A_key_immune_genes_heatmap.png
#   ./figures/Fig4B_pathway_summary_heatmap.png
# =========================================================

options(stringsAsFactors = FALSE)

library(dplyr)
library(tibble)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)

# -----------------------------
# Paths
# -----------------------------
project_dir <- getwd()
results_dir <- file.path(project_dir, "GSE200107", "results")
fig_dir <- file.path(project_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Theme helper
# -----------------------------
theme_sts <- function(base_size = 14) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# -----------------------------
# Curated genes for Fig 4A
# -----------------------------
key_genes <- c(
  "CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "IL7R", "CCR7", "LTB", "LEF1", "TCF7",
  "NKG7", "GNLY", "CTSW", "GZMB", "GZMK", "PRF1", "KLRD1", "FCGR3A",
  "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "XBP1", "SDC1",
  "LYZ", "S100A8", "S100A9", "FCN1", "LST1", "CTSS", "FCER1G", "CST3", "AIF1",
  "FCER1A", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
  "IL1B", "TNF", "CXCL8", "CCL2", "CCL3", "CCL4", "MMP14", "GSN", "ADGRE2",
  "ISG15", "IFIT1", "IFIT3", "IRF7", "STAT1"
)

# -----------------------------
# Curated pathway gene sets for Fig 4B
# -----------------------------
pathway_sets <- list(
  "T_cell_activation" = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "IL7R", "CCR7", "LTB", "LEF1", "TCF7"),
  "Cytotoxic_NK" = c("NKG7", "GNLY", "CTSW", "GZMB", "GZMK", "PRF1", "KLRD1", "FCGR3A"),
  "B_cell_plasma" = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "XBP1", "SDC1"),
  "Monocyte_inflammation" = c("LYZ", "S100A8", "S100A9", "FCN1", "LST1", "CTSS", "FCER1G", "CST3", "AIF1", "IL1B", "TNF", "CXCL8", "CCL2", "CCL3", "CCL4"),
  "Antigen_presentation" = c("FCER1A", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CST3"),
  "Interferon_response" = c("ISG15", "IFIT1", "IFIT3", "IRF7", "STAT1"),
  "Remodeling" = c("MMP14", "GSN", "ADGRE2")
)

# -----------------------------
# Read pseudobulk files
# -----------------------------
pb_files <- list.files(results_dir, pattern = "^GSE200107_.*_Post_vs_Pre\\.csv$", full.names = TRUE)
pb_files <- pb_files[!grepl("_sig\\.csv$", pb_files)]

if (length(pb_files) == 0) stop("No pseudobulk files found in: ", results_dir)

# -----------------------------
# Helpers
# -----------------------------
get_celltype_from_file <- function(f) {
  ct <- basename(f)
  ct <- gsub("^GSE200107_", "", ct)
  ct <- gsub("_Post_vs_Pre\\.csv$", "", ct)
  ct <- gsub("_", " ", ct)
  ct
}

map_ensembl_to_symbol <- function(ids) {
  ids <- gsub("\\..*$", "", ids)
  ens <- unique(ids[grepl("^ENSG", ids)])

  out <- setNames(rep(NA_character_, length(ids)), ids)
  if (length(ens) == 0) {
    out[!is.na(ids)] <- ids[!is.na(ids)]
    return(out)
  }

  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ens,
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )
  map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", , drop = FALSE]
  map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]

  idx_ens <- grepl("^ENSG", ids)
  out[!idx_ens] <- ids[!idx_ens]

  if (nrow(map) > 0) {
    out[idx_ens] <- unname(map$SYMBOL[match(ids[idx_ens], map$ENSEMBL)])
    out[idx_ens & is.na(out)] <- ids[idx_ens & is.na(out)]
  } else {
    out[idx_ens] <- ids[idx_ens]
  }
  out
}

read_pseudobulk <- function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(df)) df <- df %>% rownames_to_column("gene_id")
  df$celltype <- get_celltype_from_file(f)
  df$gene_id_clean <- gsub("\\..*$", "", as.character(df$gene_id))
  df$SYMBOL <- map_ensembl_to_symbol(df$gene_id_clean)
  df$SYMBOL <- ifelse(is.na(df$SYMBOL) | df$SYMBOL == "", df$gene_id_clean, df$SYMBOL)
  df
}

# =========================================================
# Fig 4A: key immune genes heatmap
# =========================================================
all_hits <- purrr::map_dfr(pb_files, read_pseudobulk)

# keep only curated key immune genes, and exclude technical clusters
heat_df <- all_hits %>%
  filter(SYMBOL %in% key_genes) %>%
  filter(!celltype %in% c("Rare", "Other")) %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(celltype, SYMBOL) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(celltype, SYMBOL, log2FoldChange)

heat_mat <- heat_df %>%
  pivot_wider(names_from = celltype, values_from = log2FoldChange)

heat_mat_df <- as.data.frame(heat_mat)
rownames(heat_mat_df) <- heat_mat_df$SYMBOL
heat_mat_df$SYMBOL <- NULL
heat_mat_df <- as.matrix(heat_mat_df)
heat_mat_df[is.na(heat_mat_df)] <- 0

# order by variance
heat_mat_df <- heat_mat_df[order(apply(heat_mat_df, 1, var), decreasing = TRUE), , drop = FALSE]

# optional cell ordering
preferred_cols <- c("T cell", "CD4 T", "CD8 T", "B cell", "Monocyte", "Dendritic", "NK cell", "Plasma")
current_cols <- colnames(heat_mat_df)
ordered_cols <- c(intersect(preferred_cols, current_cols), setdiff(current_cols, preferred_cols))
heat_mat_df <- heat_mat_df[, ordered_cols, drop = FALSE]

png(file.path(fig_dir, "Fig4A_key_immune_genes_heatmap.png"), width = 1700, height = 1400, res = 200)
pheatmap(
  heat_mat_df,
  scale = "row",
  color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 45
)
dev.off()

write.csv(heat_df, file.path(fig_dir, "Fig4A_key_immune_genes_table.csv"), row.names = FALSE)

# =========================================================
# Fig 4B: pathway summary across cell types
# =========================================================
# Use the same pseudobulk tables and compute pathway score as
# mean log2FC of genes in each pathway gene set.

celltype_scores <- purrr::map_dfr(pb_files, function(f) {
  df <- read_pseudobulk(f)

  # keep one row per symbol
  df <- df %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup()

  ct <- unique(df$celltype)

  out <- purrr::imap_dfr(pathway_sets, function(genes, pathway_name) {
    sub <- df %>% filter(SYMBOL %in% genes)
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      celltype = ct,
      pathway = pathway_name,
      pathway_score = mean(sub$log2FoldChange, na.rm = TRUE),
      n_genes_used = nrow(sub),
      stringsAsFactors = FALSE
    )
  })

  out
})

celltype_scores <- celltype_scores %>%
  filter(!celltype %in% c("Rare", "Other"))

pathway_mat <- celltype_scores %>%
  dplyr::select(pathway, celltype, pathway_score) %>%
  pivot_wider(names_from = celltype, values_from = pathway_score)

pathway_mat_df <- as.data.frame(pathway_mat)
rownames(pathway_mat_df) <- pathway_mat_df$pathway
pathway_mat_df$pathway <- NULL
pathway_mat_df <- as.matrix(pathway_mat_df)
pathway_mat_df[is.na(pathway_mat_df)] <- 0

desired_pathways <- c("T_cell_activation", "Cytotoxic_NK", "B_cell_plasma", "Monocyte_inflammation", "Antigen_presentation", "Interferon_response", "Remodeling")
desired_pathways <- intersect(desired_pathways, rownames(pathway_mat_df))
pathway_mat_df <- pathway_mat_df[desired_pathways, , drop = FALSE]

preferred_cols2 <- c("T cell", "CD4 T", "CD8 T", "B cell", "Monocyte", "Dendritic", "NK cell", "Plasma")
current_cols2 <- colnames(pathway_mat_df)
ordered_cols2 <- c(intersect(preferred_cols2, current_cols2), setdiff(current_cols2, preferred_cols2))
pathway_mat_df <- pathway_mat_df[, ordered_cols2, drop = FALSE]

write.csv(celltype_scores, file.path(fig_dir, "Fig4B_pathway_summary_table.csv"), row.names = FALSE)

png(file.path(fig_dir, "Fig4B_pathway_summary_heatmap.png"), width = 1600, height = 1000, res = 200)
pheatmap(
  pathway_mat_df,
  scale = "row",
  color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)
dev.off()

cat("\n✅ Fig 4A and 4B completed\n")