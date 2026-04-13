#!/usr/bin/env Rscript

# =========================================================
# Pollen allergy transcriptomics figures (consistent style)
#
# Expected directory layout:
#   <working directory>/
#     GSE206149/
#       results/
#         GSE206149_SLIT_vs_Placebo_sig.csv
#         GSE206149_SLIT_vs_Placebo_GO.csv
#         GSE206149_Placebo_vs_Healthy_sig.csv
#         GSE206149_Placebo_vs_Healthy_GO.csv
#     GSE200107/
#       results/
#         GSE200107_pseudobulk_summary.csv
#         GSE200107_<CellType>_Post_vs_Pre_sig.csv
#         GSE200107_combined_seurat_annotated_relabelled.rds
#
# Outputs:
#   <working directory>/figures/
# =========================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Packages
# -----------------------------
pkgs_cran <- c("dplyr", "tibble", "stringr", "readr", "ggplot2", "purrr")
pkgs_bioc <- c("Seurat", "pheatmap", "EnhancedVolcano", "enrichplot")

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
  library(purrr)
  library(Seurat)
  library(pheatmap)
  library(EnhancedVolcano)
  library(enrichplot)
})

# -----------------------------
# 1) Paths
# -----------------------------
project_dir <- getwd()
gse206149_dir <- file.path(project_dir, "GSE206149")
gse200107_dir <- file.path(project_dir, "GSE200107")
fig_dir <- file.path(project_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gse206149_results <- file.path(gse206149_dir, "results")
gse200107_results <- file.path(gse200107_dir, "results")

# -----------------------------
# 2) Global style
# -----------------------------
theme_sts <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    )
}

save_fig <- function(p, filename, width = 7, height = 6, dpi = 300) {
  ggsave(file.path(fig_dir, filename), plot = p, width = width, height = height, dpi = dpi)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    message("Missing file: ", path)
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# -----------------------------
# 3) Figure 1: study design
# -----------------------------
p_fig1 <- ggplot() +
  xlim(0, 10) + ylim(0, 10) +
  theme_void() +
  annotate("rect", xmin = 0.6, xmax = 4.2, ymin = 6.2, ymax = 8.8,
           fill = "#FDE2E2", color = "#C44E52", linewidth = 0.8) +
  annotate("text", x = 2.4, y = 8.25, label = "GSE206149", fontface = "bold", size = 5) +
  annotate("text", x = 2.4, y = 7.4, label = "Bulk RNA-seq", size = 4.2) +
  annotate("text", x = 2.4, y = 6.7, label = "SLIT vs Placebo", size = 4.2) +
  annotate("rect", xmin = 5.8, xmax = 9.4, ymin = 6.2, ymax = 8.8,
           fill = "#E1F0F9", color = "#4C78A8", linewidth = 0.8) +
  annotate("text", x = 7.6, y = 8.25, label = "GSE200107", fontface = "bold", size = 5) +
  annotate("text", x = 7.6, y = 7.4, label = "scRNA-seq + VDJ", size = 4.2) +
  annotate("text", x = 7.6, y = 6.7, label = "Pre vs Post SLIT", size = 4.2) +
  annotate("segment", x = 4.35, xend = 5.55, y = 7.5, yend = 7.5,
           arrow = arrow(length = unit(0.2, "inches")), linewidth = 1.0) +
  annotate("text", x = 5.0, y = 6.9, label = "Integration", size = 4.2, fontface = "bold") +
  annotate("text", x = 5.0, y = 5.8, label = "Mechanism of SLIT in pollen allergy", size = 5, fontface = "bold") +
  ggtitle("Figure 1. Study design")

save_fig(p_fig1, "Fig1_study_design.png", width = 9, height = 6)

# -----------------------------
# 4) Figure 2A: bulk volcano
# -----------------------------
# -----------------------------
# 4) Figure 2A: bulk volcano (USE FULL RESULT FILE, NOT _sig.csv)
# -----------------------------
slit_bulk_all <- safe_read_csv(file.path(gse206149_results, "GSE206149_SLIT_vs_Placebo.csv"))

if (!is.null(slit_bulk_all) && nrow(slit_bulk_all) > 0) {
  if (!"gene_id" %in% colnames(slit_bulk_all)) {
    slit_bulk_all <- slit_bulk_all %>% rownames_to_column("gene_id")
  }

  # make sure required columns exist
  if (!"pvalue" %in% colnames(slit_bulk_all) && "padj" %in% colnames(slit_bulk_all)) {
    slit_bulk_all$pvalue <- slit_bulk_all$padj
  }

  slit_bulk_all <- slit_bulk_all %>%
    mutate(
      significant = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
      pvalue_plot = ifelse(is.na(pvalue) | pvalue <= 0, NA_real_, pvalue)
    )

  p_fig2a <- ggplot(slit_bulk_all, aes(x = log2FoldChange, y = -log10(pvalue_plot))) +
    geom_point(aes(color = significant), alpha = 0.65, size = 1.4, na.rm = TRUE) +
    scale_color_manual(values = c("No" = "grey80", "Yes" = "#D62728")) +
    theme_sts() +
    labs(
      title = "Bulk RNA-seq: SLIT vs Placebo",
      x = "Log2 fold change",
      y = "-log10(p-value)",
      color = "Significant"
    )

  save_fig(p_fig2a, "Fig2A_bulk_volcano.png", width = 7.8, height = 6.2)
}

# -----------------------------
# 5) Figure 2B: GO enrichment dotplot (wrapped labels)
# -----------------------------
# -----------------------------
# 5) Figure 2B: GSEA from full bulk DE table
# -----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)

slit_bulk_all <- safe_read_csv(file.path(gse206149_results, "GSE206149_SLIT_vs_Placebo.csv"))

if (!is.null(slit_bulk_all) && nrow(slit_bulk_all) > 0) {

  # Make sure gene_id exists
  if (!"gene_id" %in% colnames(slit_bulk_all)) {
    slit_bulk_all <- slit_bulk_all %>% rownames_to_column("gene_id")
  }

  # Remove Ensembl version suffix if present
  slit_bulk_all$gene_id <- gsub("\\..*$", "", slit_bulk_all$gene_id)

  # Keep rows with usable log2FC
  slit_bulk_all <- slit_bulk_all %>%
    filter(!is.na(log2FoldChange))

  # Map ENSEMBL -> ENTREZ
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(slit_bulk_all$gene_id),
    keytype = "ENSEMBL",
    columns = c("ENTREZID")
  )

  gsea_input <- slit_bulk_all %>%
    left_join(map, by = c("gene_id" = "ENSEMBL")) %>%
    filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>%
    group_by(ENTREZID) %>%
    summarise(stat = log2FoldChange[which.max(abs(log2FoldChange))], .groups = "drop")

  # Named ranked list for GSEA
  gene_list <- gsea_input$stat
  names(gene_list) <- gsea_input$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Run GSEA
  gsea_res <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_df <- as.data.frame(gsea_res) %>%
      arrange(p.adjust) %>%
      slice_head(n = 15) %>%
      mutate(Description_wrapped = stringr::str_wrap(Description, width = 35))

    # Save table
    write.csv(gsea_df, file.path(gse206149_results, "GSE206149_SLIT_vs_Placebo_GSEA.csv"), row.names = FALSE)

    # Plot
    p_fig2b <- ggplot(gsea_df, aes(x = NES, y = reorder(Description_wrapped, NES))) +
      geom_point(aes(size = setSize, color = p.adjust)) +
      scale_color_gradient(low = "#D62728", high = "#4C78A8", trans = "reverse") +
      theme_sts() +
      theme(axis.text.y = element_text(size = 8)) +
      labs(
        title = "GSEA (SLIT vs Placebo)",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        color = "Adjusted p",
        size = "Gene set size"
      )

    save_fig(p_fig2b, "Fig2B_bulk_GSEA.png", width = 9.0, height = 6.8)
  }
}

# -----------------------------
# 6) Figure 3A: scRNA UMAP by cell type
# -----------------------------
seu_path <- file.path(gse200107_results, "GSE200107_combined_seurat_annotated_relabelled.rds")
if (file.exists(seu_path)) {
  combined <- readRDS(seu_path)
  DefaultAssay(combined) <- "RNA"
  if (!"celltype" %in% colnames(combined@meta.data) && "broad_celltype" %in% colnames(combined@meta.data)) {
    combined$celltype <- combined$broad_celltype
  }

  p_fig3a <- DimPlot(combined, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
    theme_sts() +
    labs(title = "scRNA-seq UMAP by cell type", x = NULL, y = NULL)

  save_fig(p_fig3a, "Fig3A_scRNA_UMAP_celltypes.png", width = 8.5, height = 6.8)
}

# -----------------------------
# Fig3B: cell-type composition
# -----------------------------
# Use the annotated Seurat object already loaded as `combined`
# and create a clean composition plot from all cells.

# -----------------------------
# Fig3B: cell-type composition (FIXED)
# -----------------------------
library(dplyr)
library(purrr)

meta_pb <- combined@meta.data %>%
  as.data.frame() %>%
  mutate(
    broad_celltype = if (is.list(broad_celltype)) {
      purrr::map_chr(broad_celltype, ~ as.character(.x)[1])
    } else {
      as.character(broad_celltype)
    }
  )

cell_counts <- meta_pb %>%
  filter(!is.na(broad_celltype), broad_celltype != "") %>%
  dplyr::count(broad_celltype, name = "n_cells") %>%
  mutate(percent = 100 * n_cells / sum(n_cells)) %>%
  filter(broad_celltype != "Rare")

p_fig3b <- ggplot(cell_counts, aes(x = reorder(broad_celltype, percent), y = percent)) +
  geom_col(fill = "#4C78A8", width = 0.75) +
  coord_flip() +
  theme_sts() +
  labs(
    title = "scRNA-seq cell composition",
    x = "Cell type",
    y = "Percent of cells"
  )

save_fig(p_fig3b, "Fig3B_scRNA_cell_composition.png", width = 8.0, height = 6.3)


