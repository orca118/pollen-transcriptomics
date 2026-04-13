#!/usr/bin/env Rscript

# =========================================================
# GSE200107 scRNA-seq / VDJ analysis focused on SLIT
# Main question: what changes from PRE to POST SLIT?
#
# Workflow:
#   1) Download GEO metadata + raw files
#   2) Build one Seurat object per expression sample
#   3) Merge, QC, normalize, cluster, UMAP
#   4) JoinLayers() for Seurat v5 marker testing
#   5) Cluster markers + broad immune annotation
#   6) Module scores (Th2, Treg, cytotoxic, myeloid, B cell)
#   7) Paired pseudobulk DE (Post vs Pre) by broad cell type
#   8) Optional VDJ file index export
# =========================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Packages
# -----------------------------
cran_pkgs <- c(
  "dplyr", "tibble", "stringr", "readr", "ggplot2", "patchwork", "purrr"
)
bioc_pkgs <- c(
  "GEOquery", "Matrix", "DESeq2", "EnhancedVolcano", "pheatmap",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Please install Seurat before running this script.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(GEOquery)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -----------------------------
# 1) Paths
# -----------------------------
accession <- "GSE200107"
project_dir <- file.path(getwd(), accession)
raw_dir <- file.path(project_dir, "raw")
meta_dir <- file.path(project_dir, "metadata")
results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "plots")

for (d in c(project_dir, raw_dir, meta_dir, results_dir, plots_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# 2) Helper functions
# -----------------------------
clean_text <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_patient_id <- function(title) {
  m <- stringr::str_match(title, "JCP_Patient_([0-9]+)")
  out <- m[, 2]
  ifelse(is.na(out), NA_character_, out)
}

parse_timepoint <- function(title) {
  title <- tolower(title)
  dplyr::case_when(
    stringr::str_detect(title, "_pre_") ~ "Pre",
    stringr::str_detect(title, "_post_") ~ "Post",
    TRUE ~ "Unknown"
  )
}

parse_assay <- function(title) {
  title <- tolower(title)
  dplyr::case_when(
    stringr::str_detect(title, "_expression") ~ "Expression",
    stringr::str_detect(title, "_vdj") ~ "VDJ",
    TRUE ~ "Unknown"
  )
}

safe_read_10x_mtx <- function(prefix) {
  candidates <- list(
    mtx = c(paste0(prefix, "_matrix.mtx.gz"), paste0(prefix, ".mtx.gz")),
    features = c(paste0(prefix, "_features.tsv.gz"), paste0(prefix, "_genes.tsv.gz"), paste0(prefix, ".features.tsv.gz")),
    barcodes = c(paste0(prefix, "_barcodes.tsv.gz"), paste0(prefix, ".barcodes.tsv.gz"))
  )

  find_first <- function(paths) {
    hit <- paths[file.exists(paths)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }

  mtx_file <- find_first(candidates$mtx)
  features_file <- find_first(candidates$features)
  barcodes_file <- find_first(candidates$barcodes)

  if (is.na(mtx_file)) stop("Missing matrix file for prefix: ", prefix)
  if (is.na(features_file)) stop("Missing features file for prefix: ", prefix)
  if (is.na(barcodes_file)) stop("Missing barcodes file for prefix: ", prefix)

  mat <- Matrix::readMM(mtx_file)
  feats <- data.table::fread(features_file, header = FALSE, data.table = FALSE)
  barcodes <- readr::read_tsv(barcodes_file, col_names = FALSE, show_col_types = FALSE)[[1]]

  feature_names <- if (ncol(feats) >= 2) feats[[2]] else feats[[1]]
  feature_names <- make.unique(as.character(feature_names))

  if (nrow(mat) != length(feature_names) || ncol(mat) != length(barcodes)) {
    stop(
      "Dimension mismatch for prefix ", prefix,
      " | matrix: ", nrow(mat), " x ", ncol(mat),
      " | features: ", length(feature_names),
      " | barcodes: ", length(barcodes)
    )
  }

  rownames(mat) <- feature_names
  colnames(mat) <- barcodes
  mat
}

make_seurat_from_prefix <- function(prefix, meta_row) {
  mat <- safe_read_10x_mtx(prefix)

  sample_name <- as.character(meta_row$geo_accession)
  if (is.na(sample_name) || !nzchar(sample_name)) sample_name <- as.character(meta_row$title)
  sample_name <- gsub("[^A-Za-z0-9_]+", "_", sample_name)

  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )

  obj$sample_id <- sample_name
  obj$geo_accession <- sample_name
  obj$title <- as.character(meta_row$title)
  obj$patient_id <- if (!is.null(meta_row$patient_id)) as.character(meta_row$patient_id) else NA_character_
  obj$timepoint <- if (!is.null(meta_row$timepoint)) as.character(meta_row$timepoint) else "Unknown"
  obj$assay <- if (!is.null(meta_row$assay)) as.character(meta_row$assay) else "Expression"
  obj$treatment <- "SLIT"
  obj$sample_key <- if (!is.null(meta_row$sample_key)) as.character(meta_row$sample_key) else sample_name

  obj
}

run_go <- function(sig_df, out_prefix) {
  if (nrow(sig_df) < 10) return(NULL)

  ens <- gsub("\\..*$", "", sig_df$gene_id)
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ens),
    keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  )
  entrez <- unique(na.omit(map$ENTREZID))
  if (length(entrez) < 10) return(NULL)

  ego <- enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)

  write.csv(as.data.frame(ego), file.path(results_dir, paste0(out_prefix, "_GO.csv")), row.names = FALSE)
  p <- dotplot(ego, showCategory = 20) + ggtitle(paste0(out_prefix, " GO BP"))
  ggsave(file.path(plots_dir, paste0(out_prefix, "_GO_dotplot.png")), p, width = 10, height = 7, dpi = 300)
  ego
}

safe_pseudobulk_deseq <- function(seu, celltype_name, min_cells_per_sample = 20) {
  sub <- subset(seu, subset = broad_celltype == celltype_name)
  if (ncol(sub) < 50) return(NULL)

  meta <- sub@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell_id")

  sample_counts <- meta %>% count(sample_id, name = "n_cells")
  keep_samples <- sample_counts %>% filter(n_cells >= min_cells_per_sample) %>% pull(sample_id)
  sub <- subset(sub, subset = sample_id %in% keep_samples)
  if (ncol(sub) < 50) return(NULL)

  agg <- AggregateExpression(
    sub,
    group.by = "sample_id",
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA

  pb_meta <- sub@meta.data %>%
    as.data.frame() %>%
    distinct(sample_id, patient_id, timepoint) %>%
    filter(sample_id %in% colnames(agg)) %>%
    arrange(match(sample_id, colnames(agg)))

  rownames(pb_meta) <- pb_meta$sample_id
  pb_meta <- pb_meta[colnames(agg), , drop = FALSE]

  pb_meta$timepoint <- factor(pb_meta$timepoint, levels = c("Pre", "Post"))
  pb_meta$patient_id <- factor(pb_meta$patient_id)

  if (nlevels(droplevels(pb_meta$timepoint)) < 2) return(NULL)
  if (nrow(pb_meta) < 4) return(NULL)

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(agg)),
    colData = pb_meta,
    design = ~ patient_id + timepoint
  )

  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  if (nrow(dds) < 10) return(NULL)

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("timepoint", "Post", "Pre"))
  res_df <- as.data.frame(res) %>% rownames_to_column("gene_id")
  res_df$celltype <- celltype_name

  write.csv(res_df, file.path(results_dir, paste0("GSE200107_pseudobulk_", celltype_name, "_Post_vs_Pre.csv")), row.names = FALSE)

  sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  write.csv(sig, file.path(results_dir, paste0("GSE200107_pseudobulk_", celltype_name, "_sig.csv")), row.names = FALSE)
  if (nrow(sig) >= 10) run_go(sig, paste0("GSE200107_pseudobulk_", celltype_name, "_Post_vs_Pre"))

  list(dds = dds, res = res_df, sig = sig)
}

# -----------------------------
# 3) Download GEO metadata
# -----------------------------
message("Downloading GEO metadata for ", accession)
geo_obj <- getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
eset <- geo_obj[[1]]
meta <- pData(eset) |> as.data.frame() |> rownames_to_column("geo_row")

meta$title <- as.character(meta$title)
meta$geo_accession <- as.character(meta$geo_accession)
meta$assay <- parse_assay(meta$title)
meta$timepoint <- parse_timepoint(meta$title)
meta$patient_id <- parse_patient_id(meta$title)
meta$sample_key <- clean_text(meta$title)

expression_meta <- meta %>% filter(assay == "Expression")
vdj_meta <- meta %>% filter(assay == "VDJ")

write.csv(meta, file.path(meta_dir, paste0(accession, "_sample_metadata_all.csv")), row.names = FALSE)
write.csv(expression_meta, file.path(meta_dir, paste0(accession, "_expression_metadata.csv")), row.names = FALSE)
write.csv(vdj_meta, file.path(meta_dir, paste0(accession, "_vdj_metadata.csv")), row.names = FALSE)

message("Expression samples: ", nrow(expression_meta))
message("VDJ samples: ", nrow(vdj_meta))

# -----------------------------
# 4) Download and locate raw files
# -----------------------------
message("Downloading supplementary files for ", accession)
getGEOSuppFiles(accession, makeDirectory = FALSE, baseDir = project_dir)

raw_tar <- list.files(project_dir, pattern = "GSE200107_RAW\\.tar$", full.names = TRUE)
if (length(raw_tar) == 0) raw_tar <- list.files(project_dir, pattern = "\\.tar$", full.names = TRUE)
if (length(raw_tar) == 0) stop("Could not find raw tar file in ", project_dir)
untar(raw_tar[1], exdir = raw_dir)

expr_mtx <- list.files(raw_dir, pattern = "_matrix\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE)
expr_mtx <- expr_mtx[!stringr::str_detect(expr_mtx, "contig_annotations")]
if (length(expr_mtx) == 0) stop("No expression matrix files found after extracting the tarball.")

expr_file_info <- tibble(file = expr_mtx) %>%
  mutate(
    basename = basename(file),
    prefix = stringr::str_remove(basename, "_matrix\\.mtx\\.gz$"),
    gsm = stringr::str_extract(prefix, "^GSM[0-9]+")
  ) %>%
  inner_join(expression_meta, by = c("gsm" = "geo_accession")) %>%
  mutate(
    geo_accession = gsm,
    sample_id = gsm
  )

if (nrow(expr_file_info) == 0) stop("Could not match expression matrix files to GEO accessions.")
write.csv(expr_file_info, file.path(meta_dir, paste0(accession, "_matrix_file_lookup.csv")), row.names = FALSE)
message("Matched expression matrices: ", nrow(expr_file_info))

# -----------------------------
# 5) Build Seurat objects
# -----------------------------
seu_list <- list()
qc_rows <- list()

for (i in seq_len(nrow(expr_file_info))) {
  row <- expr_file_info[i, ]
  prefix <- file.path(dirname(row$file), row$prefix)
  sample_label <- as.character(row$geo_accession)
  sample_label <- gsub("[^A-Za-z0-9_]+", "_", sample_label)

  message("Reading sample: ", sample_label, " | ", row$title)
  obj <- make_seurat_from_prefix(prefix, row)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$sample_id <- sample_label
  obj$geo_accession <- sample_label

  seu_list[[sample_label]] <- obj
  qc_rows[[sample_label]] <- data.frame(
    sample_id = sample_label,
    n_cells = ncol(obj),
    median_features = median(obj$nFeature_RNA),
    median_counts = median(obj$nCount_RNA),
    median_percent_mt = median(obj$percent.mt)
  )
}

qc_summary <- bind_rows(qc_rows)
write.csv(qc_summary, file.path(results_dir, paste0(accession, "_sample_qc_summary.csv")), row.names = FALSE)

message("Merging Seurat objects")
combined <- Reduce(function(x, y) merge(x, y), seu_list)
combined$sample_id <- factor(combined$sample_id)
combined$timepoint <- factor(combined$timepoint, levels = c("Pre", "Post", "Unknown"))
combined$patient_id <- factor(combined$patient_id)

p_qc <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave(file.path(plots_dir, paste0(accession, "_QC_violin.png")), p_qc, width = 12, height = 4, dpi = 300)

# -----------------------------
# 6) Standard Seurat processing
# -----------------------------
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, features = rownames(combined))
combined <- RunPCA(combined, features = VariableFeatures(combined))
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

saveRDS(combined, file.path(results_dir, paste0(accession, "_combined_seurat_prejoin.rds")))

p_umap_tp <- DimPlot(combined, reduction = "umap", group.by = "timepoint") + ggtitle("UMAP by timepoint")
p_umap_sample <- DimPlot(combined, reduction = "umap", group.by = "sample_id") + ggtitle("UMAP by sample")
p_umap_cluster <- DimPlot(combined, reduction = "umap", label = TRUE) + ggtitle("UMAP clusters")

ggsave(file.path(plots_dir, paste0(accession, "_UMAP_timepoint.png")), p_umap_tp, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, paste0(accession, "_UMAP_sample.png")), p_umap_sample, width = 10, height = 7, dpi = 300)
ggsave(file.path(plots_dir, paste0(accession, "_UMAP_clusters.png")), p_umap_cluster, width = 8, height = 6, dpi = 300)

# -----------------------------
# 7) JoinLayers + marker discovery
# -----------------------------
message("Joining layers for marker analysis")
combined <- JoinLayers(combined)
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "seurat_clusters"
combined <- droplevels(combined)

markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.10
)
write.csv(markers, file.path(results_dir, paste0(accession, "_cluster_markers_joined.csv")), row.names = FALSE)

if (nrow(markers) > 0) {
  top_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5)
  write.csv(top_markers, file.path(results_dir, paste0(accession, "_top_markers_per_cluster.csv")), row.names = FALSE)
}

# -----------------------------
# 8) Broad cell-type annotation
#    based on cluster-average marker profiles
# -----------------------------
marker_sets <- list(
  Tcell = c("CD3D", "CD3E", "TRAC"),
  Th2 = c("GATA3", "IL4", "IL13", "CCR4", "PTGDR2"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"),
  NK_Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "CTSW"),
  Bcell = c("MS4A1", "CD79A", "CD74", "BANK1"),
  Plasma = c("MZB1", "XBP1", "JCHAIN", "SDC1"),
  Mono = c("LYZ", "LST1", "FCN1", "S100A8", "S100A9"),
  DC = c("FCER1A", "CLEC10A", "ITGAX", "CST3"),
  Prolif = c("MKI67", "TOP2A", "STMN1")
)

avg_expr <- AverageExpression(combined, group.by = "seurat_clusters", assays = "RNA", slot = "data")$RNA
score_mat <- sapply(marker_sets, function(gs) {
  present <- intersect(gs, rownames(avg_expr))
  if (length(present) == 0) return(rep(NA_real_, ncol(avg_expr)))
  colMeans(avg_expr[present, , drop = FALSE])
})
score_mat <- as.data.frame(score_mat)
score_mat$cluster <- colnames(avg_expr)

assign_cluster_type <- function(x) {
  vals <- as.numeric(x[setdiff(names(x), "cluster")])
  names(vals) <- setdiff(names(x), "cluster")
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return("Other")
  best <- names(vals)[which.max(vals)]
  if (max(vals) < 0.05) return("Other")
  if (best == "Tcell") return("T_cell")
  if (best == "Th2") return("Th2_like")
  if (best == "Treg") return("Treg")
  if (best == "NK_Cytotoxic") return("NK_cytotoxic")
  if (best == "Bcell") return("B_cell")
  if (best == "Plasma") return("Plasma")
  if (best == "Mono") return("Monocyte")
  if (best == "DC") return("Dendritic_cell")
  if (best == "Prolif") return("Proliferating")
  "Other"
}

cluster_map <- score_mat %>% rowwise() %>% mutate(broad_celltype = assign_cluster_type(cur_data())) %>% ungroup()
write.csv(cluster_map, file.path(results_dir, paste0(accession, "_cluster_annotation_auto.csv")), row.names = FALSE)

cluster_to_type <- setNames(cluster_map$broad_celltype, cluster_map$cluster)
combined$broad_celltype <- unname(cluster_to_type[as.character(combined$seurat_clusters)])
combined$broad_celltype[is.na(combined$broad_celltype)] <- "Other"
combined$broad_celltype <- factor(combined$broad_celltype)

p_umap_type <- DimPlot(combined, reduction = "umap", group.by = "broad_celltype", label = TRUE) + ggtitle("UMAP broad cell types")
ggsave(file.path(plots_dir, paste0(accession, "_UMAP_broad_celltypes.png")), p_umap_type, width = 10, height = 7, dpi = 300)

# -----------------------------
# 9) Module scores for biology of interest
# -----------------------------
modules <- list(
  Th2_Score = c("GATA3", "IL4", "IL5", "IL13", "CCR4", "PTGDR2"),
  Treg_Score = c("FOXP3", "IL2RA", "CTLA4", "IL10", "TIGIT", "IKZF2"),
  IFN_Score = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS3"),
  Cytotoxic_Score = c("NKG7", "GNLY", "PRF1", "GZMB", "CTSW"),
  Myeloid_Score = c("LYZ", "LST1", "FCN1", "S100A8", "S100A9"),
  Bcell_Score = c("MS4A1", "CD79A", "CD74", "BANK1")
)

for (nm in names(modules)) {
  present <- intersect(modules[[nm]], rownames(combined))
  if (length(present) >= 3) {
    combined <- AddModuleScore(combined, features = list(present), name = nm)
  }
}

score_cols <- grep("Score1$", colnames(combined@meta.data), value = TRUE)
if (length(score_cols) > 0) {
  p_scores <- VlnPlot(combined, features = score_cols, group.by = "timepoint", pt.size = 0.05, ncol = 2)
  ggsave(file.path(plots_dir, paste0(accession, "_module_scores_timepoint.png")), p_scores, width = 12, height = 10, dpi = 300)
}

saveRDS(combined, file.path(results_dir, paste0(accession, "_combined_seurat_annotated.rds")))


# -----------------------------
# 10) Pseudobulk pre vs post SLIT by cell type
# -----------------------------
safe_pseudobulk_deseq <- function(seu, ct, min_cells_per_sample = 20) {

  # helper: flatten list columns safely
  as_chr1 <- function(x) {
    if (is.list(x)) {
      purrr::map_chr(x, ~ as.character(.x)[1])
    } else {
      as.character(x)
    }
  }

  meta_pb <- seu@meta.data %>%
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

  meta_ct <- meta_pb %>%
    filter(broad_celltype == ct)

  if (nrow(meta_ct) == 0) return(NULL)

  # keep only patients with both Pre and Post
  keep_patients <- meta_ct %>%
    dplyr::count(patient_id, timepoint, name = "n_cells") %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(n_timepoints = dplyr::n_distinct(timepoint), .groups = "drop") %>%
    filter(n_timepoints == 2) %>%
    pull(patient_id)

  meta_ct <- meta_ct %>%
    filter(patient_id %in% keep_patients)

  if (nrow(meta_ct) == 0) return(NULL)

  # use the RNA counts layer
  expr_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
  expr_mat <- as.matrix(expr_mat)

  common_cells <- intersect(colnames(expr_mat), rownames(meta_ct))
  if (length(common_cells) == 0) return(NULL)

  expr_mat <- expr_mat[, common_cells, drop = FALSE]
  meta_ct <- meta_ct[common_cells, , drop = FALSE]

  # create pseudobulk groups: patient + timepoint
  meta_ct$group_id <- paste(meta_ct$patient_id, meta_ct$timepoint, sep = "__")

  # require minimum number of cells per pseudobulk sample
  group_sizes <- table(meta_ct$group_id)
  keep_groups <- names(group_sizes)[group_sizes >= min_cells_per_sample]
  meta_ct <- meta_ct %>% filter(group_id %in% keep_groups)

  if (nrow(meta_ct) == 0) return(NULL)

  # rebuild after filtering
  expr_mat <- expr_mat[, rownames(meta_ct), drop = FALSE]

  group_levels <- unique(meta_ct$group_id)

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

  # need at least 2 samples per patient pair
  if (ncol(pb_mat) < 4) return(NULL)

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
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id")

  sig <- res_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

  safe_name <- gsub("[^A-Za-z0-9]+", "_", ct)

  write.csv(
    res_df,
    file.path(results_dir, paste0(accession, "_", safe_name, "_Post_vs_Pre.csv")),
    row.names = FALSE
  )

  write.csv(
    sig,
    file.path(results_dir, paste0(accession, "_", safe_name, "_Post_vs_Pre_sig.csv")),
    row.names = FALSE
  )

  return(list(res = res_df, sig = sig, dds = dds))
}

# -----------------------------
# run pseudobulk DE
# -----------------------------
library(Matrix)
library(DESeq2)

expr_mat <- GetAssayData(combined, assay = "RNA", layer = "counts")

# match cells
common_cells <- intersect(colnames(expr_mat), rownames(meta_pb))
expr_mat <- expr_mat[, common_cells, drop = FALSE]
meta_pb <- meta_pb[common_cells, , drop = FALSE]

pb_results <- list()
summary_rows <- list()

for (ct in celltypes_to_test) {
  message("Pseudobulk testing: ", ct)

  cells_ct <- rownames(meta_pb)[meta_pb$broad_celltype == ct]
  if (length(cells_ct) < 50) next

  meta_ct <- meta_pb[cells_ct, , drop = FALSE]

  # keep paired patients only
  keep_patients <- meta_ct %>%
    count(patient_id, timepoint) %>%
    count(patient_id) %>%
    filter(n == 2) %>%
    pull(patient_id)

  meta_ct <- meta_ct %>% filter(patient_id %in% keep_patients)
  cells_ct <- rownames(meta_ct)

  if (length(cells_ct) < 30) next

  # pseudobulk aggregation
  groups <- paste(meta_ct$patient_id, meta_ct$timepoint, sep = "__")
  group_levels <- unique(groups)

  pb_list <- lapply(group_levels, function(g) {
    cell_ids <- rownames(meta_ct)[groups == g]
    Matrix::rowSums(expr_mat[, cell_ids, drop = FALSE])
  })

  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- group_levels

  pb_meta <- data.frame(
    sample = group_levels,
    patient_id = sub("__.*$", "", group_levels),
    timepoint = sub("^.*__", "", group_levels),
    row.names = group_levels
  )
  pb_meta$timepoint <- factor(pb_meta$timepoint, levels = c("Pre", "Post"))

  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_mat),
    colData = pb_meta,
    design = ~ patient_id + timepoint
  )

  keep_genes <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_genes, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("timepoint", "Post", "Pre"))
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id")

  sig <- res_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

  # -----------------------------
  # SAVE FILES (IMPORTANT)
  # -----------------------------
  safe_name <- gsub("[^A-Za-z0-9]+", "_", ct)

  write.csv(
    res_df,
    file.path(results_dir, paste0(accession, "_", safe_name, "_Post_vs_Pre.csv")),
    row.names = FALSE
  )

  write.csv(
    sig,
    file.path(results_dir, paste0(accession, "_", safe_name, "_Post_vs_Pre_sig.csv")),
    row.names = FALSE
  )

  # store summary
  summary_rows[[ct]] <- data.frame(
    celltype = ct,
    n_cells = length(cells_ct),
    n_genes_tested = nrow(res_df),
    n_sig = nrow(sig),
    stringsAsFactors = FALSE
  )

  pb_results[[ct]] <- list(res = res_df, sig = sig)
}

# -----------------------------
# SAVE SUMMARY
# -----------------------------
if (length(summary_rows) > 0) {
  summary_df <- dplyr::bind_rows(summary_rows)

  write.csv(
    summary_df,
    file.path(results_dir, paste0(accession, "_pseudobulk_summary.csv")),
    row.names = FALSE
  )
}

saveRDS(pb_results, file.path(results_dir, paste0(accession, "_pseudobulk_results.rds")))

cat("\n✅ Pseudobulk analysis COMPLETE\n")

# -----------------------------
# 11) Optional VDJ file index
# -----------------------------
vdj_files <- list.files(raw_dir, pattern = "contig_annotations.*\\.csv$|contig_annotations.*\\.csv.gz$", recursive = TRUE, full.names = TRUE)
if (length(vdj_files) > 0) {
  vdj_summary <- data.frame(file = vdj_files, sample = basename(vdj_files), stringsAsFactors = FALSE)
  write.csv(vdj_summary, file.path(results_dir, paste0(accession, "_vdj_file_index.csv")), row.names = FALSE)
}

# -----------------------------
# 12) Final outputs
# -----------------------------
message("Done.")
message("Key files saved in: ", project_dir)
