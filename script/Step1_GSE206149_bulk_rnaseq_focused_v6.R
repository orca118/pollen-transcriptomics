#!/usr/bin/env Rscript

# =========================================================
# GSE206149 focused bulk RNA-seq analysis (v6)
# Primary contrast: SLIT vs Placebo
# Secondary contrast: Placebo vs Healthy
# Optional contrast: SCIT vs Placebo
#
# Fixes vs earlier version:
#   - more careful GEO metadata parsing
#   - more robust count-file detection
#   - explicit orientation checks (genes x samples)
#   - explicit sample matching and zero-count diagnostics
#   - stops early with informative messages if mapping is wrong
# =========================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Packages
# -----------------------------
pkgs_cran <- c("dplyr", "tibble", "stringr", "readr", "ggplot2", "tidyr", "purrr")
pkgs_bioc <- c(
  "GEOquery", "DESeq2", "EnhancedVolcano", "pheatmap",
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
accession <- "GSE206149"
project_dir <- file.path(getwd(), accession)
meta_dir <- file.path(project_dir, "metadata")
results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "plots")
for (d in c(project_dir, meta_dir, results_dir, plots_dir)) {
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
  NA_character_
}

infer_group_from_row <- function(x) {
  blob <- paste(x, collapse = " | ")
  blob <- tolower(blob)
  if (str_detect(blob, "healthy|non-allergic|normal control")) return("Healthy")
  if (str_detect(blob, "placebo")) return("Placebo")
  if (str_detect(blob, "scit|subcutaneous")) return("SCIT")
  if (str_detect(blob, "slit|sublingual")) return("SLIT")
  NA_character_
}

read_count_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    # keep ext as inner extension when gzipped
    inner <- tolower(sub(".*\\.([A-Za-z0-9]+)\\.gz$", "\\1", path))
    ext <- inner
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

pick_count_file <- function(root_dir) {
  files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(csv|tsv|txt)(\\.gz)?$", files, ignore.case = TRUE)]
  if (length(files) == 0) stop("No candidate count files found in: ", root_dir)

  # Prefer files with names suggesting raw counts / counts matrix.
  priority <- files[grepl("raw_counts|counts|expression|matrix|GSE206149", basename(files), ignore.case = TRUE)]
  if (length(priority) == 0) priority <- files

  # Choose first file that looks like a matrix with at least a few numeric columns.
  for (f in priority) {
    message("Checking candidate file: ", f)
    ok <- tryCatch({
      df <- read_count_table(f)
      if (ncol(df) < 3) return(FALSE)
      # At least one column after the first should be mostly numeric
      num_frac <- sapply(df[-1], function(z) mean(!is.na(suppressWarnings(as.numeric(as.character(z))))))
      any(num_frac > 0.7)
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(f)
  }
  stop("Could not identify a valid count table among candidates.")
}

make_sample_key <- function(title, geo_accession, sample_lib, subject_id) {
  vals <- c(title, geo_accession, sample_lib, subject_id)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return(NA_character_)
  clean_name(vals[1])
}

match_metadata_to_counts <- function(counts, meta) {
  cn <- colnames(counts)
  rn <- rownames(meta)

  # direct match first
  common <- intersect(cn, rn)
  if (length(common) >= 4) {
    counts <- counts[, common, drop = FALSE]
    meta <- meta[common, , drop = FALSE]
    return(list(counts = counts, meta = meta, method = "direct"))
  }

  # Try matching by a flexible sample key.
  meta$sample_key <- mapply(
    make_sample_key,
    meta$title,
    meta$geo_accession,
    if ("sample_lib" %in% names(meta)) meta$sample_lib else NA_character_,
    if ("subject_id" %in% names(meta)) meta$subject_id else NA_character_,
    USE.NAMES = FALSE
  )

  # If counts columns match sample_key, use that.
  common2 <- intersect(clean_name(cn), meta$sample_key)
  if (length(common2) >= 4) {
    map <- setNames(meta$geo_accession[match(common2, meta$sample_key)], common2)
    names(map) <- common2
    names_idx <- match(clean_name(cn), common2)
    keep <- !is.na(names_idx)
    counts <- counts[, keep, drop = FALSE]
    colnames(counts) <- map[clean_name(colnames(counts))[keep]]
    meta <- meta[colnames(counts), , drop = FALSE]
    rownames(meta) <- meta$geo_accession
    return(list(counts = counts, meta = meta, method = "sample_key"))
  }

  # Fallback: look for title substrings in count column names.
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

  stop(
    "Could not match count columns to metadata. ",
    "Please inspect colnames(counts) and meta$geo_accession / meta$title."
  )
}

orient_if_needed <- function(counts, meta) {
  # Desired: genes in rows, samples in columns.
  # If columns look like genes and rows look like samples, transpose.
  cn <- colnames(counts)
  rn <- rownames(counts)

  sample_hits_cols <- sum(clean_name(cn) %in% clean_name(meta$geo_accession), na.rm = TRUE)
  sample_hits_rows <- sum(clean_name(rn) %in% clean_name(meta$geo_accession), na.rm = TRUE)

  if (sample_hits_rows > sample_hits_cols) {
    counts <- t(counts)
    message("Transposed count matrix because rows matched sample IDs better than columns.")
  }

  # Ensure numeric matrix.
  counts <- as.matrix(counts)
  mode(counts) <- "numeric"

  # Remove rows/cols with all NA if any.
  counts[is.na(counts)] <- 0
  counts
}

run_contrast <- function(dds, contrast_vec, out_prefix) {
  res <- results(dds, contrast = contrast_vec)
  res_df <- as.data.frame(res) |> rownames_to_column("gene_id")
  write.csv(res_df, file.path(results_dir, paste0(out_prefix, ".csv")), row.names = FALSE)

  safe_png(file.path(plots_dir, paste0(out_prefix, "_volcano.png")), {
    print(
      EnhancedVolcano(
        res_df,
        lab = ifelse(is.na(res_df$padj), "", res_df$gene_id),
        x = "log2FoldChange",
        y = "pvalue",
        title = out_prefix,
        pCutoff = 0.05,
        FCcutoff = 1
      )
    )
  })

  sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  write.csv(sig, file.path(results_dir, paste0(out_prefix, "_sig.csv")), row.names = FALSE)
  sig
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

# -----------------------------
# 3) Load GEO metadata
# -----------------------------
message("Loading GEO metadata...")
geo_obj <- getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
eset <- geo_obj[[1]]
meta <- pData(eset) |> as.data.frame() |> rownames_to_column("geo_row")
meta <- collapse_characteristics(meta)

# Core labels
meta$group <- NA_character_
for (i in seq_len(nrow(meta))) {
  vals <- unlist(meta[i, ], use.names = FALSE)
  meta$group[i] <- infer_group_from_row(vals)
}

# Use more specific fields if present
if ("treatment" %in% names(meta)) {
  meta$group <- ifelse(is.na(meta$group), vapply(meta$treatment, normalize_group, character(1)), meta$group)
}
if ("treatment_raw" %in% names(meta)) {
  meta$group <- ifelse(is.na(meta$group), vapply(meta$treatment_raw, normalize_group, character(1)), meta$group)
}

meta$tissue <- case_when(
  str_detect(tolower(paste(meta$title, meta$source_name_ch1, sep = " | ")), "nasal") ~ "Nasal",
  str_detect(tolower(paste(meta$title, meta$source_name_ch1, sep = " | ")), "pbmc|blood") ~ "PBMC",
  TRUE ~ "Unknown"
)

meta$group <- factor(meta$group, levels = c("Healthy", "Placebo", "SCIT", "SLIT"))
meta$tissue <- factor(meta$tissue, levels = c("Nasal", "PBMC", "Unknown"))
rownames(meta) <- meta$geo_accession

write.csv(meta, file.path(meta_dir, paste0(accession, "_metadata_parsed.csv")), row.names = FALSE)
message("Group counts:")
print(table(meta$group, useNA = "ifany"))
print(table(meta$tissue, useNA = "ifany"))

# -----------------------------
# 4) Download and locate count file
# -----------------------------
message("Downloading supplementary files...")
getGEOSuppFiles(accession, makeDirectory = FALSE, baseDir = project_dir)
count_path <- pick_count_file(project_dir)
message("Selected count file: ", count_path)

# -----------------------------
# 5) Read and normalize count matrix
# -----------------------------
counts_df <- read_count_table(count_path)
message("Count table dimensions (raw): ", nrow(counts_df), " x ", ncol(counts_df))

# First column is typically gene identifier.
gene_col <- names(counts_df)[1]
counts_df[[gene_col]] <- gsub("\\..*$", "", as.character(counts_df[[gene_col]]))
counts <- counts_df %>% column_to_rownames(gene_col) %>% as.data.frame()

# Keep only numeric-ish columns.
counts[] <- lapply(counts, function(x) suppressWarnings(as.numeric(as.character(x))))
counts[is.na(counts)] <- 0

# Orient if needed.
counts <- orient_if_needed(counts, meta)

# -----------------------------
# 6) Match samples to metadata and sanity-check
# -----------------------------
matched <- match_metadata_to_counts(counts, meta)
counts <- matched$counts
meta <- matched$meta
message("Sample matching method: ", matched$method)
message("Matched samples: ", ncol(counts))

if (ncol(counts) < 4) {
  stop("Too few matched samples after mapping. Check sample IDs manually.")
}

# Make sure gene-by-sample shape is correct.
if (nrow(counts) < 1000 && ncol(counts) > nrow(counts)) {
  stop("Count matrix still looks transposed after matching. Please inspect the file manually.")
}

# Remove all-zero genes, and drop samples with zero library size.
counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
libsize <- colSums(counts, na.rm = TRUE)
write.csv(
  data.frame(sample = names(libsize), library_size = libsize),
  file.path(results_dir, paste0(accession, "_library_sizes.csv")),
  row.names = FALSE
)

if (all(libsize == 0)) {
  stop("All samples have 0 counts after import. The file is not the correct count matrix or it was read incorrectly.")
}

zero_samples <- names(libsize)[libsize == 0]
if (length(zero_samples) > 0) {
  message("Dropping zero-count samples: ", paste(zero_samples, collapse = ", "))
  counts <- counts[, libsize > 0, drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
}

# Final alignment checks.
stopifnot(identical(colnames(counts), rownames(meta)))
stopifnot(any(colSums(counts) > 0))

# -----------------------------
# 7) Quick QC plot
# -----------------------------
lib_df <- data.frame(sample = colnames(counts), library_size = colSums(counts))
lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample[order(lib_df$library_size)])

p_lib <- ggplot(lib_df, aes(x = sample, y = library_size)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = paste0(accession, " library sizes"), x = "Sample", y = "Counts")

ggsave(file.path(plots_dir, paste0(accession, "_library_sizes.png")), p_lib, width = 8, height = 10, dpi = 300)

# -----------------------------
# 8) DESeq2 analysis
# -----------------------------
# Keep the model simple and hypothesis-driven.
# If both tissue groups are present, adjust for tissue.
design_formula <- ~ group
if (nlevels(droplevels(meta$tissue)) > 1) {
  design_formula <- ~ tissue + group
}

message("Using design: ", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = design_formula
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
saveRDS(dds, file.path(results_dir, paste0(accession, "_dds.rds")))

vsd <- vst(dds, blind = TRUE)
write.csv(as.data.frame(assay(vsd)), file.path(results_dir, paste0(accession, "_vst_matrix.csv")))

pca_df <- plotPCA(vsd, intgroup = intersect(c("group", "tissue"), names(colData(dds))), returnData = TRUE)
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = colData(dds)$group, shape = colData(dds)$tissue)) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(title = paste0(accession, " PCA"), color = "Group", shape = "Tissue")
ggsave(file.path(plots_dir, paste0(accession, "_PCA.png")), p_pca, width = 8, height = 6, dpi = 300)

message("resultsNames(dds):")
print(resultsNames(dds))

# -----------------------------
# 9) Main contrasts
# -----------------------------
contrast_specs <- list()
if (all(c("SLIT", "Placebo") %in% levels(meta$group))) {
  contrast_specs[["SLIT_vs_Placebo"]] <- c("group", "SLIT", "Placebo")
}
if (all(c("Placebo", "Healthy") %in% levels(meta$group))) {
  contrast_specs[["Placebo_vs_Healthy"]] <- c("group", "Placebo", "Healthy")
}
if (all(c("SCIT", "Placebo") %in% levels(meta$group))) {
  contrast_specs[["SCIT_vs_Placebo"]] <- c("group", "SCIT", "Placebo")
}

sig_list <- list()
for (nm in names(contrast_specs)) {
  out_prefix <- paste0(accession, "_", nm)
  sig_list[[nm]] <- run_contrast(dds, contrast_specs[[nm]], out_prefix)
  run_go(sig_list[[nm]], out_prefix)
}

# -----------------------------
# 10) Heatmap for main contrast
# -----------------------------
main_name <- "SLIT_vs_Placebo"
main_sig_path <- file.path(results_dir, paste0(accession, "_", main_name, "_sig.csv"))
if (file.exists(main_sig_path)) {
  sig <- read.csv(main_sig_path)
  if (nrow(sig) >= 5) {
    top_genes <- head(sig$gene_id[order(sig$padj, sig$pvalue)], 30)
    top_genes <- intersect(top_genes, rownames(assay(vsd)))
    if (length(top_genes) >= 2) {
      mat <- assay(vsd)[top_genes, , drop = FALSE]
      safe_png(file.path(plots_dir, paste0(accession, "_SLIT_vs_Placebo_heatmap_top30.png")), {
        pheatmap(mat, scale = "row", show_colnames = FALSE, main = "SLIT vs Placebo")
      })
    }
  }
}

# -----------------------------
# 11) Summary table
# -----------------------------
summary_table <- tibble(
  contrast = c("SLIT_vs_Placebo", "Placebo_vs_Healthy", "SCIT_vs_Placebo"),
  file = c(
    paste0(accession, "_SLIT_vs_Placebo.csv"),
    paste0(accession, "_Placebo_vs_Healthy.csv"),
    paste0(accession, "_SCIT_vs_Placebo.csv")
  )
)
write.csv(summary_table, file.path(results_dir, paste0(accession, "_contrast_summary.csv")), row.names = FALSE)

message("Done.")
message("Outputs written to: ", project_dir)
