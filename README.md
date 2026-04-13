# Pollen Allergy SLIT Transcriptomics Project

This repository contains the analysis workflow for a multi-scale transcriptomics study of sublingual immunotherapy (SLIT) in pollen allergy.

## Overview

The project integrates:

- Bulk RNA-seq analysis for **GSE206149**
- Single-cell RNA-seq / VDJ analysis for **GSE200107**
- A restart/helper script for the scRNA workflow
- A validation script using **GSE206152**

The analysis is organized into separate scripts and output folders so the workflow is easy to rerun and review.

## Repository structure

```text
.
├── Step1_GSE206149_bulk_rnaseq_focused_v6.R
├── Step2_GSE200107_scRNA_SLIT_pre_post_v1.R
├── Step3_GSE200107_merged_restart_and_helper.R
├── Step4_gse_206149_gse_206152_validation.r
├── figures/
├── results/
├── metadata/
├── plots/
└── README.md
```

## What each script does

**Step1_GSE206149_bulk_rnaseq_focused_v6.R**  
Processes bulk RNA-seq for GSE206149, runs differential expression and pathway analysis, and writes output tables and plots into project folders.

**Step2_GSE200107_scRNA_SLIT_pre_post_v1.R**  
Processes the scRNA-seq / VDJ dataset, builds the Seurat object, performs clustering and annotation, and runs paired pseudobulk Post vs Pre SLIT analysis.

**Step3_GSE200107_merged_restart_and_helper.R**  
Reloads the saved Seurat object, joins layers for Seurat v5 compatibility, regenerates markers, applies manual annotation, and reruns paired pseudobulk analysis.

**Step4_gse_206149_gse_206152_validation.r**  
Validates the GSE206149 SLIT signature in GSE206152 by testing overlap, concordance, pathway consistency, and signature scoring.

## Expected outputs

The scripts generate:

- `metadata/` for parsed GEO metadata and sample mapping tables
- `results/` for DE tables, summary tables, RDS files, and intermediate outputs
- `plots/` for QC plots, volcano plots, heatmaps, and validation graphics
- `figures/` for final report-ready panels

## How to run

Run each script from the repository root so the relative paths resolve correctly:

```r
source("Step1_GSE206149_bulk_rnaseq_focused_v6.R")
source("Step2_GSE200107_scRNA_SLIT_pre_post_v1.R")
source("Step3_GSE200107_merged_restart_and_helper.R")
source("Step4_gse_206149_gse_206152_validation.r")
```

## Notes

- Large raw GEO downloads are not committed to the repository.
- Generated analysis artifacts can be kept if you want a fully reproducible archive, but consider excluding very large raw files.
- The figures in `figures/` are the final report-ready panels.
