# Pathway-Level Immune Signatures Associated with Sublingual Immunotherapy in Allergic Rhinitis

This repository contains the R code used for the analyses and figure generation reported in:

**Pathway-Level Immune Signatures Associated with Sublingual Immunotherapy in Allergic Rhinitis**

**Christopher Wang**  
Montgomery Blair High School, USA

*The National High School Journal of Science*, 2026  
Received: May 5, 2026  
Accepted: July 30, 2026  
Electronic access: September 30, 2026

---

## Overview

This study investigates transcriptomic patterns associated with sublingual immunotherapy (SLIT) for allergic rhinitis using complementary bulk and single-cell RNA-sequencing datasets.

The analysis integrates:

- **Bulk RNA-seq:** GSE206149
- **Single-cell RNA-seq:** GSE200107
- Differential expression analysis
- Gene set enrichment analysis (GSEA)
- Patient-level pseudobulk differential expression
- Cell-type-specific pathway-level analysis

The central finding is that SLIT-associated transcriptomic differences are more apparent at the **pathway level** than at the level of individual genes.

---

## Datasets

### GSE206149 — Bulk RNA-seq

The bulk RNA-seq dataset contains **255 samples** and was used to compare transcriptional patterns between SLIT-treated and placebo samples.

Analyses included:

- DESeq2 normalization
- Differential expression analysis
- Gene ranking by log2 fold change
- Gene Ontology Biological Process enrichment using GSEA

Only a small number of genes showed large gene-level changes, while multiple immune-related pathways showed significant enrichment.

### GSE200107 — Single-cell RNA-seq

The single-cell dataset contains **15 expression libraries from 7 patients** with paired pre- and post-treatment samples.

Major immune-cell populations analyzed included:

- CD4+ T cells
- CD8+ T cells
- B cells
- Monocytes
- Dendritic cells
- Natural killer cells
- Plasma cells

A patient-level pseudobulk framework was used to assess cell-type-specific transcriptional responses.

---

## Repository Structure

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
