# Multi-Omics & scRNA-seq Analysis of Cellular Senescence

## Overview
This repository contains R scripts, processed data tables, and publication-ready visualizations investigating molecular signatures, pathway activation, and transcriptional dynamics across cellular senescence stages (Proliferating → Early Senescent → SASP Senescent).

## Key Analyses Included
* **Phospho-Site & Pathway Activation:** Multi-pathway profiling covering DNA damage/p53, cell cycle pocket proteins, mTORC1/translational SASP, and cGAS-STING/interferon signaling.
* **Single-Cell & Trajectory Analysis:** UMAP dimensionality reduction, feature plots, dot plots of pathway markers, and Slingshot pseudotime trajectory inference.
* **Differential Expression & Effect Sizes:** Volcano plot profiling (overexpressed and underexpressed gene sets), Venn comparisons, and forest plots with 95% confidence intervals.
* **Flow Cytometry Visualization:** FACS pseudocolor density profiling.

## Repository Contents
* `scrna senesence.R` — Main R script used for data processing, statistical modeling, and figure generation.
* `00_Project_File_Inventory.csv` — Master index of generated output files.
* `10_*.png` — Pathway marker dot plots and FACS pseudocolor plots.
* `11_*.csv` & `11_*.png` — Differential expression tables (all, overexpressed, and underexpressed genes) and UMAP/volcano visualizations.
* `12_*.csv` & `12_*.png` — Forest plot effect size tables (95% CI), Venn diagrams, and Slingshot pseudotime trajectory plots.

## Requirements & R Packages
To reproduce the figures and tables in RStudio, install the following packages:
```r
install.packages(c("ggplot2", "ggrepel", "pheatmap", "dplyr", "patchwork"))
```

## How to Use
1. Clone or download this repository.
2. Open `scrna senesence.R` in **RStudio**.
3. Run the script to reproduce the `.csv` summary tables and `.png` figures
