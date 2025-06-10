Package: TranscriptoPathR
Title: A Comprehensive Shiny App for RNA-Seq Analysis
Version: 0.1.0
Authors@R: 
    person(Nuri Alpay, Temiz, email = temizna@umn.edu, role = c("aut", "cre"))
Description: 
    TransriptoPathR provides an interactive Shiny interface for complete RNA-Seq analysis,
    including data upload, normalization, quality control, differential expression,
    gene expression visualization, pathway analysis, and GSEA using CRAN and Bioconductor tools.
License: MIT
Encoding: UTF-8
LazyData: true
RoxygenNote: 7.2.3

Depends:
    R (>= 4.1.0)

Imports:
    shiny,
    shinythemes,
    ggplot2,
    dplyr,
    tidyr,
    stringr,
    readr,
    tibble,
    reshape2,
    DT,
    RColorBrewer,
    DESeq2,
    GEOquery,
    Biobase,
    clusterProfiler,
    enrichplot,
    ComplexHeatmap,
    msigdbr,
    ReactomePA,
    pathview,
    pathfindR,
    AnnotationDbi,
    org.Hs.eg.db,
    org.Mm.eg.db

Suggests:
    testthat,
    knitr,
    rmarkdown

VignetteBuilder: knitr

---

# TranscriptoPathR

**TransriptoPathR** is a comprehensive Shiny-based application for end-to-end RNA-Seq data analysis designed for bench scientists with no or minimal coding experience. It provides an interactive GUI for both novice and advanced users to perform quality control, normalization, differential expression, pathway analysis, and gene set enrichment analysis (GSEA) with no coding experience.

## ✨ Features

- **Data Input**  
  - Upload RNA-Seq raw count matrix and design metadata file (CSV or XLSX)  

- **Sample Selection**  
  - Interactive filtering of samples for downstream analysis

- **Gene Expression Visualization**  
  - Boxplots of selected genes grouped by metadata categories

- **Quality Control (QC)**  
  - PCA, sample distance heatmaps, mean-variance plots, variance histograms
- **Genes of interest Heatmap Visualization**
  - Gene expression heatmap of uploaded genes of interest grouped by metadata categories

- **Differential Expression Analysis**  
  - DESeq2-based analysis with customizable thresholds and conditions
  - Users can download DE table.
- **Heatmap**  
  - Top DE genes with hierarchical clustering and group annotations

- **Volcano and MA Plots**  
  - Visualization of differential expression results with gene labeling

- **Cross Plot**  
  - Compare DE results between two separete comparisons or experiments
  - A cross plot, Venn diagrams, heatmaps and pathway comparison is calcuated and plotted.

- **Pathway Analysis**  
  - GO, KEGG, Reactome enrichment using clusterProfiler  
  - Visualizations: dot plots, cnet, circular, emap  
  - KEGG pathway rendering with `pathfindR` + gene heatmaps
  - Pathway results can be downloaded as tables.
- **Non-overlap Pathway Analysis**  
  - Genes overlapping in the initial pathway analysis is removed and the same pathway analysis is repeated with this refined geneset.  
  - Visualizations: dot plots, cnet, circular, emap

- **GSEA (Gene Set Enrichment Analysis)**  
  - Supports MSigDB collections: Hallmark, GO, KEGG, Reactome  
  - Dot plots, enrichment plots,  and enrichment tables

- **Transcirption Factor Enrichment Analysis**
  -TRANSFAC and JASPAR, ENCODE and ChEA, TRUSST, hTFtarget and TFLink databases can be used to test for transciprtion factor enrichment
  - Dot plots and   ridgeplots
- **Cancer Gene Census Comparison**                     
  - Up and down regulated genes are compared to Cancer Gene Census
  - Venn Diagram of overlaps and table of overlapping genes
- **Dimension Reduction and Clustering Analysis**                     
  - A preliminary Principle Component Analysis to help user understand the underlying clusters and their relations within the dataset.
  - This is a starting analysis for more comprehensive consensus clustering or NMF analysis.
  - Reconstructed data heatmap, sample- sample correlation heatmap showing possible clusters, expression heatmap of genes contributing to components (rudimentary meta-gene analysis), and pathway contributing to components
- **User Session Log**                     
  - Inputs and actions taken by the user in order. Session log can be saved.
 

## 🧬 Supported Species

- Homo sapiens
- Mus musculus

## 📦 Installation

Install required CRAN and Bioconductor dependencies before installing TransriptoPathR:

```r
# Install Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
  "DESeq2", "GEOquery", "Biobase", "clusterProfiler", "enrichplot",
  "ComplexHeatmap", "msigdbr", "ReactomePA", "pathview",
  "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi"
))

# Install pathfindR (if not already)
install.packages("pathfindR")
install.packages("devtools")
library(devtools) 
# Clone or download this repo, then from root directory:
install_github("temizna/TranscriptoPathR")
```

## 🚀 Running the App

After installation, launch the app from R or RStudio:

```r
library(TranscriptoPathR)
TranscriptoPathR::run_app()
```

This will open the app in your default web browser.

## 📂 File Requirements

**Counts file**
Counts file can be any raw count file including raw counts from The Cancer Genome Atlas or any Gene Expression Omnibus database submission.
TPM, FPKM, and CPM are NOT supported.  
- Rows = genes, Columns = sample names  
- Gene identifiers: gene symbols or Ensembl IDs  
- Format: `.csv` or `.xlsx`

**Design file**  
Design file descirbes the experiment and the conditions/genotypes/treatments, etc used. 
- Rows = sample names (must match columns in count matrix)  
- Columns = metadata variables (e.g., `condition`, `batch`)  
- Format: `.csv` or `.xlsx`


## 🛠 Development

All modules are structured as reusable Shiny modules. Additional features can be added easily. Utility functions are separated into a `utils.R` file.

To ensure all required Bioconductor packages are installed:

```r
# R/zzz.R
.onLoad <- function(libname, pkgname) {
  required_bioc <- c(
    "DESeq2", "GEOquery", "Biobase", "clusterProfiler", "enrichplot",
    "ComplexHeatmap", "msigdbr", "ReactomePA", "pathview",
    "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi"
  )

  missing_pkgs <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    message("Installing missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
  }
}
```

Also ensure the following are loaded in your `run_app()` function:

```r
run_app <- function() {
  options(shiny.maxRequestSize = 250 * 1024^2)

  library(shiny)
  library(shinythemes)
  library(DESeq2)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(GEOquery)
  library(pathview)
  library(pathfindR)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(msigdbr)

  shiny::shinyApp(
    ui = TransriptoPathR::app_ui(),
    server = TransriptoPathR::app_server
  )
}
```


## 📝 License

MIT © Nuri Alpay Temiz

