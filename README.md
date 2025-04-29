# SPLIT: Spatial Purification of Layered Intracellular Transcripts

🚧 **This package is under active development.**\
For now, use the **Quick Start** guide below to get up and running quickly.\
A comprehensive tutorial and vignettes are coming soon.

------------------------------------------------------------------------

## 📦 Installation

To install SPLIT from GitHub:

``` r
# Install SPLIT
remotes::install_github("bdsc-tds/SPLIT")
```

## 🚀 Quick Start

⚠️ **IMPORTANT:**\
SPLIT currently requires **doublet-mode** RCTD results from the original [spacexr GitHub repository](https://github.com/dmcable/spacexr) or its faster [HD fork](https://github.com/jpromeror/spacexr/tree/HD), **not** from the newly released [Bioconductor version](https://www.bioconductor.org/packages/release/bioc/html/spacexr.html).\
🚧 **Compatibility with Bioconductor's spacexr is coming soon.**

If you already have your **Xenium** dataset as a Seurat object (`xe`) and **RCTD** results from **doublet-mode** decomposition in `RCTD`, you can run SPLIT purification like this:

```{r}
library(SPLIT)
library(spacexr)
library(dplyr)
library(Seurat)
library(ggplot2)

# Post-process RCTD output
RCTD <- SPLIT::run_post_process_RCTD(RCTD)

# Run SPLIT purification
res_split <- SPLIT::purify(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_purify_singlets = TRUE # optional
)

# Create a purified Seurat object
xe_purified <- CreateSeuratObject(
  counts = res_split$purified_counts,
  meta.data = res_split$cell_meta,
  assay = "Xenium"
)

# Optional: Filter, normalize and visualize
xe_purified <- subset(xe_purified, subset = nCount_Xenium > 5)
xe_purified <- xe_purified %>%
  SCTransform(assay = "Xenium") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20)
UMAPPlot(xe_purified, group.by = c("first_type"), label = T, repel = T) + theme(aspect.ratio = 1)
```
