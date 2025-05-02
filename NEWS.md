# SPLIT 0.1.0

Initial release of the SPLIT (Spatial Purification of Layered Intracellular Transcripts) package.

## Features

-   Purifies spatial transcriptomics data by removing background contamination using RCTD deconvolution results.
-   Computes local contamination scores and models background diffusion.
-   Refines cell-type annotations based on neighborhood transcriptomic homogeneity (SPLIT-shift).
-   Outputs a purified count matrix with improved cell-type specificity.
-   Scores computed by SPLIT can be used to identify highly contaminating cell types and assess cell co-localization patterns.

⚠️ **IMPORTANT:**\
SPLIT currently requires **doublet-mode** RCTD results from the original [spacexr GitHub repository](https://github.com/dmcable/spacexr) or its faster [HD fork](https://github.com/jpromeror/spacexr/tree/HD), **not** from the newly released [Bioconductor version](https://www.bioconductor.org/packages/release/bioc/html/spacexr.html).\
**Compatibility with Bioconductor's spacexr is coming soon.**
