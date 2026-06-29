# SPLIT 0.3.0

- **Residual contamination removal**: new optional flag `DO_remove_residual_contamination = TRUE`
  in `purify()` zeros out genes with no reference support in any of the cell types assigned
  to a given cell — signal that proportional redistribution alone cannot remove. The reference
  specificity threshold is automatically derived from the reference distribution, requiring no
  manual tuning. Typically removes an additional 2–5% of counts on top of standard SPLIT, with
  a measurable improvement in marker gene specificity. See the dedicated
  [vignette](https://github.com/bdsc-tds/SPLIT/blob/main/doc/Residual_Contamination_Removal.html)
  for a full comparison against standard SPLIT.
- Remove `epsilon` from the purification formula in `rctd_free_purify()`: genes with no
  reference support in the primary cell type are now explicitly zeroed rather than receiving
  a small residual value proportional to epsilon, giving cleaner output when
  `DO_remove_residual_contamination = FALSE` as well.
- Update vignette `Run_RCTD_and_SPLIT_on_Xenium.Rmd` with notes on SPLIT's annotation-method
  agnostic interface and the availability of py-RCTD as a faster deconvolution backend.
- Add new vignette `Residual_Contamination_Removal.Rmd` showcasing the new feature with a
  quantitative comparison of Raw, SPLIT, and SPLIT + residual contamination removal on a
  public Xenium breast cancer dataset.

All changes are backward compatible. `DO_remove_residual_contamination` defaults to `FALSE`,
so existing `purify()` calls are unaffected.

# SPLIT 0.2.3

- fix a silent data corruption bug in `convert_rctd_result_to_purify_input()`
  where missing `weight_first_type` and `weight_second_type` columns in
  non-post-processed RCTD objects caused deconvolution weights to be silently
  recycled, incorrectly removing most transcripts from normal cells (#23).
  The function now validates all required fields upfront and automatically
  applies `run_post_process_RCTD()` with a warning if post-processing has
  not been performed.
  
# SPLIT 0.2.2

- Add a validation check in `purify()`, `add_spatial_metric()` and `add_transcriptomics_metric()` that verifies
  `SPLIT::run_post_process_RCTD()` has been applied to the RCTD object
  before proceeding. An informative error with remediation instructions
  is now thrown if not. Non-post-processed RCTD objects will be formally
  rejected from the next release onwards (#20).
  
# SPLIT 0.2.1

- Add `reconstruct_rctd_from_rctdpy()` to bridge Python `rctd-py` (https://pypi.org/project/rctd-py/) output and the SPLIT post-processing pipeline. The function reads parquet/HDF5 files saved from `rctd-py`, assembles a `spacexr` S4 `RCTD` object, and runs `run_post_process_RCTD(lite = TRUE)`, enabling users who perform deconvolution in Python to use SPLIT without re-running RCTD in R. See `?reconstruct_rctd_from_rctdpy` for the expected file layout and a Python saving snippet.
- Add `arrow` and `rhdf5` to `Imports`.

All changes are backward compatible and do not affect results obtained with earlier versions of SPLIT.

# SPLIT 0.2.0

- **SPLIT is now annotation-method agnostic**: `SPLIT::purify()` no longer requires RCTD output and accepts any deconvolution result. All you need is a cells x cell-types weight matrix, a genes x cell-types reference matrix, and optionally a primary cell-type vector (otherwise inferred as the argmax of the weights).
- Add `SPLIT::convert_rctd_result_to_purify_input()` to extract deconvolution weights, primary cell-type vector and reference matrix from RCTD output for use with the new agnostic interface.
- **Full-transcriptome compatible**: `SPLIT::purify()` now scales to large full-transcriptome platforms such as ATERA (~18,000 genes) and runs to completion in minutes.
- Generate vignette comparing ATERA and Xenium and demonstrating SPLIT purification on both.
- Drop hard dependency on `spacexr`; RCTD-based purification remains fully supported and backward compatible via the `rctd` argument.
- Drop hard dependency on `BiocParallel`; parallelisation is no longer required for the core purification step.

All changes are backward compatible. Results obtained with earlier versions of SPLIT are not affected when using the legacy `rctd` interface.

# SPLIT 0.1.3

- Reassign residual (i.e., removed) counts to neighboring cells of the contaminating cell type via `SPLIT::reassign_residual_counts()`.
- Generate corresponding [vignette](https://github.com/bdsc-tds/SPLIT/blob/main/doc/Reassign_residual_transcripts.html).

All changes are backward compatible and do not affect the results obtained with the earlier version of SPLIT.

# SPLIT 0.1.2

- Accelerate and optimaze memory usage in `SPLIT::run_post_process_RCTD()` and `SPLIT::purify()`
- Make `SPLIT::purify()`compatible with large-scale full-transcriptome data (i.e., VisiumHD)
- Generate vignette of running SPLIT on VisiumHD

All changes are backward compatible and do not affect the results obtained with the earlier version of SPLIT.


# SPLIT 0.1.1
- Fix bugs related to `BiocParallel` compatibility on Windows.
- Accelerate and optimaze memory usage in `SPLIT::run_post_process_RCTD()` and `SPLIT::purify()`
- Fix some bugs in `SPLIT::purify()`

All changes are backward compatible and do not affect the results obtained with the earlier version of SPLIT.

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
