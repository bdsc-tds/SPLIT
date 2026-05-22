#' Reconstruct an RCTD Object from rctd-py Parquet Output
#'
#' Reads the parquet/HDF5 files saved by the Python \pkg{rctd-py} pipeline,
#' assembles a \code{spacexr} S4 \code{RCTD} object, and runs
#' \code{SPLIT::run_post_process_RCTD(lite = TRUE)}.
#'
#' @section Saving rctd-py output for compatibility:
#' The following Python snippet produces the exact file layout expected by this
#' function. Run it after fitting your rctd-py model, replacing \code{xe},
#' \code{xe_rctd_results}, and \code{sc_reference_for_xe} with your own
#' objects:
#'
#' \preformatted{
#' import h5py, numpy as np, pandas as pd
#' from pathlib import Path
#'
#' save_dir = Path(f"rctd_results_xe_{annot_level}_{reference}")
#' save_dir.mkdir(exist_ok=True)
#'
#' # Cell IDs of pixels that passed the pixel mask
#' cell_ids = xe.obs_names[xe_rctd_results.pixel_mask]
#' assert len(cell_ids) == xe_rctd_results.weights.shape[0]
#' assert len(cell_ids) == len(xe_rctd_results.spot_class)
#'
#' # cell_ids.parquet  -- barcodes of retained pixels
#' pd.DataFrame({"cell_id": cell_ids}).to_parquet(
#'     save_dir / "cell_ids.parquet", index=False)
#'
#' # weights.parquet  -- full weight matrix (cells x cell_types), no index
#' pd.DataFrame(
#'     xe_rctd_results.weights,
#'     columns=xe_rctd_results.cell_type_names
#' ).to_parquet(save_dir / "weights.parquet")
#'
#' # weights_doublet.parquet  -- doublet weights, columns MUST be w_1 / w_2
#' pd.DataFrame(
#'     xe_rctd_results.weights_doublet,
#'     columns=["w_1", "w_2"]
#' ).to_parquet(save_dir / "weights_doublet.parquet")
#'
#' # spot_results.parquet  -- per-pixel classification results
#' #   spot_class   : integer (0=reject 1=singlet 2=doublet_certain
#' #                           3=doublet_uncertain)
#' #   first_class / second_class : boolean
#' #   first_type_name / second_type_name : human-readable cell type strings
#' df = pd.DataFrame({
#'     "cell_id"         : cell_ids,
#'     "spot_class"      : xe_rctd_results.spot_class,
#'     "first_type"      : xe_rctd_results.first_type,
#'     "second_type"     : xe_rctd_results.second_type,
#'     "first_class"     : xe_rctd_results.first_class,
#'     "second_class"    : xe_rctd_results.second_class,
#'     "min_score"       : xe_rctd_results.min_score,
#'     "singlet_score"   : xe_rctd_results.singlet_score,
#'     "first_type_name" : [xe_rctd_results.cell_type_names[i]
#'                          for i in xe_rctd_results.first_type],
#'     "second_type_name": [xe_rctd_results.cell_type_names[i]
#'                          for i in xe_rctd_results.second_type],
#' })
#' df.to_parquet(save_dir / "spot_results.parquet")
#'
#' # pixel_mask.parquet  -- mask over ALL original obs_names (not just retained)
#' pd.DataFrame({
#'     "cell_id"    : xe.obs_names,
#'     "pixel_mask" : xe_rctd_results.pixel_mask
#' }).to_parquet(save_dir / "pixel_mask.parquet")
#'
#' # metadata.parquet  -- ordered list of cell type names
#' pd.DataFrame(
#'     {"cell_type_names": xe_rctd_results.cell_type_names}
#' ).to_parquet(save_dir / "metadata.parquet")
#'
#' # reference_profiles.h5  -- gene expression profiles of the reference
#' #   profiles        : float array (cell_types x genes)
#' #   cell_type_names : byte-string array
#' #   gene_names      : byte-string array
#' with h5py.File(save_dir / "reference_profiles.h5", "w") as f:
#'     f.create_dataset("profiles",
#'                      data=sc_reference_for_xe.profiles)
#'     f.create_dataset("cell_type_names",
#'                      data=np.asarray(
#'                          sc_reference_for_xe.cell_type_names, dtype="S"))
#'     f.create_dataset("gene_names",
#'                      data=np.asarray(
#'                          sc_reference_for_xe.gene_names, dtype="S"))
#'
#' # de_genes.parquet  -- differentially expressed genes (optional, not read
#' #                       by reconstruct_rctd_from_rctdpy but kept for
#' #                       provenance)
#' pd.DataFrame(
#'     {"DE_genes": sc_reference_for_xe.get_de_genes()}
#' ).to_parquet(save_dir / "de_genes.parquet")
#' }
#'
#' @section Expected files in \code{save_dir}:
#' \describe{
#'   \item{\code{cell_ids.parquet}}{Column \code{cell_id}: barcodes of cells
#'     passing the pixel mask.}
#'   \item{\code{weights.parquet}}{Cells x cell-types full weight matrix.}
#'   \item{\code{weights_doublet.parquet}}{Cells x 2 doublet weights
#'     (\code{w_1}, \code{w_2}).}
#'   \item{\code{spot_results.parquet}}{Per-cell results: \code{spot_class}
#'     (integer), \code{first_type_name}, \code{second_type_name},
#'     \code{first_class} (logical), \code{second_class} (logical),
#'     \code{min_score}, \code{singlet_score}.}
#'   \item{\code{pixel_mask.parquet}}{All original cell IDs with boolean
#'     \code{pixel_mask} column.}
#'   \item{\code{metadata.parquet}}{Column \code{cell_type_names}.}
#'   \item{\code{reference_profiles.h5}}{HDF5 with datasets \code{profiles}
#'     (cell_types x genes), \code{gene_names}, \code{cell_type_names}.}
#' }
#'
#' @param save_dir Character scalar. Path to the rctd-py output directory.
#' @param min_weight Numeric scalar. Minimum weight threshold passed to
#'   \code{SPLIT::run_post_process_RCTD()}. Default \code{0.01}.
#' @param class_df Optional \code{data.frame} with rownames equal to
#'   \code{cell_type_names} and a single column \code{class} mapping each
#'   cell type to a higher-level grouping used by SPLIT. When \code{NULL}
#'   (default) each cell type maps to itself (identity grouping).
#' @param spot_class_map Named character vector mapping integer spot-class
#'   codes (as strings) to spacexr factor levels. The default follows the
#'   standard rctd-py encoding: \code{0 = reject}, \code{1 = singlet},
#'   \code{2 = doublet_certain}, \code{3 = doublet_uncertain}. Override only
#'   if your rctd-py version uses a different encoding.
#'
#' @return A post-processed \code{RCTD} S4 object as returned by
#'   \code{SPLIT::run_post_process_RCTD()}.
#'
#' @seealso \code{\link[SPLIT]{run_post_process_RCTD}}
#'
#' @importFrom arrow read_parquet
#' @importFrom rhdf5 h5read
#' @importFrom Matrix Matrix
#' @importFrom methods new
#'
#' @export
reconstruct_rctd_from_rctdpy <- function(
    save_dir,
    min_weight     = 0.01,
    class_df       = NULL,
    spot_class_map = c(
      "0" = "reject",
      "1" = "singlet",
      "2" = "doublet_certain",
      "3" = "doublet_uncertain"
    )
) {

  # -- Input validation -------------------------------------------------------
  stopifnot(
    is.character(save_dir), length(save_dir) == 1L, dir.exists(save_dir),
    is.numeric(min_weight), length(min_weight) == 1L,
    is.character(spot_class_map), !is.null(names(spot_class_map))
  )

  # -- 1. Load parquet files --------------------------------------------------
  cell_ids_df     <- arrow::read_parquet(file.path(save_dir, "cell_ids.parquet"))
  weights_df      <- arrow::read_parquet(file.path(save_dir, "weights.parquet"))
  weights_dbl_df  <- arrow::read_parquet(file.path(save_dir, "weights_doublet.parquet"))
  spot_results_df <- arrow::read_parquet(file.path(save_dir, "spot_results.parquet"))
  pixel_mask_df   <- arrow::read_parquet(file.path(save_dir, "pixel_mask.parquet"))
  metadata_df     <- arrow::read_parquet(file.path(save_dir, "metadata.parquet"))

  cell_ids        <- cell_ids_df$cell_id
  cell_type_names <- metadata_df$cell_type_names

  # -- 2. Build results_df ----------------------------------------------------
  spot_class_levels <- c("singlet", "doublet_certain", "doublet_uncertain", "reject")

  spot_class_mapped <- spot_class_map[as.character(spot_results_df$spot_class)]
  if (anyNA(spot_class_mapped)) {
    bad <- unique(spot_results_df$spot_class[is.na(spot_class_mapped)])
    stop(
      "Unmapped spot_class integer(s): ", paste(bad, collapse = ", "),
      ". Update `spot_class_map` to cover these values."
    )
  }

  results_df <- data.frame(
    spot_class = factor(spot_class_mapped, levels = spot_class_levels),
    first_type = factor(spot_results_df$first_type_name,  levels = cell_type_names),
    second_type = factor(spot_results_df$second_type_name, levels = cell_type_names),
    first_class = factor(
      ifelse(spot_results_df$first_class, "doublet_certain", "singlet"),
      levels = spot_class_levels
    ),
    second_class = factor(
      ifelse(spot_results_df$second_class, "doublet_certain", "singlet"),
      levels = spot_class_levels
    ),
    min_score     = spot_results_df$min_score,
    singlet_score = spot_results_df$singlet_score,
    row.names     = cell_ids
  )

  # -- 3. Full weight matrix  (cells x cell_types) ----------------------------
  weights_mat <- as.matrix(weights_df)
  rownames(weights_mat) <- cell_ids
  colnames(weights_mat) <- cell_type_names

  # -- 4. Doublet weight matrix  (cells x 2) ----------------------------------
  weights_doublet_mat <- as.matrix(weights_dbl_df[, c("w_1", "w_2")])
  rownames(weights_doublet_mat) <- cell_ids
  colnames(weights_doublet_mat) <- c("first_type", "second_type")

  # -- 5. Reference profiles from HDF5  (genes x cell_types) -----------------
  h5_path        <- file.path(save_dir, "reference_profiles.h5")
  prof           <- rhdf5::h5read(h5_path, "profiles")   # cell_types x genes in Python
  gene_names_ref <- rhdf5::h5read(h5_path, "gene_names")
  ct_names_ref   <- rhdf5::h5read(h5_path, "cell_type_names")

  profiles_df <- as.data.frame(t(prof))                  # transpose -> genes x cell_types
  rownames(profiles_df) <- gene_names_ref
  colnames(profiles_df) <- ct_names_ref

  cell_type_info <- list(
    info   = list(profiles_df, ct_names_ref, length(ct_names_ref)),
    renorm = NULL
  )

  # -- 6. class_df: cell type -> higher-level class mapping -------------------
  if (is.null(class_df)) {
    class_df <- data.frame(
      class     = cell_type_names,
      row.names = cell_type_names
    )
  } else {
    missing_types <- setdiff(cell_type_names, rownames(class_df))
    if (length(missing_types) > 0L) {
      stop(
        "Supplied class_df is missing rows for: ",
        paste(missing_types, collapse = ", ")
      )
    }
  }

  # -- 7. Minimal SpatialRNA placeholder -------------------------------------
  # Raw counts are not accessed when lite = TRUE; dummy values are sufficient.
  all_cell_ids <- pixel_mask_df$cell_id
  n_cells      <- length(all_cell_ids)

  coords <- data.frame(
    x = rep(0, n_cells),
    y = rep(0, n_cells),
    row.names = all_cell_ids
  )
  dummy_counts <- Matrix::Matrix(
    0L,
    nrow     = 1L,
    ncol     = n_cells,
    dimnames = list("dummy_gene", all_cell_ids),
    sparse   = TRUE
  )
  spatialRNA_obj <- methods::new(
    "SpatialRNA",
    coords = coords,
    counts = dummy_counts,
    nUMI   = setNames(rep(1L, n_cells), all_cell_ids)
  )

  # -- 8. Assemble RCTD S4 object --------------------------------------------
  RCTD <- methods::new(
    "RCTD",
    spatialRNA         = spatialRNA_obj,
    originalSpatialRNA = spatialRNA_obj,
    reference          = methods::new("Reference"),
    results            = list(
      results_df_xe   = results_df,
      weights         = weights_mat,
      weights_doublet = weights_doublet_mat
    ),
    cell_type_info   = cell_type_info,
    config           = list(RCTDmode = "doublet"),
    internal_vars    = list(class_df = class_df),
    de_results       = list(),
    internal_vars_de = list()
  )

  # -- 9. SPLIT post-processing ----------------------------------------------
  RCTD_processed <- SPLIT::run_post_process_RCTD(RCTD, lite = TRUE, min_weight = min_weight)

  return(RCTD_processed)
}