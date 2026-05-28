#' Convert RCTD result to \code{rctd_free_purify} input format
#'
#' Converts an RCTD (Robust Cell Type Decomposition) object into a list
#' of components compatible with \code{\link{rctd_free_purify}}. This
#' function extracts per-cell deconvolution results, converts them to
#' sparse matrices, and prepares the primary cell type vector and
#' reference profiles.
#'
#' If the RCTD object has not yet been post-processed by SPLIT,
#' \code{\link{run_post_process_RCTD}} is called automatically with a
#' warning.
#'
#' @param rctd An object of class \code{RCTD} containing cell-type
#'   decomposition results, typically produced by
#'   \code{spacexr::run.RCTD()}.
#' @param ... Additional arguments passed to
#'   \code{\link{run_post_process_RCTD}} if post-processing is needed.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{primary_cell_type}{Character vector of dominant cell type
#'     per cell.}
#'   \item{deconvolution_weights}{Sparse matrix (cells x cell types)
#'     of deconvolution weights.}
#'   \item{reference}{Matrix of reference expression profiles
#'     (genes x cell types).}
#' }
#'
#' @details
#' This function is useful for compatibility testing or migrating
#' pipelines from the RCTD output structure to the input format
#' expected by \code{\link{rctd_free_purify}}. The helper functions
#' \code{\link{res_df_2_ijx}} and \code{\link{build_sparse_from_ijx}}
#' are used internally.
#'
#' @seealso
#'   \code{\link{rctd_free_purify}},
#'   \code{\link{res_df_2_ijx}},
#'   \code{\link{build_sparse_from_ijx}},
#'   \code{\link{run_post_process_RCTD}}
#'
#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#'
#' @export
convert_rctd_result_to_purify_input <- function(rctd, ...) {

  ## ---- 1. Check rctd@results is populated -------------------------------
  if (is.null(rctd@results)) {
    stop(
      "'rctd@results' is NULL. ",
      "Run spacexr::run.RCTD() first."
    )
  }

  ## ---- 2. Auto-run post-processing if not yet done ----------------------
  if (!"results_df_old" %in% names(rctd@results)) {
    warning(
      "The RCTD object has not been post-processed by SPLIT.\n",
      "Calling SPLIT::run_post_process_RCTD() automatically.\n",
      "To suppress this warning, run it explicitly first:\n",
      "  rctd <- SPLIT::run_post_process_RCTD(rctd)"
    )
    rctd <- SPLIT::run_post_process_RCTD(rctd, ...)
  }

  ## ---- 3. Validate results_df and required columns ----------------------
  if (is.null(rctd@results$results_df)) {
    stop(
      "'rctd@results$results_df' is NULL after post-processing.\n",
      "This is unexpected — please check your RCTD object."
    )
  }

  results_df   <- rctd@results$results_df
  required_cols <- c("first_type", "second_type",
                     "weight_first_type", "weight_second_type",
                     "spot_class")
  missing_cols  <- required_cols[!required_cols %in% names(results_df)]

  if (length(missing_cols) > 0L) {
    stop(
      "The following required columns are missing from ",
      "'rctd@results$results_df':\n",
      "  ", paste(missing_cols, collapse = ", "), "\n",
      "Ensure SPLIT::run_post_process_RCTD() completed ",
      "successfully."
    )
  }

  ## ---- 4. Validate reference profiles -----------------------------------
  if (is.null(rctd@cell_type_info) ||
      is.null(rctd@cell_type_info[[1L]][[1L]])) {
    stop(
      "'rctd@cell_type_info' is NULL or incomplete.\n",
      "The RCTD object may not have been run to completion."
    )
  }

  ## ---- 5. Extract primary and secondary cell types ----------------------
  primary_cell_type   <- as.vector(results_df$first_type)
  names(primary_cell_type) <- rownames(results_df)

  ## ---- 6. Build sparse weight matrix ------------------------------------
  ijx            <- res_df_2_ijx(rctd = rctd)
  weights_sparse <- build_sparse_from_ijx(ijx)

  ## ---- 7. Return --------------------------------------------------------
  list(
    primary_cell_type     = primary_cell_type,
    deconvolution_weights = weights_sparse,
    reference             = rctd@cell_type_info[[1L]][[1L]]
  )
}

#' Extract i–j–x triplets from RCTD results
#'
#' @description
#' Parses the `results_df` slot of an RCTD object to produce triplet vectors
#' (`i_vec`, `j_vec`, `x_vec`) representing the sparse cell-type weight matrix.
#' Each row (cell) is expanded to contain weights for its associated cell types.
#'
#' @param rctd An object of class `RCTD` with populated `@results` and
#'   `@results$weights` slots.
#' @param cell_types Optional character vector of cell type names. If `NULL`,
#'   inferred from the column names of `rctd@results$weights`.
#'
#' @return A **list** with the following components:
#' \describe{
#'   \item{`i_vec`}{Character vector of cell (spot) IDs.}
#'   \item{`j_vec`}{Character vector of cell-type names.}
#'   \item{`x_vec`}{Numeric vector of deconvolution weights.}
#'   \item{`col_ids`}{Character vector of all cell types (column order).}
#' }
#'
#' @details
#' Handles both singlet and doublet spots. For uncertain doublets
#' (`spot_class == "doublet_uncertain"`), it retains the full decomposition
#' weights from the RCTD results.
#'
#' @seealso [convert_rctd_result_to_purify_input()], [build_sparse_from_ijx()]
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%

res_df_2_ijx <-  function(rctd, cell_types = NULL){
  results_df <- rctd@results$results_df
  results_df$first_type <- as.vector(results_df$first_type)
  results_df$second_type <- as.vector(results_df$second_type)

  if(is.null(cell_types)){
    cell_types <- colnames(rctd@results$weights)
  }

  doublets_uncertain <- results_df %>% filter(spot_class == "doublet_uncertain") %>% rownames()
  #print(length(doublets_uncertain))

  # cells for which 1-2 cell types to consider (singlets and doublets certain):
  results_df <- results_df %>% filter(spot_class %in% c("singlet", "doublet_certain"))

  i_vec <- rep(rownames(results_df), 2)
  j_vec <- c(results_df$first_type, results_df$second_type)
  x_vec <- c(results_df$weight_first_type, results_df$weight_second_type)
  values_to_keep <- which(!is.na(j_vec)) # remove elements corresponding to secondary cell type in highly confident singlets
  i_vec <- i_vec[values_to_keep]
  j_vec <- j_vec[values_to_keep]
  x_vec <- x_vec[values_to_keep]

  # for doublets uncertain, use full decomposition

  i_vec_du <- rep(doublets_uncertain, length(cell_types))
  j_vec_du <- rep(cell_types, each = length(doublets_uncertain))
  x_vec_du <- rctd@results$weights[doublets_uncertain,] %>% as.vector()

  return(
    list(
      i_vec = c(i_vec, i_vec_du),
      j_vec = c(j_vec, j_vec_du),
      x_vec = c(x_vec, x_vec_du),
      col_ids = cell_types
    )
  )
}


#' Build sparse deconvolution matrix from i–j–x triplets
#'
#' @description
#' Constructs a sparse matrix representation of deconvolution weights from
#' i–j–x triplet lists (cell IDs, cell-type IDs, and weight values), as
#' produced by [`res_df_2_ijx()`].
#'
#' @param ijx A **list** containing elements `i_vec`, `j_vec`, `x_vec`, and
#'   `col_ids`, as returned by [`res_df_2_ijx()`].
#' @param row_order Optional **character vector** specifying the desired row
#'   ordering (cell IDs). If `NULL`, the unique order from `ijx$i_vec` is used.
#'
#' @return A **sparse matrix** (`dgCMatrix`) of deconvolution weights
#'   (cells × cell types).
#'
#' @seealso [res_df_2_ijx()], [convert_rctd_result_to_purify_input()]
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom rlang `%||%`

build_sparse_from_ijx <- function(ijx, row_order = NULL) {
  i <- match(ijx$i_vec, row_order %||% unique(ijx$i_vec))
  j <- match(ijx$j_vec, ijx$col_ids)

  sparseMatrix(
    i = i,
    j = j,
    x = ijx$x_vec,
    dims = c(length(unique(ijx$i_vec)), length(ijx$col_ids)),
    dimnames = list(unique(ijx$i_vec), ijx$col_ids)
  )
}
