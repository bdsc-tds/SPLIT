#' Purify Gene Expression Counts Using Deconvolution Results and Reference
#'
#' @description
#' This function performs signal purification of observed single-cell or spatial
#' gene expression counts based on deconvolution weights and a reference
#' cell-type expression matrix. It provides a memory-efficient, chunked
#' implementation that scales to large datasets without exceeding memory limits.
#'
#' The algorithm adjusts each cell's expression profile to reduce cross-cell-type
#' contamination according to its inferred deconvolution weights and primary
#' cell-type assignment.
#'
#' @param counts A **matrix** or `dgCMatrix` of raw expression counts
#'   (genes × cells). This typically corresponds to a single-cell or
#'   spatial transcriptomics assay.
#'
#' @param primary_cell_type A **character** or **factor** vector giving the
#'   primary (dominant) cell type for each cell (must have names matching
#'   column names of `counts`).
#'
#' @param deconvolution_weights A **matrix** or `dgCMatrix` of cell-type
#'   deconvolution weights (cells × cell types). Each row should sum to 1
#'   or will be scaled to 1 internally if `DO_require_sumup_to_one = TRUE`.
#'
#' @param reference A **matrix** of reference gene expression values
#'   (cell types × genes). This represents average or marker-based
#'   profiles for each cell type in the same gene space as `counts`.
#'
#' @param cells_to_purify Optional **character vector** specifying a subset
#'   of cells to purify. If provided, the remaining cells will still be included
#'   in the resulting object in their original (raw) form.
#'   This parameter is intended for convenience when purifying a subset of cells.
#'   However, if a large proportion of cells should remain unmodified, it is more
#'   efficient to restrict `counts` and `deconvolution_weights` to the cells that
#'   require purification, and then concatenate the purified and raw results
#'   afterward.
#'
#' @param DO_run_in_chunks Logical; whether to process cells in memory-safe
#'   chunks (`TRUE`, default) or compute all cells at once (`FALSE`).
#'
#' @param chunk_size Integer; the number of cells per block when
#'   `DO_run_in_chunks = TRUE`. Larger blocks are faster but require more RAM.
#'
#' @param DO_require_sumup_to_one Logical; if `TRUE` (default), rescale
#'   deconvolution weights per cell so that each row sums to 1 when deviations
#'   are detected.
#'
#' @details
#' For each cell, the function computes:
#' \deqn{D_{cell,gene} = \sum_{ct} w_{cell,ct} \times ref_{ct,gene}}
#' \deqn{N_{cell,gene} = w_{cell,primary(ct)} \times ref_{primary(ct),gene}}
#' The purified expression is then:
#' \deqn{purified_{gene,cell} = counts_{gene,cell} \times \frac{N_{cell,gene}}{D_{cell,gene}}}
#'
#' In chunked mode, computation is performed over blocks of cells to reduce
#' peak memory use. The result is reconstructed as a sparse matrix.
#'
#' @return A **sparse matrix** (`dgCMatrix`) of purified gene expression
#' counts (genes × cells), aligned with the input `counts` matrix.
#'
#' @section Memory and Performance:
#' - Chunked computation (`DO_run_in_chunks = TRUE`) is recommended for
#'   large datasets (tens of thousands of cells or more).
#' - Setting `DO_run_in_chunks = FALSE` is faster but may require >40 GB of RAM
#'   depending on dataset size.
#' - The function uses sparse-matrix algebra via the **Matrix** package
#'   for improved memory efficiency.
#'
#' @examples
#' \dontrun{
#' purified <- purify_rctd_free(
#'   counts = GetAssayData(xe, assay = "Xenium", layer = "counts"),
#'   primary_cell_type = new_input$primary_cell_type,
#'   deconvolution_weights = new_input$deconvolution_weights,
#'   reference = t(new_input$reference),
#'   DO_run_in_chunks = TRUE,
#'   chunk_size = 50000
#' )
#' }
#'
#' @import Matrix
#' @import SingleCellExperiment
#' @importFrom utils tail
#' @importFrom methods as
#'
#' @export # temporarily, later replace with wrapper that calls correct function

rctd_free_purify <- function(
    counts,                 # genes x cells
    primary_cell_type,      # vector
    deconvolution_weights,  # cells x weights
    reference,              # genes x cell types
    cells_to_purify = NULL, # vector
    DO_run_in_chunks = TRUE,
    chunk_size = 50000,
    DO_require_sumup_to_one = TRUE,
    DO_output_sce = TRUE
){

  # --- Core identifiers --------------------------------------------------------
  cell_names_counts <- colnames(counts)
  gene_names_counts <- rownames(counts)
  cell_names_weights <- rownames(deconvolution_weights)
  gene_names_reference <- colnames(reference)

  shared_cells <- intersect(cell_names_counts, cell_names_weights)
  shared_genes <- intersect(gene_names_counts, gene_names_reference)

  # --- Subset to shared genes/cells -------------------------------------------
  counts <- counts[shared_genes, shared_cells]
  primary_cell_type <- primary_cell_type[shared_cells]
  deconvolution_weights <- deconvolution_weights[shared_cells, ]

  cell_types_weights <- sort(colnames(deconvolution_weights))
  cell_types_reference <- sort(rownames(reference))

  if (!identical(cell_types_weights, cell_types_reference)) {
    stop("Reference and deconvolution cell types do not match!")
  } else {
    cell_types <- cell_types_weights
    deconvolution_weights <- deconvolution_weights[, cell_types]
  }

  reference <- reference[cell_types, shared_genes]

  counts <- Matrix::Matrix(counts, sparse = TRUE)
  deconvolution_weights <- Matrix::Matrix(deconvolution_weights, sparse = TRUE)

  # --- Normalize weights ------------------------------------------------------
  row_sum_weights <- Matrix::rowSums(deconvolution_weights)
  if (min(row_sum_weights) < 1 || max(row_sum_weights) > 1) {
    warning("Some deconvolution weights do not sum to 1.")
    if (DO_require_sumup_to_one) {
      message("Rescaling deconvolution weights per cell.")
      deconvolution_weights <- deconvolution_weights / row_sum_weights
    }
  }

  # --- Filter cells and genes -------------------------------------------------
  if (length(cell_names_counts) < length(shared_cells)) {
    warning("Some cells missing deconvolution results; removed from output.")
  }

  if (!is.null(cells_to_purify)) {
    cells_to_keep_raw <- setdiff(shared_cells, cells_to_purify)

    if (length(cells_to_keep_raw) > 0) {
      message("Setting deconvolution weights to identity for ",
              length(cells_to_keep_raw), " non-purified cells (sparse update).")

      # Prepare indices for nonzero (1) entries
      keep_types <- primary_cell_type[cells_to_keep_raw]
      valid_idx <- keep_types %in% colnames(deconvolution_weights)

      if (any(!valid_idx)) {
        warning(sum(!valid_idx), " primary cell types not found in deconvolution matrix; skipped.")
      }

      valid_cells <- cells_to_keep_raw[valid_idx]
      valid_types <- keep_types[valid_idx]

      # --- Sparse update strategy ---
      # 1. Zero out only relevant rows (set entire row to 0)
      #    This is safe for dgCMatrix and avoids dense coercion.
      deconvolution_weights[valid_cells, ] <- 0

      # 2. Assign 1s efficiently for (cell, type) pairs
      ones_matrix <- Matrix::sparseMatrix(
        i = match(valid_cells, rownames(deconvolution_weights)),
        j = match(valid_types, colnames(deconvolution_weights)),
        x = 1,
        dims = dim(deconvolution_weights),
        dimnames = dimnames(deconvolution_weights)
      )

      # 3. Combine (replace 0s with new 1s)
      deconvolution_weights <- pmax(deconvolution_weights, ones_matrix)
      rm(ones_matrix); gc(FALSE)
    }
  }

  if (length(gene_names_counts) < length(shared_genes)) {
    warning("Some genes not present in the reference; removed from output.")
  }

  # --- Precompute helper vectors ---------------------------------------------
  primary_type_weight <- deconvolution_weights[cbind(shared_cells, primary_cell_type)]
  n_celltypes_per_cell <- Matrix::rowSums(deconvolution_weights[shared_cells, ] > 0)
  epsilon <- 1e-10

  # --- Chunked or full computation -------------------------------------------
  if (DO_run_in_chunks) {

    block_indices <- split(seq_along(shared_cells),
                           ceiling(seq_along(shared_cells) / chunk_size))
    n_blocks <- length(block_indices)

    row_idx_list <- vector("list", n_blocks)
    col_idx_list <- vector("list", n_blocks)
    values_list <- vector("list", n_blocks)

    for (b in seq_len(n_blocks)) {
      idx <- block_indices[[b]]
      cells_in_block <- shared_cells[idx]

      message("Processing block ", b, "/", n_blocks,
              " (", idx[1], "-", tail(idx, 1), ")")

      # --- Compute block numerator & denominator -----------------------------
      block_denominator <- as.matrix(deconvolution_weights[cells_in_block, , drop = FALSE]) %*%
        reference + epsilon

      block_numerator <- primary_type_weight[idx] *
        reference[primary_cell_type[idx], , drop = FALSE] +
        epsilon / n_celltypes_per_cell[idx]

      # --- Correct counts ----------------------------------------------------
      block_corrected <- t(block_numerator / block_denominator) *
        counts[, idx, drop = FALSE]

      # --- Extract nonzero entries ------------------------------------------
      nonzero <- which(block_corrected != 0, arr.ind = TRUE)
      if (nrow(nonzero) > 0) {
        row_idx_list[[b]] <- nonzero[, 1]
        col_idx_list[[b]] <- idx[nonzero[, 2]]
        values_list[[b]] <- block_corrected[nonzero]
      }

      rm(block_denominator, block_numerator, block_corrected, nonzero)
      gc(FALSE)
    }

    # --- Combine sparse results ---------------------------------------------
    row_idx <- unlist(row_idx_list, use.names = FALSE)
    col_idx <- unlist(col_idx_list, use.names = FALSE)
    values <- unlist(values_list, use.names = FALSE)

    purified_counts <- Matrix::sparseMatrix(
      i = row_idx, j = col_idx, x = values,
      dims = dim(counts),
      dimnames = dimnames(counts)
    )

  } else {
    denominator <- deconvolution_weights[shared_cells, cell_types] %*% reference + epsilon
    numerator <- primary_type_weight * reference[primary_cell_type[shared_cells], ] +
      epsilon / n_celltypes_per_cell
    purified_counts <- t(numerator / denominator) * counts
    colnames(purified_counts) <- shared_cells
  }


  # --- assemble final output as SingleCellExperiment -------------------------

  # Compute per-cell metadata

  # Define purification status: "purified" vs "raw"

  n_cell_types <- Matrix::rowSums(deconvolution_weights>0)

  purification_status <- ifelse(
    n_cell_types > 1,
    "purified",
    "raw"
  )

  col_metadata <- data.frame(
    primary_cell_type = primary_cell_type[shared_cells],
    w1_primary = as.numeric(primary_type_weight),
    cell_id = shared_cells,
    first_type = primary_cell_type[shared_cells],
    purification_status = purification_status,
    n_cell_types = n_cell_types
  )

  if(!DO_output_sce){
    return(list(purified_counts = purified_counts, cell_meta = col_metadata))
  } else {

    # Build SCE object
    sce_out <- SingleCellExperiment(
      assays = list(
        counts = counts,               # original counts
        purified_counts = purified_counts     # purified (corrected) matrix
      ),
      colData = col_metadata
    )

    # Store reference info for reproducibility
    metadata(sce_out)$reference_cell_types <- rownames(reference)
    metadata(sce_out)$deconvolution_summary <- list(
      n_cells = length(shared_cells),
      n_genes = length(shared_genes),
      run_in_chunks = DO_run_in_chunks,
      chunk_size = chunk_size,
      rescaled_weights = DO_require_sumup_to_one
    )

    message("Returning SingleCellExperiment with purified counts assay.")
    return(sce_out)
  }
}
