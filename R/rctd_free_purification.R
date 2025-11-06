rctd_free_purify <- function(
  counts, # genes x cells
  primary_cell_type, # vect
  deconvolution_weights, # cells x weights
  reference, # genes x cell types
  cells_to_purify = NULL # vect
){

  cell_types_from_weights <- colnames(deconvolution_weights) %>% sort()
  cell_types_from_reference <- colnames(reference) %>% sort()

  if(!identical(cell_types_from_weights, cell_types_from_reference)){
    stop("Reference and deconvolution cell types do not match!\n")
  } else {
    ct <- cell_types_from_weights
  }

  # Make sure weights sum up to 1
  rs <- Matrix::rowSums()
  if(min(rs) < 1 | max(rs) > 1){
    warning("Deconvolution does not sum up to 1 in some cells!\n")
  }

  cells_counts <- colnames(counts)
  genes_counts <- rownames(counts)

  cells_deconvolution <- rownams(deconvolution_weights)
  genes_reference <- rownames(reference)

  common_cells <- intersect(cells_counts, cells_deconvolution)
  common_genes <- intersect(genes_counts, genes_reference)

  if(length(cells_counts) < length(common_cells)){
    warning("Some cells do not have decomposition results, they will be removed from the resulting corrected count matrix!\n")
  }

  if(length(genes_counts) < length(common_genes)){
    warning("Some genes are not present in the reference and will be removed from the resulting corrected count matrix!\n")
  }

  W1 <- deconvolution_weights[cbind(common_cells, primary_cell_type[common_cells])]

  Denominator <- deconvolution_weights[common_cells, ct] %*% Matrix::t(reference[common_genes, ct])  #
  Nominator <- ""

}

# for testing (to compare with current version (v0.1.3))
convert_rctd_result_to_purify_input <- function(
    rctd
){
  results_df <- rctd@results$results_df
  primary_cell_type <- results_df$first_type  ## need
  names(primary_cell_type) <- results_df %>% rownames()

  secondary_cell_type <- results_df$second_type
  names(secondary_cell_type) <- results_df %>% rownames()

  ijx <- res_df_2_ijx(rctd = rctd)
  weights_sparse <- build_sparse_from_ijx(ijx) ## need

  return(
    list(
      primary_cell_type = primary_cell_type, # vect
      deconvolution_weights = weights_sparse, # cells x weights
      reference = rctd@cell_type_info[[1]][[1]], # genes x cell types

    )
  )
}


res_df_2_ijx <-  function(rctd, cell_types = NULL){
  results_df <- rctd@results$results_df # to remove
  if(is.null(cell_types)){
    cell_types <- colnames(rctd@results$weights)
  }

  doublets_uncertain <- results_df %>% filter(spot_class == "doublet_uncertain") %>% rownames()
  print(length(doublets_uncertain))

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


library(Matrix)
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

