# Functions below were adapted and modified from:

#get_decomposed_data() # https://github.com/dmcable/spacexr/blob/0a0861e3d1e16014a20e9b743d0e19d3b42231f3/R/postProcessing.R#L81C1-L81C20
#get_decomposed_data_full_doublet() # https://github.com/dmcable/spacexr/blob/0a0861e3d1e16014a20e9b743d0e19d3b42231f3/R/postProcessing.R#L44
#decompose_doublet_fast() # https://github.com/dmcable/spacexr/blob/0a0861e3d1e16014a20e9b743d0e19d3b42231f3/R/RCTD_helper.R#L219C1-L219C23

#' Decompose doublet into 2 profiles
#'
#' @param bead (counts) profile to decompose
#' @param weights cell-type weights (expected to sum up to 1)
#' @param gene_list vector of genes (should be (a subset of) common genes between reference data and query data)
#' @param cell_type_info slot of \link[spacexr]{run.RCDT} output containing reference profile matrix
#' @param type1 name of the main cell-type (expected `first_type`)
#' @param cell_types name(s) of all cell types present in the cell (expects `second_type` in case of decomposition of `doublet_certain` and all present cell types in case of `doublet_uncertain`)
#'
#' @return description list of the profile of the first and the second profiles,
#'in case of decomposition of `doublet_uncertain` (i.e., when `length(cell_types)>1`), the second profile corresponds to the complement (bead - first_profile)
#' @export

decompose_doublet <- function(
    bead, weights, gene_list, cell_type_info, type1, cell_types
){
  bead     <- bead[gene_list]

  N_genes  <- length(gene_list)
  epsilon  <- 1e-10

  if(!(type1 %in% cell_types)){
    cell_types <- c(type1, cell_types)
  }

  cell_types   <- unique(cell_types)
  N_cell_types <- length(cell_types)

  denom       <- rowSums(sweep(as.matrix(cell_type_info[[1]][gene_list,cell_types]), 2, weights[cell_types], FUN = "*")) + epsilon
  posterior_1 <- (weights[type1] * cell_type_info[[1]][gene_list,type1] + epsilon/N_cell_types) / denom
  expect_1    <- posterior_1 * bead
  expect_2    <- bead - expect_1
  variance    <- expect_1 * (1 - posterior_1)

  return(list(expect_1 = expect_1, expect_2 = expect_2, variance = variance))
}


#' Purify Counts Using RCTD Output
#'
#' This function purifies a query count matrix using the RCTD output, handling both certain and uncertain doublets, as well as optionally processing singlets.
#'
#' @param counts A matrix of gene expression counts, where rows represent genes and columns represent cell IDs.
#' @param results_df A data frame from the `results_df` slot of \link[spacexr]{run.RCTD} output, containing cell-type assignments and classifications.
#' @param ct_weights A matrix of cell-type weights for each spot from the RCTD output.
#' @param cell_type_info A list containing cell-type information used for purification.
#' @param DO_purify_singlets Logical; if `TRUE`, singlets will also be purified.
#' @param n_workers Integer; the number of parallel workers to use. If `NULL`, it defaults to the number of available cores minus one.
#' @param chunk_size Integer; the number of barcodes processed in each batch for parallelization. Default is 10,000.
#'
#' @return A list containing:
#' \describe{
#'   \item{purified_counts}{A matrix of purified gene expression counts.}
#'   \item{cell_meta}{A data frame with metadata for each cell, including purification status.}
#' }
#'
#' @import BiocParallel
#' @export

purify_counts_with_rctd <- function(counts, results_df, ct_weights, cell_type_info, DO_purify_singlets, n_workers = NULL, chunk_size = 10000) {

  is.certain <- c("doublet_certain")
  if(DO_purify_singlets){
    is.certain <- c(is.certain, "singlet")
  }

  doublets_certain    <- results_df[results_df$spot_class %in% is.certain,]
  doublets_uncertain  <- results_df[results_df$spot_class %in% c("doublet_uncertain"),]

  if (is.null(n_workers)) {
    n_workers <- multicoreWorkers() - 1
  }
  param <- MulticoreParam(workers = n_workers)

  gene_list <- intersect(rownames(counts), rownames(cell_type_info[[1]]))
  print(length(gene_list))
  # Function to decompose certain doublets
  decompose_certain <- function(counts, barcode, results_df, ct_weights, gene_list, cell_type_info) {
    tryCatch({

      bead <- counts[gene_list, barcode]

      type1 <- as.vector(results_df[barcode, "first_type"])
      type2 <- as.vector(results_df[barcode, "second_type"])

      w1 <- as.vector(results_df[barcode, "weight_first_type"])
      w2 <- as.vector(results_df[barcode, "weight_second_type"])

      if(is.na(type2)){ # highly confident singlet -> Do Not Purify!
        return(list(barcode = barcode, res = bead))
      }

      if(is.na(w2)){
        w1 <- 1
        w2 <- 0
      }

      wgts <- c(w1, w2)
      wgts <- wgts / sum(wgts)

      names(wgts) <- c(type1, type2)

      doub_res <- decompose_doublet(
        bead,
        wgts,
        gene_list,
        cell_type_info,
        type1,
        c(type1, type2)
      )
      return(list(barcode = barcode, res = doub_res$expect_1))
    }, error = function(e) {
      message(sprintf("Error processing barcode %s: %s", barcode, e$message))
      return(NULL)
    })
  }

  # Function to decompose uncertain doublets
  decompose_uncertain <- function(counts, barcode, results_df, ct_weights, gene_list, cell_type_info) {
    tryCatch({

      bead <- counts[gene_list, barcode]

      type1 <- as.vector(results_df[barcode, "first_type"])
      type2 <- as.vector(results_df[barcode, "second_type"])

      if(is.na(type2)){ # highly confident singlet -> Do Not Purify!
        return(list(barcode = barcode, res = bead))
      }

      wgts <- ct_weights[barcode,]
      wgts <- wgts / sum(wgts)
      types <- names(wgts)

      doub_res <- decompose_doublet(
        bead,
        wgts,
        gene_list,
        cell_type_info,
        type1,
        types
      )
      return(list(barcode = barcode, res = doub_res$expect_1))
    }, error = function(e) {
      message(sprintf("Error processing barcode %s: %s", barcode, e$message))
      return(NULL)
    })
  }

  # Helper function to process chunks
  process_chunks <- function(barcodes, decompose_func) {
    results_list <- list()
    for (i in seq(1, length(barcodes), by = chunk_size)) {
      cat(round(i*100/length(barcodes)))
      chunk_barcodes <- barcodes[i:min(i + chunk_size - 1, length(barcodes))]
      chunk_results <- bplapply(
        chunk_barcodes,
        function(barcode) {
          decompose_func(counts, barcode, results_df, ct_weights, gene_list, cell_type_info)
        },
        BPPARAM = param
      )
      results_list <- c(results_list, chunk_results)
    }
    results_list <- Filter(Negate(is.null), results_list)
    return(results_list)
  }

  # Process certain doublets
  cat("Processing certain doublets...\n")
  certain_results <- process_chunks(rownames(doublets_certain), decompose_certain)
  res_certain_mtrx <- matrix(NA, nrow = length(gene_list), ncol = nrow(doublets_certain), dimnames = list(gene_list, rownames(doublets_certain)))
  for (res in certain_results) {
    res_certain_mtrx[, res$barcode] <- res$res
  }

  # Process uncertain doublets
  cat("Processing uncertain doublets...\n")
  uncertain_results <- process_chunks(rownames(doublets_uncertain), decompose_uncertain)
  res_uncertain_mtrx <- matrix(NA, nrow = length(gene_list), ncol = nrow(doublets_uncertain), dimnames = list(gene_list, rownames(doublets_uncertain)))
  for (res in uncertain_results) {
    res_uncertain_mtrx[, res$barcode] <- res$res
  }

  # Process singlets
  cat("Processing singlets...\n")
  if(!DO_purify_singlets){

  }
  # Combine results
  all_DGE <- cbind(res_certain_mtrx, res_uncertain_mtrx)
  cell_ids <- c(colnames(res_certain_mtrx), colnames(res_uncertain_mtrx))


  if(!DO_purify_singlets){
    singlets <- results_df[results_df$spot_class == "singlet",]
    res_singlet_mtrx <- counts[gene_list, rownames(singlets)]
    all_DGE <- cbind(all_DGE, res_singlet_mtrx)
    cell_ids <- c(cell_ids, colnames(res_singlet_mtrx))
  }

  colnames(all_DGE) <- cell_ids
  rownames(all_DGE) <- gene_list

  #all_DGE <- all_DGE[,colnames(counts)]
  cell_meta <- results_df[colnames(all_DGE),]

  cell_meta$purification_status <- "purified"
  if(!DO_purify_singlets){
    cell_meta$purification_status[cell_meta$spot_class == "singlet"] <- "raw"
  }
  cell_meta$purification_status[is.na(cell_meta$second_type)] <- "raw" # Very confident singlets are not purified

  return(list(purified_counts = all_DGE, cell_meta = cell_meta))
}

#' Balance raw and purified data by merging high quality raw data and otherwise purified data
#'
#' Merges raw and purified data into one dataset by keeping raw counts
#' for high quality cells and replacing contaminated cells with their purified profiles. Reject cells are removed
#'
#'
#' @param xe_raw raw seurat object
#' @param xe_purified purified seurat object (post \code{purify_counts_with_rctd})
#' @param threshold value below which cell is considered as contaminated and is replaced by purified profile.  For the moment, it's a single value, but should accept cell-type-specific vector later on #TODO
#' @param score_name name of the param to threshold on
#'
#' @export

balance_raw_and_purified_data <- function(
    xe_raw,
    xe_purified,
    threshold = .5,
    score_name = "sp_second_type_neighborhood_weight"
){
  if(!score_name %in% colnames(xe_raw@meta.data)){
    stop("score", score_name, "is not available in `xe_raw`, please compute it first!")
  }

  if(!"first_type" %in% colnames(xe_raw@meta.data)){
    stop("`first_type` is not available in `xe_raw`, please compute it first!")
  }

  if(!"spot_class" %in% colnames(xe_raw@meta.data)){
    stop("`spot_class` is not available in `xe_raw`, please compute it first!")
  }

  cells_to_replace_with_purified <-
    xe_raw@meta.data %>%
    filter(spot_class != "reject" & .data[[score_name]] < threshold) %>%
    rownames()

  cells_to_remove <-
    xe_raw@meta.data %>%
    filter(spot_class == "reject") %>%
    rownames()

  cells_to_keep_raw <- setdiff(colnames(xe_raw), c(cells_to_replace_with_purified, cells_to_remove))

}

