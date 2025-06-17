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
#  expect_2    <- bead - expect_1
  variance    <- expect_1 * (1 - posterior_1)

  return(list(expect_1 = expect_1,
             # expect_2 = expect_2,
              variance = variance))
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

  doublets_certain    <- results_df[results_df$spot_class %in% is.certain,] %>% rownames()
  doublets_uncertain  <- results_df[results_df$spot_class %in% c("doublet_uncertain"),] %>% rownames()

  if (is.null(n_workers)) {
    n_workers <- max(1, BiocParallel::multicoreWorkers() - 1)
  }

  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers)
  }

  gene_list <- intersect(rownames(counts), rownames(cell_type_info[[1]]))
  print(length(gene_list))
  # Function to decompose certain doublets
  decompose_certain <- function(bead, results_df_bead, ct_weights_bead, gene_list, cell_type_info) {
    tryCatch({

      bead <- bead[gene_list, ]
      barcode <- rownames(results_df_bead)[1]
      type1 <- as.vector(results_df_bead[barcode, "first_type"])
      type2 <- as.vector(results_df_bead[barcode, "second_type"])

      w1 <- as.vector(results_df_bead[barcode, "weight_first_type"])
      w2 <- as.vector(results_df_bead[barcode, "weight_second_type"])

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
  decompose_uncertain <- function(bead, results_df_bead, ct_weights_bead, gene_list, cell_type_info) {
    tryCatch({

      bead <- bead[gene_list,]
      barcode <- rownames(results_df_bead)[1]
      type1 <- as.vector(results_df_bead[barcode, "first_type"])
      type2 <- as.vector(results_df_bead[barcode, "second_type"])

      if(is.na(type2)){ # highly confident singlet -> Do Not Purify!
        return(list(barcode = barcode, res = bead))
      }

      wgts <- ct_weights_bead
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
      cat(round(i * 100 / length(barcodes)), "%\n")

      # Select only relevant barcodes for this chunk
      chunk_barcodes <- barcodes[i:min(i + chunk_size - 1, length(barcodes))]



      # Process in parallel, passing only the subsetted counts
      chunk_results <- bplapply(
        chunk_barcodes,
        function(barcode) {
          decompose_func(counts[, barcode, drop = FALSE], results_df[barcode,], ct_weights[barcode,], gene_list, cell_type_info)
        },
        BPPARAM = BPPARAM
      )

      results_list <- c(results_list, chunk_results)
    }

    results_list <- Filter(Negate(is.null), results_list)
    return(results_list)
  }


  # Process certain doublets
  cat("Processing certain doublets...\n")
  certain_results <- process_chunks(doublets_certain, decompose_certain)
  res_certain_mtrx <- matrix(NA, nrow = length(gene_list), ncol = length(doublets_certain), dimnames = list(gene_list, doublets_certain))
  for (res in certain_results) {
    res_certain_mtrx[, res$barcode] <- res$res
  }

  # Process uncertain doublets
  cat("Processing uncertain doublets...\n")
  uncertain_results <- process_chunks(doublets_uncertain, decompose_uncertain)
  res_uncertain_mtrx <- matrix(NA, nrow = length(gene_list), ncol = length(doublets_uncertain), dimnames = list(gene_list, doublets_uncertain))
  for (res in uncertain_results) {
    res_uncertain_mtrx[, res$barcode] <- res$res
  }

  # Combine results
  cat("Combaning doublets results ...\n")
  purified <- cbind(res_certain_mtrx, res_uncertain_mtrx)
  cell_ids <- c(colnames(res_certain_mtrx), colnames(res_uncertain_mtrx))

  # Process singlets
  if(!DO_purify_singlets){
    cat("Processing singlets...\n")
    singlets <- results_df[results_df$spot_class == "singlet",]
    res_singlet_mtrx <- counts[gene_list, rownames(singlets)]
    purified <- cbind(purified, res_singlet_mtrx)
    cell_ids <- c(cell_ids, colnames(res_singlet_mtrx))
  }

  colnames(purified) <- cell_ids
  rownames(purified) <- gene_list

  #purified <- purified[,colnames(counts)]
  cell_meta <- results_df[colnames(purified),]

  cell_meta$purification_status <- "purified"
  if(!DO_purify_singlets){
    cell_meta$purification_status[cell_meta$spot_class == "singlet"] <- "raw"
  }
  cell_meta$purification_status[is.na(cell_meta$second_type)] <- "raw" # Very confident singlets are not purified

  cell_meta$cell_id <- rownames(cell_meta)

  return(list(purified_counts = purified, cell_meta = cell_meta))
}


#' Purify Data Using RCTD Output
#'
#' This function purifies a query count matrix using the RCTD output, handling both certain and uncertain doublets, as well as optionally processing singlets.
#'
#' @param counts A matrix of gene expression counts, where rows represent genes and columns represent cell IDs.
#' @param rctd RCTD output
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

purify <- function(counts, rctd, DO_purify_singlets, n_workers = NULL, chunk_size = 10000) {

  results_df <- rctd@results$results_df

  common_cells <- intersect(colnames(counts), rownames(results_df))
  results_df <- results_df[common_cells, ]

  cell_type_info <- rctd@cell_type_info[[1]]
  ct_weights <- rctd@results$weights
  ct_weights <- ct_weights[common_cells, colnames(cell_type_info[[1]])]

  counts <- counts[,common_cells]

  return(purify_counts_with_rctd(
    counts = counts,
    results_df = results_df,
    ct_weights = ct_weights,
    cell_type_info = cell_type_info,
    DO_purify_singlets = DO_purify_singlets,
    n_workers = n_workers,
    chunk_size = chunk_size
  )
  )
}












