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
  #bead     <- bead[gene_list]

  N_genes  <- length(gene_list)
  epsilon  <- 1e-10

  if(!(type1 %in% cell_types)){
    cell_types <- c(type1, cell_types)
  }

  cell_types   <- unique(cell_types)
  N_cell_types <- length(cell_types)

  denom       <- rowSums(sweep(as.matrix(cell_type_info[[1]][gene_list,cell_types]), 2, weights[cell_types], FUN = "*")) + epsilon
  posterior_1 <- (weights[type1] * cell_type_info[[1]][gene_list,type1] + epsilon/N_cell_types) / denom
  vec <- posterior_1 * bead  # vector of length G, mostly 0
  nz_idx <- which(as.numeric(vec) != 0)
  vec_sparse <- Matrix::sparseVector(i = nz_idx, x = vec[nz_idx], length = length(vec))
  return(vec_sparse)
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
#' @param DO_parallel Logical; if `TRUE`, biocparallel will be used to accelerate computation. For the moment, non-parallel version perform equally well
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

purify_counts_with_rctd <- function(counts, results_df, ct_weights, cell_type_info, DO_purify_singlets, DO_parallel = FALSE, n_workers = NULL, chunk_size = 10000) {

  sparse <- inherits(counts, "sparseMatrix")

  is.certain <- c("doublet_certain")
  if(DO_purify_singlets){
    is.certain <- c(is.certain, "singlet")
  }

  doublets_certain    <- results_df[results_df$spot_class %in% is.certain,] %>% rownames()
  doublets_uncertain  <- results_df[results_df$spot_class %in% c("doublet_uncertain"),] %>% rownames()

  if(DO_parallel){
    if (is.null(n_workers)) {
      n_workers <- min(4, BiocParallel::multicoreWorkers() - 1)
    }
    if (.Platform$OS.type == "windows") {
      BPPARAM <- BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
    } else {
      BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers)
    }
  } else {
    BPPARAM <- NULL
  }

  gene_list <- rownames(counts)
  cat("N_genes = ", length(gene_list))

  # Function to decompose certain doublets
  decompose_certain <- function(bead, results_df_bead, ct_weights_bead, gene_list, cell_type_info) {
    tryCatch({
      #bead <- bead[gene_list,]
      barcode <- rownames(results_df_bead)[1]
      type1 <- as.vector(results_df_bead[barcode, "first_type"])
      type2 <- as.vector(results_df_bead[barcode, "second_type"])

      w1 <- as.vector(results_df_bead[barcode, "weight_first_type"])
      w2 <- as.vector(results_df_bead[barcode, "weight_second_type"])

      if(is.na(type2)){ # highly confident singlet -> Do Not Purify!

        ###
        nz_idx <- which(as.numeric(bead) != 0)
        vec_sparse <- Matrix::sparseVector(i = nz_idx, x = bead[nz_idx], length = length(bead))
        return(list(barcode = barcode, res = vec_sparse))
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
      return(list(barcode = barcode, res = doub_res))
    }, error = function(e) {
      message(sprintf("Error processing barcode %s: %s", barcode, e$message))
      return(NULL)
    })
  }

  # Function to decompose uncertain doublets
  decompose_uncertain <- function(bead, results_df_bead, ct_weights_bead, gene_list, cell_type_info) {
    tryCatch({

      #bead <- bead[gene_list,]
      barcode <- rownames(results_df_bead)[1]
      type1 <- as.vector(results_df_bead[barcode, "first_type"])
      type2 <- as.vector(results_df_bead[barcode, "second_type"])

      if(is.na(type2)){ # highly confident singlet -> Do Not Purify!
        nz_idx <- which(as.numeric(bead) != 0)
        vec_sparse <- Matrix::sparseVector(i = nz_idx, x = bead[nz_idx], length = length(bead))
        return(list(barcode = barcode, res = vec_sparse))
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
      return(list(barcode = barcode, res = doub_res))
    }, error = function(e) {
      message(sprintf("Error processing barcode %s: %s", barcode, e$message))
      return(NULL)
    })
  }

  # Helper function to process chunks
  process_chunks <- function(
    barcodes, decompose_func,
    parallel = FALSE,
    BPPARAM = NULL
  ) {

    if (length(barcodes) == 0) {
      warning("No barcodes provided. Returning empty result.")
      return(list())
    }

    results_list <- list()

    for (i in seq(1, length(barcodes), by = chunk_size)) {
      cat(round(i * 100 / length(barcodes)), "%\n")

      # Select chunk of barcodes
      chunk_barcodes <- barcodes[i:min(i + chunk_size - 1, length(barcodes))]

      # Subset once per chunk to avoid copying large matrix per worker
      counts_chunk     <- counts[gene_list, chunk_barcodes, drop = FALSE]
      counts_chunk     <- as(counts_chunk, "dgCMatrix")
      results_chunk_df <- results_df[chunk_barcodes, , drop = FALSE]
      weights_chunk    <- ct_weights[chunk_barcodes, , drop = FALSE]

      # Build list of per-barcode input for parallel loop
      inputs <- lapply(seq_along(chunk_barcodes), function(j) {
        list(
          count  = counts_chunk[, j, drop = FALSE],
          result = results_chunk_df[j, , drop = FALSE],
          weight = weights_chunk[j, ]
        )
      })

      # Define lightweight worker function
      worker_fun <- function(input) {
        decompose_func(
          input$count,
          input$result,
          input$weight,
          gene_list,
          cell_type_info
        )
      }

      # Run
      chunk_results <- if (parallel) {
        bplapply(inputs, worker_fun, BPPARAM = BPPARAM)
      } else {
        lapply(inputs, worker_fun)
      }

      results_list <- c(results_list, chunk_results)
    }

    results_list <- Filter(Negate(is.null), results_list)
    return(results_list)
  }


  # Preallocate sparse matrix for result list
  build_sparse_result_matrix <- function(results_list, gene_list) {
    message("Building sparse matrix from sparseVectors... \n")

    n_genes <- length(gene_list)
    n_cells <- length(results_list)

    message("Computing N nz... \n")
    nnz_total <- sum(sapply(results_list, function(res) if (!is.null(res$res)) length(res$res@i) else 0))
    message(" DONE\n")

    i_vec <- integer(nnz_total)
    j_vec <- integer(nnz_total)
    x_vec <- numeric(nnz_total)
    col_ids <- character(n_cells)

    counter <- 1L

    for (j in seq_along(results_list)) {
      res <- results_list[[j]]
      if (is.null(res) || is.null(res$res)) next

      sv <- res$res
      if (!inherits(sv, "dsparseVector")) {
        stop(sprintf("Expected dsparseVector, got %s", class(sv)))
      }

      if (length(sv) != n_genes) {
        stop(sprintf("Vector at column %d has %d genes, expected %d", j, length(sv), n_genes))
      }
      if (any(sv@i > n_genes)) {
        warning(sprintf("Invalid sv@i at column %d: max index %d, length %d",
                        j, max(sv@i), length(sv)))
      }

      n_nz <- length(sv@i)
      if (n_nz > 0) {
        idx <- counter:(counter + n_nz - 1)
        i_vec[idx] <- sv@i
        j_vec[idx] <- rep.int(j, n_nz)
        x_vec[idx] <- sv@x
        counter <- counter + n_nz
      }

      col_ids[j] <- res$barcode

      if (j %% 10000 == 0 || j == n_cells) {
        cat(sprintf("Processed %d / %d\n", j, n_cells))
      }
    }

    # Trim to actual used size
    i_vec <- i_vec[1:(counter - 1)]
    j_vec <- j_vec[1:(counter - 1)]
    x_vec <- x_vec[1:(counter - 1)]

    if (length(i_vec) == 0L) {
      warning("No non-zero entries found. Returning empty sparse matrix.")
      return(Matrix::Matrix(0, nrow = n_genes, ncol = n_cells,
                            dimnames = list(gene_list, col_ids), sparse = TRUE))
    }

    sparse_mat <- Matrix::sparseMatrix(
      i = i_vec,
      j = j_vec,
      x = x_vec,
      dims = c(n_genes, n_cells),
      dimnames = list(gene_list, col_ids)
    )

    return(sparse_mat)
  }



  # Process certain doublets
  all_doublet_results <- list()
  nz <- 1

  tak <- Sys.time()
  cat("Processing certain doublets...\n")
  cat(length(doublets_certain), "\n")

  for (r in process_chunks(doublets_certain, decompose_certain, parallel = DO_parallel, BPPARAM = BPPARAM)) {
    all_doublet_results[[nz]] <- r
    nz <- nz + 1
  }

  cat("Processing uncertain doublets...\n")
  cat(length(doublets_uncertain), "\n")

  for (r in process_chunks(doublets_uncertain, decompose_uncertain, parallel = DO_parallel, BPPARAM = BPPARAM)) {
    all_doublet_results[[nz]] <- r
    nz <- nz + 1
  }

  tik <- Sys.time()
  #cat("Purification completed in ", tik-tak)
  #cat("object.size(all_doublet_results): ")
  #cat(object.size(all_doublet_results))

  # Combine results
  tak <- Sys.time()
  cat("Combaning doublets results ...\n")
  purified <- build_sparse_result_matrix(all_doublet_results, gene_list)
  cell_ids <- colnames(purified)

  tik <- Sys.time()
  #cat("Combining results completed in ", tik-tak)
  #cat("object.size(purified): ")
  #cat(object.size(purified))

  # Process singlets
  if(!DO_purify_singlets){
    cat("Processing singlets...\n")
    singlets <- results_df[results_df$spot_class == "singlet", , drop = FALSE]
    res_singlet_mtrx <- Matrix::Matrix(counts[gene_list, rownames(singlets), drop = FALSE], sparse = sparse)
    purified <- cbind(purified, res_singlet_mtrx)

    cell_ids <- c(cell_ids, colnames(res_singlet_mtrx))
  }

  colnames(purified) <- cell_ids
  rownames(purified) <- gene_list

  #purified <- purified[,colnames(counts)]
  cell_meta <- results_df[colnames(purified),]

  if(nrow(cell_meta) > 0){
    cell_meta$purification_status <- "purified"
    if(!DO_purify_singlets){
      cell_meta$purification_status[cell_meta$spot_class == "singlet"] <- "raw"
    }
    cell_meta$purification_status[is.na(cell_meta$second_type)] <- "raw" # Very confident singlets are not purified

    cell_meta$cell_id <- rownames(cell_meta)
  } else {
    warning("No cells were purified!")
  }
  return(list(purified_counts = purified, cell_meta = cell_meta))
}


#' Purify Data Using RCTD Output
#'
#' This function purifies a query count matrix using the RCTD output, handling both certain and uncertain doublets, as well as optionally processing singlets.
#'
#' @param counts A matrix of gene expression counts, where rows represent genes and columns represent cell IDs.
#' @param rctd RCTD output
#' @param DO_purify_singlets Logical; if `TRUE`, singlets will also be purified.
#' @param DO_parallel Logical; if `TRUE`, biocparallel will be used to accelerate computation. For the moment, non-parallel version perform equally well
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

purify <- function(counts, rctd, DO_purify_singlets, gene_list = NULL, DO_parallel = FALSE, n_workers = NULL, chunk_size = 10000) {

  results_df <- rctd@results$results_df

  common_cells <- intersect(colnames(counts), rownames(results_df))
  results_df <- results_df[common_cells, ]

  cell_type_info <- rctd@cell_type_info[[1]]
  ct_weights <- rctd@results$weights
  ct_weights <- ct_weights[common_cells, colnames(cell_type_info[[1]])]

  if(is.null(gene_list)){
    gene_list <- intersect(rownames(counts), rownames(cell_type_info[[1]]))
  } else{
    gene_list <- gene_list[gene_list %in% intersect(rownames(counts), rownames(cell_type_info[[1]]))]
  }

  counts <- counts[gene_list, common_cells]
  print(dim(counts))

  return(purify_counts_with_rctd(
    counts = counts,
    results_df = results_df,
    ct_weights = ct_weights,
    cell_type_info = cell_type_info,
    DO_purify_singlets = DO_purify_singlets,
    DO_parallel = DO_parallel,
    n_workers = n_workers,
    chunk_size = chunk_size
  )
  )
}












