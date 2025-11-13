#' Split Purified Count Data
#'
#' This function splits purified count data into two components based on the first and second cell type annotations. It ensures that negative values in the residual split profile are replaced with zeros.
#'
#' @param counts A matrix of raw count data where rows represent genes and columns represent cells.
#' @param rctd An object containing RCTD results, including spatial and cell type information.
#' @param DO_purify_singlets A logical value indicating whether to purify singlet cells.
#' @param n_workers An optional integer specifying the number of workers for parallel processing. Defaults to `NULL`.
#' @param chunk_size An integer defining the chunk size for processing. Defaults to `10000`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `purified_counts`: A matrix of purified count data, with two components labeled `_1` and `_2` for each cell.
#'     \item `cell_meta`: A data frame containing metadata for the purified cells, including cell type annotations and decomposition order.
#'   }
#'
#' @note The function modifies negative residual values in `purified_2`, replacing them with zeros.
#' @export


split_cells <- function(counts, rctd, DO_purify_singlets, n_workers = NULL, chunk_size = 10000){

  purified <- purify(counts = counts, rctd = rctd, DO_purify_singlets = DO_purify_singlets, n_workers = n_workers, chunk_size = chunk_size)

  purified_1 <- purified$purified_counts
  cell_meta  <- purified$cell_meta
  cell_meta$cell_id <- rownames(cell_meta)

  purified_2 <- counts[rownames(purified_1), colnames(purified_1)] - purified_1

  if(sum(purified_2<0) > 0){
    warning("Residual split profile has ", sum(purified_2<0), "negative values that are replaced with 0s")
    purified_2[purified_2<0] <- 0
  }

  cell_meta_1 <- cell_meta
  cell_meta_1$cell_type <- cell_meta_1$first_type

  colnames(purified_1) <- paste0(colnames(purified_1), "_1")
  rownames(cell_meta_1) <- paste0(rownames(cell_meta_1), "_1")
  cell_meta_1$decomposition_order <- "first"

  cell_meta_2 <- cell_meta
  cell_meta_2$cell_type <- cell_meta_2$second_type
  cell_meta_2$purification_status[cell_meta_2$purification_status == "raw"] <- "null"
  colnames(purified_2) <- paste0(colnames(purified_2), "_2")
  rownames(cell_meta_2) <- paste0(rownames(cell_meta_2), "_2")
  cell_meta_2$decomposition_order <- "second"


  return(list(
    purified_counts = cbind(purified_1, purified_2),
    cell_meta = rbind(cell_meta_1, cell_meta_2)
  ))
}


#' Split Purified Count Data
#' @export


balance_split <- function(
    xe_raw,
    xe_purified,
    spot_class_key = "spot_class",
    DO_purify_singlets = TRUE,
    DO_split_singlets = TRUE,
    DO_split_doublets_uncertain = FALSE, # keep first only
    default_assay = "Xenium"
){

  if(!spot_class_key %in% colnames(xe_raw@meta.data)){
    stop("spot_class_key", spot_class_key, "is not available in `xe_raw`, please compute it first!")
  }

  if(!"first_type" %in% colnames(xe_raw@meta.data)){
    stop("`first_type` is not available in `xe_raw`, please compute it first!")
  }

  cells_to_remove <-
    xe_raw@meta.data %>%
    filter(.data[[spot_class_key]] == "reject") %>%
    rownames()

  if(DO_purify_singlets){
    cells_to_keep_raw <-
      xe_raw@meta.data %>%
      filter(is.na(second_type)) %>%
      rownames()
    cat(cells_to_keep_raw[1:10], "\n")
    cat(length(cells_to_keep_raw))
  } else {
    cells_to_keep_raw <-
      xe_raw@meta.data %>%
      filter(.data[[spot_class_key]] == "singlet") %>%
      rownames()
    cat(cells_to_keep_raw[1:10], "\n")
    cat(length(cells_to_keep_raw))
  }
  cells_to_keep_raw <- setdiff(cells_to_keep_raw, cells_to_remove)

  cells_to_replace_with_purified <- setdiff(colnames(xe_raw), c(cells_to_keep_raw, cells_to_remove))
  cat(length(cells_to_replace_with_purified), " - all cells_to_replace_with_purified")
  cells_to_replace_with_purified <- cells_to_replace_with_purified[cells_to_replace_with_purified %in% colnames(xe_purified)]
  cat(length(cells_to_replace_with_purified), " - in colnames(xe_purified)")

  xe_raw$purification_status <- "raw"
  xe_raw@meta.data[cells_to_replace_with_purified,"purification_status"] <- xe_purified@meta.data[cells_to_replace_with_purified,"purification_status"]
  xe_raw@meta.data[cells_to_remove,"purification_status"] <- "removed"

  common_genes <- intersect(rownames(xe_raw), rownames(xe_purified))

  message("raw... \n")
  raw <- GetAssayData(xe_raw, assay = default_assay, layer = "counts")[common_genes, cells_to_keep_raw]
  cell_meta_raw <- xe_raw@meta.data[cells_to_keep_raw,]
  cell_meta_raw$decomposition_order <- "raw"
  cell_meta_raw$cell_id <- rownames(cell_meta_raw)
  cell_meta_raw$cell_type <- cell_meta_raw$first_type

  message("pur_1 ... \n")
  cat(cells_to_replace_with_purified[1:5], "\n")
  cat(length(cells_to_replace_with_purified))
  purified_1 <- GetAssayData(xe_purified, assay = default_assay, layer = "counts")[common_genes, cells_to_replace_with_purified]
  cell_meta  <- xe_purified@meta.data[colnames(purified_1), ]
  cell_meta$cell_id <- rownames(cell_meta)

  message("pur_2 ... \n")
  purified_2 <- GetAssayData(xe_raw, assay = default_assay, layer = "counts")[rownames(purified_1), colnames(purified_1)] - purified_1
  if(sum(purified_2<0) > 0){
    warning(purified_2[purified_2<0] %>% as.vector() %>% summary())
    warning("Residual split profile has ", sum(purified_2<0), "negative values that are replaced with 0s")
    purified_2[purified_2<0] <- 0
  }

  cell_meta_1 <- cell_meta
  cell_meta_1$cell_type <- cell_meta_1$first_type
  colnames(purified_1) <- paste0(colnames(purified_1), "_1")
  rownames(cell_meta_1) <- paste0(rownames(cell_meta_1), "_1")
  cell_meta_1$decomposition_order <- "first"

  cell_meta_2 <- cell_meta
  cell_meta_2$cell_type <- cell_meta_2$second_type
  cell_meta_2$purification_status[cell_meta_2$purification_status == "raw"] <- "null"
  colnames(purified_2) <- paste0(colnames(purified_2), "_2")
  rownames(cell_meta_2) <- paste0(rownames(cell_meta_2), "_2")
  cell_meta_2$decomposition_order <- "second"


  if(!DO_split_singlets){
    if(DO_purify_singlets){
      # remove second profile of singlets
      cells_to_keep_second_profile <- cell_meta_2 %>% filter(spot_class != "singlet") %>% rownames()

      purified_2 <- purified_2[,cells_to_keep_second_profile]
      cell_meta_2 <- cell_meta_2[cells_to_keep_second_profile,]
    }
  } else {
    if(!DO_purify_singlets){
      warning("Cannot split singlets as they are raw")
    }
  }

  if(!DO_split_doublets_uncertain){
    cells_to_keep_second_profile <- cell_meta_2 %>% filter(spot_class != "doublet_uncertain") %>% rownames()

    purified_2 <- purified_2[,cells_to_keep_second_profile]
    cell_meta_2 <- cell_meta_2[cells_to_keep_second_profile,]
  }

  common_cols <- intersect(colnames(cell_meta_raw), colnames(cell_meta_1))

  count_matrix <- cbind(raw, cbind(purified_1, purified_2))
  meta_data = rbind(cell_meta_raw[,common_cols], rbind(cell_meta_1[,common_cols], cell_meta_2[,common_cols]))

  xe_balanced <- CreateSeuratObject(counts = count_matrix, assay = default_assay, meta.data = meta_data)
  return(xe_balanced)
}


