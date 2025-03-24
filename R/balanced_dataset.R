#' @keywords internal
swap_expr <- rlang::expr({
  cells_to_swap_label <- meta_data %>%
    filter(!first_type_neighborhood_agreement & !first_type_class_neighborhood_agreement) %>%
    filter(second_type_class == first_type_class_neighborhood) %>%
    rownames()

  meta_data <- meta_data %>%
    mutate(across(c(first_type, second_type, weight_first_type, weight_second_type),
                  .fns = list(before_swap = identity))) %>%  # Create backup in one step
    mutate(
      swap = rownames(.) %in% cells_to_swap_label,  # Logical column for swapping
      first_type = if_else(swap, second_type_before_swap, first_type),
      second_type = if_else(swap, first_type_before_swap, second_type),
      weight_first_type = if_else(swap, weight_second_type_before_swap, weight_first_type),
      weight_second_type = if_else(swap, weight_first_type_before_swap, weight_second_type)
    )

  # Replace corrected counts for swapped cells
  cells_to_swap_correcred_profile <- cells_to_swap_label[cells_to_swap_label %in% colnames(xe_purified)] #intersect(cells_to_swap_label, cells_to_replace_with_purified)
  count_matrix[, cells_to_swap_correcred_profile] <-
    GetAssayData(xe_raw, assay = default_assay, layer = "counts")[common_genes, cells_to_swap_correcred_profile] -
    GetAssayData(xe_purified, assay = default_assay, layer = "counts")[common_genes, cells_to_swap_correcred_profile]
  count_matrix[count_matrix < 0] <- 0
})

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
#' @param DO_swap_lables A logical indicating whether to swap first and second cell types for cells which label does not agree with its transcriptomic neighborhod lable
#'
#' @export

balance_raw_and_purified_data_by_score <- function(
    xe_raw,
    xe_purified,
    threshold = .15,
    score_name = c("neighborhood_weights_second_type", "second_type_neighbors_N",  "second_type_neighbors_no_reject_N"),
    spot_class_key = "spot_class",
    DO_swap_lables = FALSE,
    default_assay = "Xenium"
){
  score_name <- score_name[1]
  if(!score_name %in% colnames(xe_raw@meta.data)){
    stop("score", score_name, "is not available in `xe_raw`, please compute it first!")
  }

  if(!"first_type" %in% colnames(xe_raw@meta.data)){
    stop("`first_type` is not available in `xe_raw`, please compute it first!")
  }

  if(!spot_class_key %in% colnames(xe_raw@meta.data)){
    stop("`spot_class` is not available in `xe_raw`, please compute it first!")
  }

  cells_to_replace_with_purified <-
    xe_raw@meta.data %>%
    filter((.data[[spot_class_key]] != "reject" & .data[[score_name]] > threshold) |
             .data[[spot_class_key]] == "doublet_uncertain") %>%
    rownames()

  cells_to_remove <-
    xe_raw@meta.data %>%
    filter(.data[[spot_class_key]] == "reject") %>%
    rownames()

  cells_to_keep_raw <- setdiff(colnames(xe_raw), c(cells_to_replace_with_purified, cells_to_remove))
  cells_to_replace_with_purified <- setdiff(cells_to_replace_with_purified, cells_to_remove)
  # restrict to cells present in purified data
  cells_to_replace_with_purified <- cells_to_replace_with_purified[cells_to_replace_with_purified %in% colnames(xe_purified)]

  xe_raw$purification_status <- "raw"
  xe_raw@meta.data[cells_to_replace_with_purified,"purification_status"] <- xe_purified@meta.data[cells_to_replace_with_purified,"purification_status"]
  xe_raw@meta.data[cells_to_remove,"purification_status"] <- "removed"

  common_genes <- intersect(rownames(xe_raw), rownames(xe_purified))
  count_matrix <- cbind(GetAssayData(xe_raw, assay = default_assay, layer = "counts")[common_genes, cells_to_keep_raw],
                        GetAssayData(xe_purified, assay = default_assay, layer = "counts")[common_genes, cells_to_replace_with_purified])

  meta_data = xe_raw@meta.data[colnames(count_matrix),]

  if(DO_swap_lables){
    eval(swap_expr)  # Expands the swap logic in-place
  } else {
    meta_data$swap <- FALSE
  }

  xe_balanced <- CreateSeuratObject(counts = count_matrix, assay = default_assay, meta.data = meta_data)

  return(xe_balanced)
}

#' Balance raw and purified data using spot class categories
#'
#' Merges raw and purified data into one dataset by retaining raw counts
#' for singlet cells and replacing other cell categories with their purified profiles.
#' Reject cells are removed from the dataset.
#'
#' @param xe_raw A raw Seurat object containing the original spatial transcriptomics data.
#' @param xe_purified A purified Seurat object, typically post \code{purify_counts_with_rctd}, containing corrected cell profiles.
#' @param spot_class_key A character string specifying the metadata column name in \code{xe_raw} that indicates the spot classification (default: `"spot_class"`).
#' @param DO_swap_lables A logical indicating whether to swap first and second cell types for cells which label does not agree with its transcriptomic neighborhod lable
#'
#' @return A balanced Seurat object with a merged count matrix where:
#' \itemize{
#'   \item Singlet cells retain their raw counts.
#'   \item Cells classified as non-singlet (but not rejected) are replaced with purified counts.
#'   \item Rejected cells are removed.
#' }
#'
#' @details The function categorizes cells into three groups:
#' \itemize{
#'   \item Cells classified as `"singlet"` retain their raw data.
#'   \item Cells classified as anything other than `"singlet"` (excluding `"reject"`) are replaced with their purified profiles.
#'   \item Cells classified as `"reject"` are completely removed from the dataset.
#' }
#'
#' @export

balance_raw_and_purified_data_by_spot_class <- function(
    xe_raw,
    xe_purified,
    spot_class_key = "spot_class",
    DO_swap_lables = FALSE,
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

  cells_to_replace_with_purified <-
    xe_raw@meta.data %>%
    filter(.data[[spot_class_key]] != "singlet") %>%
    rownames()
  cells_to_replace_with_purified <- setdiff(cells_to_replace_with_purified, cells_to_remove)

  cells_to_keep_raw <-
    xe_raw@meta.data %>%
    filter(.data[[spot_class_key]] == "singlet") %>%
    rownames()


  xe_raw$purification_status <- "raw"
  xe_raw@meta.data[cells_to_replace_with_purified,"purification_status"] <- xe_purified@meta.data[cells_to_replace_with_purified,"purification_status"]
  xe_raw@meta.data[cells_to_remove,"purification_status"] <- "removed"

  common_genes <- intersect(rownames(xe_raw), rownames(xe_purified))
  count_matrix <- cbind(GetAssayData(xe_raw, assay = default_assay, layer = "counts")[common_genes, cells_to_keep_raw],
                        GetAssayData(xe_purified, assay = default_assay, layer = "counts")[common_genes, cells_to_replace_with_purified])

  meta_data = xe_raw@meta.data[colnames(count_matrix),]

  if(DO_swap_lables){
    eval(swap_expr)  # Expands the swap logic in-place
  } else {
    meta_data$swap <- FALSE
  }

  xe_balanced <- CreateSeuratObject(counts = count_matrix, assay = default_assay, meta.data = meta_data)
  return(xe_balanced)
}

