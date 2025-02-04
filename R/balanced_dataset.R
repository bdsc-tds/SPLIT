
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

balance_raw_and_purified_data_by_score <- function(
    xe_raw,
    xe_purified,
    threshold = .15,
    score_name = "neighborhood_weights_second_type",
    spot_class_key = "spot_class"
){
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

  xe_raw$purification_status <- "raw"
  xe_raw@meta.data[cells_to_replace_with_purified,"purification_status"] <- xe_purified@meta.data[cells_to_replace_with_purified,"purification_status"]
  xe_raw@meta.data[cells_to_remove,"purification_status"] <- "removed"

  common_genes <- intersect(rownames(xe_raw), rownames(xe_purified))
  count_matrix <- cbind(GetAssayData(xe_raw, layer = "counts")[common_genes, cells_to_keep_raw],
                        GetAssayData(xe_purified, layer = "counts")[common_genes, cells_to_replace_with_purified])

  xe_balanced <- CreateSeuratObject(counts = count_matrix, meta.data = xe_raw@meta.data[colnames(count_matrix),])

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
    spot_class_key = "spot_class"
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
  count_matrix <- cbind(GetAssayData(xe_raw, layer = "counts")[common_genes, cells_to_keep_raw],
                        GetAssayData(xe_purified, layer = "counts")[common_genes, cells_to_replace_with_purified])

  xe_balanced <- CreateSeuratObject(counts = count_matrix, meta.data = xe_raw@meta.data[colnames(count_matrix),])
  return(xe_balanced)
}
