#' Extract Singlet Cell Data
#'
#' This function filters an `xe` object to retain only the cells classified as "singlet"
#' based on the classification stored in the RCTD results.
#'
#' @param xe A Seurat object containing spatial transcriptomics data.
#' @param rctd An object containing RCTD results.
#' @param spot_class_key A character string specifying the column name in `rctd@results$results_df`
#'   that contains spot classification information. Default is `"spot_class"`.
#'
#' @return A subset of `xe` containing only singlet-classified cells.
#'
#' @import spacexr
#' @import Seurat
#' @import dplyr
#' @export
#'
get_singlet_data <- function(
    xe,
    rctd,
    spot_class_key = "spot_class"
){
  results_df <- rctd@results$results_df
  singlets_ids <- results_df %>% filter(.data[[spot_class_key]] == "singlet") %>% rownames()
  xe <- subset(xe, cells = singlets_ids)
  return(xe)
}
