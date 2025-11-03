#' Compute Transcriptomic and Spatial Swapping Scores (not used, for exploratory purpose only)
#'
#' @description
#' Computes per-cell transcriptomic and spatial "swapping" scores that quantify
#' inconsistencies between the RCTD cell-type assignment and the cell’s local
#' transcriptomic and spatial neighborhood composition. The function combines
#' results from RCTD deconvolution with both transcriptomic and spatial
#' neighborhood summaries to identify potential cell-type misassignments or
#' mixed signals.
#'
#' @param rctd A output of \link[spacexr]{run.RCDT}
#' containing a `results$results_df` slot with per-spot deconvolution metrics.
#' @param tr_neigh_df A `data.frame` containing transcriptomic neighborhood
#' information per cell. Must include a `cell_id` column and neighborhood-based
#' annotations such as `first_type_neighborhood`,
#' `first_type_neighborhood_agreement`, etc.
#' @param sp_neigh_df A `data.frame` containing spatial neighborhood
#' information per cell. Must include a `cell_id` column and neighborhood-derived
#' features such as `first_type_neighbors_N`, `second_type_neighbors_N`,
#' and neighborhood weights.
#'
#' @details
#' This function:
#' 1. Merges transcriptomic and spatial neighborhood data with the RCTD results table.
#' 2. Computes two separate scores:
#'    - **Transcriptomic swapping score**: measures inconsistency between the
#'      assigned (`first_type`, `second_type`) and neighborhood-predicted
#'      cell types or classes.
#'    - **Spatial swapping score**: compares weighted neighbor support for the
#'      first and second assigned types.
#' 3. Combines these into a `total_swapping_score` and includes additional summary
#'    statistics (`max_swapping_score`, `min_swapping_score`).
#'
#' The scores can be used to visualize or filter potentially misassigned or
#' ambiguous cells in spatial transcriptomics datasets.
#'
#' @return
#' A `data.frame` with one row per cell and the following additional columns:
#' \describe{
#'   \item{transcriptomic_swapping_score}{Numeric score summarizing transcriptomic-level inconsistencies.}
#'   \item{spatial_swapping_score}{Difference between neighborhood weights of first and second types.}
#'   \item{total_swapping_score}{Sum of transcriptomic and spatial swapping scores.}
#'   \item{max_swapping_score, min_swapping_score}{Summary extrema across both score types.}
#' }
#'
#' @examples
#' \dontrun{
#' df_int <- compute_swapping_score(
#'   rctd = RCTD,
#'   tr_neigh_df = tr_neigh_df,
#'   sp_neigh_df = sp_neigh_df
#' )
#' }
#'
#' @export

compute_swapping_score <- function(rctd, tr_neigh_df, sp_neigh_df){

  df <- RCTD@results$results_df
  rownames(tr_neigh_df) <- tr_neigh_df$cell_id
  rownames(sp_neigh_df) <- sp_neigh_df$cell_id

  df <- cbind(df, tr_neigh_df[rownames(df),])
  df <- cbind(df, sp_neigh_df[rownames(df),])

  col_of_interest <- c("cell_id", "spot_class", "first_type", "second_type", "first_class", "second_class",
                       "singlet_score", "w1_larger_w2", "weight_first_type",
                       "annot_min_singlet_score", "annot_max_weight", "annot_max_doublet_weight",
                       "first_type_neighborhood", "first_type_neighborhood_agreement", "first_type_class_neighborhood", "first_type_class_neighborhood_agreement",
                       "second_type_neighborhood", "second_type_neighborhood_agreement", "second_type_class_neighborhood", "second_type_class_neighborhood_agreement",
                       "annotated_neighbors_N", "first_type_neighbors_N", "second_type_neighbors_N", "same_second_type_neighbors_N", "neighborhood_weights_first_type", "neighborhood_weights_second_type", "sum_nCount_neighborhood_spilling_type")

  df_int <- df[, col_of_interest]

  # First, identify candidates for swapping:
  # 1 - they have secondary cell type,
  # 2 - `first_type_neighborhood_agreement == FALSE`
  # 3 -

  # Swap if second_type = first_type_neighborhood & second_type_neighbors_N == 0

  # Swapping is only possible for cells that have a secondary cell type, if there is no, make score either NA to -inf

  ## Transcriptomic scores
  # !is.na(second_type) -> +1
  # argmax first type tr neigh == second_type -> +1 or +2
  # agrmax  first type tr neigh == first_type -> -1 or -2

  df_int <- df_int %>%
    rowwise() %>%
    mutate(
      transcriptomic_swapping_score =
        sum(
          2*(second_type == first_type_neighborhood),
          (second_class == first_type_class_neighborhood),
          (first_type == second_type_neighborhood),
          (first_class == second_type_class_neighborhood),
          (second_type == annot_min_singlet_score),
          (second_type == annot_max_weight),
          (second_type == annot_max_doublet_weight),
          na.rm = TRUE
        )/8 -
        sum(
          2*(first_type == first_type_neighborhood),
          2*(second_type != first_type_neighborhood),
          (second_class == first_type_class_neighborhood),
          (second_type == second_type_neighborhood),
          (second_class == second_type_class_neighborhood),
          (first_type == annot_min_singlet_score),
          (first_type == annot_max_weight),
          (first_type == annot_max_doublet_weight),
          na.rm = TRUE
        )/10,

      spatial_swapping_score =
        neighborhood_weights_first_type - neighborhood_weights_second_type,
      # sum(
      #   (second_type_neighbors_N == 0),
      #   (first_type_neighbors_N > 0),
      #   first_type_neighbors_N,
      #   na.rm = TRUE
      # ) -
      # sum(
      #   (second_type_neighbors_N > 0),
      #   second_type_neighbors_N,
      #   na.rm = TRUE
      # ),

      total_swapping_score = transcriptomic_swapping_score + spatial_swapping_score,
      max_swapping_score = max(transcriptomic_swapping_score, spatial_swapping_score),
      min_swapping_score = max(transcriptomic_swapping_score, spatial_swapping_score)
    ) %>%
    ungroup()


  # argmax first class ft neigh == second class -> +1
  # argmax first class ft neigh == first class -> -1

  # argmax second type tr neigh == second_type -> -1 or -2
  # agrmax  second type tr neigh == first_type -> +1 or +2

  # argmax second class ft neigh == second class -> -1
  # argmax second class ft neigh == first class -> +1

  ## Spatial scores
  # n second type neighbors < 1 -> +1 or +2
  # n second class neighbors < 1 -> +1 or +2

  # n first type neighbors > 1 -> +1

  # n

  # return a vector of indications or sub-scores for each cell or directly a summed up score (better the first one, so that we can visuallize and investigate)
  df_int <- as.data.frame(df_int)
  rownames(df_int) <- df_int$cell_id
  return(df_int)

}

