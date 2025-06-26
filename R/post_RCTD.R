#' Updates `score_mat` of the `Run.RCTD` output  (Deprecated)
#'
#' Adds `singlet_scores` as a diagonal and removes cell types with a low weight in full cell-type decomposition.
#'
#' This function modifies the `score_mat` by adding `singlet_scores` as a diagonal matrix and removes cell types
#' from the candidate list if their corresponding weight is below a user-defined threshold. The updated
#' `score_mat` and `singlet_scores` are saved in the `rctd` object.
#'
#' @param rctd An object resulting from \link[spacexr]{Run.RCTD}.
#' @param min_weight A threshold (numeric) to keep cell types as candidates. Cell types with a weight below this
#'        threshold are removed from the `score_mat`. Default is 0.01, which is the same as in `Run.RCTD`.
#' @param verbose Logical. If `TRUE`, the function will print messages about removed low-weight cell types.
#'        Default is `FALSE`.
#'
#' @return An updated `rctd` object with modified `score_mat` and `singlet_scores`.

update_score_mat_RCTD <- function(
    rctd,
    min_weight = .01, # same as in run.RCTD, use higher values to remove non-relevant candidates
    verbose = FALSE,
    BPPARAM = bpparam(),
    n_workers = NULL
){
  # Deprecated: This function is no longer used and will be removed in a future release.
  intermediate_function <- function(...) {
    warning("`update_score_mat_RCTD()` is deprecated and no longer used. It will be removed in a future version.", call. = FALSE)
    invisible(NULL)
  }

  score_mat        <- rctd@results$score_mat
  weights          <- rctd@results$weights
  singlet_scores   <- rctd@results$singlet_scores

  # update score_mat by adding `singlet_scores` as a diagonal and then removing
  # candidate cell types that do not have sufficient weight

  cell_types <- colnames(weights)

  if (is.null(n_workers)) {
    n_workers <- min(4, BiocParallel::multicoreWorkers() - 1)
  }


  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers)
  }


  result_list <- bplapply(seq_along(score_mat), function(i) {
    smat <- as.matrix(score_mat[[i]])
    diag(smat) <- singlet_scores[[i]]

    new_candidates <- cell_types[weights[i, ] > min_weight]
    cur_ct <- colnames(smat)
    keep_ct <- intersect(cur_ct, new_candidates)

    smat_sub <- smat[keep_ct, keep_ct, drop = FALSE]
    sscore_sub <- singlet_scores[[i]][keep_ct]

    list(score_mat = smat_sub, singlet_scores = sscore_sub)
  }, BPPARAM = BPPARAM)

  score_mat_xe <- lapply(result_list, function(x) x$score_mat)
  singlet_scores_xe <- lapply(result_list, function(x) x$singlet_scores)

  rctd@results$score_mat_xe <- score_mat_xe
  rctd@results$singlet_scores_xe <- singlet_scores_xe
  return(rctd)
}

#' Correct First Type in Very Confident Singlets
#'
#' This function updates the `first_type` and `spot_class` fields in the results of an RCTD object for cells classified as highly confident singlets. It ensures accurate cell-type annotations for singlets with only one candidate cell type and updates related fields accordingly.
#'
#' @param rctd An `RCTD` object containing spatial transcriptomics annotation results. The object must have the `singlet_scores_xe` field in `rctd@results`. Run `update_score_mat()` prior to using this function if the field is missing.
#' @param min_weight Minimum weight threshold for the score matrix update (default: 0.05).
#'
#' @return An updated `RCTD` object with the following fields modified in `rctd@results$results_df_xe`:
#' \itemize{
#'   \item \code{first_type}: Updated for confident singlets to reflect correct annotation.
#'   \item \code{second_type}: Set to \code{NA} for confident singlets.
#'   \item \code{spot_class}: Updated to "singlet" for highly confident cells.
#'   \item \code{max_doublet_weight}: Added to capture the maximum weight for doublets.
#'   \item \code{n_candidates}: Number of candidate cell types for each spot.
#'   \item \code{rctd_weights_entropy}: Entropy of rctd weights in full cell-type decomposition.
#'
#' }
#'
#' @details
#' The function identifies cells with only one candidate cell type (i.e., highly confident singlets). It corrects the `first_type` assignment, removes the `second_type` (assignes to NA), and updates the `spot_class` to "singlet." This ensures proper classification of confident cells while leaving other cells unchanged.
#'
#' @examples
#' \dontrun{
#' rctd <- update_score_mat(rctd)
#' rctd <- correct_singlets(rctd)
#' }
#'
#' @import dplyr


correct_singlets <- function(
    rctd,
    min_weight = 0.01
){
  #if(!("singlet_scores_xe" %in% names(rctd@results)))
  # stop("No `singlet_scores_xe` field in rctd@results, run `update_score_mat()` first!")

  if("scond_type" %in% colnames(rctd@results$results_df)){
    rctd@results$results_df$second_type <- rctd@results$results_df$scond_type
  }
  df <- rctd@results$results_df

  # select cells with only one candidate cell type
  len_mat <-  apply(rctd@results$weights > min_weight, 1, FUN = sum) %>% as.numeric() # number of cell type candidates (dont use singlet score, as that one has not been updated)
  confident_singlet_idx <- which(len_mat == 1)
  no_cell_type_idx <- which(len_mat == 0)
  rejects_idx <- which(rctd@results$results_df$spot_class == "reject")

  #names(rctd@results$singlet_scores_xe) <- rctd@results$results_df %>% rownames()

  #get their singlet scores
  argmax_weight <- apply(rctd@results$weights,
                         1,
                         function(x){
                           r <- x %>% which.max() %>% names()
                           if(is.null(r))
                             r <- NA
                           return(r)
                         }, simplify = T) %>% unlist()

  # update first type to make sure highly confined cells have correct annotation
  first_type_updated <- rctd@results$results_df$first_type
  first_type_updated[confident_singlet_idx] <- argmax_weight[confident_singlet_idx]
  first_type_updated[no_cell_type_idx] <- NA_character_

  # replace second type with NA as now RCTD assigns it to a random cell type (usually the first of the available cell types)
  second_type_updated <- rctd@results$results_df$second_type
  second_type_updated[confident_singlet_idx] <- NA_character_
  second_type_updated[no_cell_type_idx] <- NA_character_

  # update spot class of highly confident cells to singlets
  spot_class_upd <- rctd@results$results_df$spot_class
  spot_class_upd[setdiff(confident_singlet_idx, rejects_idx)] <- "singlet"

  spot_class_upd[no_cell_type_idx] <- "reject"
  spot_class_upd <- ordered(spot_class_upd, levels = c("reject", "doublet_uncertain", "doublet_certain", "singlet"))
  max_doublet_weight <- apply(rctd@results$weights_doublet, 1, max) %>% unname()

  # compute entropy of rctd weight score
  weights <- rctd@results$weights
  weights[weights<0] <- 0
  weights <- spacexr::normalize_weights(weights)
  weights_entr <- apply(weights, 1, entropy::entropy) %>% unname()

  df <- rctd@results$results_df %>%
    mutate(
      first_type = first_type_updated,
      second_type = second_type_updated,
      max_doublet_weight = max_doublet_weight,
      spot_class = spot_class_upd,
      n_candidates = len_mat,
      rctd_weights_entropy = weights_entr
    )

  rctd@results$results_df_xe <- df
  return(rctd)
}


#' Update scores in RCTD results data frame
#'
#' This function updates existing scores in `rctd@results$results_df_xe`, computes new scores, and adds flags to facilitate further analysis. The function integrates additional metadata, refines singlet scores for cell type annotations, computes score differences, and checks for relationships between cell types.
#'
#' @param rctd An `RCTD` object containing spatial transcriptomics data and results.
#'
#' @return An updated `RCTD` object with modified and new fields in `rctd@results$results_df_xe`, including:
#' \itemize{
#'   \item \code{singlet_score_first}: Singlet score for the cell type assigned to `first_type`.
#'   \item \code{singlet_score_second}: Singlet score for the cell type assigned to `second_type` (if applicable).
#'   \item \code{delta_singlet_score_first_second}: Difference between \code{singlet_score_first} and \code{singlet_score_second}.
#'   \item \code{score_diff}: Difference between \code{singlet_score_first} and \code{min_score}.
#'   \item \code{delta_singlet_score}: Difference between the smallest and second smallest singlet scores.
#'   \item \code{delta_singlet_score_class}: Similar to \code{delta_singlet_score}, but excludes scores from the same class as \code{first_type}.
#'   \item \code{weight_first_type}: Weight for the cell type assigned to `first_type`.
#'   \item \code{weight_second_type}: Weight for the cell type assigned to `second_type` (if applicable).
#'   \item \code{same_class}: Logical flag indicating whether `first_type` and `second_type` belong to the same class.
#' }
#'
#' @details
#' The function refines singlet scores for `first_type` and `second_type`, computes score differences, and introduces new metrics such as `delta_singlet_score` and `delta_singlet_score_class`. These metrics help distinguish between highly confident cell types and ambiguous cases. Additionally, the function evaluates whether `first_type` and `second_type` belong to the same class, using an internal class mapping.
#'

update_scores_RCTD <- function(rctd, lite = TRUE){

  df <- rctd@results$results_df_xe

  if(!lite){
    # Ensure first_type and second_type are available in the data frame
    df <- df %>% mutate(
      singlet_score_first = sapply(1:nrow(df), function(i) {
        ft <- df$first_type[i] %>% as.vector()
        return(rctd@results$singlet_scores[[i]][ft] %>% unname())
      }),
      singlet_score_second = sapply(1:nrow(df), function(i) {
        st <- df$second_type[i] %>% as.vector()
        if (is.na(st)) return(NA)
        return(rctd@results$singlet_scores[[i]][st] %>% unname())
      }),
      delta_singlet_score_first_second = singlet_score_second - singlet_score_first
    )

    # Update score_diff to be consistent with the new singlet_score
    df <- df %>% mutate(
      score_diff = df$singlet_score_first - df$min_score,  # Use this one
      score_diff_old = df$singlet_score - df$min_score,  # Keep this to track
      delta_singlet_score_original_first_class = singlet_score_first - singlet_score  # Difference between singlet_score and singlet_score_first
    )
  }

  # Add weight_first_type, weight_second_type
  df <- df %>% mutate(
    weight_first_type = rctd@results$weights_doublet[,"first_type"] %>% unname(),
    weight_second_type = rctd@results$weights_doublet[,"second_type"] %>% unname()
  )

  if(!lite){
    # Calculate delta_singlet_score: difference between first and second smallest singlet scores
    sorted_singlet_scores <- sapply(rctd@results$singlet_scores, FUN = function(x) sort(x, decreasing = F))
    delta_singlet_score <- sapply(sorted_singlet_scores, FUN = function(x) {
      if (length(x) == 1) return(Inf)
      return(x[2] - x[1])
    })
    df$delta_singlet_score <- delta_singlet_score

    # Calculate delta_singlet_score_class (ignore same class elements when computing delta score)
    sorted_singlet_scores_class <- sapply(sorted_singlet_scores, function(x) {
      class_vec <- rctd@internal_vars$class_df[names(x), "class"]
      mask <- class_vec != class_vec[1]  # Keep elements where their class != class[1]
      mask[1] <- TRUE  # Keep the first element
      return(x[mask])
    })

    delta_singlet_score_class <- sapply(sorted_singlet_scores_class, FUN = function(x) {
      if (length(x) == 1) return(0)
      return(x[2] - x[1])
    })

    df$delta_singlet_score_class <- delta_singlet_score_class
  }

  # Check whether first_type and second_type come from the same class
  df <- df %>%
    mutate(
      first_type_class = factor(
        rctd@internal_vars$class_df[df$first_type %>% as.vector(), "class"],
        levels = unique(rctd@internal_vars$class_df$class)
      ),
      second_type_class = factor(
        rctd@internal_vars$class_df[df$second_type %>% as.vector(), "class"],
        levels = unique(rctd@internal_vars$class_df$class)
      ),
      same_class =  first_type_class == second_type_class
    )

  rctd@results$results_df_xe <- df
  return(rctd)
}


#' Normalize `score_diff` to Compensate for Feature Count and Update `spot_class`
#'
#' This function normalizes `score_diff` in the RCTD results to account for the positive correlation between the number of features (`nFeature`) and the goodness-of-fit metrics (`singlet_score` and `min_score`). The normalization helps to reduce the bias introduced by feature count and adjusts the `spot_class` accordingly.
#'
#' @param rctd An `RCTD` object containing spatial transcriptomics data and results.
#' @param nFeature_doublet_threshold Numeric. Threshold for `score_diff_normalized` to classify a spot as a singlet. Default is 0.5.
#' @param nFeature a vector of nFeatures for `rctd@results$results_df_xe` rows, if not provided, will be computed from `rctd@spatialRNA@counts`
#' @param nCount a vector of nCount for `rctd@results$results_df_xe` rows, if not provided, will be computed from `rctd@spatialRNA@counts`
#'
#' @return An updated `RCTD` object with the following fields added or updated in `rctd@results$results_df_xe`:
#' \itemize{
#'   \item \code{nCount}: The count of reads or molecules for each spot, retrieved from \code{xe}.
#'   \item \code{nFeature}: The number of features (e.g., genes) detected for each spot, retrieved from \code{xe}.
#'   \item \code{score_diff_normalized}: The normalized \code{score_diff}, computed as \code{score_diff / nFeature}.
#'   \item \code{is_singlet_in_normalized_thresh}: Logical flag indicating whether a spot's normalized score falls below the doublet threshold.
#'   \item \code{spot_class_normalized}: Updated `spot_class` that classifies spots as "singlet" if they pass the normalized threshold and were not previously rejected.
#' }
#'
#' @details
#' The function compensates for the effect of feature count (`nFeature`) on `score_diff` by normalizing it. Spots are reclassified as singlets if their normalized score (`score_diff_normalized`) is below the threshold specified by \code{nFeature_doublet_threshold}, and they are not rejected.


normalize_score_diff_by_nFeature <- function(
    rctd,
    nFeature_doublet_threshold = 0.5,
    nFeature = NULL,
    nCount = NULL
) {

  if(is.null(nFeature) || is.null(nCount)){
    # nCount, nFeature
    rctd@results$results_df_xe <- rctd@results$results_df_xe %>% mutate(
      nCount = Matrix::colSums(rctd@spatialRNA@counts) %>% unname(),
      nFeature = Matrix::colSums(rctd@spatialRNA@counts > 0) %>% unname()
    )
  } else {
    rctd@results$results_df_xe <- rctd@results$results_df_xe %>% mutate(
      nCount = nCount,
      nFeature = nFeature
    )
  }

  # Normalization
  rctd@results$results_df_xe <- rctd@results$results_df_xe %>%
    mutate(
      score_diff_normalized = score_diff / nFeature,
      is_singlet_in_normalized_thresh = score_diff_normalized < nFeature_doublet_threshold,
      spot_class_normalized = if_else(is_singlet_in_normalized_thresh & spot_class != "reject",
                                      "singlet", spot_class) %>% ordered(levels = levels(spot_class))
    )
  return(rctd)
}

#' Compute Alternative Annotations
#'
#' This function computes three alternative annotations for each cell based on the results in the `rctd` object:
#'
#' - `annot_min_singlet_score`: The cell type with the minimum singlet score.
#' - `annot_max_weight`: The cell type with the maximum weight.
#' - `annot_max_doublet_weight`: The cell type with the maximum doublet weight, choosing the first type if the first weight is greater than the second, or the second type otherwise.
#'
#' These annotations are added to the `results_df_xe` slot of the `rctd` object.
#'
#' @param rctd An object of class `RCTD` containing the results from the spatial transcriptomics analysis.
#'
#' @return The input `rctd` object with the updated `results_df_xe` slot containing the three alternative annotations.
#'

compute_alternative_annotations <- function(rctd){

  # Compute alternative annotations
  annot_min_singlet_score <- sapply(rctd@results$singlet_scores, function(x) {
    names(x)[which.min(x)]  # directly find the name with min singlet score
  }) %>% unname()

  annot_max_weight <- apply(rctd@results$weights, 1, function(x) {
    names(x)[which.max(x)] # find the name with max weight
  }) %>% unname()

  annot_max_doublet_weight <- ifelse(
    rctd@results$results_df_xe$weight_first_type > rctd@results$results_df_xe$weight_second_type |
      is.na(rctd@results$results_df_xe$weight_second_type),
    rctd@results$results_df_xe$first_type %>% as.vector(),
    rctd@results$results_df_xe$second_type %>% as.vector()
  ) %>% unname()

  rctd@results$results_df_xe <- rctd@results$results_df_xe %>%
    mutate(
      annot_min_singlet_score = annot_min_singlet_score,
      annot_max_weight = annot_max_weight,
      annot_max_doublet_weight = annot_max_doublet_weight,
      w1_larger_w2 = first_type == annot_max_doublet_weight
    )
  return(rctd)
}

#' Computes Normalized Entropy of Annotations and Stores in `rctd@results$results_df_xe`
#'
#' This function computes the normalized entropy of the annotations for each cell type in the `rctd@results$results_df_xe` dataframe.
#' The entropy is computed based on the relative frequencies of each annotation field and normalized by the maximum possible entropy.
#'
#' The following entropy values are computed:
#' - `entropy_first_type`: Normalized entropy for the `first_type` annotation, considering annotations such as `annot_min_singlet_score`, `annot_max_weight`, and `annot_max_doublet_weight`.
#' - `entropy_second_type`: Normalized entropy for the `second_type` annotation, computed in the same way as `entropy_first_type`.
#'
#' Entropy values indicate the diversity of the annotations for each cell type. A high entropy value suggests more diverse annotations, while a low value indicates more certainty in the annotation.
#'
#' @param rctd An object of class `RCTD` containing the results of the spatial transcriptomics analysis.
#'
#' @return The input `rctd` object with two additional columns in `results_df_xe`:
#' - `entropy_first_type`: The normalized entropy for the `first_type` annotation.
#' - `entropy_second_type`: The normalized entropy for the `second_type` annotation.


compute_annotation_entropy <- function(rctd){
  annotation_fields <- grep("annot_", colnames(rctd@results$results_df_xe), value = T)

  normalized_entropy <- function(x){
    vec <- x %>% factor() %>% forcats::fct_infreq() %>% as.numeric()
    freqs <- table(vec) / length(vec)
    entropy_value <- entropy::entropy(freqs, unit = "log2")
    # Calculate maximum entropy
    max_entropy <- log2(length(vec))

    if(is.na(entropy_value)) # missing argmax doublet weight because the first type has low doublet weight,but very high overall weight such that ther is only one candidate cell type (instability in RCTD decomposition)
      entropy_value <- 0

    return(entropy_value / max_entropy)
  }

  entropy_first_type <- rctd@results$results_df_xe %>% select(all_of(c("first_type", annotation_fields))) %>%
    apply(1, FUN = normalized_entropy)
  entropy_second_type <- rctd@results$results_df_xe %>% select(all_of(c("second_type", annotation_fields))) %>%
    apply(1, FUN = normalized_entropy)

  rctd@results$results_df_xe$entropy_first_type <- entropy_first_type
  rctd@results$results_df_xe$entropy_second_type <- entropy_second_type
  return(rctd)
}

#' Computes Annotation Confidence Based on Normalized Entropy of Original and Alternative Annotations
#'
#' This function computes the annotation confidence for each cell based on the normalized entropy of the original and alternative annotations (`entropy_first_type` and `entropy_second_type`), as well as the cell's weight (`weight_first_type`).
#' The confidence for each cell type is computed as follows:
#' - `confidence_first_type = (1 - entropy_first_type) * weight_first_type`
#' - `confidence_second_type = (1 - entropy_second_type) * (1 - weight_first_type)`
#'
#' These confidence values reflect the reliability of the annotations, with higher confidence indicating greater certainty in the cell's type assignment. The entropy values capture the diversity of the annotations, and the weights are used to adjust the confidence based on the strength of the cell-type assignment.
#'
#' @param rctd An object of class `RCTD` containing the results of the spatial transcriptomics analysis.
#'
#' @return The input `rctd` object with two additional columns in `results_df_xe`:
#' - `confidence_first_type`: The computed confidence for the `first_type` annotation.
#' - `confidence_second_type`: The computed confidence for the `second_type` annotation.


compute_annotation_confidence <- function(rctd){
  rctd@results$results_df_xe <- rctd@results$results_df_xe %>%
    mutate(
      confidence_first_type = (1 - entropy_first_type) * weight_first_type,
      confidence_second_type = (1 - entropy_second_type) * (1-weight_first_type)
    )
  return(rctd)
}


#' Runs the post-processing pipeline for RCTD.
#'
#' This function sequentially applies a series of processing steps to the RCTD results,
#' including score matrix updates, singlet correction, score updates, normalization,
#' alternative annotation computation, entropy calculation, and annotation confidence computation.
#'
#' The following processing steps are performed in order:
#' - Corrects singlets using `correct_singlets()`
#' - Updates scores with `update_scores_RCTD()`
#' - Normalizes the score difference by the number of features using `normalize_score_diff_by_nFeature()`
#' - Computes alternative annotations using `compute_alternative_annotations()`
#' - Computes normalized entropy for each annotation using `compute_annotation_entropy()`
#' - Computes annotation confidence using `compute_annotation_confidence()`
#'
#' The final results are stored in the `rctd@results$results_df` slot, and the previous results are saved in
#' `rctd@results$results_df_old` for tracking purposes. The intermediate results in `results_df_xe` are cleared.
#'
#' @param rctd RCTD object containing the results.
#' @param min_weight Minimum weight threshold for the score matrix update (default: 0.05).
#' @param nFeature_doublet_threshold Threshold for doublet classification when normalizing the score difference
#'                                  by the number of features (default: 0.5).
#' @param nFeature Optional numeric vector providing the number of features for each cell,
#'                 overrides calculation based on the data if provided (default: NULL).
#' @param nCount Optional numeric vector providing the count of features for each cell,
#'                overrides calculation based on the data if provided (default: NULL).
#' @param lite Logical; if \code{TRUE}, skips computation of many aux scores that are not used downstream and were meant for the exploration stage
#'
#' @return Updated RCTD object with processed results.
#' @export

run_post_process_RCTD <- function(
    rctd,
    min_weight = 0.05,
    nFeature_doublet_threshold = 0.5,
    nFeature = NULL,
    nCount = NULL,
    lite = TRUE
){

  message("Correcting singlets ...")
  rctd <- correct_singlets(rctd = rctd, min_weight = min_weight)

  message("Updating scores ...")
  rctd <- update_scores_RCTD(rctd = rctd, lite = lite)

  message("Add coordinates to results ...")
  rctd@results$results_df_xe$x <- rctd@spatialRNA@coords[rownames(rctd@results$results_df_xe), "x"]
  rctd@results$results_df_xe$y <- rctd@spatialRNA@coords[rownames(rctd@results$results_df_xe), "y"]

  if(!lite){
    message("Normalizing score_diff by nFeature ...")
    rctd <- normalize_score_diff_by_nFeature(
      rctd = rctd,
      nFeature_doublet_threshold = nFeature_doublet_threshold,
      nFeature = nFeature,
      nCount = nCount
    )
  }

  message("Computing alternative annotations ...")
  rctd <- compute_alternative_annotations(rctd)

  #message("Computing annotation entropy ...")
  #rctd <- compute_annotation_entropy(rctd)

  #message("Computing annotation confidence ...")
  #rtrd <- compute_annotation_confidence(rctd)

  message("Replacing results_df ...")
  # store the old results to keep track and replace them with the new ones
  rctd@results$results_df_old <- rctd@results$results_df
  rctd@results$results_df <- rctd@results$results_df_xe
  rctd@results$results_df_xe <- NULL

  return(rctd)
}

