##  ==== 1. add_rctd_to_neighborhood ====

#' Add RCTD Results to Neighborhood Graph
#'
#' Integrates RCTD (Robust Cell Type Decomposition) results into a
#' neighbourhood graph, ensuring cell alignment and reshaping results
#' for downstream use.
#'
#' @param graph A list representing a neighbourhood graph, typically
#'   created by \code{compute_neighborhood}. Must contain:
#'   \describe{
#'     \item{nn_idx}{Matrix of nearest-neighbour indices.}
#'     \item{cell_id}{Character vector of cell identifiers.}
#'   }
#' @param rctd An RCTD object containing cell decomposition results in
#'   \code{rctd@results$results_df}.
#'
#' @return The input \code{graph} list extended with one matrix per
#'   column of \code{rctd@results$results_df}, where rows correspond
#'   to cells and columns to neighbours.
#'
#' @seealso \code{\link{add_spatial_metric}},
#'   \code{\link{add_transcriptomics_metric}}
#'
#' @export
add_rctd_to_neighborhood <- function(graph, rctd) {

  ## ---- Validate required graph fields -----------------------------------
  required    <- c("nn_idx", "cell_id")
  missing_req <- required[!required %in% names(graph)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from 'graph':\n",
      "  ", paste(missing_req, collapse = ", ")
    )
  }

  ## ---- Validate RCTD object ---------------------------------------------
  if (!methods::is(rctd, "RCTD") ||
      is.null(rctd@results$results_df)) {
    stop(
      "'rctd' must be an RCTD object with a non-NULL ",
      "rctd@results$results_df."
    )
  }

  ## ---- Align results to graph cells -------------------------------------
  missing_cells <- setdiff(graph$cell_id,
                           rownames(rctd@results$results_df))
  if (length(missing_cells) > 0L) {
    warning(
      length(missing_cells), " cell(s) in graph$cell_id not found ",
      "in rctd@results$results_df. These will produce NA rows."
    )
  }

  results_df           <- rctd@results$results_df[graph$cell_id, ]
  rownames(results_df) <- graph$cell_id

  ## ---- Reshape to neighbour matrix layout ------------------------------
  results_df      <- results_df[as.vector(graph$nn_idx), ]
  num_neighbors   <- ncol(graph$nn_idx)

  result_list <- lapply(colnames(results_df), function(col) {
    matrix(results_df[[col]], ncol = num_neighbors)
  })
  names(result_list) <- colnames(results_df)

  return(c(graph, result_list))
}


## ==== 2. add_infiltration_metrics_to_neighborhood ====

#' Add Neighborhood Infiltration Metrics to Spatial Neighborhood
#'
#' Calculates infiltration metrics for cell neighbourhoods within a
#' spatial context. Metrics relate to the number and type of neighbours
#' annotated to a cell's first or second type, as well as infiltration
#' metrics based on cell classes.
#'
#' Missing optional fields (\code{first_type_class},
#' \code{second_type_class}) produce a \code{\link{warning}} and the
#' affected metrics are skipped rather than causing an error.
#'
#' @param neighborhood A list containing spatial neighbourhood data.
#'   Required fields:
#'   \describe{
#'     \item{nn_idx}{Matrix of nearest-neighbour indices
#'       (cells x neighbours).}
#'     \item{first_type}{Matrix of first cell-type annotations.}
#'     \item{second_type}{Matrix of second cell-type annotations.}
#'     \item{spot_class}{Matrix of RCTD spot classifications.}
#'   }
#'   Optional fields (skipped with a warning if absent):
#'   \describe{
#'     \item{first_type_class}{Matrix of first cell-type class
#'       annotations.}
#'     \item{second_type_class}{Matrix of second cell-type class
#'       annotations.}
#'   }
#'
#' @return The input \code{neighborhood} list augmented with infiltration
#'   metric fields. Fields depending on missing optional inputs are
#'   omitted rather than causing an error.
#'
#' @importFrom methods is
#'
#' @export
add_infiltration_metrics_to_neighborhood <- function(neighborhood) {

  ## ---- 0. Required field validation -------------------------------------
  required    <- c("nn_idx", "first_type", "second_type", "spot_class")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", "), "\n",
      "Ensure add_rctd_to_neighborhood() was called first."
    )
  }

  ## ---- 0b. Optional field availability ----------------------------------
  has_first_class  <- "first_type_class"  %in% names(neighborhood)
  has_second_class <- "second_type_class" %in% names(neighborhood)

  if (!has_first_class) {
    warning(
      "'first_type_class' not found in neighborhood. ",
      "Metrics 'second_type_class_neighbors*' and ",
      "'first_type_class_neighbors*' will not be computed."
    )
  }
  if (!has_second_class && has_first_class) {
    warning(
      "'second_type_class' not found in neighborhood. ",
      "Metric 'second_type_class_neighbors*' will not be computed."
    )
  }

  ## ---- 1. Total neighbour counts ----------------------------------------
  neighborhood$total_neighbors_N <-
    rowSums(!is.na(neighborhood$nn_idx)) - 1L
  neighborhood$annotated_neighbors_N <-
    rowSums(!is.na(neighborhood$first_type)) - 1L

  ## ---- 2. Singlet neighbours --------------------------------------------
  neighborhood$total_singlets_neighbors_N <- rowSums(
    neighborhood$spot_class == "singlet", na.rm = TRUE
  )

  ## ---- 3. No-reject first_type ------------------------------------------
  neighborhood$first_type_no_reject <- neighborhood$first_type
  neighborhood$first_type_no_reject[
    neighborhood$spot_class == "reject"
  ] <- NA_character_

  ## ---- 4. Combined focal+neighbour matrices -----------------------------
  neighborhood$first_type_fist_element_second_type <- cbind(
    neighborhood$second_type[, 1L, drop = TRUE],
    neighborhood$first_type[, -1L, drop = FALSE]
  )
  neighborhood$first_type_fist_element_second_type_no_reject <- cbind(
    neighborhood$second_type[, 1L, drop = TRUE],
    neighborhood$first_type_no_reject[, -1L, drop = FALSE]
  )

  ## ---- 5. Internal helper -----------------------------------------------
  .get_neighbors <- function(mat) {
    apply(mat, 1L, function(x) which(x[-1L] == x[1L]) + 1L)
  }

  ## ---- 6. Second-type neighbour metrics ---------------------------------
  neighborhood$second_type_neighbors <- .get_neighbors(
    neighborhood$first_type_fist_element_second_type
  )
  neighborhood$second_type_neighbors_N <- vapply(
    neighborhood$second_type_neighbors, length, integer(1L)
  )

  neighborhood$second_type_neighbors_no_reject <- .get_neighbors(
    neighborhood$first_type_fist_element_second_type_no_reject
  )
  neighborhood$second_type_neighbors_no_reject_N <- vapply(
    neighborhood$second_type_neighbors_no_reject, length, integer(1L)
  )

  n_cells <- nrow(neighborhood$first_type_fist_element_second_type)
  neighborhood$second_type_singlets_neighbors <- lapply(
    seq_len(n_cells),
    function(i) {
      x <- neighborhood$first_type_fist_element_second_type[i, ]
      y <- neighborhood$spot_class[i, ]
      which((x[-1L] == x[1L]) & (y[-1L] == "singlet")) + 1L
    }
  )
  neighborhood$second_type_singlets_neighbors_N <- vapply(
    neighborhood$second_type_singlets_neighbors, length, integer(1L)
  )

  ## ---- 7. Second-type CLASS metrics (optional) --------------------------
  if (has_first_class && has_second_class) {
    neighborhood$first_type_class_fist_element_second_type_class <-
      cbind(
        neighborhood$second_type_class[, 1L, drop = TRUE],
        neighborhood$first_type_class[, -1L, drop = FALSE]
      )
    neighborhood$second_type_class_neighbors <- .get_neighbors(
      neighborhood$first_type_class_fist_element_second_type_class
    )
    neighborhood$second_type_class_neighbors_N <- vapply(
      neighborhood$second_type_class_neighbors, length, integer(1L)
    )
  }

  ## ---- 8. First-type neighbour metrics ----------------------------------
  neighborhood$first_type_neighbors <- .get_neighbors(
    neighborhood$first_type
  )
  neighborhood$first_type_neighbors_N <- vapply(
    neighborhood$first_type_neighbors, length, integer(1L)
  )

  neighborhood$first_type_singlets_neighbors <- lapply(
    seq_len(nrow(neighborhood$first_type)),
    function(i) {
      x <- neighborhood$first_type[i, ]
      y <- neighborhood$spot_class[i, ]
      which((x[-1L] == x[1L]) & (y[-1L] == "singlet")) + 1L
    }
  )
  neighborhood$first_type_singlets_neighbors_N <- vapply(
    neighborhood$first_type_singlets_neighbors, length, integer(1L)
  )

  ## ---- 9. First-type CLASS metrics (optional) ---------------------------
  if (has_first_class) {
    neighborhood$first_type_class_neighbors <- .get_neighbors(
      neighborhood$first_type_class
    )
    neighborhood$first_type_class_neighbors_N <- vapply(
      neighborhood$first_type_class_neighbors, length, integer(1L)
    )
  }

  ## ---- 10. Same second-type neighbour metrics ---------------------------
  neighborhood$same_second_type_neighbors <- .get_neighbors(
    neighborhood$second_type
  )
  neighborhood$same_second_type_neighbors_N <- vapply(
    neighborhood$same_second_type_neighbors, length, integer(1L)
  )

  return(neighborhood)
}


##  ==== 3. add_neighborhood_weight_composition ====

#' Add Neighborhood Weight Composition
#'
#' Computes the spatial composition of a neighbourhood by calculating
#' the weights of different cell types, then adds a
#' \code{neighborhood_weight_composition} matrix to the input object.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{nn_idx}, \code{first_type}, \code{second_type},
#'   \code{weight_first_type}, \code{weight_second_type}.
#'
#' @return The input \code{neighborhood} list with an additional element
#'   \code{neighborhood_weight_composition}: a matrix where each row is
#'   the cell-type weight vector for that cell's neighbourhood.
#'
#' @export
add_neighborhood_weight_composition <- function(neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("nn_idx", "first_type", "second_type",
                   "weight_first_type", "weight_second_type")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", ")
    )
  }

  ## ---- Unique cell types ------------------------------------------------
  cell_types <- sort(unique(c(
    as.vector(neighborhood$first_type),
    as.vector(neighborhood$second_type)
  )))
  cell_types <- cell_types[!is.na(cell_types)]

  ## ---- Per-neighbourhood composition ------------------------------------
  .compute_composition <- function(nbhd_idx) {
    composition      <- rep(0, length(cell_types))
    names(composition) <- cell_types

    ## Exclude unannotated cells and the focal cell (index 1)
    annotated <- nbhd_idx[!is.na(
      neighborhood$first_type[nbhd_idx, 1L, drop = TRUE]
    )]
    annotated <- annotated[-1L]

    if (length(annotated) == 0L) return(composition)

    for (i in annotated) {
      ft <- neighborhood$first_type[i,  1L, drop = TRUE]
      st <- neighborhood$second_type[i,  1L, drop = TRUE]
      w1 <- neighborhood$weight_first_type[i,  1L, drop = TRUE]
      w2 <- neighborhood$weight_second_type[i, 1L, drop = TRUE]

      if (!is.na(st)) {
        composition[ft] <- composition[ft] + w1
        composition[st] <- composition[st] + w2
      } else {
        composition[ft] <- composition[ft] + 1
      }
    }
    total <- sum(composition)
    if (total > 0) composition / total else composition
  }

  comp_mat <- t(apply(
    neighborhood$nn_idx, 1L, .compute_composition
  ))

  neighborhood$neighborhood_weight_composition <- comp_mat
  return(neighborhood)
}


##  ==== 4. add_cell_types_neighborhood_weights ====

#' Add Neighborhood Weights for Cell Types
#'
#' Calculates the weights of the first and second cell types for each
#' neighbourhood based on the neighbourhood weight composition.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{neighborhood_weight_composition}, \code{first_type},
#'   \code{second_type}.
#'
#' @return The input \code{neighborhood} list with additional elements:
#'   \code{neighborhood_weights_first_type} and
#'   \code{neighborhood_weights_second_type}.
#'
#' @export
add_cell_types_neighborhood_weights <- function(neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("neighborhood_weight_composition",
                   "first_type", "second_type")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", "), "\n",
      "Run add_neighborhood_weight_composition() first."
    )
  }

  n_cells <- nrow(neighborhood$neighborhood_weight_composition)

  weights_mat <- vapply(
    seq_len(n_cells),
    function(i) {
      ft <- neighborhood$first_type[i,  1L, drop = TRUE]
      st <- neighborhood$second_type[i, 1L, drop = TRUE]

      if (!is.na(st)) {
        c(
          neighborhood$neighborhood_weight_composition[i, ft],
          neighborhood$neighborhood_weight_composition[i, st]
        )
      } else if (!is.na(ft)) {
        c(
          neighborhood$neighborhood_weight_composition[i, ft],
          0
        )
      } else {
        c(NA_real_, NA_real_)
      }
    },
    numeric(2L)
  )

  weights_mat <- t(weights_mat)
  colnames(weights_mat) <- c("first_type", "second_type")

  neighborhood$neighborhood_weights_first_type  <-
    weights_mat[, "first_type"]
  neighborhood$neighborhood_weights_second_type <-
    weights_mat[, "second_type"]

  return(neighborhood)
}


## ==== 5. add_individual_w1_of_second_type_in_neighborhood ====

#' Add Weights for First Type of Second Type Neighbors
#'
#' Calculates the weight of the first cell type (\code{w1}) for each
#' second-type neighbour in a neighbourhood.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{second_type_neighbors},
#'   \code{weight_first_type}.
#'
#' @return The input \code{neighborhood} list with additional elements:
#'   \code{w1_second_type_in_neighborhood} (list of weight vectors) and
#'   \code{sum_w1_second_type_in_neighborhood} (numeric vector of sums).
#'
#' @export
add_individual_w1_of_second_type_in_neighborhood <- function(
    neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("second_type_neighbors", "weight_first_type")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", "), "\n",
      "Run add_infiltration_metrics_to_neighborhood() first."
    )
  }

  n_cells <- length(neighborhood$second_type_neighbors)

  neighborhood$w1_second_type_in_neighborhood <- lapply(
    seq_len(n_cells),
    function(i) {
      neighborhood$weight_first_type[
        i, neighborhood$second_type_neighbors[[i]],
        drop = TRUE
      ]
    }
  )

  neighborhood$sum_w1_second_type_in_neighborhood <- vapply(
    neighborhood$w1_second_type_in_neighborhood,
    sum,
    numeric(1L)
  )

  return(neighborhood)
}


## ==== 6. add_individual_w2_of_second_type_in_neighborhood ====

#' Add Weights of Second Cell Type's Second Annotation in Neighborhood
#'
#' Calculates the weights associated with the second type's second
#' annotation (\code{w2}) for neighbours of the same second cell type.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{same_second_type_neighbors},
#'   \code{weight_second_type}.
#'
#' @return The input \code{neighborhood} list with additional elements:
#'   \code{w2_second_type_in_neighborhood} (list of weight vectors) and
#'   \code{sum_w2_second_type_in_neighborhood} (numeric vector of sums).
#'
#' @export
add_individual_w2_of_second_type_in_neighborhood <- function(
    neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("same_second_type_neighbors", "weight_second_type")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", "), "\n",
      "Run add_infiltration_metrics_to_neighborhood() first."
    )
  }

  n_cells <- length(neighborhood$same_second_type_neighbors)

  neighborhood$w2_second_type_in_neighborhood <- lapply(
    seq_len(n_cells),
    function(i) {
      neighborhood$weight_second_type[
        i, neighborhood$same_second_type_neighbors[[i]],
        drop = TRUE
      ]
    }
  )

  neighborhood$sum_w2_second_type_in_neighborhood <- vapply(
    neighborhood$w2_second_type_in_neighborhood,
    sum,
    numeric(1L)
  )

  return(neighborhood)
}


##  ==== 7. add_max_weight_of_spilling_type_in_neighborhood ====

#' Add Maximum Weight of Spilling Cell Type in the Neighborhood
#'
#' Computes the maximum weight of a spilling cell type (the largest
#' value among \code{w1} and \code{w2} weights for neighbours of the
#' same second cell type) within a spatial neighbourhood.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{same_second_type_neighbors},
#'   \code{w1_second_type_in_neighborhood},
#'   \code{w2_second_type_in_neighborhood}.
#'
#' @return The input \code{neighborhood} list with an additional element
#'   \code{max_weight_of_spilling_type_in_neighborhood}: a numeric
#'   vector of maximum spilling-type weights per cell.
#'
#' @export
add_max_weight_of_spilling_type_in_neighborhood <- function(
    neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required <- c(
    "same_second_type_neighbors",
    "w1_second_type_in_neighborhood",
    "w2_second_type_in_neighborhood"
  )
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", "), "\n",
      "Run add_individual_w1/w2_of_second_type_in_neighborhood() ",
      "first."
    )
  }

  n_cells <- length(neighborhood$same_second_type_neighbors)

  neighborhood$max_weight_of_spilling_type_in_neighborhood <- vapply(
    seq_len(n_cells),
    function(i) {
      max(
        0,
        neighborhood$w1_second_type_in_neighborhood[[i]],
        neighborhood$w2_second_type_in_neighborhood[[i]]
      )
    },
    numeric(1L)
  )

  return(neighborhood)
}


## ==== 8. add_neigborhood_weight_composition_of_spilling_cell_type ====

#' Add Neighborhood Weight Composition of Spilling Cell Type
#'
#' Wrapper that sequentially calls the five helpers needed to compute
#' spilling-type weight composition metrics.
#'
#' @param neighborhood A list containing spatial neighbourhood data with
#'   all fields required by the underlying helper functions.
#'
#' @return The modified \code{neighborhood} list augmented with all
#'   weight composition metrics.
#'
#' @seealso
#'   \code{\link{add_neighborhood_weight_composition}},
#'   \code{\link{add_cell_types_neighborhood_weights}},
#'   \code{\link{add_individual_w1_of_second_type_in_neighborhood}},
#'   \code{\link{add_individual_w2_of_second_type_in_neighborhood}},
#'   \code{\link{add_max_weight_of_spilling_type_in_neighborhood}}
#'
#' @export
add_neigborhood_weight_composition_of_spilling_cell_type <- function(
    neighborhood) {

  neighborhood <- add_neighborhood_weight_composition(neighborhood)
  neighborhood <- add_cell_types_neighborhood_weights(neighborhood)
  neighborhood <- add_individual_w1_of_second_type_in_neighborhood(
    neighborhood
  )
  neighborhood <- add_individual_w2_of_second_type_in_neighborhood(
    neighborhood
  )
  neighborhood <- add_max_weight_of_spilling_type_in_neighborhood(
    neighborhood
  )
  return(neighborhood)
}


##  ==== 9. add_neigborhood_nCount_of_spilling_cell_type ====

#' Add Neighborhood nCount Metrics for Spilling Cell Types
#'
#' Calculates the sum of \code{nCount} values for all neighbours and
#' for spilling-type neighbours specifically.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{nn_idx}, \code{second_type_neighbors}, \code{nCount}.
#'
#' @return The input \code{neighborhood} list with additional elements:
#'   \code{sum_nCount_neighborhood} and
#'   \code{sum_nCount_neighborhood_spilling_type}.
#'
#' @export
add_neigborhood_nCount_of_spilling_cell_type <- function(neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("nn_idx", "second_type_neighbors", "nCount")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", ")
    )
  }

  n_cells <- length(neighborhood$second_type_neighbors)

  neighborhood$sum_nCount_neighborhood <- vapply(
    seq_len(n_cells),
    function(i) sum(neighborhood$nCount[i, -1L]),
    numeric(1L)
  )

  neighborhood$sum_nCount_neighborhood_spilling_type <- vapply(
    seq_len(n_cells),
    function(i) {
      idx <- neighborhood$second_type_neighbors[[i]]
      if (length(idx) == 0L) return(0)
      sum(neighborhood$nCount[i, idx])
    },
    numeric(1L)
  )

  return(neighborhood)
}


## ==== 10. add_spatial_metric ====

#' Add Spatial Metrics to Neighborhood
#'
#' End-to-end wrapper that integrates RCTD results into a spatial
#' neighbourhood and computes infiltration, weight composition, and
#' nCount metrics.
#'
#' @param spatial_neighborhood A list representing the spatial
#'   neighbourhood data.
#' @param rctd An RCTD object containing cell type decomposition
#'   results.
#'
#' @return The \code{spatial_neighborhood} list enriched with
#'   infiltration metrics, weight compositions, and nCount metrics.
#'
#' @seealso
#'   \code{\link{add_rctd_to_neighborhood}},
#'   \code{\link{add_infiltration_metrics_to_neighborhood}},
#'   \code{\link{add_neigborhood_weight_composition_of_spilling_cell_type}},
#'   \code{\link{add_neigborhood_nCount_of_spilling_cell_type}}
#'
#' @export
add_spatial_metric <- function(spatial_neighborhood, rctd) {

  if (is.null(rctd@results)) {
    stop(
      "'rctd@results' is NULL. Run spacexr::run.RCTD() first."
    )
  }
  if (!"results_df_old" %in% names(rctd@results)) {
    stop(
      "The RCTD object has not been post-processed by SPLIT.\n",
      "Please run:\n",
      "  rctd <- SPLIT::run_post_process_RCTD(rctd)"
    )
  }

  spatial_neighborhood <- add_rctd_to_neighborhood(
    graph = spatial_neighborhood,
    rctd  = rctd
  )
  spatial_neighborhood <- add_infiltration_metrics_to_neighborhood(
    neighborhood = spatial_neighborhood
  )
  spatial_neighborhood <-
    add_neigborhood_weight_composition_of_spilling_cell_type(
      neighborhood = spatial_neighborhood
    )
  if(FALSE){
    # currently does not work
    spatial_neighborhood <- add_neigborhood_nCount_of_spilling_cell_type(
      neighborhood = spatial_neighborhood
    )
  }
  return(spatial_neighborhood)
}


## ==== 11. add_annotation_from_neighbors ====

#' Annotate Neighborhood with Dominant Cell Types from Neighbors
#'
#' Determines the most frequent cell type (or class) annotation among
#' neighbours for each cell.
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{first_type}, \code{second_type}. Optional:
#'   \code{first_type_class}, \code{second_type_class}.
#'
#' @return The updated \code{neighborhood} list with added columns:
#'   \code{first_type_neighborhood},
#'   \code{first_type_neighborhood_agreement},
#'   \code{second_type_neighborhood},
#'   \code{second_type_neighborhood_agreement}, and optionally
#'   \code{first_type_class_neighborhood},
#'   \code{second_type_class_neighborhood} and their agreement fields.
#'
#' @export
add_annotation_from_neighbors <- function(neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("first_type", "second_type")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", ")
    )
  }

  ## ---- Internal helper --------------------------------------------------
  .find_most_frequent <- function(values) {
    counts <- table(values[-1L], useNA = "ifany")
    names(counts)[which.max(counts)][1L]
  }

  ## ---- First type -------------------------------------------------------
  neighborhood$first_type_neighborhood <- apply(
    neighborhood$first_type, 1L, .find_most_frequent
  )
  neighborhood$first_type_neighborhood_agreement <-
    neighborhood$first_type[, 1L, drop = TRUE] ==
    neighborhood$first_type_neighborhood

  ## ---- First type class (optional) --------------------------------------
  if ("first_type_class" %in% names(neighborhood)) {
    neighborhood$first_type_class_neighborhood <- apply(
      neighborhood$first_type_class, 1L, .find_most_frequent
    )
    neighborhood$first_type_class_neighborhood_agreement <-
      neighborhood$first_type_class[, 1L, drop = TRUE] ==
      neighborhood$first_type_class_neighborhood
  } else {
    warning(
      "'first_type_class' not found in neighborhood. ",
      "'first_type_class_neighborhood' was not computed."
    )
  }

  ## ---- Second type ------------------------------------------------------
  neighborhood$second_type_neighborhood <- apply(
    neighborhood$second_type, 1L, .find_most_frequent
  )
  neighborhood$second_type_neighborhood_agreement <-
    neighborhood$second_type[, 1L, drop = TRUE] ==
    neighborhood$second_type_neighborhood

  ## ---- Second type class (optional) -------------------------------------
  if ("second_type_class" %in% names(neighborhood)) {
    neighborhood$second_type_class_neighborhood <- apply(
      neighborhood$second_type_class, 1L, .find_most_frequent
    )
    neighborhood$second_type_class_neighborhood_agreement <-
      neighborhood$second_type_class[, 1L, drop = TRUE] ==
      neighborhood$second_type_class_neighborhood
  } else {
    warning(
      "'second_type_class' not found in neighborhood. ",
      "'second_type_class_neighborhood' was not computed."
    )
  }

  return(neighborhood)
}


## ==== 12. add_neighborhood_annotation_certainty ====

#' Annotate Neighborhood with Certainty Scores
#'
#' Calculates a certainty score for each neighbourhood based on the
#' entropy of its cell type composition. Score ranges from 0 (diverse)
#' to 1 (dominated by one cell type).
#'
#' @param neighborhood A list containing neighbourhood data. Required
#'   fields: \code{first_type}, \code{second_type}. Optional:
#'   \code{first_type_class}, \code{second_type_class}.
#'
#' @return The updated \code{neighborhood} list with certainty score
#'   fields appended.
#'
#' @details
#' Certainty is computed as:
#' \deqn{1 - H / \log(N)}
#' where \eqn{H} is the entropy of the neighbourhood cell type
#' distribution and \eqn{N} is the number of neighbours.
#'
#' @importFrom entropy entropy
#'
#' @export
add_neighborhood_annotation_certainty <- function(neighborhood) {

  ## ---- Validate required fields -----------------------------------------
  required    <- c("first_type", "second_type")
  missing_req <- required[!required %in% names(neighborhood)]
  if (length(missing_req) > 0L) {
    stop(
      "The following required fields are missing from ",
      "'neighborhood':\n",
      "  ", paste(missing_req, collapse = ", ")
    )
  }

  ## ---- Internal helper --------------------------------------------------
  .certainty <- function(x) {
    y  <- x[-1L]
    fr <- as.numeric(table(y))
    1 - entropy::entropy(fr) / log(length(y))
  }

  ## ---- First type -------------------------------------------------------
  neighborhood$first_type_neighborhood_certainty <- apply(
    neighborhood$first_type, 1L, .certainty
  )

  ## ---- First type class (optional) --------------------------------------
  if ("first_type_class" %in% names(neighborhood)) {
    neighborhood$first_type_class_neighborhood_certainty <- apply(
      neighborhood$first_type_class, 1L, .certainty
    )
  } else {
    warning(
      "'first_type_class' not found in neighborhood. ",
      "'first_type_class_neighborhood_certainty' was not computed."
    )
  }

  ## ---- Second type ------------------------------------------------------
  neighborhood$second_type_neighborhood_certainty <- apply(
    neighborhood$second_type, 1L, .certainty
  )

  ## ---- Second type class (optional) -------------------------------------
  if ("second_type_class" %in% names(neighborhood)) {
    neighborhood$second_type_class_neighborhood_certainty <- apply(
      neighborhood$second_type_class, 1L, .certainty
    )
  } else {
    warning(
      "'second_type_class' not found in neighborhood. ",
      "'second_type_class_neighborhood_certainty' was not computed."
    )
  }

  return(neighborhood)
}


## ==== 13. add_transcriptomics_metric ====

#' Add Transcriptomics-Based Metrics to Neighborhood Data
#'
#' Enriches a transcriptomics neighbourhood by integrating RCTD results,
#' propagating annotations from neighbours, and computing certainty
#' scores.
#'
#' @param transcriptomics_neighborhood A list representing the
#'   transcriptomics neighbourhood.
#' @param rctd An RCTD object containing cell type decomposition
#'   results.
#'
#' @return The updated \code{transcriptomics_neighborhood} list with
#'   RCTD-derived cell type information, propagated annotations, and
#'   certainty scores.
#'
#' @seealso
#'   \code{\link{add_rctd_to_neighborhood}},
#'   \code{\link{add_annotation_from_neighbors}},
#'   \code{\link{add_neighborhood_annotation_certainty}}
#'
#' @export
add_transcriptomics_metric <- function(transcriptomics_neighborhood,
                                       rctd) {

  if (is.null(rctd@results)) {
    stop(
      "'rctd@results' is NULL. Run spacexr::run.RCTD() first."
    )
  }
  if (!"results_df_old" %in% names(rctd@results)) {
    stop(
      "The RCTD object has not been post-processed by SPLIT.\n",
      "Please run:\n",
      "  rctd <- SPLIT::run_post_process_RCTD(rctd)"
    )
  }

  transcriptomics_neighborhood <- add_rctd_to_neighborhood(
    graph = transcriptomics_neighborhood,
    rctd  = rctd
  )
  transcriptomics_neighborhood <- add_annotation_from_neighbors(
    neighborhood = transcriptomics_neighborhood
  )
  transcriptomics_neighborhood <- add_neighborhood_annotation_certainty(
    neighborhood = transcriptomics_neighborhood
  )
  return(transcriptomics_neighborhood)
}


## ==== 14. neighborhood_analysis_to_metadata ====

#' Convert Neighborhood Analysis to Metadata
#'
#' Extracts scalar and vector elements from a neighbourhood list and
#' returns them as a \code{data.frame} indexed by \code{cell_id}.
#' Matrix and list elements are excluded; for matrices only the first
#' column is retained.
#'
#' @param neighborhood A list containing neighbourhood data including
#'   a \code{cell_id} field.
#'
#' @return A \code{data.frame} with one row per cell, indexed by
#'   \code{cell_id}.
#'
#' @export
neighborhood_analysis_to_metadata <- function(neighborhood) {

  if (!"cell_id" %in% names(neighborhood)) {
    stop("'neighborhood' must contain a 'cell_id' field.")
  }

  var_class <- vapply(
    neighborhood,
    function(x) class(x)[1L],
    character(1L)
  )

  ## Keep only scalar/vector elements (not matrix, list, igraph)
  keep <- names(var_class)[
    !var_class %in% c("matrix", "list", "igraph")
  ]

  neighborhood_df <- as.data.frame(
    neighborhood[keep],
    row.names = neighborhood$cell_id
  )

  return(neighborhood_df)
}
