#' Add RCTD Results to Neighborhood Graph
#'
#' This function integrates RCTD (Robust Cell Type Decomposition) results into a neighborhood graph.
#' It ensures that the RCTD results align with the graph's cells and formats them for downstream use.
#'
#' @param graph A list representing a neighborhood graph, typically created by \code{compute_neighborhood}.
#'   Must contain \code{nn_idx} (nearest neighbor indices) and \code{cell_id} (cell identifiers).
#' @param rctd An RCTD object containing cell decomposition results. The results are expected in
#'   \code{rctd@results$results_df}.
#'
#' @return A list containing the original neighborhood graph and a matrix for each column in
#'   \code{rctd@results$results_df}, where rows correspond to cells and columns correspond to neighbors.
#'
#' @details
#' The function extracts the decomposition results from the RCTD object and maps them to the
#' neighborhood graph structure, ensuring consistency in cell alignment. It then organizes
#' the results into matrices where each matrix corresponds to a specific result column from the RCTD results.
#'
#' @examples
#' \dontrun{
#' # Assuming `graph` is a neighborhood graph and `rctd` is an RCTD object
#' updated_graph <- add_rctd_res_to_neighborhood(graph = graph, rctd = rctd)
#' }
#'
#' @export

add_rctd_to_neighborhood <- function(
    graph,
    rctd
){
  # Ensure RCTD results align with graph cells
  results_df <- rctd@results$results_df[graph$cell_id, ] # Match cells
  rownames(results_df) <- graph$cell_id

  # Reshape results to align with nearest neighbor indices
  results_df <- results_df[as.vector(graph$nn_idx), ]
  num_neighbors <- ncol(graph$nn_idx)

  # Prepare results as a list of matrices
  result_list <- lapply(colnames(results_df), function(column_name) {
    matrix(results_df[[column_name]], ncol = num_neighbors)
  })

  # Combine graph and results
  names(result_list) <- colnames(results_df)
  return(c(graph, result_list))
}


#' Add Neighborhood Infiltration Metrics to Spatial Neighborhood
#'
#' This function calculates various infiltration metrics for cell neighborhoods
#' within a spatial context. It computes metrics related to the number and type
#' of neighbors annotated to a cell's first or second type, as well as infiltration
#' metrics based on cell classes.
#'
#' @param neighborhood A list containing spatial neighborhood data. This should
#'   include matrices for `first_type`, `second_type`, `nn.idx`, and `spot_class`.
#'
#' @return A modified version of the input `neighborhood` list with additional
#'   metrics. The added fields include:
#'   \item{`total_neighbors_N`}{Total number of neighbors for each cell.}
#'   \item{`total_singlets_neighbors_N`}{Number of singlet neighbors for each cell.}
#'   \item{`second_type_neighbors`}{Indices of neighbors annotated to the second cell type.}
#'   \item{`second_type_neighbors_N`}{Count of neighbors annotated to the second cell type.}
#'   \item{`second_type_singlets_neighbors`}{Indices of singlet neighbors annotated to the second cell type.}
#'   \item{`second_type_singlets_neighbors_N`}{Count of singlet neighbors annotated to the second cell type.}
#'   \item{`first_type_neighbors`}{Indices of neighbors annotated to the first cell type.}
#'   \item{`first_type_neighbors_N`}{Count of neighbors annotated to the first cell type.}
#'   \item{`first_type_singlets_neighbors`}{Indices of singlet neighbors annotated to the first cell type.}
#'   \item{`first_type_singlets_neighbors_N`}{Count of singlet neighbors annotated to the first cell type.}
#'   \item{`second_type_class_neighbors`}{Indices of neighbors annotated to the second cell type class.}
#'   \item{`second_type_class_neighbors_N`}{Count of neighbors annotated to the second cell type class.}
#'   \item{`first_type_class_neighbors`}{Indices of neighbors annotated to the first cell type class.}
#'   \item{`first_type_class_neighbors_N`}{Count of neighbors annotated to the first cell type class.}
#'   \item{`same_second_type_neighbors`}{Indices of neighbors annotated to the same second cell type.}
#'   \item{`same_second_type_neighbors_N`}{Count of neighbors annotated to the same second cell type.}
#'
#'
#' @export
#'
add_infiltration_metrics_to_neighborhood <- function(neighborhood) {


  # Compute the total number of neighbors
  neighborhood$total_neighbors_N <- rowSums(!is.na(neighborhood$nn_idx)) - 1
  neighborhood$annotated_neighbors_N <- rowSums(!is.na(neighborhood$first_type)) - 1

  # Compute total number of singlet neighbors
  neighborhood$total_singlets_neighbors_N <- rowSums(neighborhood$spot_class == "singlet", na.rm = TRUE)

  ## Generate vector of annotation of cell's neighbors ignoring reject cells
  neighborhood$first_type_no_reject <- neighborhood$first_type
  neighborhood$first_type_no_reject[neighborhood$spot_class == "reject"] <- "NA"

  # Create first_type_fist_element_second_type (with first element as second_type)
  # to facilitate identification of the number of the second cell type in the neighborhood
  neighborhood$first_type_fist_element_second_type <- cbind(neighborhood$second_type[, 1], neighborhood$first_type[, -1])
  neighborhood$first_type_fist_element_second_type_no_reject <- cbind(neighborhood$second_type[, 1], neighborhood$first_type_no_reject[, -1])


  # Helper function to calculate neighbors for a given type
  get_neighbors <- function(neighborhood_col) {
    apply(neighborhood_col,
          1,
          FUN = function(x) {
            which(x[-1] == x[1]) + 1
          })
  }

  # SECOND type neighbors: indices where neighbors are annotated to the residual second type
  neighborhood$second_type_neighbors <- get_neighbors(neighborhood$first_type_fist_element_second_type)
  neighborhood$second_type_neighbors_N <- sapply(neighborhood$second_type_neighbors, length)

  neighborhood$second_type_neighbors_no_reject <- get_neighbors(neighborhood$first_type_fist_element_second_type_no_reject)
  neighborhood$second_type_neighbors_no_reject_N <- sapply(neighborhood$second_type_neighbors_no_reject, length)

  # SECOND type singlet neighbors: neighbors annotated to the residual cell type and are singlets
  neighborhood$second_type_singlets_neighbors <- sapply(
    seq(nrow(neighborhood$first_type_fist_element_second_type)),
    FUN = function(i){
      x = neighborhood$first_type_fist_element_second_type[i,]
      y = neighborhood$spot_class[i,]
      which((x[-1] == x[1]) & (y[-1] == "singlet")) + 1
    })

  neighborhood$second_type_singlets_neighbors_N <- sapply(neighborhood$second_type_singlets_neighbors, length)

  # Compute infiltration of second type CLASS in the neighborhood
  neighborhood$first_type_class_fist_element_second_type_class <- cbind(neighborhood$second_type_class[, 1], neighborhood$first_type_class[, -1])
  neighborhood$second_type_class_neighbors <- get_neighbors(neighborhood$first_type_class_fist_element_second_type_class)
  neighborhood$second_type_class_neighbors_N <- sapply(neighborhood$second_type_class_neighbors, length)

  # FIRST type neighbors: indices where neighbors are annotated to the main first type
  neighborhood$first_type_neighbors <- get_neighbors(neighborhood$first_type)
  neighborhood$first_type_neighbors_N <- sapply(neighborhood$first_type_neighbors, length)

  # FIRST type singlet neighbors: neighbors annotated to the main cell type and are singlets
  neighborhood$first_type_singlets_neighbors <- sapply(
    seq(nrow(neighborhood$first_type)),
    FUN = function(i){
      x = neighborhood$first_type[i,]
      y = neighborhood$spot_class[i,]
      which((x[-1] == x[1]) & (y[-1] == "singlet")) + 1
    })
  neighborhood$first_type_singlets_neighbors_N <- sapply(neighborhood$first_type_singlets_neighbors, length)

  # Compute infiltration of first type CLASS in the neighborhood
  neighborhood$first_type_class_neighbors <- get_neighbors(neighborhood$first_type_class)
  neighborhood$first_type_class_neighbors_N <- sapply(neighborhood$first_type_class_neighbors, length)

  # SAME SECOND type neighbors: indices where neighbors are annotated to the same second type
  neighborhood$same_second_type_neighbors   <- get_neighbors(neighborhood$second_type)
  neighborhood$same_second_type_neighbors_N <- sapply(neighborhood$same_second_type_neighbors, length)

  return(neighborhood)
}


#' Add Neighborhood Weight Composition
#'
#' This function computes the spatial composition of a neighborhood by calculating the weights of different cell types.
#' It adds the computed composition as a new column (`neighborhood_weight_composition`) to the input `neighborhood` object.
#'
#' @param neighborhood A list containing information about the neighborhood, including `nn_idx` (indices of neighbors),
#' `first_type` and `second_type` (cell types of the neighbors), and `weight_first_type` and `weight_second_type` (weights of the corresponding cell types).
#'
#' @return A modified `neighborhood` object with an added column `neighborhood_weight_composition`, which is a matrix where each row corresponds
#' to the vector of cell-type weights for a neighborhood.
#'
#' @importFrom dplyr %>%
#' @export
#'
add_neighborhood_weight_composition <- function(
    neighborhood
){

  # Neighborhood composition
  compute_neighborhood_weight_composition <- function(
    nbhd
  ){
    composition <- rep(0, length(cell_types))
    names(composition) <- cell_types
    nbhd <- nbhd[!is.na(neighborhood$first_type[nbhd, 1])] # exclude cells not annotated by RCTD
    nbhd <- nbhd[-1]

    if(length(nbhd) == 0)
      return(composition)

    for(i in nbhd){
      ft <- neighborhood$first_type[i, 1]
      st <- neighborhood$second_type[i, 1]
      w1 <- neighborhood$weight_first_type[i, 1]
      w2 <- neighborhood$weight_second_type[i, 1]
      if(!is.na(st)){
        composition[ft] <- composition[ft] + w1
        composition[st] <- composition[st] + w2
      } else {
        composition[ft] <- composition[ft] + 1
      }
    }
    #composition <- unname(composition)
    return(composition/sum(composition))
  }

  # Extract unique cell types
  cell_types <- neighborhood$first_type %>% as.vector() %>% unique() %>% sort()
  # Store `neighborhood_weight_composition` -- the global weight composition of the neighborhood
  neighborhood$neighborhood_weight_composition <- apply(neighborhood$nn_idx, 1, compute_neighborhood_weight_composition, simplify = T) %>% t()

  return(neighborhood)
}


#' Add Neighborhood Weights for Cell Types
#'
#' This function calculates the weights of the first and second cell types for each neighborhood based on the neighborhood composition.
#' It assigns a weight of 0 to the second cell type if it is `NA`.
#'
#' @param neighborhood A list containing the neighborhood data, including:
#'   - `neighborhood_weight_composition`: A matrix where each row corresponds to a neighborhood, and columns represent cell types with their respective weights.
#'   - `first_type`: A matrix or vector indicating the first cell type for each neighborhood.
#'   - `second_type`: A matrix or vector indicating the second cell type for each neighborhood (can contain `NA` values).
#'
#' @return The input `neighborhood` object with an additional element `neighborhood_weights`,
#'   a matrix of dimensions `nrow(neighborhood$neighborhood_weight_composition) x 2`.
#'   Each row contains the weights of the `first_type` and `second_type` for the corresponding neighborhood.
#'
#' @export
add_cell_types_neighborhood_weights <- function(
    neighborhood
){

  neighborhood_weights <- sapply(
    1:nrow(neighborhood$neighborhood_weight_composition),
    FUN = function(i){
      ft <- neighborhood$first_type[i,1]
      st <- neighborhood$second_type[i,1]

      if(!is.na(st)){
        neighborhood_weights <- c(neighborhood$neighborhood_weight_composition[i,ft], neighborhood$neighborhood_weight_composition[i, st])
      } else{
        if(!is.na(ft)){
          neighborhood_weights <- c(neighborhood$neighborhood_weight_composition[i,ft], 0)
        } else {
          neighborhood_weights <- c(NA, NA)
        }
      }
    }
  )
  neighborhood_weights <- neighborhood_weights %>% t()
  colnames(neighborhood_weights) <- c("first_type", "second_type")

  neighborhood$neighborhood_weights_first_type  <- neighborhood_weights[, "first_type"]
  neighborhood$neighborhood_weights_second_type <- neighborhood_weights[, "second_type"]
  return(neighborhood)
}


#' Add Weights for First Type of Second Type Neighbors
#'
#' This function calculates the weight of the first cell type (`w1`) for each second-type neighbor in a neighborhood.
#' It adds this information to the `neighborhood` object as a new element.
#'
#' @param neighborhood A list containing the neighborhood data, including:
#'   - `second_type_neighbors`: A list where each element contains indices of neighbors corresponding to the second type.
#'   - `weight_first_type`: A matrix where rows correspond to neighborhoods, and columns represent weights of the first type for each neighbor.
#'
#' @return The input `neighborhood` object with an additional element,
#'   - `w1_second_type_in_neighborhood`, which is a vector of weights for the first cell type
#'   corresponding to each second-type neighbor.
#'   - `sum_w1_second_type_in_neighborhood` : A vector with sum of `w1_second_type_in_neighborhood`
#'
#' @export

add_individual_w1_of_second_type_in_neighborhood <- function(
    neighborhood
){
  neighborhood$w1_second_type_in_neighborhood <- sapply(1:length(neighborhood$second_type_neighbors), function(i){
    neighborhood$weight_first_type[i, neighborhood$second_type_neighbors[[i]]]
  })

  neighborhood$sum_w1_second_type_in_neighborhood <- lapply(neighborhood$w1_second_type_in_neighborhood, sum) %>% unlist()

  return(neighborhood)
}

# Functions to retrieve from `xenium_analysis_pipeline` project

# add scores: max w1 as second type
#' Add Weights of Second Cell Type's Second Annotation in the Neighborhood
#'
#' This function calculates the weights associated with the second type's second annotation
#' (w2) for neighbors of the same second cell type within a spatial neighborhood.
#'
#' @param neighborhood A list containing spatial neighborhood data. This list must include:
#'   - `same_second_type_neighbors`: Indices of neighbors annotated to the same second cell type.
#'   - `weight_second_type`: A matrix of weights for the second cell type annotations.
#'
#' @return The modified `neighborhood` list with an added field:
#'   - `w2_second_type_in_neighborhood`: A list where each entry contains the `w2` weights
#'     of neighbors annotated to the same second cell type for each cell.
#'   - `sum_w2_second_type_in_neighborhood` : A vector with sum of `w2_second_type_in_neighborhood`
#'
#' @export
add_individual_w2_of_second_type_in_neighborhood <- function(
    neighborhood
){
  neighborhood$w2_second_type_in_neighborhood <- sapply(1:length(neighborhood$same_second_type_neighbors), function(i){
    neighborhood$weight_second_type[i, neighborhood$same_second_type_neighbors[[i]]]
  })

  neighborhood$sum_w2_second_type_in_neighborhood <- lapply(neighborhood$w2_second_type_in_neighborhood, sum) %>% unlist()

  return(neighborhood)
}

#' Add Maximum Weight of Spilling Cell Type in the Neighborhood
#'
#' This function computes the maximum weight of a spilling cell type
#' (i.e., the largest value among `w1` and `w2` weights for neighbors of the same second cell type)
#' within a spatial neighborhood.
#'
#' @param neighborhood A list containing spatial neighborhood data. This list must include:
#'   - `same_second_type_neighbors`: Indices of neighbors annotated to the same second cell type.
#'   - `w1_second_type_in_neighborhood`: A list of `w1` weights for the same second cell type neighbors.
#'   - `w2_second_type_in_neighborhood`: A list of `w2` weights for the same second cell type neighbors.
#'
#' @return The modified `neighborhood` list with an added field:
#'   - `max_weight_of_spilling_type_in_neighborhood`: A numeric vector containing the maximum
#'     weight of the spilling cell type for each cell in the neighborhood.
#'
#' @export
add_max_weight_of_spilling_type_in_neighborhood <- function(
    neighborhood
){
  N_cells <- length(neighborhood$same_second_type_neighbors)
  neighborhood$max_weight_of_spilling_type_in_neighborhood <- sapply(1:N_cells, function(i){
    max(0, neighborhood$w1_second_type_in_neighborhood[[i]], neighborhood$w2_second_type_in_neighborhood[[i]])
  })
  return(neighborhood)
}


#' Add Neighborhood Weight Composition of Spilling Cell Type
#'
#' This function computes and adds several metrics related to the weight composition of spilling
#' cell types within a spatial neighborhood. It integrates multiple helper functions to calculate
#' neighborhood composition, cell type weights, and the maximum weight of spilling cell types.
#'
#' @param neighborhood A list containing spatial neighborhood data. The `neighborhood` list
#'   should include the required fields for the following computations:
#'   - `add_neighborhood_weight_composition`: Calculates neighborhood cell-type weight composition.
#'   - `add_cell_types_neighborhood_weights`: Extracts weights for the first and second cell types.
#'   - `add_individual_w1_of_second_type_in_neighborhood`: Extracts individual `w1` weights of
#'     the second type neighbors.
#'   - `add_individual_w2_of_second_type_in_neighborhood`: Extracts individual `w2` weights of
#'     the second type neighbors.
#'   - `add_max_weight_of_spilling_type_in_neighborhood`: Computes the maximum weight of spilling
#'     cell types in the neighborhood.
#'
#' @return The modified `neighborhood` list, augmented with additional metrics:
#'   - `neighborhood_weight_composition`: Cell-type weight composition in the neighborhood.
#'   - `neighborhood_weights`: Weights of the first and second cell types.
#'   - `w1_second_type_in_neighborhood`: List of `w1` weights for the second cell type neighbors.
#'   - `w2_second_type_in_neighborhood`: List of `w2` weights for the second cell type neighbors.
#'   - `max_weight_of_spilling_type_in_neighborhood`: Maximum weight of spilling cell types in
#'     the neighborhood.
#'
#' @seealso
#' - \code{\link{add_neighborhood_weight_composition}}
#' - \code{\link{add_cell_types_neighborhood_weights}}
#' - \code{\link{add_individual_w1_of_second_type_in_neighborhood}}
#' - \code{\link{add_individual_w2_of_second_type_in_neighborhood}}
#' - \code{\link{add_max_weight_of_spilling_type_in_neighborhood}}
#'
#' @export

add_neigborhood_weight_composition_of_spilling_cell_type <- function(
    neighborhood
){
  neighborhood <- add_neighborhood_weight_composition(neighborhood = neighborhood)
  neighborhood <- add_cell_types_neighborhood_weights(neighborhood = neighborhood)
  neighborhood <- add_individual_w1_of_second_type_in_neighborhood(neighborhood = neighborhood)
  neighborhood <- add_individual_w2_of_second_type_in_neighborhood(neighborhood = neighborhood)
  neighborhood <- add_max_weight_of_spilling_type_in_neighborhood(neighborhood = neighborhood)
  return(neighborhood)
}

#' Add Neighborhood nCount Metrics for Spilling Cell Types
#'
#' This function calculates the sum of `nCount` values for cells in a spatial neighborhood,
#' as well as the sum of `nCount` values specifically for cells of the spilling cell type
#' within each neighborhood.
#'
#' @param neighborhood A list containing spatial neighborhood data. The list should include:
#'   \itemize{
#'     \item{\code{nn_idx} A matrix of indices representing neighborhood relationships. The first column represents the focal cell, and subsequent columns represent its neighbors.}
#'     \item{\code{second_type_neighbors} A list where each element contains indices of neighbors annotated to the second cell type for the corresponding focal cell.}
#'     \item{\code{nCount} A list of numeric vectors where each vector contains `nCount` values for cells in the dataset.}
#'   }
#'
#' @return A modified version of the input \code{neighborhood} list with the following additional fields:
#'   \itemize{
#'     \item{\code{sum_nCount_neighborhood} A numeric vector containing the sum of `nCount` values for all neighbors of each focal cell.}
#'     \item{\code{sum_nCount_neighborhood_spilling_type} A numeric vector containing the sum of `nCount` values for neighbors of each focal cell that belong to the spilling cell type.}
#'   }
#'
#' @export

add_neigborhood_nCount_of_spilling_cell_type <- function(
    neighborhood
){
  N_cells <- length(neighborhood$second_type_neighbors)
  neighborhood$sum_nCount_neighborhood <- sapply(1:N_cells, function(i){
    sum(neighborhood$nCount[i, -1])
  })

  neighborhood$sum_nCount_neighborhood_spilling_type <- sapply(1:N_cells, function(i){
    sum(neighborhood$nCount[i, neighborhood$second_type_neighbors[[i]]])
  })

  return(neighborhood)
}


#' Add Spatial Metrics to Neighborhood
#'
#' This function enhances a spatial neighborhood structure by integrating RCTD results and adding various spatial metrics.
#' It incorporates infiltration metrics, weight compositions, and neighborhood count metrics, facilitating in-depth spatial analysis.
#'
#' @param spatial_neighborhood A list representing the spatial neighborhood data. It should include fields required for neighborhood and infiltration analysis.
#' @param rctd An RCTD object containing cell type decomposition results.
#'
#' @return A modified spatial neighborhood list enriched with additional metrics, including:
#' \itemize{
#'   \item Infiltration metrics.
#'   \item Weight composition of spilling cell types.
#'   \item Sum of \code{nCount} values for neighborhood and spilling cell types.
#' }
#'
#' @details This function sequentially:
#' \enumerate{
#'   \item Integrates RCTD results into the spatial neighborhood.
#'   \item Calculates infiltration metrics based on cell type and class.
#'   \item Computes the weight composition for spilling cell types.
#'   \item Adds neighborhood \code{nCount} metrics for spilling cell types.
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{add_rctd_to_neighborhood}}: Integrates RCTD results into the neighborhood.
#'   \item \code{\link{add_infiltration_metrics_to_neighborhood}}: Computes infiltration metrics.
#'   \item \code{\link{add_neigborhood_weight_composition_of_spilling_cell_type}}: Adds weight composition for spilling cell types.
#'   \item \code{\link{add_neigborhood_nCount_of_spilling_cell_type}}: Computes neighborhood \code{nCount} metrics.
#' }
#'
#' @export
#'


add_spatial_metric <- function(spatial_neighborhood, rctd){

  spatial_neighborhood <- add_rctd_to_neighborhood(graph = spatial_neighborhood, rctd = rctd)
  spatial_neighborhood <- add_infiltration_metrics_to_neighborhood(neighborhood = spatial_neighborhood)
  spatial_neighborhood <- add_neigborhood_weight_composition_of_spilling_cell_type(neighborhood = spatial_neighborhood)
  spatial_neighborhood <- add_neigborhood_nCount_of_spilling_cell_type(neighborhood = spatial_neighborhood)
  return(spatial_neighborhood)
}


#' Annotate Transcriptomics Neighborhood with Dominant Cell Types
#'
#' This function annotates each neighborhood in a transcriptomics dataset by determining
#' the most frequent cell type (or class) for each neighborhood category (`first_type` and `second_type`).
#'
#' @param neighborhood A data frame or list containing neighborhood-level transcriptomics data.
#'        It must contain `first_type`, `second_type`, and optionally `first_type_class` and `second_type_class`.
#'
#' @return The updated `neighborhood` object with the following added columns:
#'   \itemize{
#'     \item `first_type_neighborhood`: The most frequent first type in the neighborhood.
#'     \item `first_type_class_neighborhood`: The most frequent first type class in the neighborhood (if available).
#'     \item `second_type_neighborhood`: The most frequent second type in the neighborhood.
#'     \item `second_type_class_neighborhood`: The most frequent second type class in the neighborhood (if available).
#'   }
#'
#' @details The function applies a helper function `find_most_frequent()` that determines
#' the most frequently occurring value in a given set of neighborhood annotations. If the `first_type_class`
#' or `second_type_class` column is missing, a warning is issued, and the respective computation is skipped.
#'
#' @export
#'
add_annotation_from_neighbors <- function(neighborhood) {

  # Helper function to determine the most frequent annotation in a neighborhood
  find_most_frequent <- function(values) {
    value_counts <- values[-1]
    most_frequent <- table(value_counts, useNA = "ifany") %>% which.max %>% names %>% .[1]
    return(most_frequent)
  }

  # Annotate first type neighborhood
  neighborhood$first_type_neighborhood <- apply(
    neighborhood$first_type,
    1,
    find_most_frequent
  )

  neighborhood$first_type_neighborhood_agreement <- neighborhood$first_type[,1] == neighborhood$first_type_neighborhood

  # Annotate first type class neighborhood if available
  if ("first_type_class" %in% names(neighborhood)) {
    neighborhood$first_type_class_neighborhood <- apply(
      neighborhood$first_type_class,
      1,
      find_most_frequent
    )

    neighborhood$first_type_class_neighborhood_agreement <- neighborhood$first_type_class[,1] == neighborhood$first_type_class_neighborhood

  } else {
    warning("`first_type_class` does not exist in the neighborhood object, `first_type_class_neighborhood` was not computed")
  }

  # Annotate second type neighborhood
  neighborhood$second_type_neighborhood <- apply(
    neighborhood$second_type,
    1,
    find_most_frequent
  )

  neighborhood$second_type_neighborhood_agreement <- neighborhood$second_type[,1] == neighborhood$second_type_neighborhood


  # Annotate second type class neighborhood if available
  if ("second_type_class" %in% names(neighborhood)) {
    neighborhood$second_type_class_neighborhood <- apply(
      neighborhood$second_type_class,
      1,
      find_most_frequent
    )

    neighborhood$second_type_class_neighborhood_agreement <- neighborhood$second_type_class[,1] == neighborhood$second_type_class_neighborhood

  } else {
    warning("`second_type_class` does not exist in the neighborhood object, `second_type_class_neighborhood` was not computed")
  }

  return(neighborhood)
}


#' Annotate Neighborhood with Certainty Scores
#'
#' This function calculates a certainty score for each neighborhood based on
#' the entropy of its cell type composition. The certainty score ranges from 0
#' (high entropy, diverse cell types) to 1 (low entropy, dominated by a single cell type).
#'
#' @param neighborhood
#'   A data frame or list containing neighborhood-level cell type data. It must
#'   contain `first_type` and `second_type`, and optionally
#'   `first_type_class` and `second_type_class`.
#'
#' @return
#'   The updated `neighborhood` object with the following additional columns:
#'   \itemize{
#'     \item `first_type_neighborhood_certainty`: Certainty score for the first type.
#'     \item `first_type_class_neighborhood_certainty`: Certainty score for the first type class (if available).
#'     \item `second_type_neighborhood_certainty`: Certainty score for the second type.
#'     \item `second_type_class_neighborhood_certainty`: Certainty score for the second type class (if available).
#'   }
#'
#' @details
#'   Certainty is computed as:
#'   \deqn{1 - \frac{H}{\log(N)}}{
#'   1 - entropy / log(N)
#'   }
#'   where \eqn{H} is the entropy of the neighborhood's cell type distribution,
#'   and \eqn{N} is the total number of unique cell types. A higher score
#'   indicates a more homogenous neighborhood.
#'
#'   If `first_type_class` or `second_type_class` is missing, a warning is issued
#'   and the respective certainty score is not computed.
#'
#' @importFrom entropy entropy
#' @export
#'

add_neighborhood_annotation_certainty <- function(
    neighborhood
) {

  # Helper function to compute normalized certainty score based on entropy
  neighborhood_normalized_certainty <- function(x) {
    y <- x[-1]
    fr <- table(y) %>% as.numeric()
    res <- 1 - entropy::entropy(fr) / log(length(y))
    # More accurate denominator: log(min(length(y), N_cell_types))
    return(res)
  }

  # Compute certainty scores for first type neighborhoods
  neighborhood$first_type_neighborhood_certainty <- apply(
    neighborhood$first_type,
    1,
    neighborhood_normalized_certainty
  )

  # Compute certainty scores for first type class neighborhoods if available
  if ("first_type_class" %in% names(neighborhood)) {
    neighborhood$first_type_class_neighborhood_certainty <- apply(
      neighborhood$first_type_class,
      1,
      neighborhood_normalized_certainty
    )
  } else {
    warning("`first_type_class` does not exist in the neighborhood object, `first_type_class_neighborhood_certainty` was not computed")
  }

  # Compute certainty scores for second type neighborhoods
  neighborhood$second_type_neighborhood_certainty <- apply(
    neighborhood$second_type,
    1,
    FUN = neighborhood_normalized_certainty
  )

  # Compute certainty scores for second type class neighborhoods if available
  if ("second_type_class" %in% names(neighborhood)) {
    neighborhood$second_type_class_neighborhood_certainty <- apply(
      neighborhood$second_type_class,
      1,
      FUN = neighborhood_normalized_certainty
    )
  } else {
    warning("`second_type_class` does not exist in the neighborhood object, `second_type_class_neighborhood_certainty` was not computed")
  }

  return(neighborhood)
}



#' Add Transcriptomics-Based Metrics to Neighborhood Data
#'
#' This function enriches a transcriptomics neighborhood dataset by incorporating
#' additional metrics from an `rctd` object, propagating annotations from neighboring cells,
#' and computing neighborhood annotation certainty scores.
#'
#' @param transcriptomics_neighborhood
#'   A data frame or list representing the transcriptomics neighborhood,
#'   containing spatial relationships between cells.
#' @param rctd
#'   An object containing **RCTD** (Robust Cell Type Decomposition) results,
#'   which provide cell type compositions inferred from spatial transcriptomics data.
#'
#' @return
#'   The updated `transcriptomics_neighborhood` object with additional transcriptomics-based metrics:
#'   \itemize{
#'     \item **RCTD-derived cell type information** integrated into the neighborhood.
#'     \item **Propagated annotations from neighboring cells.**
#'     \item **Certainty scores** quantifying the confidence of neighborhood-level annotations.
#'   }
#'
#' @details
#'   The function performs the following steps:
#'   \enumerate{
#'     \item **Integrates RCTD data** into the transcriptomics neighborhood using `add_rctd_to_neighborhood()`.
#'     \item **Propagates annotations from neighboring cells** using `add_annotation_from_neighbors()`.
#'     \item **Computes neighborhood annotation certainty** using `add_neighborhood_annotation_certainty()`.
#'   }
#'
#' @seealso
#'   \code{\link{add_rctd_to_neighborhood}},
#'   \code{\link{add_annotation_from_neighbors}},
#'   \code{\link{add_neighborhood_annotation_certainty}}
#'
#' @export
#'

add_transcriptomics_metric <- function(transcriptomics_neighborhood, rctd) {

  # Add RCTD data to the transcriptomics neighborhood
  transcriptomics_neighborhood <- add_rctd_to_neighborhood(
    graph = transcriptomics_neighborhood,
    rctd = rctd
  )

  # Propagate annotations from neighboring cells
  transcriptomics_neighborhood <- add_annotation_from_neighbors(
    neighborhood = transcriptomics_neighborhood
  )

  # Compute neighborhood annotation certainty
  transcriptomics_neighborhood <- add_neighborhood_annotation_certainty(
    neighborhood = transcriptomics_neighborhood
  )

  return(transcriptomics_neighborhood)
}



#' Convert Neighborhood Analysis to Metadata
#'
#' This function extracts specific elements from a neighborhood list and
#' converts them into a data frame format suitable for further analysis or
#' integration with metadata. It processes variables that are either vectors
#' or matrices, retaining relevant columns and applying necessary transformations.
#'
#' @param neighborhood A list containing neighborhood data, with various
#'   variables, including vectors and matrices.
#'
#' @return A data frame with the following columns:
#'   - Variables from the input `neighborhood` list that are either vectors
#'     or matrices, with matrices reduced to their first column.
#'   - The data frame is indexed by the `cell_id` column from the input list.
#'
#' @details
#'   - The function checks for the class of each element in the `neighborhood`
#'     list and processes vectors and matrices.
#'   - For matrices, only the first column is kept.
#'   - A data frame is created from the selected variables, using `cell_id`
#'     as the row names.
#'
#' @export
#'

neighborhood_analysis_to_metadata <- function(
    neighborhood
){
  var_names <- names(neighborhood)
  var_class <- lapply(neighborhood, function(x){class(x)[1]}) %>% unlist()
  is_vct    <- var_class[which(!var_class %in% c("matrix", "list", "igraph"))] %>% names()

  variables_to_keep <- c(is_vct)
  neighborhood_df <- neighborhood[variables_to_keep]
  neighborhood_df <- as.data.frame(neighborhood_df, row.names = neighborhood_df$cell_id)

  return(neighborhood_df)

}
