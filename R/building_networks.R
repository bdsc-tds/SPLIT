#' Compute Neighborhood Graph with Optional Pruning
#'
#' This function computes or retrieves a k-nearest neighbors (KNN) graph from a given object,
#' using a specified dimensional reduction or features. Optionally, it prunes edges based on
#' a specified distance threshold.
#'
#' @param obj An object (e.g., Seurat object) containing data to compute the KNN graph.
#' @param neighbors_name A character string specifying the name of an existing neighbors object
#'   in `obj`. If `NULL`, neighbors are computed.
#' @param k_knn An integer specifying the number of nearest neighbors to use. Default is 20.
#' @param reduction A character string specifying the dimensional reduction to use
#'   (e.g., "pca", "spatial"). Default is "pca".
#' @param dims A numeric vector specifying the dimensions to use for the reduction. Default is `1:50`.
#' @param features A character vector of feature names to use for neighbor computation. Default is `NULL`.
#' @param graph_name A character string specifying the name of the graph to store in `obj`. Default is `"transcriptomics_knn"`.
#' @param DO_prune A logical value indicating whether to prune edges in the graph based on a distance threshold. Default is `FALSE`.
#' @param rad_pruning A numeric value specifying the maximum distance for retaining edges. Edges with distances greater than this value will be pruned. Default is `Inf`.
#' @param ... Additional parameters passed to `Seurat::FindNeighbors`.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{`graph`}{An igraph object representing the (possibly pruned) KNN graph.}
#'   \item{`nn_idx`}{A matrix of nearest-neighbor indices.}
#'   \item{`nn_dist`}{A matrix of distances to the nearest neighbors.}
#'   \item{`cell_id`}{A vector of cell IDs corresponding to the rows of `nn_idx` and `nn_dist`.}
#' }
#'
#' @details
#' If `DO_prune` is `TRUE`, the function removes edges from the graph where the
#' distance exceeds the `rad_pruning` threshold. A warning is issued if `rad_pruning` is
#' larger than the maximum distance in the graph, as no edges will be pruned in that case.
#'
#' @examples
#' \dontrun{
#' # Example usage with a Seurat object
#' result <- compute_neighborhood(obj = seurat_obj, reduction = "pca", k_knn = 15)
#' plot(result$graph)
#'
#' # Example with pruning
#' result <- compute_neighborhood(
#'   obj = seurat_obj,
#'   reduction = "pca",
#'   k_knn = 15,
#'   DO_prune = TRUE,
#'   rad_pruning = 0.2
#' )
#' plot(result$graph)
#' }
#'
#' @export
#'

compute_neighborhood <- function(
    obj,
    neighbors_name = NULL,
    k_knn = 20,
    reduction = "pca",
    dims = 1:50,
    features = NULL,
    graph_name = "transcriptomics_knn",
    DO_prune = FALSE,
    rad_pruning = Inf,
    ...
) {
  # Check if neighbors exist or compute them
  knn_neighbors <- if (!is.null(neighbors_name)) {
    if (!(neighbors_name %in% obj@neighbors)) {
      stop(paste("Neighbors", neighbors_name, "do not exist. Compute one or set `neighbors_name` to `NULL` to re-compute neighbors."))
    }
    obj@neighbors[[neighbors_name]]
  } else {
    obj <- Seurat::FindNeighbors(
      obj,
      reduction = reduction,
      dims = dims,
      features = features,
      k = k_knn,
      return.neighbor = TRUE,
      graph.name = graph_name,
      ...
    )
    obj@neighbors[[graph_name]]
  }

  # Extract adjacency list and clean up NAs
  adjacency_knn <- apply(knn_neighbors@nn.idx[, -1, drop = FALSE], 1, function(x) x[!is.na(x)])

  # Create graph from adjacency list
  graph  <- igraph::graph_from_adj_list(adjacency_knn)
  n_edge <- igraph::ecount(graph)

  max_dist <- max(knn_neighbors@nn.dist)

  # Pruning
  if(DO_prune){
    if(rad_pruning > max_dist){
      warning("No puning as `rad_pruning` is larger than any distance in the graph ")
    } else {

      # prune edges with large distance (remove distant neighbors)
      knn_neighbors@nn.idx[knn_neighbors@nn.dist > rad_pruning]  <- NA
      knn_neighbors@nn.dist[knn_neighbors@nn.dist > rad_pruning] <- NA

      adjacency_knn <- apply(knn_neighbors@nn.idx[, -1, drop = FALSE], 1, function(x) x[!is.na(x)])

      graph         <- igraph::graph_from_adj_list(adjacency_knn)
      n_edge_prune  <- igraph::ecount(graph)
      delta_n_edge  <- n_edge-n_edge_prune
      message(paste("N =", delta_n_edge, "(", round(100*delta_n_edge/n_edge), "%) edges were pruned"))
    }
  }

  # Prepare results
  result <- list(
    graph = graph,
    nn_idx = knn_neighbors@nn.idx,
    nn_dist = knn_neighbors@nn.dist,
    cell_id = knn_neighbors@cell.names
  )

  return(result)
}

#' Build Spatial Network
#'
#' A wrapper function around `compute_neighborhood` to create a spatial KNN graph
#' using spatial reduction data.
#'
#' @param obj An object (e.g., Seurat object) containing spatial data.
#' @param ... Additional parameters passed to `compute_neighborhood`.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{`graph`}{An igraph object representing the KNN graph.}
#'   \item{`nn_idx`}{A matrix of nearest-neighbor indices.}
#'   \item{`nn_dist`}{A matrix of distances to the nearest neighbors.}
#'   \item{`cell_id`}{A vector of cell IDs corresponding to the rows of `nn_idx` and `nn_dist`.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage with a Seurat object containing spatial data
#' spatial_result <- build_spatial_network(obj = seurat_obj)
#' plot(spatial_result$graph)
#' }
#'
#' @export

build_spatial_network <- function(obj, reduction = "spatial", dims = 1:2, graph_name = "spatial_knn", DO_prune = T, rad_pruning = 30, ...) {
  compute_neighborhood(obj = obj, reduction = reduction, dims = dims, DO_prune = DO_prune, rad_pruning = rad_pruning, ...)
}

#' Build Transcriptomics Network
#'
#' A wrapper function around `compute_neighborhood` to create a transcriptomics KNN graph
#' using transcriptomics reduction data.
#'
#' @param obj An object (e.g., Seurat object) containing transcriptomics data.
#' @param ... Additional parameters passed to `compute_neighborhood`.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{`graph`}{An igraph object representing the KNN graph.}
#'   \item{`nn_idx`}{A matrix of nearest-neighbor indices.}
#'   \item{`nn_dist`}{A matrix of distances to the nearest neighbors.}
#'   \item{`cell_id`}{A vector of cell IDs corresponding to the rows of `nn_idx` and `nn_dist`.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage with a Seurat object containing transcriptomics data
#' transcriptomics_result <- build_transcriptomics_network(obj = seurat_obj, k_knn = 15)
#' plot(spatial_result$graph)
#' }
#' @export
build_transcriptomics_network <- function(obj, reduction = "pca", graph_name = "transcriptomics_knn", ...){
  compute_neighborhood(obj = obj, reduction = reduction, graph_name = graph_name, ...)
}

