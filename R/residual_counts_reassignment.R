#' Build a residual transcript reassignment operator
#'
#' Constructs a sparse reassignment operator matrix that redistributes
#' residual transcripts from cells with excess transcripts to their spatial
#' neighbors, based on uniform weights or optional per-cell weights (`nCount`).
#'
#' Each row corresponds to a "sending" cell and each column to a "receiving" cell.
#' The resulting sparse matrix can be multiplied with residual counts to
#' redistribute them across neighbors.
#'
#' @param sp_nw A spatial network object, expected to contain:
#'   - `cell_id`: character vector of all cell identifiers.
#'   - `nn_idx`: matrix of nearest-neighbor indices (rows aligned to `cell_id`).
#'   - `second_type_neighbors_no_reject`: list of integer indices of valid neighbors per cell.
#'
#' @param cells_with_residual_transcripts Character vector of cell IDs that
#'   contribute residual transcripts (senders).
#'
#' @param nCount Optional numeric vector of per-cell weights (e.g., total transcript counts).
#'   If named, names must correspond to `sp_nw$cell_id`. If provided, weights are used
#'   to proportionally distribute transcripts among neighbors. If `NULL`, uniform
#'   redistribution is used.
#'
#' @param self_keep Numeric scalar (default = `0`). Fraction of residual
#'   transcripts that each sender retains for itself (between `0` and `1`).
#'
#' @param weight_power Numeric (default = `1`). Exponent applied to `nCount` weights
#'   before normalization. Values >1 emphasize large counts, values <1 reduce weight
#'   differences.
#'
#' @param eps Small numeric (default = `0`). Added to `nCount` before weighting to
#'   avoid division by zero.
#'
#' @return A sparse column-stochastic matrix of class `dgCMatrix` with dimensions
#'   `length(sp_nw$cell_id) x length(sp_nw$cell_id)`. Entry `[i,j]` gives the
#'   fraction of transcripts reassigned from cell `i` to cell `j`.
#'
#' @details
#' - Rows correspond to sender cells, columns to receiver cells.
#' - By default (`self_keep=0` and `nCount=NULL`), residuals are uniformly
#'   distributed among neighbors.
#' - If `nCount` is provided, neighbors receive proportions proportional to
#'   their `nCount^weight_power`.
#'
#' @importFrom Matrix sparseMatrix
#' @export

build_reassigment_operator <- function(
    sp_nw,
    cells_with_residual_transcripts,
    nCount = NULL,
    self_keep = 0,
    weight_power = 1,
    eps = 0
){
  stopifnot(requireNamespace("Matrix"))
  idx2cellid <- sp_nw$cell_id
  cellid2idx <- setNames(seq_along(idx2cellid), idx2cellid)
  n <- length(idx2cellid)

  # Normalize/align nCount to match sp_nw$cell_id
  if (!is.null(nCount)) {
    if (!is.null(names(nCount))) {
      nc <- rep(NA_real_, n); names(nc) <- idx2cellid
      nc[names(nCount)] <- as.numeric(nCount)
      nCount <- nc
    } else {
      stopifnot(length(nCount) == n)
      names(nCount) <- idx2cellid
    }
    nCount <- (pmax(nCount, 0) + eps)^weight_power
  }

  i_idx <- integer(0)
  j_idx <- integer(0)
  x_vals <- numeric(0)

  for (cid in cells_with_residual_transcripts) {
    idi <- cellid2idx[[cid]]
    recv_idx <- sp_nw$nn_idx[idi, sp_nw$second_type_neighbors_no_reject[[idi]]]
    k <- length(recv_idx)
    if (k == 0L) next

    # optional self retention
    if (self_keep > 0) {
      i_idx <- c(i_idx, idi)
      j_idx <- c(j_idx, idi)
      x_vals <- c(x_vals, self_keep)
    }

    # distribute the remainder
    if (!is.null(nCount)) {
      w_raw <- nCount[recv_idx]
      w_raw[!is.finite(w_raw)] <- 0
      s <- sum(w_raw)
      if (s > 0) {
        w <- w_raw / s
      } else {
        w <- rep(1/k, k)
      }
    } else {
      w <- rep(1/k, k)
    }

    share <- (1 - self_keep)
    i_idx  <- c(i_idx, rep.int(idi, k))
    j_idx  <- c(j_idx, recv_idx)
    x_vals <- c(x_vals, w * share)
  }

  Matrix::sparseMatrix(
    i = i_idx, j = j_idx, x = x_vals,
    dims = c(n, n),
    dimnames = list(idx2cellid, idx2cellid)
  )
}


#' Reassign residual transcript counts after purification
#'
#' This function redistributes residual transcripts (raw minus corrected counts)
#' from purified cells to their neighbors in a spatial network. Redistribution
#' can be uniform across neighbors or weighted by the sending cell's total count.
#'
#' @param raw_counts A numeric matrix of raw transcript counts (genes x cells).
#' @param corrected_counts A numeric matrix of purified transcript counts
#'   (genes x cells). Must have the same dimensions and dimnames as `raw_counts`.
#' @param spatial_network An adjacency-like object that is a output of \code{\link{build_spatial_network}} and \code{\link{add_spatial_metric}}
#' @param purification_status A named vector indicating purification status
#'   for each cell (names = cell IDs, values = e.g., `"purified"` or other).
#' @param mode Character, redistribution mode. One of:
#'   * `"uniform"` (default): redistribute evenly across neighbors
#'   * `"count_proportinal"`: redistribute proportionally to `rowSums(raw_counts)`
#' @param ... Additional arguments passed to `build_reassigment_operator`.
#'
#' @return A numeric matrix of corrected counts (genes x cells), after
#'   redistributing residual transcripts.
#'
#' @details
#' The function:
#' 1. Computes residual counts (`raw - corrected`).
#' 2. Identifies cells with residual transcripts (`purification_status == "purified"`).
#' 3. Builds a reassignment operator using the spatial network.
#' 4. Redistributes residuals according to the chosen `mode`.
#'
#' @examples
#' # corrected <- reassign_residual_counts(raw_counts, purified_counts, sp_nw, purification_status)
#'
#' @export
reassign_residual_counts <- function(
    raw_counts,
    corrected_counts,
    spatial_network,
    purification_status,
    mode = c("uniform", "count_proportinal"),
    ...
){
  # ---- Argument checks ----
  stopifnot(is.matrix(raw_counts), is.matrix(corrected_counts))
  if (!all(dim(raw_counts) == dim(corrected_counts))) {
    stop("`raw_counts` and `corrected_counts` must have the same dimensions.")
  }
  if (!all(rownames(raw_counts) == rownames(corrected_counts))) {
    stop("Row names (genes) of `raw_counts` and `corrected_counts` must match.")
  }
  if (!all(colnames(raw_counts) == colnames(corrected_counts))) {
    stop("Column names (cells) of `raw_counts` and `corrected_counts` must match.")
  }
  if (is.null(names(purification_status))) {
    stop("`purification_status` must be a named vector (names = cell IDs).")
  }

  cells <- colnames(corrected_counts)
  genes <- rownames(corrected_counts)

  # ---- Residual calculation ----
  residual_counts <- raw_counts[genes, cells, drop = FALSE] -
    corrected_counts[genes, cells, drop = FALSE]

  cells_with_residual_transcripts <- names(
    purification_status[purification_status == "purified"]
  )

  # ---- Mode selection ----
  mode <- match.arg(mode)
  if (mode == "uniform") {
    nCount <- NULL
  } else if (mode == "count_proportinal") {
    nCount <- rowSums(raw_counts)
  }

  # ---- Build reassignment operator ----
  reassignment_operator <- build_reassigment_operator(
    sp_nw = spatial_network,
    cells_with_residual_transcripts = cells_with_residual_transcripts,
    nCount = nCount,
    ...
  )

  # ---- Redistribution ----
  reassigned_counts <- residual_counts %*%
    reassignment_operator[colnames(residual_counts), colnames(residual_counts)]

  corrected_counts <- corrected_counts + reassigned_counts

  corrected_counts
}
