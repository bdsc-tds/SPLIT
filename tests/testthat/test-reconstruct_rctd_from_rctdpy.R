# test-reconstruct_rctd_from_rctdpy.R
#
# Tests for reconstruct_rctd_from_rctdpy().
#
# These tests use a small self-contained synthetic dataset written to a
# temporary directory so no real rctd-py output is needed.  The SPLIT
# post-processing step (run_post_process_RCTD) is mocked so the tests remain
# fast and do not require a fitted model.

library(testthat)
library(arrow)
library(rhdf5)
library(Matrix)

# ---------------------------------------------------------------------------
# Helpers: build a minimal synthetic rctd-py output directory
# ---------------------------------------------------------------------------

.make_synthetic_rctdpy_dir <- function(
    dir,
    n_pixels    = 20L,
    n_all_cells = 30L,
    cell_types  = c("TypeA", "TypeB", "TypeC"),
    n_genes     = 10L
) {
  n_ct <- length(cell_types)

  cell_ids  <- paste0("cell_", seq_len(n_pixels))
  all_cells <- paste0("cell_", seq_len(n_all_cells))
  gene_names <- paste0("gene_", seq_len(n_genes))

  # cell_ids.parquet
  arrow::write_parquet(
    data.frame(cell_id = cell_ids),
    file.path(dir, "cell_ids.parquet")
  )

  # weights.parquet  (cells x cell_types, rows sum to 1)
  w <- matrix(1 / n_ct, nrow = n_pixels, ncol = n_ct)
  arrow::write_parquet(
    as.data.frame(w) |> setNames(cell_types),
    file.path(dir, "weights.parquet")
  )

  # weights_doublet.parquet
  arrow::write_parquet(
    data.frame(w_1 = rep(0.7, n_pixels), w_2 = rep(0.3, n_pixels)),
    file.path(dir, "weights_doublet.parquet")
  )

  # spot_results.parquet
  # Use all four spot_class codes to exercise the full mapping
  spot_classes <- rep_len(c(0L, 1L, 2L, 3L), n_pixels)
  arrow::write_parquet(
    data.frame(
      cell_id          = cell_ids,
      spot_class       = spot_classes,
      first_type       = rep(0L, n_pixels),
      second_type      = rep(1L, n_pixels),
      first_class      = rep(FALSE, n_pixels),
      second_class     = rep(FALSE, n_pixels),
      min_score        = runif(n_pixels, 100, 200),
      singlet_score    = runif(n_pixels, 100, 200),
      first_type_name  = rep(cell_types[1], n_pixels),
      second_type_name = rep(cell_types[2], n_pixels)
    ),
    file.path(dir, "spot_results.parquet")
  )

  # pixel_mask.parquet
  mask <- c(rep(TRUE, n_pixels), rep(FALSE, n_all_cells - n_pixels))
  arrow::write_parquet(
    data.frame(cell_id = all_cells, pixel_mask = mask),
    file.path(dir, "pixel_mask.parquet")
  )

  # metadata.parquet
  arrow::write_parquet(
    data.frame(cell_type_names = cell_types),
    file.path(dir, "metadata.parquet")
  )

  # reference_profiles.h5  (cell_types x genes)
  h5_path <- file.path(dir, "reference_profiles.h5")
  profiles <- matrix(runif(n_ct * n_genes), nrow = n_ct, ncol = n_genes)
  rhdf5::h5createFile(h5_path)
  rhdf5::h5write(profiles,    h5_path, "profiles")
  rhdf5::h5write(
    as.character(cell_types), h5_path, "cell_type_names"
  )
  rhdf5::h5write(
    as.character(gene_names), h5_path, "gene_names"
  )

  # de_genes.parquet (optional, not read by the function)
  arrow::write_parquet(
    data.frame(DE_genes = gene_names[seq_len(5)]),
    file.path(dir, "de_genes.parquet")
  )

  invisible(dir)
}

# Mock SPLIT::run_post_process_RCTD so tests do not require a real model
.mock_run_post_process <- function(rctd, lite, min_weight) rctd

# ---------------------------------------------------------------------------
# Fixtures shared across tests
# ---------------------------------------------------------------------------

setup_dir <- local({
  d <- tempfile("rctdpy_test_")
  dir.create(d)
  .make_synthetic_rctdpy_dir(d)
  d
})

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("function errors on non-existent save_dir", {
  expect_error(
    reconstruct_rctd_from_rctdpy("/this/path/does/not/exist"),
    regexp = NULL   # any error is sufficient
  )
})

test_that("function errors on save_dir that is not a string", {
  expect_error(reconstruct_rctd_from_rctdpy(123L))
  expect_error(reconstruct_rctd_from_rctdpy(NULL))
})

test_that("function errors on non-numeric min_weight", {
  expect_error(
    reconstruct_rctd_from_rctdpy(setup_dir, min_weight = "high")
  )
})

test_that("function errors on unmapped spot_class integer", {
  # Introduce an integer code (99) not present in spot_class_map
  bad_map <- c("0" = "reject", "1" = "singlet")   # missing 2 and 3
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    expect_error(
      reconstruct_rctd_from_rctdpy(setup_dir, spot_class_map = bad_map),
      regexp = "Unmapped spot_class integer"
    )
  )
})

test_that("function errors when supplied class_df is missing cell types", {
  incomplete_class_df <- data.frame(
    class     = "GroupA",
    row.names = "TypeA"   # TypeB and TypeC are missing
  )
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    expect_error(
      reconstruct_rctd_from_rctdpy(setup_dir, class_df = incomplete_class_df),
      regexp = "class_df is missing rows for"
    )
  )
})

test_that("returned object is an RCTD S4 instance", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result <- reconstruct_rctd_from_rctdpy(setup_dir)
      expect_s4_class(result, "RCTD")
    }
  )
})

test_that("results$results_df_xe has correct dimensions and rownames", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result    <- reconstruct_rctd_from_rctdpy(setup_dir)
      rdf       <- result@results$results_df_xe
      cell_ids  <- arrow::read_parquet(
        file.path(setup_dir, "cell_ids.parquet")
      )$cell_id

      expect_equal(nrow(rdf), length(cell_ids))
      expect_equal(rownames(rdf), cell_ids)
      expect_true(all(
        c("spot_class", "first_type", "second_type",
          "min_score", "singlet_score") %in% colnames(rdf)
      ))
    }
  )
})

test_that("spot_class factor has correct levels and no NAs", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result <- reconstruct_rctd_from_rctdpy(setup_dir)
      sc     <- result@results$results_df_xe$spot_class

      expect_equal(
        levels(sc),
        c("singlet", "doublet_certain", "doublet_uncertain", "reject")
      )
      expect_false(anyNA(sc))
    }
  )
})

test_that("weights matrix has correct dimensions and column names", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result      <- reconstruct_rctd_from_rctdpy(setup_dir)
      wmat        <- result@results$weights
      cell_types  <- c("TypeA", "TypeB", "TypeC")
      cell_ids    <- arrow::read_parquet(
        file.path(setup_dir, "cell_ids.parquet")
      )$cell_id

      expect_equal(nrow(wmat), length(cell_ids))
      expect_equal(ncol(wmat), length(cell_types))
      expect_equal(colnames(wmat), cell_types)
      expect_equal(rownames(wmat), cell_ids)
    }
  )
})

test_that("weights_doublet matrix has columns first_type and second_type", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result <- reconstruct_rctd_from_rctdpy(setup_dir)
      wd     <- result@results$weights_doublet

      expect_equal(ncol(wd), 2L)
      expect_equal(colnames(wd), c("first_type", "second_type"))
    }
  )
})

test_that("cell_type_info profiles are transposed to genes x cell_types", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result  <- reconstruct_rctd_from_rctdpy(setup_dir)
      prof    <- result@cell_type_info$info[[1]]

      expect_equal(nrow(prof), 10L)   # n_genes
      expect_equal(ncol(prof), 3L)    # n_cell_types
      expect_equal(colnames(prof), c("TypeA", "TypeB", "TypeC"))
      expect_true(all(startsWith(rownames(prof), "gene_")))
    }
  )
})

test_that("identity class_df is built when class_df = NULL", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result   <- reconstruct_rctd_from_rctdpy(setup_dir)
      cdf      <- result@internal_vars$class_df
      ct_names <- c("TypeA", "TypeB", "TypeC")

      expect_equal(rownames(cdf), ct_names)
      expect_equal(cdf$class,     ct_names)
    }
  )
})

test_that("custom class_df is passed through correctly", {
  custom_class_df <- data.frame(
    class     = c("GroupA", "GroupA", "GroupB"),
    row.names = c("TypeA",  "TypeB",  "TypeC")
  )
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result <- reconstruct_rctd_from_rctdpy(
        setup_dir, class_df = custom_class_df
      )
      expect_equal(
        result@internal_vars$class_df$class,
        c("GroupA", "GroupA", "GroupB")
      )
    }
  )
})

test_that("SpatialRNA placeholder has correct number of cells", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result      <- reconstruct_rctd_from_rctdpy(setup_dir)
      all_cell_ids <- arrow::read_parquet(
        file.path(setup_dir, "pixel_mask.parquet")
      )$cell_id

      expect_equal(
        ncol(result@spatialRNA@counts),
        length(all_cell_ids)
      )
      expect_equal(
        colnames(result@spatialRNA@counts),
        all_cell_ids
      )
    }
  )
})

test_that("config slot has RCTDmode = 'doublet'", {
  with_mocked_bindings(
    run_post_process_RCTD = .mock_run_post_process,
    .package = "SPLIT",
    {
      result <- reconstruct_rctd_from_rctdpy(setup_dir)
      expect_equal(result@config$RCTDmode, "doublet")
    }
  )
})
