make_mock_net_result <- function() {
  list(
    go_term_df = data.frame(
      go_term = c("immune_response", "cell_migration", "apoptotic_process",
                   "signal_transduction", "inflammatory_response"),
      cluster = c(1, 1, 2, 2, 1),
      degree = c(5, 3, 4, 2, 6),
      regulation = c("up", "down", "up", "down", "up"),
      avg_log2FC = c(1.2, -0.8, 0.9, -1.1, 1.5),
      leiden_resolution = rep(0.03, 5),
      stringsAsFactors = FALSE
    )
  )
}

test_that("write_cluster_results returns correct cluster terms", {
  net_result <- make_mock_net_result()

  result <- write_cluster_results(net_result, "immune_response", outdir = NULL)

  expect_s3_class(result, "data.frame")
  expect_true(all(result$cluster == 1))
  expect_equal(nrow(result), 3)
})

test_that("write_cluster_results returns NULL for missing GO term", {
  net_result <- make_mock_net_result()

  expect_warning(
    result <- write_cluster_results(net_result, "nonexistent_term", outdir = NULL),
    "not found in network"
  )

  expect_null(result)
})

test_that("write_cluster_results writes CSV when outdir is set", {
  net_result <- make_mock_net_result()
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_cluster_out/")

  result <- write_cluster_results(
    net_result, "apoptotic_process",
    outdir = outdir, label = "test_label"
  )

  expect_true(file.exists(paste0(outdir, "network-cluster-test_label.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_cluster_results uses default label from cluster id", {
  net_result <- make_mock_net_result()
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_cluster_default/")

  result <- write_cluster_results(
    net_result, "apoptotic_process",
    outdir = outdir
  )

  expect_true(file.exists(paste0(outdir, "network-cluster-cluster-2.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_cluster_results handles hyphen-to-underscore in go_term", {
  net_result <- list(
    go_term_df = data.frame(
      go_term = c("cell_cell_adhesion", "other_term"),
      cluster = c(1, 2),
      degree = c(3, 4),
      regulation = c("up", "down"),
      avg_log2FC = c(0.5, -0.5),
      leiden_resolution = c(0.03, 0.03),
      stringsAsFactors = FALSE
    )
  )

  result <- write_cluster_results(net_result, "cell-cell-adhesion", outdir = NULL)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$go_term, "cell_cell_adhesion")
})

test_that("write_cluster_results writes inside outdir without a trailing slash", {
  outdir <- file.path(tempdir(), "test_cluster_noslash")
  write_cluster_results(make_mock_net_result(), "immune_response",
                        outdir = outdir, label = "x")

  expect_true(file.exists(file.path(outdir, "network-cluster-x.csv")))
  expect_false(file.exists(paste0(outdir, "network-cluster-x.csv")))
  unlink(outdir, recursive = TRUE)
})
