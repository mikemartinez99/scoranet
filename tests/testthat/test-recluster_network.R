# Network of 8 loosely linked modules, so resolution changes the cluster count.
make_recluster_net <- function() {
  set.seed(3)
  collection <- unlist(lapply(1:8, function(m) {
    pool <- paste0("M", m, "_", 1:8)
    nxt <- paste0("M", m %% 8 + 1, "_", 1:2)
    sets <- lapply(1:4, function(i) c(sample(pool, 5), nxt))
    stats::setNames(sets, paste0("m", m, "_term_", 1:4))
  }), recursive = FALSE)
  results_set <- data.frame(avg_log2FC = rep(c(1, -1), 16), row.names = names(collection))
  suppressMessages(generate_overlap_network(
    results_set, collection, unique(unlist(collection)), collection_name = "rc"
  ))
}

test_that("recluster_network changes the clustering and keeps the layout", {
  net <- make_recluster_net()
  low  <- suppressMessages(recluster_network(net, resolution = 0.2))
  high <- suppressMessages(recluster_network(net, resolution = 3))

  expect_s3_class(high, "scoranet_network")
  expect_equal(high$leiden_resolution, 3)
  expect_true(all(high$go_term_df$leiden_resolution == 3))
  expect_gt(nrow(high$cluster_stats), nrow(low$cluster_stats))
  expect_equal(nrow(high$cluster_themes), nrow(high$cluster_stats))
  expect_equal(igraph::V(high$network)$cluster, high$go_term_df$cluster)
  expect_identical(high$layout, net$layout)
  expect_identical(high$overlap_matrix, net$overlap_matrix)
})

test_that("recluster_network at the same settings reproduces the original", {
  net <- make_recluster_net()
  again <- suppressMessages(recluster_network(net, resolution = net$leiden_resolution))
  expect_identical(again$go_term_df$cluster, net$go_term_df$cluster)
})

test_that("recluster_network leaves the session seed alone", {
  net <- make_recluster_net()
  set.seed(42); before <- stats::runif(1)
  set.seed(42); invisible(suppressMessages(recluster_network(net, 0.5))); after <- stats::runif(1)
  expect_identical(before, after)
})

test_that("recluster_network writes CSV tables when outdir is set", {
  net <- make_recluster_net()
  outdir <- file.path(tempdir(), "test_recluster")
  suppressMessages(recluster_network(net, 0.5, outdir = outdir))
  expect_setequal(list.files(outdir), c(
    "rc_network_summary.csv", "rc-cluster-stats.csv",
    "rc-cluster-membership.csv", "rc-cluster-themes.csv"
  ))
  unlink(outdir, recursive = TRUE)
})
