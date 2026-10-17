# Small network with two modules and direction scores.
make_plot_net <- function(with_z = TRUE, ...) {
  collection <- list(
    up_term    = paste0("U", 1:6),
    up_term2   = paste0("U", 2:7),
    up_term3   = paste0("U", 3:8),
    down_term  = paste0("D", 1:6),
    down_term2 = paste0("D", 2:7),
    down_term3 = paste0("D", 3:8)
  )
  results_set <- data.frame(avg_log2FC = c(1, 1, 1, -1, -1, -1),
                            row.names = names(collection))
  if (with_z) {
    results_set$median_log2FC <- results_set$avg_log2FC
    results_set$n_up   <- c(6, 6, 6, 0, 0, 0)
    results_set$n_down <- c(0, 0, 0, 6, 6, 6)
    results_set$zscore <- c(6, 6, 6, -6, -6, -6) / sqrt(6)
  }
  suppressMessages(generate_overlap_network(
    results_set, collection, unique(unlist(collection)),
    collection_name = "A_vs_B", ...
  ))
}

test_that("plot_network returns one ggplot per type", {
  net <- make_plot_net()
  for (t in c("leiden", "regulation", "zscore", "themes")) {
    expect_s3_class(plot_network(net, t), "ggplot")
  }
  expect_match(plot_network(net)$labels$title, "A_vs_B \\(Leiden\\)")
  expect_equal(plot_network(net, "leiden", title = "Custom")$labels$title, "Custom")
})

test_that("plot_network returns a named list for several types", {
  net <- make_plot_net()
  two <- plot_network(net, c("leiden", "themes"))
  expect_named(two, c("leiden", "themes"))
  expect_named(plot_network(net, "all"), c("regulation", "leiden", "zscore", "themes"))
})

test_that("plot_network handles networks without direction scores", {
  net <- make_plot_net(with_z = FALSE)
  expect_error(plot_network(net, "zscore"), "direction scores")
  all <- suppressMessages(plot_network(net, "all"))
  expect_named(all, c("regulation", "leiden", "themes"))
})

test_that("plot_network saves each PNG with its CSV tables", {
  net <- make_plot_net()
  outdir <- file.path(tempdir(), "test_plot_network_reg")
  suppressMessages(plot_network(net, "regulation", outdir = outdir, file_stem = "x"))

  expect_setequal(list.files(outdir),
                  c("x-regulation.png", "x_network_summary.csv", "x-cluster-membership.csv"))
  unlink(outdir, recursive = TRUE)

  outdir <- file.path(tempdir(), "test_plot_network_themes")
  suppressMessages(plot_network(net, "themes", outdir = outdir))
  expect_setequal(list.files(outdir), c(
    "A_vs_B-themes.png", "A_vs_B_network_summary.csv", "A_vs_B-cluster-stats.csv",
    "A_vs_B-cluster-membership.csv", "A_vs_B-cluster-themes.csv"
  ))
  unlink(outdir, recursive = TRUE)
})

test_that("plot_network rejects objects that are not network results", {
  expect_error(plot_network(list(go_term_df = data.frame())), "generate_overlap_network")
})

test_that("print and plot methods work for scoranet_network", {
  net <- make_plot_net()
  expect_output(print(net), "scoranet_network")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot(net, type = "regulation"), "ggplot")
})
