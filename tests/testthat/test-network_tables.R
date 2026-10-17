make_table_net <- function(outdir = NULL) {
  collection <- list(
    a1 = paste0("A", 1:6), a2 = paste0("A", 2:7), a3 = paste0("A", 3:8),
    b1 = paste0("B", 1:6), b2 = paste0("B", 2:7)
  )
  results_set <- data.frame(avg_log2FC = c(1, 1, 1, -1, -1), row.names = names(collection))
  suppressMessages(generate_overlap_network(
    results_set, collection, unique(unlist(collection)),
    collection_name = "tables", outdir = outdir
  ))
}

test_that("network_tables returns each table with the expected columns", {
  net <- make_table_net()

  expect_named(network_tables(net, "summary"), c("metric", "count"))
  expect_equal(network_tables(net, "summary")$count,
               c(igraph::vcount(net$network), igraph::ecount(net$network),
                 nrow(net$cluster_stats)))
  expect_named(network_tables(net, "cluster_stats"), c("cluster", "avg_degree", "n_terms"))

  membership <- network_tables(net, "membership")
  expect_true(all(c(names(net$go_term_df), "genes") %in% names(membership)))
  expect_false(is.unsorted(membership$cluster))
  expect_equal(network_tables(net, "themes"), net$cluster_themes)

  all <- network_tables(net)
  expect_named(all, c("summary", "cluster_stats", "membership", "themes"))
  expect_named(network_tables(net, c("themes", "summary")), c("themes", "summary"))
})

test_that("CSV files match network_tables()", {
  outdir <- file.path(tempdir(), "test_network_tables")
  net <- make_table_net(outdir)

  csv <- utils::read.csv(file.path(outdir, "tables-cluster-membership.csv"))
  tab <- network_tables(net, "membership")
  expect_equal(csv$go_term, tab$go_term)
  expect_equal(csv$cluster, tab$cluster)
  expect_equal(csv$genes, tab$genes)

  csv <- utils::read.csv(file.path(outdir, "tables_network_summary.csv"))
  expect_equal(csv$count, network_tables(net, "summary")$count)

  unlink(outdir, recursive = TRUE)
})
