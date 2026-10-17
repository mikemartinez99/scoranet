make_test_inputs <- function(n_pathways = 5) {
  genes <- paste0("Gene", seq_len(20))
  pathway_names <- paste0("pathway_", LETTERS[seq_len(n_pathways)])

  collection <- lapply(seq_len(n_pathways), function(i) {
    sample(genes, size = sample(5:10, 1))
  })
  names(collection) <- pathway_names

  fc_values <- rnorm(n_pathways, mean = 0, sd = 1)
  results_set <- data.frame(
    avg_log2FC = fc_values,
    row.names = pathway_names
  )

  list(
    results_set = results_set,
    collection = collection,
    de_genes = genes
  )
}

test_that("generate_overlap_network returns expected list structure", {
  set.seed(42)
  inp <- make_test_inputs(8)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_network",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("network", "clustering", "overlap_matrix",
                          "go_term_df", "cluster_stats", "leiden_resolution",
                          "objective_function", "overlap_threshold", "seed",
                          "term_genes", "cluster_themes", "plots"))
})

test_that("generate_overlap_network produces igraph object", {
  set.seed(42)
  inp <- make_test_inputs(8)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  expect_s3_class(result$network, "igraph")
})

test_that("overlap_matrix is square with correct dimensions", {
  set.seed(42)
  inp <- make_test_inputs(6)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  mat <- result$overlap_matrix
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), ncol(mat))
  expect_true(all(diag(mat) == 1))
})

test_that("overlap_matrix values are between 0 and 1", {
  set.seed(42)
  inp <- make_test_inputs(6)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  expect_true(all(result$overlap_matrix >= 0))
  expect_true(all(result$overlap_matrix <= 1))
})

test_that("go_term_df has expected columns", {
  set.seed(42)
  inp <- make_test_inputs(8)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  expected_cols <- c("go_term", "cluster", "degree", "regulation",
                     "avg_log2FC", "leiden_resolution")
  expect_true(all(expected_cols %in% colnames(result$go_term_df)))
})

test_that("regulation labels are only 'up' or 'down'", {
  set.seed(42)
  inp <- make_test_inputs(8)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  regulations <- result$go_term_df$regulation
  expect_true(all(regulations %in% c("up", "down")))
})

test_that("cluster_stats has expected columns", {
  set.seed(42)
  inp <- make_test_inputs(8)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  expect_true(all(c("cluster", "avg_degree", "n_terms") %in%
                    colnames(result$cluster_stats)))
})

test_that("leiden_resolution is preserved in output", {
  set.seed(42)
  inp <- make_test_inputs(8)

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    leiden_resolution = 0.05,
    outdir = NULL,
    plot_to_screen = FALSE
  )

  expect_equal(result$leiden_resolution, 0.05)
})

test_that("generate_overlap_network writes files when outdir is set", {
  set.seed(42)
  inp <- make_test_inputs(8)
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_network_out/")

  result <- generate_overlap_network(
    results_set = inp$results_set,
    collection = inp$collection,
    de_genes = inp$de_genes,
    collection_name = "test_net",
    file_stem = "test_stem",
    outdir = outdir,
    plot_to_screen = FALSE
  )

  expect_true(file.exists(paste0(outdir, "test_stem_network_summary.csv")))
  expect_true(file.exists(paste0(outdir, "test_stem-cluster-stats.csv")))
  expect_true(file.exists(paste0(outdir, "test_stem-plain.png")))
  expect_true(file.exists(paste0(outdir, "test_stem-clustered-leiden.png")))

  unlink(outdir, recursive = TRUE)
})

test_that("generate_overlap_network warns with fewer than 2 pathways", {
  results_set <- data.frame(avg_log2FC = 1.0, row.names = "only_pathway")
  collection <- list(only_pathway = c("G1", "G2"))
  de_genes <- c("G1", "G2")

  expect_warning(
    try(
      generate_overlap_network(
        results_set = results_set,
        collection = collection,
        de_genes = de_genes,
        collection_name = "sparse_net",
        outdir = NULL,
        plot_to_screen = FALSE
      ),
      silent = TRUE
    ),
    "Not enough pathways"
  )
})

test_that("hyphens in names are replaced with underscores", {
  results_set <- data.frame(
    avg_log2FC = c(1.0, -0.5, 0.3),
    row.names = c("path-A", "path-B", "path-C")
  )
  collection <- list(
    `path-A` = c("G1", "G2", "G3"),
    `path-B` = c("G2", "G3", "G4"),
    `path-C` = c("G1", "G3", "G5")
  )
  de_genes <- c("G1", "G2", "G3", "G4", "G5")

  result <- generate_overlap_network(
    results_set = results_set,
    collection = collection,
    de_genes = de_genes,
    collection_name = "hyphen_test",
    outdir = NULL,
    plot_to_screen = FALSE
  )

  node_names <- igraph::V(result$network)$name
  expect_false(any(grepl("-", node_names)))
})

test_that("generate_overlap_network calls up, down, and mixed from the z-score", {
  collection <- list(
    up_term    = paste0("U", 1:6),
    up_term2   = paste0("U", 2:7),
    down_term  = paste0("D", 1:6),
    down_term2 = paste0("D", 2:7),
    mixed_term = c("U1", "U2", "D1", "D2"),
    mixed2     = c("U2", "U3", "D2", "D3")
  )
  results_set <- data.frame(
    avg_log2FC    = c(1, 1, -1, -1, 0.1, -0.1),
    median_log2FC = c(1, 1, -1, -1, 0, 0),
    n_up          = c(6, 6, 0, 0, 2, 2),
    n_down        = c(0, 0, 6, 6, 2, 2),
    zscore        = c(6, 6, -6, -6, 0, 0) / sqrt(6),
    row.names     = names(collection)
  )

  res <- generate_overlap_network(
    results_set, collection, de_genes = unique(unlist(collection)),
    outdir = NULL, plot_to_screen = FALSE, leiden_resolution = 0.5
  )

  reg <- stats::setNames(res$go_term_df$regulation, res$go_term_df$go_term)
  expect_equal(unname(reg[c("up_term", "down_term", "mixed_term")]),
               c("up", "down", "mixed"))
  expect_true(all(c("median_log2FC", "n_up", "n_down", "zscore") %in%
                    colnames(res$go_term_df)))
  expect_s3_class(res$plots$zscore, "ggplot")
})

test_that("generate_overlap_network falls back to fold-change sign without z-scores", {
  set.seed(42)
  inp <- make_test_inputs(8)

  res <- generate_overlap_network(
    inp$results_set, inp$collection, inp$de_genes,
    outdir = NULL, plot_to_screen = FALSE
  )

  expect_true(all(res$go_term_df$regulation %in% c("up", "down")))
  expect_true(all(is.na(res$go_term_df$zscore)))
  expect_null(res$plots$zscore)
  expect_s3_class(res$plots$plain, "ggplot")
})

test_that("generate_overlap_network clustering is reproducible and leaves the session seed alone", {
  inp <- make_test_inputs(10)
  run <- function() generate_overlap_network(
    inp$results_set, inp$collection, inp$de_genes,
    outdir = NULL, plot_to_screen = FALSE
  )$go_term_df$cluster

  set.seed(1); a <- run()
  set.seed(999); b <- run()
  expect_identical(a, b)

  set.seed(42); before <- stats::runif(1)
  set.seed(42); invisible(run()); after <- stats::runif(1)
  expect_identical(before, after)
})

test_that("generate_overlap_network returns NULL when no terms overlap", {
  collection <- list(a = c("G1", "G2"), b = c("G3", "G4"), c = c("G5", "G6"))
  results_set <- data.frame(avg_log2FC = c(1, -1, 1), row.names = names(collection))

  res <- suppressMessages(generate_overlap_network(
    results_set, collection, unlist(collection),
    outdir = NULL, plot_to_screen = FALSE
  ))

  expect_null(res)
})
