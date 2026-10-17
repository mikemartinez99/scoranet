test_that("build_network_inputs returns correct structure with valid input", {
  ego_df <- data.frame(
    Description = c("immune response", "cell migration", "apoptotic process"),
    p.adjust = c(0.01, 0.03, 0.04),
    geneID = c("GeneA/GeneB/GeneC", "GeneB/GeneD", "GeneA/GeneE"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE"),
    avg_log2FC = c(1.5, -0.8, 2.0, 0.3, -1.2),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_type(result, "list")
  expect_named(result, c("collection", "results_set", "de_genes"))
  expect_type(result$collection, "list")
  expect_s3_class(result$results_set, "data.frame")
  expect_type(result$de_genes, "character")
})

test_that("build_network_inputs creates correct collection from geneID", {
  ego_df <- data.frame(
    Description = c("pathway_A", "pathway_B"),
    p.adjust = c(0.01, 0.02),
    geneID = c("G1/G2/G3", "G2/G4"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("G1", "G2", "G3", "G4"),
    avg_log2FC = c(1.0, -1.0, 0.5, 0.8),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_equal(length(result$collection), 2)
  expect_equal(result$collection[["pathway_A"]], c("G1", "G2", "G3"))
  expect_equal(result$collection[["pathway_B"]], c("G2", "G4"))
})

test_that("build_network_inputs computes correct mean fold-change per pathway", {
  ego_df <- data.frame(
    Description = c("pathway_X"),
    p.adjust = c(0.01),
    geneID = c("G1/G2"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("G1", "G2", "G3"),
    avg_log2FC = c(2.0, 4.0, 10.0),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_equal(result$results_set["pathway_X", "avg_log2FC"], 3.0)
})

test_that("build_network_inputs returns NULL when no significant terms", {
  ego_df <- data.frame(
    Description = c("pathway_A"),
    p.adjust = c(0.1),
    geneID = c("G1/G2"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("G1", "G2"),
    avg_log2FC = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_null(result)
})

test_that("build_network_inputs returns NULL for NULL ego input", {
  de_df <- data.frame(
    gene = c("G1"),
    avg_log2FC = c(1.0),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(NULL, de_df)

  expect_null(result)
})

test_that("build_network_inputs removes duplicate descriptions", {
  ego_df <- data.frame(
    Description = c("pathway_A", "pathway_A", "pathway_B"),
    p.adjust = c(0.01, 0.02, 0.03),
    geneID = c("G1/G2", "G1/G3", "G4"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("G1", "G2", "G3", "G4"),
    avg_log2FC = c(1.0, 2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_equal(length(result$collection), 2)
  expect_true("pathway_A" %in% names(result$collection))
  expect_true("pathway_B" %in% names(result$collection))
})

test_that("build_network_inputs returns 0 FC when no DE genes overlap pathway", {
  ego_df <- data.frame(
    Description = c("pathway_A"),
    p.adjust = c(0.01),
    geneID = c("X1/X2"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("G1", "G2"),
    avg_log2FC = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_equal(result$results_set["pathway_A", "avg_log2FC"], 0)
})

test_that("build_network_inputs captures all unique DE genes", {
  ego_df <- data.frame(
    Description = c("pathway_A"),
    p.adjust = c(0.01),
    geneID = c("G1"),
    stringsAsFactors = FALSE
  )

  de_df <- data.frame(
    gene = c("G1", "G2", "G2", "G3"),
    avg_log2FC = c(1.0, 2.0, 2.0, 3.0),
    stringsAsFactors = FALSE
  )

  result <- build_network_inputs(ego_df, de_df)

  expect_equal(sort(result$de_genes), c("G1", "G2", "G3"))
})

test_that("build_network_inputs computes mean, median, counts, and z-score", {
  ego_df <- data.frame(
    Description = c("mostly_up", "balanced"),
    p.adjust = c(0.01, 0.02),
    geneID = c("G1/G2/G3/G4/G5/G6", "G1/G7"),
    stringsAsFactors = FALSE
  )
  # mostly_up: five genes at +0.6 and one outlier at -4
  de_df <- data.frame(
    gene = paste0("G", 1:7),
    avg_log2FC = c(rep(0.6, 5), -4, -0.6),
    stringsAsFactors = FALSE
  )

  rs <- build_network_inputs(ego_df, de_df)$results_set

  expect_named(rs, c("avg_log2FC", "median_log2FC", "n_up", "n_down", "zscore"))
  expect_equal(rs["mostly_up", "avg_log2FC"], mean(c(rep(0.6, 5), -4)))
  expect_lt(rs["mostly_up", "avg_log2FC"], 0)       # the mean is pulled negative
  expect_equal(rs["mostly_up", "median_log2FC"], 0.6)
  expect_equal(rs["mostly_up", "n_up"], 5)
  expect_equal(rs["mostly_up", "n_down"], 1)
  expect_equal(rs["mostly_up", "zscore"], (5 - 1) / sqrt(6))
  expect_equal(rs["balanced", "zscore"], 0)
})
