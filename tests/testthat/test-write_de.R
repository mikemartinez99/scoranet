# Build one cluster's element in the run_cluster_de() output format.
mock_de_cluster <- function(degs, tested = degs) {
  list(degs = degs, tested = tested)
}

test_that("write_de writes DE CSV files", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_out/")

  de_list <- list(
    ClusterA = mock_de_cluster(
      degs = data.frame(
        gene = c("G1", "G2"),
        avg_log2FC = c(1.0, -0.5),
        p_val_adj = c(0.01, 0.03),
        stringsAsFactors = FALSE
      ),
      tested = data.frame(
        gene = c("G1", "G2", "G9"),
        avg_log2FC = c(1.0, -0.5, 0.1),
        p_val_adj = c(0.01, 0.03, 0.8),
        stringsAsFactors = FALSE
      )
    ),
    ClusterB = mock_de_cluster(
      degs = data.frame(
        gene = c("G3"),
        avg_log2FC = c(2.0),
        p_val_adj = c(0.001),
        stringsAsFactors = FALSE
      )
    )
  )

  write_de(de_list, outdir)

  expect_true(file.exists(paste0(outdir, "ClusterA_DE_Results.csv")))
  expect_true(file.exists(paste0(outdir, "ClusterB_DE_Results.csv")))
  expect_true(file.exists(paste0(outdir, "ClusterA_tested_genes.csv")))
  expect_true(file.exists(paste0(outdir, "ClusterB_tested_genes.csv")))

  expect_equal(nrow(read.csv(paste0(outdir, "ClusterA_DE_Results.csv"))), 2)
  expect_equal(nrow(read.csv(paste0(outdir, "ClusterA_tested_genes.csv"))), 3)

  unlink(outdir, recursive = TRUE)
})

test_that("write_de skips empty DE clusters", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_empty/")

  de_list <- list(
    ClusterA = mock_de_cluster(
      degs = data.frame(
        gene = character(0),
        avg_log2FC = numeric(0),
        stringsAsFactors = FALSE
      ),
      tested = data.frame(
        gene = c("G1", "G2"),
        avg_log2FC = c(0.1, -0.2),
        stringsAsFactors = FALSE
      )
    )
  )

  expect_message(
    write_de(de_list, outdir),
    "no DEGs"
  )

  expect_false(file.exists(paste0(outdir, "ClusterA_DE_Results.csv")))
  # Tested genes are still written when a cluster has no DEGs
  expect_true(file.exists(paste0(outdir, "ClusterA_tested_genes.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_de sanitizes cluster names with slashes and spaces", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_sanitize/")

  de_list <- list(
    `Cell Type/Sub 1` = mock_de_cluster(
      degs = data.frame(
        gene = "G1",
        avg_log2FC = 1.0,
        stringsAsFactors = FALSE
      )
    )
  )

  write_de(de_list, outdir)

  expect_true(file.exists(paste0(outdir, "Cell_Type_Sub_1_DE_Results.csv")))
  expect_true(file.exists(paste0(outdir, "Cell_Type_Sub_1_tested_genes.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_de returns invisible NULL", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_return/")

  de_list <- list(
    A = mock_de_cluster(
      degs = data.frame(gene = "G1", avg_log2FC = 1.0, stringsAsFactors = FALSE)
    )
  )

  result <- write_de(de_list, outdir)

  expect_null(result)

  unlink(outdir, recursive = TRUE)
})

test_that("write_de writes inside outSubDir without a trailing slash", {
  outdir <- file.path(tempdir(), "test_de_noslash")
  de_list <- list(ClusterA = mock_de_cluster(
    degs = data.frame(gene = "G1", avg_log2FC = 1, p_val_adj = 0.01)
  ))

  write_de(de_list, outdir)

  expect_true(file.exists(file.path(outdir, "ClusterA_DE_Results.csv")))
  expect_false(file.exists(paste0(outdir, "ClusterA_DE_Results.csv")))
  unlink(outdir, recursive = TRUE)
})
