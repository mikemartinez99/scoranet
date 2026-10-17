test_that("write_go handles GO without split", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_go_nosplit/")

  mock_ego <- data.frame(
    Description = "immune response",
    GeneRatio = "5/100",
    p.adjust = 0.01,
    geneID = "G1/G2/G3/G4/G5",
    stringsAsFactors = FALSE
  )

  go_list <- list(ClusterA = mock_ego)

  write_go(go_list, outdir)

  expect_true(file.exists(paste0(outdir, "ClusterA_GO_ORA_Results.csv")))

  written <- read.csv(paste0(outdir, "ClusterA_GO_ORA_Results.csv"))
  expect_true("geneRatio_num" %in% colnames(written))
  expect_true("geneRatio_den" %in% colnames(written))
  expect_equal(written$geneRatio_num, 5)
  expect_equal(written$geneRatio_den, 100)

  unlink(outdir, recursive = TRUE)
})

test_that("write_go handles GO with split direction", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_go_split/")

  mock_up <- data.frame(
    Description = "pathway_up",
    GeneRatio = "3/50",
    p.adjust = 0.02,
    geneID = "G1/G2/G3",
    stringsAsFactors = FALSE
  )

  mock_down <- data.frame(
    Description = "pathway_down",
    GeneRatio = "2/50",
    p.adjust = 0.04,
    geneID = "G4/G5",
    stringsAsFactors = FALSE
  )

  go_list <- list(ClusterA = list(up = mock_up, down = mock_down))

  write_go(go_list, outdir, split_direction = TRUE)

  expect_true(file.exists(paste0(outdir, "ClusterA_GO_ORA_up_Results.csv")))
  expect_true(file.exists(paste0(outdir, "ClusterA_GO_ORA_down_Results.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_go errors clearly when split_direction = TRUE gets unsplit results", {
  mock_ego <- data.frame(
    Description = "immune response",
    GeneRatio = "5/100",
    p.adjust = 0.01,
    stringsAsFactors = FALSE
  )

  expect_error(
    write_go(list(ClusterA = mock_ego), tempdir(), split_direction = TRUE),
    "up"
  )
})

test_that("write_go sanitizes cluster names with slashes and spaces", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_go_sanitize/")

  go_list <- list(
    `Cell Type/Sub 1` = data.frame(
      Description = "immune response",
      GeneRatio = "5/100",
      stringsAsFactors = FALSE
    )
  )

  write_go(go_list, outdir)

  expect_true(file.exists(paste0(outdir, "Cell_Type_Sub_1_GO_ORA_Results.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_go returns invisible NULL", {
  go_list <- list(A = data.frame(Description = "x", GeneRatio = "1/10"))
  outdir <- paste0(tempdir(), "/test_go_return/")

  expect_null(write_go(go_list, outdir))

  unlink(outdir, recursive = TRUE)
})

test_that("write_go errors when split results are written with split_direction = FALSE", {
  mock_ego <- data.frame(Description = "x", GeneRatio = "1/10")
  go_list <- list(ClusterA = list(up = mock_ego, down = mock_ego))

  expect_error(
    write_go(go_list, tempdir(), split_direction = FALSE),
    "split_direction = TRUE"
  )
})

test_that("write_go writes inside outSubDir without a trailing slash", {
  outdir <- file.path(tempdir(), "test_go_noslash")
  go_list <- list(ClusterA = data.frame(
    Description = "immune response", GeneRatio = "5/100", p.adjust = 0.01,
    geneID = "G1/G2/G3/G4/G5"
  ))

  write_go(go_list, outdir)

  expect_true(file.exists(file.path(outdir, "ClusterA_GO_ORA_Results.csv")))
  expect_false(file.exists(paste0(outdir, "ClusterA_GO_ORA_Results.csv")))
  unlink(outdir, recursive = TRUE)
})
