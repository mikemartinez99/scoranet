# Mock run_cluster_de() output: ECM genes up, ribosomal genes down.
mock_ora_input <- function() {
  up   <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "FN1", "SPARC",
            "LUM", "DCN", "FBN1", "ELN", "POSTN", "COL6A1", "COL6A2", "COL6A3")
  down <- c("RPL3", "RPL4", "RPL5", "RPL7", "RPL8", "RPL10", "RPL11", "RPL13",
            "RPS2", "RPS3", "RPS4X", "RPS5", "RPS6", "RPS8", "RPS9")
  degs <- data.frame(
    gene = c(up, down),
    avg_log2FC = c(rep(1.5, length(up)), rep(-1.5, length(down))),
    stringsAsFactors = FALSE
  )
  list(
    Fib = list(degs = degs, tested = degs),
    Empty = list(degs = degs[0, ], tested = degs)
  )
}

test_that("run_go_ora returns one enrichResult per cluster by default", {
  skip_if_not_installed("org.Hs.eg.db")
  outdir <- paste0(tempdir(), "/test_ora_combined/")

  res <- suppressMessages(suppressWarnings(
    run_go_ora(mock_ora_input(), org.Hs.eg.db::org.Hs.eg.db, outputDir = outdir)
  ))

  expect_s4_class(res$Fib, "enrichResult")
  expect_null(res$Empty)

  qc <- read.csv(paste0(outdir, "GO_ORA_QC_summary.csv"))
  expect_false("direction" %in% colnames(qc))
  expect_equal(qc$set, c("Fib", "Empty"))

  unlink(outdir, recursive = TRUE)
})

test_that("run_go_ora with split_direction = TRUE returns up and down results", {
  skip_if_not_installed("org.Hs.eg.db")
  outdir <- paste0(tempdir(), "/test_ora_split/")

  res <- suppressMessages(suppressWarnings(
    run_go_ora(mock_ora_input(), org.Hs.eg.db::org.Hs.eg.db,
               split_direction = TRUE, outputDir = outdir)
  ))

  expect_named(res$Fib, c("up", "down"))
  expect_s4_class(res$Fib$up, "enrichResult")
  expect_s4_class(res$Fib$down, "enrichResult")
  expect_true(all(res$Fib$up@result$direction == "up"))
  expect_true(all(res$Fib$down@result$direction == "down"))
  expect_null(res$Empty)

  # Up DEGs are the ECM genes, so a collagen term should lead the up results
  expect_match(as.data.frame(res$Fib$up)$Description[1], "collagen|extracellular")

  qc <- read.csv(paste0(outdir, "GO_ORA_QC_summary.csv"))
  expect_equal(qc$direction, c("up", "down", "up", "down"))

  unlink(outdir, recursive = TRUE)
})

test_that("split run_go_ora output must be written with write_go(split_direction = TRUE)", {
  skip_if_not_installed("org.Hs.eg.db")
  qcdir <- paste0(tempdir(), "/test_ora_split_qc/")
  outdir <- paste0(tempdir(), "/test_ora_split_write/")

  res <- suppressMessages(suppressWarnings(
    run_go_ora(mock_ora_input(), org.Hs.eg.db::org.Hs.eg.db,
               split_direction = TRUE, outputDir = qcdir)
  ))

  expect_error(write_go(res, outdir, split_direction = FALSE), "split_direction = TRUE")

  suppressMessages(write_go(res, outdir, split_direction = TRUE))
  expect_true(file.exists(paste0(outdir, "Fib_GO_ORA_up_Results.csv")))
  expect_true(file.exists(paste0(outdir, "Fib_GO_ORA_down_Results.csv")))

  unlink(c(qcdir, outdir), recursive = TRUE)
})

test_that("run_go_ora with simplify = TRUE removes redundant terms", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("GOSemSim")
  outdir <- paste0(tempdir(), "/test_ora_simplify/")
  input <- mock_ora_input()

  full <- suppressMessages(suppressWarnings(
    run_go_ora(input, org.Hs.eg.db::org.Hs.eg.db, outputDir = outdir)
  ))
  simp <- suppressMessages(suppressWarnings(
    run_go_ora(input, org.Hs.eg.db::org.Hs.eg.db, simplify = TRUE,
               simplify_cutoff = 0.5, outputDir = outdir)
  ))

  n_full <- nrow(as.data.frame(full$Fib))
  n_simp <- nrow(as.data.frame(simp$Fib))
  expect_lt(n_simp, n_full)

  qc <- read.csv(paste0(outdir, "GO_ORA_QC_summary.csv"))
  expect_true("after_simplify" %in% colnames(qc))
  expect_equal(qc$retained[qc$set == "Fib"], n_full)
  expect_equal(qc$after_simplify[qc$set == "Fib"], n_simp)

  unlink(outdir, recursive = TRUE)
})

test_that("run_go_ora errors when simplify = TRUE is used with ont = 'ALL'", {
  skip_if_not_installed("org.Hs.eg.db")
  expect_error(
    run_go_ora(mock_ora_input(), org.Hs.eg.db::org.Hs.eg.db, ont = "ALL",
               simplify = TRUE, outputDir = tempdir()),
    "ALL"
  )
})
