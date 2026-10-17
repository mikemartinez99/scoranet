# Mock GO table: overlapping terms built from one set of genes.
mock_go_table <- function(genes, prefix) {
  data.frame(
    Description = paste0(prefix, "_term_", 1:4),
    p.adjust = 0.001,
    geneID = c(
      paste(genes[1:6], collapse = "/"),
      paste(genes[3:8], collapse = "/"),
      paste(genes[5:10], collapse = "/"),
      paste(genes[c(1:3, 8:10)], collapse = "/")
    ),
    stringsAsFactors = FALSE
  )
}

# Mock run_cluster_de() output with 10 up and 10 down DEGs.
mock_net_de <- function() {
  degs <- data.frame(
    gene = c(paste0("Up", 1:10), paste0("Dn", 1:10)),
    avg_log2FC = c(rep(1.5, 10), rep(-1.5, 10)),
    stringsAsFactors = FALSE
  )
  list(`Fib A/1` = list(degs = degs, tested = degs))
}

run_quiet <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  suppressMessages(suppressWarnings(run_comparison_networks(...)))
}

test_that("run_comparison_networks builds one network per cluster by default", {
  outdir <- paste0(tempdir(), "/test_net_combined/")
  go_list <- list(`Fib A/1` = mock_go_table(paste0("Up", 1:10), "comb"))

  res <- run_quiet(go_list, mock_net_de(), comp_label = "A_vs_B", base_outdir = outdir)

  expect_named(res, "Fib A_1")
  expect_s3_class(res[["Fib A_1"]]$network, "igraph")
  expect_s3_class(res[["Fib A_1"]]$plots$plain, "ggplot")
  expect_s3_class(res[["Fib A_1"]]$plots$leiden, "ggplot")
  expect_true(file.exists(paste0(outdir, "A_vs_B/networks/Fib_A_1-plain.png")))
  expect_true(file.exists(paste0(outdir, "A_vs_B/networks/Fib_A_1-clustered-leiden.png")))
  expect_s3_class(res[["Fib A_1"]]$plots$zscore, "ggplot")
  expect_true(file.exists(paste0(outdir, "A_vs_B/networks/Fib_A_1-zscore.png")))

  membership <- read.csv(paste0(outdir, "A_vs_B/networks/Fib_A_1-cluster-membership.csv"))
  expect_true(all(c("go_term", "cluster", "regulation", "zscore", "genes") %in% names(membership)))
  expect_equal(nrow(membership), nrow(res[["Fib A_1"]]$go_term_df))
  expect_false(is.unsorted(membership$cluster))
  expect_true(all(nzchar(membership$genes)))

  themes <- read.csv(paste0(outdir, "A_vs_B/networks/Fib_A_1-cluster-themes.csv"))
  expect_equal(nrow(themes), length(unique(membership$cluster)))

  unlink(outdir, recursive = TRUE)
})

test_that("run_comparison_networks with split_direction = TRUE builds up and down networks", {
  outdir <- paste0(tempdir(), "/test_net_split/")
  go_list <- list(`Fib A/1` = list(
    up   = mock_go_table(paste0("Up", 1:10), "up"),
    down = mock_go_table(paste0("Dn", 1:10), "down")
  ))

  res <- run_quiet(go_list, mock_net_de(), comp_label = "A_vs_B",
                   base_outdir = outdir, split_direction = TRUE)

  expect_named(res, c("Fib A_1_up", "Fib A_1_down"))
  expect_true(all(res[["Fib A_1_up"]]$go_term_df$regulation == "up"))
  expect_true(all(res[["Fib A_1_down"]]$go_term_df$regulation == "down"))

  netdir <- paste0(outdir, "A_vs_B/networks/")
  expect_true(file.exists(paste0(netdir, "Fib_A_1_up-plain.png")))
  expect_true(file.exists(paste0(netdir, "Fib_A_1_down-clustered-leiden.png")))

  unlink(outdir, recursive = TRUE)
})

test_that("run_comparison_networks skips an empty direction but keeps the other", {
  outdir <- paste0(tempdir(), "/test_net_split_empty/")
  go_list <- list(`Fib A/1` = list(
    up   = mock_go_table(paste0("Up", 1:10), "up"),
    down = NULL
  ))

  res <- run_quiet(go_list, mock_net_de(), comp_label = "A_vs_B",
                   base_outdir = outdir, split_direction = TRUE)

  expect_named(res, "Fib A_1_up")

  unlink(outdir, recursive = TRUE)
})

test_that("run_comparison_networks errors when split_direction does not match go_list", {
  outdir <- paste0(tempdir(), "/test_net_mismatch/")
  combined <- list(`Fib A/1` = mock_go_table(paste0("Up", 1:10), "comb"))
  split <- list(`Fib A/1` = list(
    up   = mock_go_table(paste0("Up", 1:10), "up"),
    down = mock_go_table(paste0("Dn", 1:10), "down")
  ))

  expect_error(
    run_quiet(split, mock_net_de(), comp_label = "A_vs_B", base_outdir = outdir),
    "split_direction = TRUE"
  )
  expect_error(
    run_quiet(combined, mock_net_de(), comp_label = "A_vs_B",
              base_outdir = outdir, split_direction = TRUE),
    "up"
  )

  unlink(outdir, recursive = TRUE)
})

test_that("run_comparison_networks errors on old-format de_list", {
  outdir <- paste0(tempdir(), "/test_net_oldde/")
  go_list <- list(`Fib A/1` = mock_go_table(paste0("Up", 1:10), "comb"))
  old_de <- list(`Fib A/1` = mock_net_de()[[1]]$degs)

  expect_error(
    run_quiet(go_list, old_de, comp_label = "A_vs_B", base_outdir = outdir),
    "run_cluster_de"
  )

  unlink(outdir, recursive = TRUE)
})
