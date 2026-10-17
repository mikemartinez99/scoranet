# Mock comparison: 8 modules of 4 overlapping GO terms each (up, down, or
# mixed genes), for two cell types.
mock_diag_inputs <- function() {
  set.seed(7)
  up <- paste0("Up", 1:40)
  dn <- paste0("Dn", 1:40)
  mods <- lapply(1:8, function(m) {
    pool <- switch(m %% 3 + 1, up[m * 4 + 1:8], dn[m * 4 + 1:8], c(up[m + 1:4], dn[m + 1:4]))
    data.frame(
      Description = paste0("m", m, "_term_", 1:4), p.adjust = 0.001,
      geneID = sapply(1:4, function(i) paste(sample(pool, 5), collapse = "/"))
    )
  })
  degs <- data.frame(
    gene = c(up, dn), avg_log2FC = c(runif(40, 0.5, 3), -runif(40, 0.5, 3)),
    comparison = "Car_vs_Bleo"
  )
  go_tab <- do.call(rbind, mods)
  list(
    de = list(`Fib A/1` = list(degs = degs, tested = degs),
              Myo = list(degs = degs[1:30, ], tested = degs)),
    go = list(`Fib A/1` = go_tab, Myo = go_tab[1:16, ])
  )
}

mock_diag_networks <- function(inp, outdir, ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  suppressMessages(run_comparison_networks(
    inp$go, inp$de, comp_label = "Car_vs_Bleo", base_outdir = outdir,
    overlap_threshold = 0.2, ...
  ))
}

test_that("plot_leiden_sweep returns a ggplot with complete sweep data", {
  inp <- mock_diag_inputs()
  outdir <- paste0(tempdir(), "/test_diag_sweep/")
  nets <- mock_diag_networks(inp, outdir)

  p <- plot_leiden_sweep(nets[["Fib A_1"]], resolutions = c(0.5, 1, 2), n_seeds = 3)

  expect_s3_class(p, "ggplot")
  expect_named(p$data, c("resolution", "metric", "mean", "lo", "hi"))
  expect_equal(nrow(p$data), 3 * 3)
  expect_false(anyNA(p$data$mean))
  ari <- p$data$mean[p$data$metric == "Stability (ARI)"]
  expect_true(all(ari >= -1 & ari <= 1))
  expect_error(plot_leiden_sweep(nets[["Fib A_1"]], n_seeds = 1), "at least 2")

  unlink(outdir, recursive = TRUE)
})

test_that("plot_overlap_threshold counts kept and dropped terms", {
  inp <- mock_diag_inputs()
  outdir <- paste0(tempdir(), "/test_diag_threshold/")
  net <- mock_diag_networks(inp, outdir)[["Fib A_1"]]

  p <- plot_overlap_threshold(net, thresholds = c(0.1, 0.5, 0.9))

  expect_s3_class(p, "patchwork")
  expect_s3_class(p[[1]], "ggplot")
  sweep <- p[[2]]$data
  terms <- sweep[sweep$panel == "Terms", ]
  totals <- tapply(terms$value, terms$threshold, sum)
  expect_true(all(totals == nrow(net$overlap_matrix)))
  kept <- terms$value[terms$series == "Kept"]
  expect_true(all(diff(kept) <= 0))   # stricter thresholds never keep more terms

  unlink(outdir, recursive = TRUE)
})

test_that("plot_driver_genes shows at most top_n genes per cluster", {
  inp <- mock_diag_inputs()
  outdir <- paste0(tempdir(), "/test_diag_drivers/")
  net <- mock_diag_networks(inp, outdir)[["Fib A_1"]]

  p <- plot_driver_genes(net, de_df = inp$de[["Fib A/1"]]$degs, top_n = 3)

  expect_s3_class(p, "ggplot")
  expect_true(all(table(p$data$cluster) <= 3))
  expect_true(all(p$data$share > 0 & p$data$share <= 1))
  expect_false(anyNA(p$data$log2FC))

  grey <- plot_driver_genes(net, top_n = 3)
  expect_true(all(is.na(grey$data$log2FC)))

  unlink(outdir, recursive = TRUE)
})

test_that("plot_driver_genes sizes dots across the observed share range", {
  inp <- mock_diag_inputs()
  outdir <- paste0(tempdir(), "/test_diag_driver_size/")
  net <- mock_diag_networks(inp, outdir)[["Fib A_1"]]

  p <- plot_driver_genes(net, top_n = 3)
  size <- ggplot2::ggplot_build(p)$data[[1]]$size
  if (diff(range(p$data$share)) > 0) {
    expect_equal(range(size), c(1.5, 7))
  } else {
    expect_true(all(size == 5))
  }

  unlink(outdir, recursive = TRUE)
})

test_that("plot_comparison_overview counts DEGs, terms, and clusters", {
  inp <- mock_diag_inputs()
  outdir <- paste0(tempdir(), "/test_diag_overview/")
  nets <- mock_diag_networks(inp, outdir)

  p <- plot_comparison_overview(inp$de, inp$go, nets)
  d <- p$data

  expect_s3_class(p, "ggplot")
  myo_up <- d$value[d$celltype == "Myo" & d$metric == "DEGs" & d$direction == "up"]
  expect_equal(myo_up, 30)
  expect_equal(d$value[d$celltype == "Myo" & d$metric == "Enriched GO terms"], 16)
  expect_equal(d$value[d$celltype == "Fib A/1" & d$metric == "Network clusters"],
               length(unique(nets[["Fib A_1"]]$go_term_df$cluster)))

  unlink(outdir, recursive = TRUE)
})

test_that("run_network_diagnostics builds and saves every figure, including split networks", {
  inp <- mock_diag_inputs()
  outdir <- paste0(tempdir(), "/test_diag_all/")
  split_go <- lapply(inp$go, function(g) list(up = g, down = g))
  nets <- mock_diag_networks(list(de = inp$de, go = split_go), outdir, split_direction = TRUE)
  diag_dir <- paste0(outdir, "diagnostics/")

  res <- suppressMessages(run_network_diagnostics(
    nets, inp$de, split_go, outdir = diag_dir,
    resolutions = c(0.5, 1), n_seeds = 2, thresholds = c(0.1, 0.3)
  ))

  expect_s3_class(res$overview, "ggplot")
  expect_true(all(c("up", "down") %in% res$overview$data$direction))
  expect_named(res$networks, names(nets))
  first <- res$networks[[1]]
  expect_named(first, c("leiden_sweep", "overlap_threshold", "driver_genes"))
  expect_false(anyNA(first$driver_genes$data$log2FC))   # split names map back to de_list
  expect_true(file.exists(paste0(diag_dir, "comparison-overview.png")))
  expect_true(file.exists(paste0(diag_dir, .safe_filename(names(nets)[1]), "-leiden-sweep.png")))

  unlink(outdir, recursive = TRUE)
})

test_that("plot_driver_genes explains an outdated network result", {
  old_result <- list(go_term_df = data.frame(go_term = "a", cluster = 1))
  expect_error(plot_driver_genes(old_result), "older version")
})

test_that("run_network_diagnostics saves figures inside outdir without a trailing slash", {
  inp <- mock_diag_inputs()
  outdir <- file.path(tempdir(), "test_diag_noslash")
  nets <- mock_diag_networks(inp, paste0(outdir, "/nets/"))
  diag_dir <- file.path(outdir, "diagnostics")

  suppressMessages(run_network_diagnostics(
    nets, inp$de, inp$go, outdir = diag_dir,
    resolutions = c(0.5, 1), n_seeds = 2, thresholds = c(0.1, 0.3)
  ))

  expect_true(file.exists(file.path(diag_dir, "comparison-overview.png")))
  expect_false(file.exists(paste0(diag_dir, "comparison-overview.png")))
  unlink(outdir, recursive = TRUE)
})
