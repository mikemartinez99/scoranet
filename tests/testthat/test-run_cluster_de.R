# Simulated Seurat object: two clusters of 60 cells, half "Car" and half
# "Bleo". Genes 1-20 are strongly up in Car within ClusterA. A third cluster
# has only 4 cells, too few for FindMarkers.
mock_seurat <- function(split_layers = FALSE) {
  set.seed(1)
  counts <- matrix(
    stats::rpois(200 * 124, 2), nrow = 200,
    dimnames = list(paste0("G", 1:200), paste0("c", 1:124))
  )
  counts[1:20, 1:30] <- counts[1:20, 1:30] + 10

  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts))
  obj$celltype <- c(rep("ClusterA", 60), rep("ClusterB", 60), rep("Tiny", 4))
  obj$sample <- c(
    rep("Car", 30), rep("Bleo", 30),
    rep(c("Car", "Bleo"), 30),
    "Car", rep("Bleo", 3)
  )

  if (split_layers) {
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
  }

  Seurat::NormalizeData(obj, verbose = FALSE)
}

run_de_quiet <- function(obj, ...) {
  suppressMessages(suppressWarnings(run_cluster_de(
    obj, cluster_col = "celltype", ident_col = "sample",
    ident1 = "Car", ident2 = "Bleo", ...
  )))
}

test_that("run_cluster_de returns degs and tested for each cluster", {
  res <- run_de_quiet(mock_seurat(), clusters = c("ClusterA", "ClusterB"))

  expect_named(res, c("ClusterA", "ClusterB"))
  expect_named(res$ClusterA, c("degs", "tested"))
  expect_gt(nrow(res$ClusterA$degs), 0)
  expect_true(all(res$ClusterA$degs$gene %in% res$ClusterA$tested$gene))
  spiked <- res$ClusterA$degs[res$ClusterA$degs$gene %in% paste0("G", 1:20), ]
  expect_gte(nrow(spiked), 15)
  expect_true(all(spiked$res == "Upregulated in Car"))
})

test_that("run_cluster_de passes seurat_logfc_threshold to FindMarkers", {
  obj <- mock_seurat()

  loose  <- run_de_quiet(obj, clusters = "ClusterA", seurat_logfc_threshold = 0)
  strict <- run_de_quiet(obj, clusters = "ClusterA", seurat_logfc_threshold = 1)

  expect_gt(nrow(loose$ClusterA$tested), nrow(strict$ClusterA$tested))
  expect_true(all(abs(strict$ClusterA$tested$avg_log2FC) >= 1))
})

test_that("run_cluster_de passes seurat_min_pct to FindMarkers", {
  obj <- mock_seurat()

  loose  <- run_de_quiet(obj, clusters = "ClusterA", seurat_min_pct = 0,
                         seurat_logfc_threshold = 0)
  strict <- run_de_quiet(obj, clusters = "ClusterA", seurat_min_pct = 0.99,
                         seurat_logfc_threshold = 0)

  expect_gt(nrow(loose$ClusterA$tested), nrow(strict$ClusterA$tested))
  expect_true(all(pmax(strict$ClusterA$tested$pct.1,
                       strict$ClusterA$tested$pct.2) >= 0.99))
})

test_that("run_cluster_de skips a cluster that FindMarkers cannot test", {
  obj <- mock_seurat()

  expect_warning(
    res <- suppressMessages(run_cluster_de(
      obj, clusters = c("ClusterA", "Tiny"), cluster_col = "celltype",
      ident_col = "sample", ident1 = "Car", ident2 = "Bleo"
    )),
    "DE failed"
  )
  expect_null(res$Tiny$degs)
  expect_null(res$Tiny$tested)
  expect_gt(nrow(res$ClusterA$degs), 0)
})

test_that(".seurat_assay_version reports the assay class", {
  obj <- mock_seurat()
  expected <- if (inherits(obj[["RNA"]], "Assay5")) "v5" else "v3"

  expect_equal(.seurat_assay_version(obj, "RNA"), expected)
})

test_that("run_cluster_de joins split Seurat v5 layers before DE", {
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0", "needs Seurat v5")
  obj <- mock_seurat(split_layers = TRUE)

  expect_message(
    res <- suppressWarnings(run_cluster_de(
      obj, clusters = "ClusterA", cluster_col = "celltype",
      ident_col = "sample", ident1 = "Car", ident2 = "Bleo"
    )),
    "joining"
  )
  expect_gt(nrow(res$ClusterA$tested), 0)
  expect_gt(nrow(res$ClusterA$degs), 0)
})
