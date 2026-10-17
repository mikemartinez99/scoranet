#' Run Cluster-Specific Differential Expression Analysis
#'
#' Performs differential expression between two sample identities within each
#' specified cell-type cluster using \code{Seurat::FindMarkers}. Results are
#' filtered by adjusted p-value and absolute log2 fold-change thresholds, and
#' annotated with a regulation-direction label.
#'
#' @param obj A \code{Seurat} object containing the RNA assay and metadata
#'   columns specified by \code{cluster_col} and \code{ident_col}.
#' @param clusters Character vector of cluster names (values in the
#'   \code{cluster_col} metadata column) to iterate over.
#' @param cluster_col Character string naming the metadata column that holds
#'   cell-type cluster assignments. Default \code{"final_celltype"}.
#' @param ident_col Character string naming the metadata column to set as the
#'   active identity for the DE comparison (e.g., \code{"sample"},
#'   \code{"condition"}). Default \code{"sample"}.
#' @param ident1 Character string. The first identity level (numerator in the
#'   log fold-change calculation). Positive fold-changes indicate genes
#'   up-regulated in \code{ident1}.
#' @param ident2 Character string. The second identity level (denominator).
#' @param seurat_min_pct Numeric. Seurat \code{FindMarkers} pre-filter:
#'   minimum fraction of cells in either group that must express a gene for
#'   it to be tested (passed as \code{FindMarkers(min.pct = )}). Like
#'   \code{seurat_logfc_threshold}, it controls which genes appear in
#'   \code{$tested}. Set explicitly because the \code{FindMarkers} default
#'   differs between Seurat 4 (0.1) and Seurat 5 (0.01). Default \code{0.1},
#'   the Seurat 4 default.
#' @param seurat_logfc_threshold Numeric. Seurat \code{FindMarkers}
#'   pre-filter: minimum absolute log2 fold-change for a gene to be tested at
#'   all (passed as \code{FindMarkers(logfc.threshold = )}). It controls which
#'   genes appear in \code{$tested}, and so the universe used by
#'   \code{run_go_ora(use_tested_universe = TRUE)}. Set explicitly because the
#'   \code{FindMarkers} default differs between Seurat 4 (0.25) and Seurat 5
#'   (0.1). Default \code{0.25}, the Seurat 4 default. This is not the DEG cutoff; use
#'   \code{logfc_thresh} for that.
#' @param padj_thresh Numeric. Maximum adjusted p-value for a gene to be
#'   retained in the filtered output. Default \code{0.05}.
#' @param logfc_thresh Numeric. DEG cutoff: minimum absolute log2 fold-change
#'   for a tested gene to be kept in \code{$degs}. Default \code{0.5}. Applied
#'   after \code{FindMarkers}; see \code{seurat_logfc_threshold} for the
#'   pre-filter that decides which genes are tested.
#'
#' @return A named list (one element per cluster). Each element is a list
#'   with two data frames:
#'   \describe{
#'     \item{degs}{Genes passing \code{padj_thresh} and \code{logfc_thresh}.
#'       Contains the standard \code{FindMarkers} columns plus \code{gene}
#'       (gene symbol), \code{cluster} (the cluster label), \code{comparison}
#'       (\code{"<ident1>_vs_<ident2>"}), and \code{res} (direction label:
#'       \code{"Upregulated in <ident1>"} or
#'       \code{"Upregulated in <ident2>"}).}
#'     \item{tested}{Unfiltered \code{FindMarkers} output: every gene tested
#'       in the cluster, with the same \code{gene}, \code{cluster}, and
#'       \code{comparison} columns.}
#'   }
#'   If \code{FindMarkers} fails for a cluster (for example, too few cells in
#'   one identity), a warning is issued and both elements are \code{NULL}.
#'
#' @details The function subsets the Seurat object to each cluster, sets the
#'   active identity to \code{ident_col}, and runs
#'   \code{Seurat::FindMarkers(assay = "RNA")}.
#'
#'   If the \code{RNA} assay is a Seurat 5 assay (class \code{Assay5}), its
#'   layers are joined with \code{JoinLayers} first, because
#'   \code{FindMarkers} cannot run on layers split by sample (for example,
#'   after \code{merge()}). Seurat 3/4 assays are used as they are.
#'
#' @importFrom Seurat FindMarkers Idents
#'
#' @export
run_cluster_de <- function(
    obj,
    clusters,
    cluster_col = "final_celltype",
    ident_col = "sample",
    ident1,
    ident2,
    seurat_min_pct = 0.1,
    seurat_logfc_threshold = 0.25,
    padj_thresh = 0.05,
    logfc_thresh = 0.5
) {

  if (.seurat_assay_version(obj, "RNA") == "v5") {
    cli::cli_alert_info("Seurat v5 assay detected; joining {.field RNA} layers.")
    # Looked up at runtime: JoinLayers only exists in SeuratObject >= 5.0,
    # and this branch only runs when a v5 assay (and so SeuratObject 5) is present.
    join_layers <- getExportedValue("SeuratObject", "JoinLayers")
    obj[["RNA"]] <- join_layers(obj[["RNA"]])
  }

  results_list <- lapply(clusters, function(clust) {

    cli::cli_progress_step("Running DE for {.val {clust}} ({ident1} vs {ident2})")

    sub <- obj[, obj@meta.data[[cluster_col]] == clust]

    Seurat::Idents(sub) <- ident_col

    markers <- tryCatch(
      Seurat::FindMarkers(
        sub,
        assay = "RNA",
        ident.1 = ident1,
        ident.2 = ident2,
        min.pct = seurat_min_pct,
        logfc.threshold = seurat_logfc_threshold,
        verbose = FALSE
      ),
      error = function(e) e
    )

    if (inherits(markers, "error")) {
      cli::cli_progress_done(result = "failed")
      cli::cli_warn(c(
        "DE failed for {.val {clust}}; skipping.",
        "x" = "{conditionMessage(markers)}"
      ))
      return(list(degs = NULL, tested = NULL))
    }

    markers$gene <- rownames(markers)
    markers$cluster <- clust
    markers$comparison <- paste0(ident1, "_vs_", ident2)

    markers.filt <- markers[
      markers$p_val_adj <= padj_thresh & abs(markers$avg_log2FC) > logfc_thresh,
    ]

    markers.filt$res <- ifelse(
      markers.filt$avg_log2FC > 0,
      paste0("Upregulated in ", ident1),
      paste0("Upregulated in ", ident2)
    )

    n_up <- sum(markers.filt$avg_log2FC > 0)
    n_down <- nrow(markers.filt) - n_up
    cli::cli_progress_done()
    cli::cli_bullets(c(
      "*" = "Tested genes: {nrow(markers)}",
      "*" = "DEGs: {nrow(markers.filt)}",
      "*" = "Up in {ident1}: {n_up}",
      "*" = "Up in {ident2}: {n_down}"
    ))
    cli::cli_text("")

    return(list(
        degs = markers.filt,
        tested = markers))
  })

  names(results_list) <- clusters
  return(results_list)
}
