#' Run Gene Ontology Over-Representation Analysis on Multiple DEG Sets
#'
#' Performs GO enrichment analysis for each cluster-specific DEG list using
#' \code{clusterProfiler::enrichGO}. Produces a QC summary table of tested
#' versus retained gene sets per cluster.
#'
#' @param de_list The output of \code{\link{run_cluster_de}}: a named list
#'   whose elements each contain a \code{$degs} data frame (with a \code{gene}
#'   column of DEG symbols) and a \code{$tested} data frame (with a
#'   \code{gene} column of every gene tested in that cluster).
#' @param org_db An \code{OrgDb} annotation object (e.g.,
#'   \code{org.Mm.eg.db::org.Mm.eg.db} for mouse,
#'   \code{org.Hs.eg.db::org.Hs.eg.db} for human).
#' @param universe Character vector of background gene symbols used as the
#'   enrichment universe for every cluster (e.g., \code{rownames(seurat_obj)}).
#'   If \code{NULL} (default), \code{enrichGO} uses all genes annotated in
#'   \code{org_db}. Ignored when \code{use_tested_universe = TRUE}.
#' @param use_tested_universe Logical. If \code{TRUE}, each cluster's
#'   universe is the \code{gene} column of \code{de_list[[cluster]]$tested}
#'   (the genes \code{FindMarkers} actually tested in that cluster), and
#'   \code{universe} is ignored. Default \code{FALSE}.
#' @param keyType Character string specifying the gene identifier type.
#'   Default \code{"SYMBOL"}.
#' @param qvalueCutoff Numeric. Threshold on the Storey q-value (the
#'   \code{qvalue} column, an estimate of the false discovery rate from the
#'   \pkg{qvalue} method) for retaining enriched terms. This is separate from
#'   the Benjamini-Hochberg adjusted p-value in \code{p.adjust}, which
#'   \code{enrichGO} filters at 0.05. Default \code{0.1}.
#' @param ont Character string specifying the GO sub-ontology to test:
#'   \code{"BP"} (Biological Process), \code{"MF"} (Molecular Function),
#'   \code{"CC"} (Cellular Component), or \code{"ALL"}. Default \code{"BP"}.
#' @param split_direction Logical. If \code{TRUE}, each cluster's DEGs are
#'   split by the sign of \code{avg_log2FC} and ORA is run separately on the
#'   up-regulated (\code{avg_log2FC > 0}) and down-regulated genes, using the
#'   same universe for both. Results from a split run must be written with
#'   \code{write_go(..., split_direction = TRUE)}. Default \code{FALSE}.
#' @param simplify Logical. If \code{TRUE}, redundant GO terms are removed from
#'   each result with \code{clusterProfiler::simplify()}, which groups terms by
#'   semantic similarity in the GO hierarchy (Wang measure) and keeps the term
#'   with the lowest adjusted p-value from each group. Applied to the
#'   significant terms only, and to the up and down results separately when
#'   \code{split_direction = TRUE}. Not available with \code{ont = "ALL"}.
#'   Default \code{FALSE}, because the overlap networks in
#'   \code{\link{run_comparison_networks}} already group redundant terms, and
#'   the number of terms in a network cluster is itself informative.
#' @param simplify_cutoff Numeric \code{[0, 1]}. Semantic-similarity cutoff
#'   passed to \code{clusterProfiler::simplify()}; terms more similar than
#'   this are treated as redundant. Higher values remove fewer terms. Only used
#'   when \code{simplify = TRUE}. Default \code{0.7}.
#' @param outputDir Character path to the directory where the QC summary CSV
#'   is written, with or without a trailing \code{/}. Created if it does not
#'   exist. If \code{NULL} (default), no file is written; the QC summary is
#'   still printed and attached to the result as the \code{"qc"} attribute.
#' @param outname Character filename for the QC summary CSV. Default
#'   \code{"GO_ORA_QC_summary.csv"}.
#'
#' @return A named list (one element per cluster). Elements are \code{NULL}
#'   for clusters with no DEGs. Otherwise:
#'   \itemize{
#'     \item \code{split_direction = FALSE}: a single \code{enrichResult}.
#'     \item \code{split_direction = TRUE}: a list with \code{$up} and
#'       \code{$down} \code{enrichResult} objects; either is \code{NULL} when
#'       that direction has no DEGs.
#'   }
#'   Each result table gains a \code{cluster} column, and a \code{direction}
#'   column when split. With \code{simplify = TRUE}, each result holds only
#'   the significant, non-redundant terms. A QC summary CSV is also written to
#'   \code{<outputDir>/<outname>}, with one row per cluster, or one row per
#'   cluster and direction when split. Its \code{tested} and \code{retained}
#'   columns count terms before simplifying; an \code{after_simplify} column
#'   is added when \code{simplify = TRUE}.
#'
#' @details For each cluster, the function:
#' \enumerate{
#'   \item Extracts unique gene symbols from the DE data frame, split into up
#'     and down sets when \code{split_direction = TRUE}.
#'   \item Calls \code{clusterProfiler::enrichGO} with \code{minGSSize = 5},
#'     \code{maxGSSize = 300}, BH adjustment, and \code{readable = TRUE}.
#'   \item Annotates the result with the cluster name.
#'   \item Optionally removes redundant terms with
#'     \code{clusterProfiler::simplify()}.
#'   \item Prints a per-cluster summary of tested/retained gene sets.
#' }
#'
#' @export
run_go_ora <- function(de_list,
                       org_db,
                       universe = NULL,
                       use_tested_universe = FALSE,
                       keyType = "SYMBOL",
                       qvalueCutoff = 0.1,
                       ont = "BP",
                       split_direction = FALSE,
                       simplify = FALSE,
                       simplify_cutoff = 0.7,
                       outputDir = NULL,
                       outname = "GO_ORA_QC_summary.csv") {

  if (use_tested_universe && !is.null(universe)) {
    cli::cli_alert_info(
      "{.arg use_tested_universe} is {.val TRUE}; ignoring the supplied {.arg universe}."
    )
  }

  if (simplify && toupper(ont) == "ALL") {
    cli::cli_abort(c(
      "{.arg simplify = TRUE} is not available with {.code ont = \"ALL\"}.",
      "i" = "Use {.code ont = \"BP\"}, {.val MF}, or {.val CC}, or set {.code simplify = FALSE}."
    ))
  }

  # Term counts before simplifying, keyed by cluster and direction, for the QC table.
  pre_counts <- list()
  .key <- function(clust, direction = NULL) paste(c(clust, direction), collapse = " | ")

  .enrich <- function(genes, clust_universe, clust, direction = NULL) {
    label <- if (is.null(direction)) "" else paste0(" (", direction, ")")
    if (length(genes) == 0) {
      cli::cli_alert_warning("Skipping {.val {clust}}{label} (no DEGs)")
      return(NULL)
    }

    ego <- clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org_db,
      keyType       = keyType,
      universe      = clust_universe,
      ont           = ont,
      minGSSize     = 5,
      maxGSSize     = 300,
      pAdjustMethod = "BH",
      qvalueCutoff  = qvalueCutoff,
      readable      = TRUE
    )

    if (!is.null(ego)) {
      ego@result$cluster <- clust
      if (!is.null(direction)) ego@result$direction <- direction

      n_sig <- nrow(as.data.frame(ego))
      pre_counts[[.key(clust, direction)]] <<- c(tested = nrow(ego@result), retained = n_sig)

      if (simplify && n_sig > 0) {
        ego <- clusterProfiler::simplify(
          ego, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min
        )
        cli::cli_alert_info(
          "Simplified {.val {clust}}{label}: {n_sig} -> {nrow(as.data.frame(ego))} terms"
        )
      }
    }

    ego
  }

  ora_results <- lapply(names(de_list), function(clust) {

    degs <- de_list[[clust]]$degs
    tested <- de_list[[clust]]$tested

    if (is.null(degs) || nrow(degs) == 0) {
      cli::cli_alert_warning("Skipping {.val {clust}} (no DEGs)")
      return(NULL)
    }

    clust_universe <- universe
    if (use_tested_universe) {
      if (is.null(tested) || is.null(tested$gene) || length(tested$gene) == 0) {
        cli::cli_abort(c(
          "No tested genes found for {.val {clust}}.",
          "i" = "{.arg use_tested_universe = TRUE} needs {.code de_list[[\"{clust}\"]]$tested$gene}, as returned by {.fn run_cluster_de}."
        ))
      }
      clust_universe <- unique(tested$gene)
    }

    if (!split_direction) {
      return(.enrich(unique(degs$gene), clust_universe, clust))
    }

    if (is.null(degs$avg_log2FC)) {
      cli::cli_abort(c(
        "{.arg split_direction = TRUE} needs an {.field avg_log2FC} column in {.code de_list[[\"{clust}\"]]$degs}.",
        "i" = "Use the output of {.fn run_cluster_de}."
      ))
    }

    list(
      up   = .enrich(unique(degs$gene[degs$avg_log2FC > 0]), clust_universe, clust, "up"),
      down = .enrich(unique(degs$gene[degs$avg_log2FC < 0]), clust_universe, clust, "down")
    )
  })

  names(ora_results) <- names(de_list)

  .qc_row <- function(x, set, direction = NULL) {
    counts <- pre_counts[[.key(set, direction)]]
    row <- if (is.null(x) || is.null(counts) || counts[["tested"]] == 0) {
      label <- if (is.null(direction)) "" else paste0(" (", direction, ")")
      cli::cli_alert_warning("No results for {.val {set}}{label}")
      out <- data.frame(set = set, tested = 0, retained = 0, prop = NA)
      if (simplify) out$after_simplify <- 0
      out
    } else {
      total <- counts[["tested"]]
      filt <- counts[["retained"]]
      prop <- round(filt / total, 3)
      items <- c(
        "Genesets tested"      = total,
        "After filtering"      = filt,
        "Retention proportion" = prop
      )
      out <- data.frame(set = set, tested = total, retained = filt, prop = prop)
      if (simplify) {
        out$after_simplify <- nrow(as.data.frame(x))
        items <- c(items, "After simplify" = out$after_simplify)
      }
      cli::cli_dl(items)
      out
    }
    if (!is.null(direction)) row <- cbind(row[1], direction = direction, row[-1])
    row
  }

  qc_rows <- lapply(names(ora_results), function(i) {
    x <- ora_results[[i]]
    if (!split_direction) {
      cli::cli_h2("GO set: {i}")
      return(.qc_row(x, i))
    }
    do.call(rbind, lapply(c("up", "down"), function(d) {
      cli::cli_h2("GO set: {i} ({d})")
      .qc_row(x[[d]], i, d)
    }))
  })
  qc_df <- do.call(rbind, qc_rows)

  attr(ora_results, "qc") <- qc_df
  if (!is.null(outputDir)) {
    utils::write.csv(qc_df, file = paste0(.as_dir(outputDir), outname))
  }
  return(ora_results)
}
