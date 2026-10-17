# Plot types: default title suffix, PNG suffix, PNG size, and the tables
# (network_tables() names) saved alongside the PNG.
.plot_specs <- list(
  regulation = list(title = " (pathway overlap)", file = "-regulation.png",
                    width = 10, height = 10,
                    tables = c("summary", "membership")),
  leiden     = list(title = " (Leiden)", file = "-clustered-leiden.png",
                    width = 12.5, height = 10,
                    tables = c("summary", "cluster_stats", "membership")),
  zscore     = list(title = " (direction z-score)", file = "-zscore.png",
                    width = 11, height = 10,
                    tables = c("summary", "membership")),
  themes     = list(title = " (cluster themes)", file = "-themes.png",
                    width = 15, height = 11,
                    tables = c("summary", "cluster_stats", "membership", "themes"))
)

# Table names and their CSV file suffixes.
.table_files <- c(
  summary       = "_network_summary.csv",
  cluster_stats = "-cluster-stats.csv",
  membership    = "-cluster-membership.csv",
  themes        = "-cluster-themes.csv"
)

#' Plot a GO-Term Overlap Network
#'
#' Draws one or more views of a network built by
#' \code{\link{generate_overlap_network}} (or \code{\link{run_comparison_networks}}
#' / \code{\link{recluster_network}}). All types except \code{themes} share the
#' layout stored in the network result, so nodes sit in the same place in each.
#'
#' @param net_result A \code{scoranet_network} object from
#'   \code{\link{generate_overlap_network}}, \code{\link{run_comparison_networks}}
#'   or \code{\link{recluster_network}}.
#' @param type Character. One or more of:
#'   \describe{
#'     \item{\code{"leiden"} (default)}{Nodes filled by Leiden cluster, with
#'       the cluster legend on the left.}
#'     \item{\code{"regulation"}}{Nodes filled by direction: up, down, or
#'       mixed.}
#'     \item{\code{"zscore"}}{Nodes filled by a continuous blue-white-red
#'       gradient of the direction z-score. Needs z-scores in the network
#'       (see \code{\link{build_network_inputs}}).}
#'     \item{\code{"themes"}}{The Leiden plot annotated with cluster themes: a
#'       shaded hull around each cluster with a callout naming its theme, and a
#'       layout that keeps each cluster together (so node positions differ from
#'       the other plots).}
#'     \item{\code{"all"}}{All four.}
#'   }
#'   In \code{regulation}, \code{leiden} and \code{themes}, node shape shows
#'   direction (circle = up in \code{ident1}, square = up in \code{ident2},
#'   diamond = mixed). Edge width shows Jaccard overlap; node size shows the
#'   number of DE genes in the term.
#' @param label_by Character. Text used to label each cluster in the
#'   \code{themes} plot:
#'   \describe{
#'     \item{\code{"ancestor"} (default)}{Up to \code{theme_max_ancestors}
#'       higher-level GO terms that the cluster's terms fall under in the GO
#'       hierarchy (via GO.db), with how many of the cluster's terms they
#'       cover, e.g. \code{"cytokine production / response to cytokine
#'       (22/39)"}. Falls back to the hub term when no ancestor qualifies or
#'       the results have no GO IDs.}
#'     \item{\code{"keywords"}}{The most distinctive words in the cluster's
#'       term names.}
#'     \item{\code{"phrase"}}{The top term name, as used by the wordclouds.}
#'     \item{\code{"hub"}}{The cluster's most connected term.}
#'   }
#' @param min_terms Integer. In the \code{themes} plot, clusters with fewer
#'   terms get a hull but no callout label, to limit clutter. They are still
#'   listed colour-coded in the cluster legend. Default \code{3}.
#' @param title Character plot title. If \code{NULL} (default), the title is
#'   the network's \code{collection_name} plus the plot type, e.g.
#'   \code{"<name> (Leiden)"}. A title given here is used for every type.
#' @param outdir Character path to an output directory. If non-\code{NULL},
#'   each plot is saved as a PNG along with the CSV tables it shows (see
#'   Side effects). The directory is created if it does not exist. Default
#'   \code{NULL} (nothing is written).
#' @param file_stem Character filename prefix. Defaults to the
#'   \code{file_stem} stored in \code{net_result}.
#'
#' @return One \code{ggplot} object when one type is drawn, otherwise a named
#'   list of \code{ggplot} objects. \code{zscore} is skipped (with a message)
#'   when the network has no z-scores and other types are requested; asking
#'   for \code{zscore} alone is then an error.
#'
#' @section Side effects:
#' When \code{outdir} is set, each plot is saved as
#' \code{<file_stem>-regulation.png}, \code{-clustered-leiden.png},
#' \code{-zscore.png} or \code{-themes.png}, and the
#' \code{\link{network_tables}} behind it are saved as CSVs in the same folder:
#' \tabular{ll}{
#'   \code{regulation} \tab summary, membership \cr
#'   \code{leiden} \tab summary, cluster_stats, membership \cr
#'   \code{zscore} \tab summary, membership \cr
#'   \code{themes} \tab summary, cluster_stats, membership, themes
#' }
#' Each CSV is written once per call.
#'
#' @seealso \code{\link{generate_overlap_network}},
#'   \code{\link{recluster_network}}, \code{\link{network_tables}}
#'
#' @export
plot_network <- function(net_result,
                         type = "leiden",
                         label_by = c("ancestor", "keywords", "phrase", "hub"),
                         min_terms = 3,
                         title = NULL,
                         outdir = NULL,
                         file_stem = NULL) {
  .check_net_result(net_result)
  type <- match.arg(type, c(names(.plot_specs), "all"), several.ok = TRUE)
  if ("all" %in% type) type <- names(.plot_specs)
  type <- unique(type)
  label_by <- match.arg(label_by)

  network <- net_result$network
  if ("zscore" %in% type && all(is.na(igraph::V(network)$zscore))) {
    if (length(type) == 1) {
      cli::cli_abort(c(
        "This network has no direction z-scores, so there is no {.val zscore} plot.",
        "i" = "Build it from {.fn build_network_inputs} output to get z-scores."
      ))
    }
    cli::cli_alert_info("Skipped the {.val zscore} plot: the network has no direction z-scores.")
    type <- setdiff(type, "zscore")
  }

  reg_labels <- c(paste0("Upregulated in ", net_result$ident1),
                  paste0("Upregulated in ", net_result$ident2), "Mixed")
  .title <- function(t) .or_default(title, paste0(net_result$collection_name, .plot_specs[[t]]$title))

  plots <- lapply(stats::setNames(type, type), function(t) {
    if (t == "themes") {
      .themes_plot(
        network, .grouped_layout(network, net_result$seed), reg_labels,
        title = .title(t), cluster_palette = net_result$cluster_palette,
        themes = net_result$cluster_themes, label_by = label_by,
        min_terms = min_terms
      )
    } else {
      .network_plot(
        network, net_result$layout, reg_labels,
        fill_by = if (t == "leiden") "cluster" else t,
        title = .title(t), cluster_palette = net_result$cluster_palette
      )
    }
  })

  if (!is.null(outdir)) {
    outdir <- .as_dir(outdir)
    stem <- .or_default(file_stem, net_result$file_stem)
    for (t in type) {
      spec <- .plot_specs[[t]]
      .save_network_plot(plots[[t]], paste0(outdir, stem, spec$file),
                         width = spec$width, height = spec$height)
    }
    tables <- unique(unlist(lapply(.plot_specs[type], `[[`, "tables")))
    .write_network_tables(net_result, outdir, stem, tables)
  }

  if (length(plots) == 1) plots[[1]] else plots
}

#' Tables Behind a GO-Term Overlap Network
#'
#' Returns the tables that \code{\link{generate_overlap_network}},
#' \code{\link{recluster_network}} and \code{\link{plot_network}} write as CSV
#' files, so they can be used in R without reading the files back.
#'
#' @param net_result A \code{scoranet_network} object from
#'   \code{\link{generate_overlap_network}}, \code{\link{run_comparison_networks}}
#'   or \code{\link{recluster_network}}.
#' @param which Character. One or more of:
#'   \describe{
#'     \item{\code{"summary"}}{Columns \code{metric} (\code{nodes},
#'       \code{edges}, \code{clusters}) and \code{count}. Saved as
#'       \code{<stem>_network_summary.csv}.}
#'     \item{\code{"cluster_stats"}}{Columns \code{cluster},
#'       \code{avg_degree} and \code{n_terms}. Saved as
#'       \code{<stem>-cluster-stats.csv}.}
#'     \item{\code{"membership"}}{Every GO term with its Leiden cluster: the
#'       columns of \code{go_term_df} plus \code{genes}, the term's DE genes
#'       separated by \code{"/"}. Sorted by cluster, then by degree. Saved as
#'       \code{<stem>-cluster-membership.csv}.}
#'     \item{\code{"themes"}}{\code{cluster_themes}: one row per cluster with
#'       its size and candidate labels. Saved as
#'       \code{<stem>-cluster-themes.csv}.}
#'     \item{\code{"all"} (default)}{All four.}
#'   }
#'
#' @return A data frame when one table is requested, otherwise a named list of
#'   data frames.
#'
#' @seealso \code{\link{generate_overlap_network}}, \code{\link{plot_network}}
#'
#' @export
network_tables <- function(net_result, which = "all") {
  .check_net_result(net_result)
  which <- match.arg(which, c(names(.table_files), "all"), several.ok = TRUE)
  if ("all" %in% which) which <- names(.table_files)
  which <- unique(which)

  build <- list(
    summary = function() {
      data.frame(
        metric = c("nodes", "edges", "clusters"),
        count  = c(igraph::vcount(net_result$network),
                   igraph::ecount(net_result$network),
                   nrow(net_result$cluster_stats)),
        stringsAsFactors = FALSE
      )
    },
    cluster_stats = function() {
      net_result$cluster_stats[, c("cluster", "avg_degree", "n_terms")]
    },
    membership = function() {
      membership <- net_result$go_term_df
      membership$genes <- vapply(net_result$term_genes[membership$go_term],
                                 function(g) paste(g, collapse = "/"), character(1))
      membership[order(membership$cluster, -membership$degree), ]
    },
    themes = function() net_result$cluster_themes
  )

  tables <- lapply(stats::setNames(which, which), function(w) build[[w]]())
  if (length(tables) == 1) tables[[1]] else tables
}

#' Re-cluster a GO-Term Overlap Network at a New Resolution
#'
#' Runs Leiden clustering again on the network stored in a network result,
#' without rebuilding the overlap matrix. Use it after
#' \code{\link{plot_leiden_sweep}} suggests a better resolution. Cluster
#' membership, cluster statistics and cluster themes are recalculated; the
#' network, overlap matrix and plot layout are kept, so plots before and after
#' re-clustering have the same node positions.
#'
#' To change \code{overlap_threshold}, call
#' \code{\link{generate_overlap_network}} again instead: terms that become
#' isolated are dropped when the network is built, and their statistics are
#' not kept.
#'
#' @param net_result A \code{scoranet_network} object from
#'   \code{\link{generate_overlap_network}} or
#'   \code{\link{run_comparison_networks}}.
#' @param resolution Numeric. The new Leiden resolution.
#' @param objective_function Character, \code{"modularity"} or \code{"CPM"}.
#'   Defaults to the objective stored in \code{net_result}.
#' @param seed Integer seed for Leiden clustering. Defaults to the seed stored
#'   in \code{net_result}. \code{NULL} disables seeding.
#' @param outdir Character path to an output directory. If non-\code{NULL},
#'   all four \code{\link{network_tables}} are written there as CSVs. Default
#'   \code{NULL}.
#' @param file_stem Character filename prefix. Defaults to the
#'   \code{file_stem} stored in \code{net_result}.
#'
#' @return The updated \code{scoranet_network} object. \code{clustering},
#'   \code{go_term_df}, \code{cluster_stats}, \code{cluster_themes},
#'   \code{cluster_palette}, \code{leiden_resolution},
#'   \code{objective_function}, \code{seed}, and the vertex attributes
#'   \code{cluster} and \code{color} are replaced.
#'
#' @seealso \code{\link{plot_leiden_sweep}}, \code{\link{plot_network}}
#'
#' @export
recluster_network <- function(net_result,
                              resolution,
                              objective_function = net_result$objective_function,
                              seed = net_result$seed,
                              outdir = NULL,
                              file_stem = NULL) {
  .check_net_result(net_result)
  objective_function <- match.arg(objective_function, c("modularity", "CPM"))

  net_result <- .apply_clustering(net_result, resolution, objective_function, seed)

  if (!is.null(outdir)) {
    .write_network_tables(net_result, outdir, .or_default(file_stem, net_result$file_stem))
  }
  net_result
}

#' @export
print.scoranet_network <- function(x, ...) {
  cat("<scoranet_network> ", x$collection_name, "\n",
      "  Nodes: ", igraph::vcount(x$network),
      "  Edges: ", igraph::ecount(x$network),
      "  Clusters: ", nrow(x$cluster_stats), "\n",
      "  Leiden resolution: ", x$leiden_resolution, " (", x$objective_function, ")",
      "  Overlap threshold: ", x$overlap_threshold, "\n",
      "  Plot with plot_network(); tables with network_tables().\n",
      sep = "")
  invisible(x)
}

#' @export
plot.scoranet_network <- function(x, y, ...) {
  p <- plot_network(x, ...)
  if (inherits(p, "ggplot")) print(p) else lapply(p, print)
  invisible(p)
}

# Write network_tables() as CSVs: `which` picks the tables.
.write_network_tables <- function(net_result, outdir, file_stem, which = "all") {
  outdir <- .as_dir(outdir)
  tables <- network_tables(net_result, which)
  if (is.data.frame(tables)) {
    tables <- stats::setNames(list(tables), setdiff(which, "all"))
  }
  for (nm in names(tables)) {
    file <- paste0(outdir, file_stem, .table_files[[nm]])
    utils::write.csv(tables[[nm]], file = file, row.names = FALSE)
    cli::cli_alert_success(
      "Saved {nm} CSV: {.file {normalizePath(file, winslash = '/', mustWork = FALSE)}}"
    )
  }
  invisible(file.path(outdir, paste0(file_stem, .table_files[names(tables)])))
}

# `outdir` with a trailing slash, created if missing.
.as_dir <- function(outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (grepl("/$", outdir)) outdir else paste0(outdir, "/")
}

# Stop unless `x` is a network result from this version of the package.
.check_net_result <- function(x) {
  needed <- c("network", "layout", "go_term_df", "cluster_themes")
  if (!is.list(x) || !all(needed %in% names(x)) || !inherits(x$network, "igraph")) {
    cli::cli_abort(c(
      "{.arg net_result} must be a network result from {.fn generate_overlap_network}.",
      "i" = "Results from an older build of scoranet have no stored layout; rebuild them."
    ))
  }
  invisible(x)
}
