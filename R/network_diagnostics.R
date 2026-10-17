#' Leiden Resolution Sweep for a GO-Term Network
#'
#' Re-clusters a network across a range of Leiden resolutions, with several
#' seeds at each, to help choose \code{leiden_resolution}. Three metrics are
#' plotted against resolution:
#' \describe{
#'   \item{Clusters}{Number of Leiden clusters (mean and range across seeds).}
#'   \item{Modularity}{Standard (resolution 1) weighted modularity of the
#'     partition: how well separated the clusters are (mean and range).}
#'   \item{Stability (ARI)}{Mean adjusted Rand index between every pair of
#'     seeds. Values near 1 mean the partition does not depend on the seed.}
#' }
#' A good resolution sits on a plateau where the cluster count is steady and
#' stability is high. The resolution used for \code{net_result} is marked with
#' a dashed line.
#'
#' @param net_result One network result from \code{\link{run_comparison_networks}}
#'   or \code{\link{generate_overlap_network}}, or an \code{igraph} network.
#' @param resolutions Numeric vector of resolutions to test. Default
#'   \code{seq(0.2, 3, by = 0.2)}, suited to the modularity objective; use much
#'   smaller values for \code{"CPM"}.
#' @param n_seeds Integer. Seeds per resolution (at least 2). Default \code{10}.
#' @param objective_function Character, \code{"modularity"} or \code{"CPM"}.
#'   If \code{NULL} (default), the objective stored in \code{net_result} is used,
#'   or \code{"modularity"} for a bare \code{igraph} network.
#' @param seed Integer. First seed; seeds \code{seed, seed + 1, ...} are used.
#'   Default \code{1234}.
#' @param title Character plot title. Default \code{"Leiden resolution sweep"}.
#'
#' @return A \code{ggplot} object. The sweep results are in its \code{data}
#'   element (columns \code{resolution}, \code{metric}, \code{mean}, \code{lo},
#'   \code{hi}).
#'
#' @seealso \code{\link{run_network_diagnostics}}
#'
#' @importFrom igraph membership modularity compare
#'
#' @export
plot_leiden_sweep <- function(net_result,
                              resolutions = seq(0.2, 3, by = 0.2),
                              n_seeds = 10,
                              objective_function = NULL,
                              seed = 1234,
                              title = "Leiden resolution sweep") {

  network <- .diag_network(net_result)
  if (n_seeds < 2) {
    cli::cli_abort("{.arg n_seeds} must be at least 2 to measure stability.")
  }
  if (igraph::vcount(network) < 3) {
    cli::cli_abort("The network has fewer than 3 nodes; a resolution sweep is not informative.")
  }

  if (is.null(objective_function)) {
    objective_function <- .or_default(
      if (is.list(net_result)) net_result$objective_function, "modularity"
    )
  }
  objective_function <- match.arg(objective_function, c("modularity", "CPM"))
  current <- if (is.list(net_result) && identical(net_result$objective_function, objective_function)) {
    net_result$leiden_resolution
  }

  seeds <- seed + seq_len(n_seeds) - 1
  weights <- igraph::E(network)$weight
  seed_pairs <- utils::combn(n_seeds, 2)

  rows <- lapply(resolutions, function(r) {
    mems <- lapply(seeds, function(s) {
      igraph::membership(.run_leiden(network, r, objective_function, s))
    })
    n_clust <- vapply(mems, function(m) length(unique(m)), numeric(1))
    modul <- vapply(mems, function(m) igraph::modularity(network, m, weights = weights), numeric(1))
    ari <- apply(seed_pairs, 2, function(p) {
      a <- as.integer(mems[[p[1]]])
      b <- as.integer(mems[[p[2]]])
      # ARI is 0/0 when both partitions are a single cluster; identical
      # partitions agree perfectly, so score them as 1.
      if (identical(a, b)) return(1)
      igraph::compare(a, b, method = "adjusted.rand")
    })
    data.frame(
      resolution = r,
      metric = c("Clusters", "Modularity", "Stability (ARI)"),
      mean = c(mean(n_clust), mean(modul), mean(ari)),
      lo   = c(min(n_clust), min(modul), min(ari)),
      hi   = c(max(n_clust), max(modul), max(ari))
    )
  })
  sweep <- do.call(rbind, rows)
  sweep$metric <- factor(sweep$metric, levels = c("Clusters", "Modularity", "Stability (ARI)"))

  p <- ggplot2::ggplot(sweep, ggplot2::aes(x = .data$resolution, y = .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                         fill = "grey80", alpha = 0.6) +
    ggplot2::geom_line(colour = "grey25") +
    ggplot2::geom_point(colour = "grey25", size = 1.5) +
    ggplot2::facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = paste0(objective_function, " objective, ", n_seeds,
                        " seeds per resolution (line = mean, band = range)"),
      x = "Leiden resolution", y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95"))

  if (!is.null(current)) {
    p <- p + ggplot2::geom_vline(xintercept = current, linetype = "dashed", colour = "#D55E00")
  }
  p
}

#' Overlap Threshold Diagnostic for a GO-Term Network
#'
#' Shows how the Jaccard \code{overlap_threshold} shapes a network. Terms with
#' no overlap at or above the threshold are dropped as isolates, so a high
#' threshold can silently remove many terms. Two panels:
#' \enumerate{
#'   \item A histogram of all non-zero pairwise Jaccard overlaps between terms,
#'     with the threshold used for \code{net_result} marked.
#'   \item Terms kept, terms dropped as isolates, and edges, across a range of
#'     thresholds.
#' }
#'
#' @param net_result One network result from \code{\link{run_comparison_networks}}
#'   or \code{\link{generate_overlap_network}} (it must contain
#'   \code{overlap_matrix}).
#' @param thresholds Numeric vector of thresholds to evaluate. Default
#'   \code{seq(0.05, 0.5, by = 0.05)}.
#' @param title Character plot title. Default \code{"Overlap threshold diagnostic"}.
#'
#' @return A \code{patchwork} object combining two \code{ggplot}s. Edit the
#'   whole figure with \code{&} (e.g. \code{p & ggplot2::theme_minimal()}) or one
#'   panel with \code{p[[1]]} / \code{p[[2]]}.
#'
#' @seealso \code{\link{run_network_diagnostics}}
#'
#' @export
plot_overlap_threshold <- function(net_result,
                                   thresholds = seq(0.05, 0.5, by = 0.05),
                                   title = "Overlap threshold diagnostic") {

  mat <- net_result$overlap_matrix
  if (is.null(mat) || nrow(mat) < 2) {
    cli::cli_abort("{.arg net_result} needs an {.field overlap_matrix} with at least 2 terms.")
  }
  current <- net_result$overlap_threshold
  vals <- mat[upper.tri(mat)]
  n_terms <- nrow(mat)

  hist_df <- data.frame(jaccard = vals[vals > 0])
  p_hist <- ggplot2::ggplot(hist_df, ggplot2::aes(x = .data$jaccard)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey60", colour = "white") +
    ggplot2::labs(
      x = "Pairwise Jaccard overlap", y = "Term pairs",
      subtitle = paste0(length(hist_df$jaccard), " of ", length(vals),
                        " term pairs share at least one gene")
    ) +
    ggplot2::theme_bw()
  if (!is.null(current)) {
    p_hist <- p_hist + ggplot2::geom_vline(xintercept = current, linetype = "dashed",
                                           colour = "#D55E00")
  }

  sweep <- do.call(rbind, lapply(thresholds, function(t) {
    adj <- mat >= t
    diag(adj) <- FALSE
    kept <- sum(rowSums(adj) > 0)
    data.frame(
      threshold = t,
      panel = c("Terms", "Terms", "Edges"),
      series = c("Kept", "Dropped (isolates)", "Edges"),
      value = c(kept, n_terms - kept, sum(adj[upper.tri(adj)]))
    )
  }))
  sweep$panel <- factor(sweep$panel, levels = c("Terms", "Edges"))
  sweep$series <- factor(sweep$series, levels = c("Kept", "Dropped (isolates)", "Edges"))

  p_sweep <- ggplot2::ggplot(sweep, ggplot2::aes(x = .data$threshold, y = .data$value,
                                                 colour = .data$series)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("Kept" = "#009E73", "Dropped (isolates)" = "#D55E00", "Edges" = "grey30")
    ) +
    ggplot2::labs(x = "Overlap threshold", y = "Count") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95"),
                   legend.position = "bottom")
  if (!is.null(current)) {
    p_sweep <- p_sweep + ggplot2::geom_vline(xintercept = current, linetype = "dashed",
                                             colour = "#D55E00")
  }

  patchwork::wrap_plots(p_hist, p_sweep, ncol = 1) +
    patchwork::plot_annotation(title = title)
}

#' Driver Genes for Each Network Cluster
#'
#' For each Leiden cluster, finds the DE genes that appear in the most GO terms
#' of that cluster: the genes the cluster's theme is built on. Shown as a dot
#' plot with clusters as columns and genes as rows. Dot size is the share of
#' the cluster's terms that contain the gene; colour is the gene's log2 fold
#' change.
#'
#' A gene's \emph{share} in cluster \eqn{c} is
#' \eqn{n_{terms}(g, c) / N_c}: the number of the cluster's GO terms whose DE
#' genes include \eqn{g}, divided by the number of terms in the cluster. Genes
#' are ranked by \eqn{n_{terms}}, with ties broken by larger
#' \eqn{|\log_2 FC|}, and the top \code{top_n} are shown. This measures how
#' widely a gene is annotated across the cluster, not how specific it is to
#' that cluster: broadly annotated genes can drive several clusters, and in a
#' two-term cluster every shared gene has a share of 100\%. Dot sizes span the
#' range of shares shown (not 0--100\%), so read exact values from the legend.
#'
#' @param net_result One network result from \code{\link{run_comparison_networks}}
#'   or \code{\link{generate_overlap_network}} (it must contain
#'   \code{term_genes} and \code{go_term_df}).
#' @param de_df Optional data frame of DE results with columns \code{gene} and
#'   \code{avg_log2FC} (e.g. \code{de_results[["Adipocytes"]]$degs}), used to
#'   colour genes by fold change. If \code{NULL}, dots are grey.
#' @param top_n Integer. Genes shown per cluster. Default \code{5}.
#' @param title Character plot title. Default \code{"Driver genes per cluster"}.
#'
#' @return A \code{ggplot} object. Its \code{data} element holds, for every
#'   shown gene and cluster, \code{n_terms} (cluster terms containing the gene),
#'   \code{share}, and \code{log2FC}.
#'
#' @seealso \code{\link{run_network_diagnostics}}
#'
#' @export
plot_driver_genes <- function(net_result, de_df = NULL, top_n = 5,
                              title = "Driver genes per cluster") {

  term_genes <- net_result$term_genes
  gtd <- net_result$go_term_df
  if (is.null(term_genes)) {
    cli::cli_abort(c(
      "{.arg net_result} has no {.field term_genes} element.",
      "i" = "It was made by an older version of the package. Run {.code devtools::load_all()}, then re-run {.fn run_comparison_networks}."
    ))
  }
  if (is.null(gtd) || nrow(gtd) == 0) {
    cli::cli_abort("{.arg net_result} has no network nodes, so there are no driver genes to show.")
  }

  # Network gene names have "-" replaced by "_"; map back to the DE gene names.
  fc_lookup <- NULL
  name_lookup <- NULL
  if (!is.null(de_df)) {
    clean <- gsub("-", "_", de_df$gene)
    fc_lookup <- stats::setNames(de_df$avg_log2FC, clean)
    name_lookup <- stats::setNames(de_df$gene, clean)
  }

  rows <- lapply(sort(unique(gtd$cluster)), function(cl) {
    terms <- gtd$go_term[gtd$cluster == cl]
    genes <- unlist(lapply(term_genes[terms], unique), use.names = FALSE)
    if (length(genes) == 0) return(NULL)
    counts <- sort(table(genes), decreasing = TRUE)
    df <- data.frame(gene = names(counts), n_terms = as.integer(counts),
                     stringsAsFactors = FALSE)
    df$log2FC <- if (!is.null(fc_lookup)) unname(fc_lookup[df$gene]) else NA_real_
    df <- df[order(-df$n_terms, -abs(ifelse(is.na(df$log2FC), 0, df$log2FC))), ]
    df <- utils::head(df, top_n)
    df$cluster <- cl
    df$share <- df$n_terms / length(terms)
    df
  })
  drivers <- do.call(rbind, rows)
  if (!is.null(name_lookup)) {
    mapped <- unname(name_lookup[drivers$gene])
    drivers$gene <- ifelse(is.na(mapped), drivers$gene, mapped)
  }

  # Order genes by the first cluster they drive, so the plot reads diagonally.
  first <- tapply(drivers$cluster, drivers$gene, min)
  gene_order <- names(sort(first, decreasing = TRUE))
  drivers$gene <- factor(drivers$gene, levels = unique(gene_order))
  drivers$cluster_label <- factor(paste("Cluster", drivers$cluster),
                                  levels = paste("Cluster", sort(unique(drivers$cluster))))

  # Size spans the observed range of shares, not 0 to 100%, so close values
  # (e.g. 88% vs 100%) are still told apart. Equal shares get one size.
  share_range <- range(drivers$share)
  share_labels <- function(x) paste0(round(100 * x), "%")
  size_scale <- if (diff(share_range) < 1e-9) {
    ggplot2::scale_size(name = "Share of\ncluster terms", range = c(5, 5),
                        breaks = share_range[1], labels = share_labels)
  } else {
    ggplot2::scale_size(name = "Share of\ncluster terms", range = c(1.5, 7),
                        limits = share_range,
                        breaks = unique(seq(share_range[1], share_range[2], length.out = 4)),
                        labels = share_labels)
  }

  p <- ggplot2::ggplot(drivers, ggplot2::aes(x = .data$cluster_label, y = .data$gene,
                                             size = .data$share)) +
    size_scale +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (all(is.na(drivers$log2FC))) {
    p + ggplot2::geom_point(fill = "grey60", shape = 21, colour = "grey25")
  } else {
    lim <- max(abs(drivers$log2FC), na.rm = TRUE)
    p +
      ggplot2::geom_point(ggplot2::aes(fill = .data$log2FC), shape = 21, colour = "grey25") +
      ggplot2::scale_fill_gradient2(name = "log2FC", low = "#278EEA", mid = "white",
                                    high = "#E67A6E", midpoint = 0, limits = c(-lim, lim))
  }
}

#' Overview of DE, ORA, and Network Results Across Cell Types
#'
#' One figure summarising a comparison: for each cell type, the number of DEGs
#' (split by direction), enriched GO terms, and network clusters. Shows at a
#' glance where the biology is concentrated.
#'
#' @param de_list Output of \code{\link{run_cluster_de}}.
#' @param go_list Optional output of \code{\link{run_go_ora}} (combined or
#'   split by direction).
#' @param net_list Optional output of \code{\link{run_comparison_networks}}.
#' @param title Character plot title. Default \code{"Comparison overview"}.
#'
#' @return A \code{ggplot} object. Its \code{data} element holds the counts
#'   (columns \code{celltype}, \code{metric}, \code{direction}, \code{value}).
#'
#' @seealso \code{\link{run_network_diagnostics}}
#'
#' @export
plot_comparison_overview <- function(de_list, go_list = NULL, net_list = NULL,
                                     title = "Comparison overview") {

  idents <- .diag_idents(de_list)
  rows <- list()
  add <- function(celltype, metric, direction, value) {
    rows[[length(rows) + 1]] <<- data.frame(celltype = celltype, metric = metric,
                                            direction = direction, value = value)
  }
  n_sig <- function(x) if (is.null(x)) 0 else nrow(as.data.frame(x))

  for (ct in names(de_list)) {
    degs <- de_list[[ct]]$degs
    fc <- if (is.null(degs)) numeric(0) else degs$avg_log2FC
    add(ct, "DEGs", "up", sum(fc > 0))
    add(ct, "DEGs", "down", sum(fc < 0))

    if (!is.null(go_list)) {
      ego <- go_list[[ct]]
      if (.diag_is_split(ego)) {
        add(ct, "Enriched GO terms", "up", n_sig(ego$up))
        add(ct, "Enriched GO terms", "down", n_sig(ego$down))
      } else {
        add(ct, "Enriched GO terms", "all", n_sig(ego))
      }
    }

    if (!is.null(net_list)) {
      key <- gsub("/", "_", ct)
      n_clust <- function(net) if (is.null(net)) 0 else length(unique(net$go_term_df$cluster))
      if (!is.null(net_list[[key]])) {
        add(ct, "Network clusters", "all", n_clust(net_list[[key]]))
      } else {
        add(ct, "Network clusters", "up", n_clust(net_list[[paste0(key, "_up")]]))
        add(ct, "Network clusters", "down", n_clust(net_list[[paste0(key, "_down")]]))
      }
    }
  }

  overview <- do.call(rbind, rows)
  overview$metric <- factor(overview$metric,
                            levels = c("DEGs", "Enriched GO terms", "Network clusters"))
  overview$direction <- factor(overview$direction, levels = c("up", "down", "all"))
  overview$celltype <- factor(overview$celltype, levels = rev(names(de_list)))

  # Side-by-side bars; reversed dodge group puts "up" above "down" in each row.
  overview$dodge_group <- factor(overview$direction, levels = rev(levels(overview$direction)))
  dodge <- ggplot2::position_dodge2(width = 0.8, preserve = "single", padding = 0.1)

  ggplot2::ggplot(overview, ggplot2::aes(x = .data$value, y = .data$celltype,
                                         fill = .data$direction,
                                         group = .data$dodge_group)) +
    ggplot2::geom_col(width = 0.8, position = dodge) +
    ggplot2::geom_text(ggplot2::aes(label = .data$value), position = dodge,
                       hjust = -0.2, size = 3) +
    ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(x))),
                                expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c(up = "#E67A6E", down = "#278EEA", all = "grey60"),
      labels = c(up = paste0("Upregulated in ", idents[1]),
                 down = paste0("Upregulated in ", idents[2]),
                 all = "Combined"),
      drop = TRUE
    ) +
    ggplot2::labs(title = title, x = "Count", y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95"),
                   legend.position = "bottom")
}

#' Run All Network Diagnostic Figures
#'
#' Convenience wrapper that builds every diagnostic figure for a comparison:
#' the overview across cell types (\code{\link{plot_comparison_overview}}) and,
#' for each network, the Leiden resolution sweep
#' (\code{\link{plot_leiden_sweep}}), the overlap threshold diagnostic
#' (\code{\link{plot_overlap_threshold}}), and the driver-gene plot
#' (\code{\link{plot_driver_genes}}).
#'
#' @param net_list Output of \code{\link{run_comparison_networks}}.
#' @param de_list Optional output of \code{\link{run_cluster_de}}. Needed for
#'   the overview and to colour driver genes by fold change.
#' @param go_list Optional output of \code{\link{run_go_ora}}, used in the
#'   overview.
#' @param outdir Character path. If not \code{NULL}, every figure is saved
#'   there as a PNG: \code{comparison-overview.png} and, per network,
#'   \code{<network>-leiden-sweep.png}, \code{<network>-overlap-threshold.png},
#'   and \code{<network>-driver-genes.png}. Default \code{NULL}.
#' @param resolutions,n_seeds Passed to \code{\link{plot_leiden_sweep}}.
#' @param thresholds Passed to \code{\link{plot_overlap_threshold}}.
#' @param top_n Passed to \code{\link{plot_driver_genes}}.
#'
#' @return A list with \code{overview} (a \code{ggplot}, or \code{NULL} without
#'   \code{de_list}) and \code{networks}: a named list with one element per
#'   network, each a list of \code{leiden_sweep}, \code{overlap_threshold}, and
#'   \code{driver_genes}. Figures that could not be built (for example, a
#'   network too small for a sweep) are \code{NULL}, with a warning message.
#'
#' @export
run_network_diagnostics <- function(net_list,
                                    de_list = NULL,
                                    go_list = NULL,
                                    outdir = NULL,
                                    resolutions = seq(0.2, 3, by = 0.2),
                                    n_seeds = 10,
                                    thresholds = seq(0.05, 0.5, by = 0.05),
                                    top_n = 5) {

  if (!is.null(outdir) && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  .try_plot <- function(label, expr) {
    tryCatch(expr, error = function(e) {
      cli::cli_alert_warning("Skipped {label}: {conditionMessage(e)}")
      NULL
    })
  }
  .save <- function(plot, name, width, height) {
    if (!is.null(outdir) && !is.null(plot)) {
      .save_network_plot(plot, paste0(outdir, name), width = width, height = height)
    }
  }

  overview <- NULL
  if (!is.null(de_list)) {
    cli::cli_h2("Comparison overview")
    overview <- .try_plot("comparison overview",
                          plot_comparison_overview(de_list, go_list, net_list))
    .save(overview, "comparison-overview.png", width = 11,
          height = 2 + 0.6 * length(de_list))
  }

  # Map network names (cluster names with "/" -> "_", plus "_up"/"_down" when
  # split) back to the run_cluster_de() entries.
  de_by_key <- if (!is.null(de_list)) stats::setNames(de_list, gsub("/", "_", names(de_list)))

  networks <- lapply(names(net_list), function(nm) {
    cli::cli_h2("Diagnostics: {nm}")
    net <- net_list[[nm]]
    key <- sub("_(up|down)$", "", nm)
    de_entry <- .or_default(de_by_key[[nm]], de_by_key[[key]])
    stem <- .safe_filename(nm)

    figs <- list(
      leiden_sweep = .try_plot(paste(nm, "Leiden sweep"),
        plot_leiden_sweep(net, resolutions = resolutions, n_seeds = n_seeds,
                          title = paste0(nm, ": Leiden resolution sweep"))),
      overlap_threshold = .try_plot(paste(nm, "overlap threshold"),
        plot_overlap_threshold(net, thresholds = thresholds,
                               title = paste0(nm, ": overlap threshold diagnostic"))),
      driver_genes = .try_plot(paste(nm, "driver genes"),
        plot_driver_genes(net, de_df = de_entry$degs, top_n = top_n,
                          title = paste0(nm, ": driver genes per cluster")))
    )

    .save(figs$leiden_sweep, paste0(stem, "-leiden-sweep.png"), width = 7, height = 8)
    .save(figs$overlap_threshold, paste0(stem, "-overlap-threshold.png"), width = 9, height = 8)
    n_clust <- length(unique(net$go_term_df$cluster))
    .save(figs$driver_genes, paste0(stem, "-driver-genes.png"),
          width = 3 + 0.5 * n_clust, height = 2 + 0.25 * n_clust * top_n)
    figs
  })
  names(networks) <- names(net_list)

  invisible(list(overview = overview, networks = networks))
}

# Network from a network result or a bare igraph object.
.diag_network <- function(net_result) {
  if (inherits(net_result, "igraph")) return(net_result)
  if (is.list(net_result) && inherits(net_result$network, "igraph")) return(net_result$network)
  cli::cli_abort("{.arg net_result} must be a network result or an {.cls igraph} object.")
}

# TRUE for a split run_go_ora() element (a list with $up / $down).
.diag_is_split <- function(x) {
  is.list(x) && !is.data.frame(x) && any(c("up", "down") %in% names(x))
}

# Numerator / denominator identities from run_cluster_de() output.
.diag_idents <- function(de_list) {
  for (entry in de_list) {
    for (df in list(entry$tested, entry$degs)) {
      if (!is.null(df$comparison) && length(df$comparison) > 0) {
        parts <- strsplit(df$comparison[1], "_vs_", fixed = TRUE)[[1]]
        if (length(parts) == 2) return(parts)
      }
    }
  }
  c("ident1", "ident2")
}
