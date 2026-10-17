#' Run GO-ORA Overlap Networks for All Clusters in a Comparison
#'
#' Iterates over every cluster in a named list of GO enrichment results,
#' builds a Jaccard-overlap network of enriched GO terms, applies Leiden
#' community detection, and saves network plots and summary tables. Clusters
#' that lack enrichment results or have fewer than two significant pathways
#' are skipped with a message. With \code{split_direction = TRUE}, a separate
#' network is built for the up- and down-regulated results of each cluster.
#'
#' @param go_list A named list of GO results from \code{\link{run_go_ora}}.
#'   With \code{split_direction = FALSE}, each element is a single
#'   \code{enrichResult}; with \code{split_direction = TRUE}, each element is
#'   a list with \code{$up} and \code{$down} \code{enrichResult} objects.
#'   Names must match the names in \code{de_list}.
#' @param de_list The output of \code{\link{run_cluster_de}}: a named list
#'   whose elements each contain a \code{$degs} data frame with at least the
#'   columns \code{gene} (gene symbol) and \code{avg_log2FC}.
#' @param comp_label Character string identifying the comparison (e.g.,
#'   \code{"CD206_Car_vs_Bleo"}). Used to build output sub-directory paths and
#'   as a prefix in plot titles.
#' @param base_outdir Character path to the root output directory. A
#'   sub-directory \code{<comp_label>/networks/} is created automatically.
#' @param overlap_threshold Numeric \code{[0, 1]}. Minimum Jaccard overlap
#'   between two GO term gene sets required to draw an edge. Default
#'   \code{0.1}.
#' @param leiden_resolution Numeric. Resolution for Leiden clustering. With
#'   the default modularity objective, values around 1 are standard; higher
#'   values give more, smaller clusters. Use \code{\link{plot_leiden_sweep}} to
#'   choose a value. Default \code{1}.
#' @param objective_function Character. Leiden objective passed to
#'   \code{igraph::cluster_leiden}: \code{"modularity"} (default; resolution
#'   is comparable across networks) or \code{"CPM"} (constant Potts model;
#'   resolution acts as an edge-density cutoff and needs much smaller values,
#'   e.g. 0.03).
#' @param seed Integer seed for Leiden clustering and the plot layout, so
#'   results are reproducible. It does not change the random state of your
#'   session. \code{NULL} disables seeding. Default \code{1234}.
#' @param split_direction Logical. Set to \code{TRUE} when \code{go_list} comes
#'   from \code{run_go_ora(split_direction = TRUE)}. Each direction is then
#'   built from its own results and only the DEGs in that direction
#'   (\code{avg_log2FC > 0} for up, \code{< 0} for down). Must match the
#'   structure of \code{go_list}; a mismatch is an error. Default \code{FALSE}.
#' @param zscore_threshold Numeric. Passed to
#'   \code{\link{generate_overlap_network}}: nodes with a direction z-score of
#'   at least this value are called up, at most minus this value down, and
#'   mixed otherwise. Default \code{1}.
#' @param theme_label Character. Text used to label each cluster in the
#'   themes plot:
#'   \describe{
#'     \item{\code{"ancestor"} (default)}{Up to \code{theme_max_ancestors}
#'       higher-level GO terms that the cluster's terms fall under in the GO
#'       hierarchy (via GO.db), with how many of the cluster's terms they
#'       cover, e.g. \code{"cytokine production / response to cytokine
#'       (22/39)"}. Each is chosen to cover many of the remaining terms while
#'       staying specific. Falls back to the hub term when no
#'       ancestor qualifies or the results have no GO IDs.}
#'     \item{\code{"keywords"}}{The most distinctive words in the cluster's
#'       term names.}
#'     \item{\code{"phrase"}}{The top term name, as used by the wordclouds.}
#'     \item{\code{"hub"}}{The cluster's most connected term.}
#'   }
#' @param theme_ancestor_min_share Numeric \code{(0, 1]}. Minimum share of a
#'   cluster's terms an ancestor must cover to be used as a label. Default
#'   \code{0.15}.
#' @param theme_max_ancestors Integer. Maximum number of ancestors in one
#'   cluster label. Default \code{2}.
#' @param theme_ancestor_max_size Integer. GO terms with more descendants than
#'   this are considered too generic to be a label (for reference: "cell
#'   cycle" has about 800, "response to stress" about 1,300, "cellular
#'   process" about 16,000). Default \code{1000}.
#' @param theme_min_terms Integer. Clusters with fewer terms get a hull but no
#'   callout label in the themes plot, to limit clutter. They are still listed
#'   colour-coded in the cluster legend. Default \code{3}.
#' @param theme_n_words Integer. Number of keywords per cluster label when
#'   \code{theme_label = "keywords"}. Default \code{3}.
#'
#' @return A flat named list of network result lists, one per network that
#'   could be built. Each element is the list returned by
#'   \code{\link{generate_overlap_network}}. Names are the cluster names (with
#'   \code{"/"} replaced by \code{"_"}), plus a \code{"_up"} or \code{"_down"}
#'   suffix when \code{split_direction = TRUE}.
#'
#' @section Output files:
#' Files are written to \code{<base_outdir>/<comp_label>/networks/} with the
#' file stem \code{<cluster>} (or \code{<cluster>_up} / \code{<cluster>_down}
#' when split), including \code{<stem>-cluster-membership.csv} (which GO
#' terms belong to which Leiden cluster, with their DE genes) and
#' \code{<stem>-cluster-themes.csv} (one row per cluster with its labels).
#' See \code{\link{generate_overlap_network}} for the full file list.
#'
#' @seealso \code{\link{build_network_inputs}},
#'   \code{\link{generate_overlap_network}}
#'
#' @export
run_comparison_networks <- function(
    go_list,
    de_list,
    comp_label,
    base_outdir,
    overlap_threshold = 0.1,
    leiden_resolution = 1,
    objective_function = c("modularity", "CPM"),
    seed = 1234,
    split_direction = FALSE,
    zscore_threshold = 1,
    theme_label = c("ancestor", "keywords", "phrase", "hub"),
    theme_min_terms = 3,
    theme_n_words = 3,
    theme_ancestor_min_share = 0.15,
    theme_ancestor_max_size = 1000,
    theme_max_ancestors = 2
) {

  net_outdir <- paste0(base_outdir, comp_label, "/networks/")
  if (!dir.exists(net_outdir)) dir.create(net_outdir, recursive = TRUE)

  .is_split <- function(x) {
    is.list(x) && !is.data.frame(x) && any(c("up", "down") %in% names(x))
  }

  # Numerator / denominator identities, from the comparison column that
  # run_cluster_de() adds ("<ident1>_vs_<ident2>").
  .idents <- function(de_entry) {
    for (df in list(de_entry$tested, de_entry$degs)) {
      if (!is.null(df$comparison) && length(df$comparison) > 0) {
        parts <- strsplit(df$comparison[1], "_vs_", fixed = TRUE)[[1]]
        if (length(parts) == 2) return(parts)
      }
    }
    NULL
  }

  .build <- function(ego, de_df, clust, direction = NULL, idents = NULL) {
    label <- if (is.null(direction)) "" else paste0(" (", direction, ")")
    file_stem <- .safe_filename(clust)
    if (!is.null(direction)) file_stem <- paste0(file_stem, "_", direction)

    if (is.null(ego) || is.null(de_df) || nrow(de_df) == 0) {
      cli::cli_alert_warning("Skipping network for {.val {clust}}{label}")
      return(NULL)
    }

    inputs <- build_network_inputs(ego, de_df)

    if (is.null(inputs) || length(inputs$collection) < 2) {
      cli::cli_alert_warning("Not enough pathways for {.val {clust}}{label}")
      return(NULL)
    }

    generate_overlap_network(
      results_set       = inputs$results_set,
      collection        = inputs$collection,
      de_genes          = inputs$de_genes,
      collection_name   = paste0(comp_label, "_", file_stem),
      file_stem         = file_stem,
      overlap_threshold = overlap_threshold,
      leiden_resolution = leiden_resolution,
      objective_function = objective_function,
      seed              = seed,
      outdir            = net_outdir,
      ident1            = idents[1],
      ident2            = idents[2],
      zscore_threshold  = zscore_threshold,
      theme_label       = theme_label,
      theme_min_terms   = theme_min_terms,
      theme_n_words     = theme_n_words,
      theme_ancestor_min_share = theme_ancestor_min_share,
      theme_ancestor_max_size = theme_ancestor_max_size,
      theme_max_ancestors = theme_max_ancestors
    )
  }

  results <- list()

  for (clust in names(go_list)) {

    ego <- go_list[[clust]]
    de_entry <- de_list[[clust]]

    if (is.data.frame(de_entry)) {
      cli::cli_abort(c(
        "{.code de_list[[\"{clust}\"]]} is a data frame, not {.fn run_cluster_de} output.",
        "i" = "Pass the list returned by {.fn run_cluster_de}, whose elements contain {.field degs}."
      ))
    }

    if (!split_direction && .is_split(ego)) {
      cli::cli_abort(c(
        "{.code go_list[[\"{clust}\"]]} contains {.field up}/{.field down} results from {.code run_go_ora(split_direction = TRUE)}.",
        "i" = "Use {.code run_comparison_networks(..., split_direction = TRUE)}."
      ))
    }
    if (split_direction && !is.null(ego) && !.is_split(ego)) {
      cli::cli_abort(c(
        "{.arg split_direction = TRUE} needs {.code go_list[[\"{clust}\"]]} to be a list with {.field up} and {.field down} results.",
        "i" = "{.fn run_go_ora} returns one combined result per cluster unless run with {.code split_direction = TRUE}."
      ))
    }

    degs <- de_entry$degs
    idents <- .idents(de_entry)
    clust_key <- gsub("/", "_", clust)

    if (!split_direction) {
      res <- .build(ego, degs, clust, idents = idents)
      if (!is.null(res)) results[[clust_key]] <- res
      next
    }

    for (direction in c("up", "down")) {
      dir_degs <- if (is.null(degs)) NULL else if (direction == "up") {
        degs[degs$avg_log2FC > 0, , drop = FALSE]
      } else {
        degs[degs$avg_log2FC < 0, , drop = FALSE]
      }
      res <- .build(ego[[direction]], dir_degs, clust, direction, idents)
      if (!is.null(res)) results[[paste0(clust_key, "_", direction)]] <- res
    }
  }

  return(results)
}

#' Build Network Input Data from GO Enrichment and DE Results
#'
#' Extracts significant GO terms (adjusted p-value < 0.05) from an
#' \code{enrichResult} object and pairs them with fold-change information from
#' a differential-expression table. The output is a list containing:
#' (1) a named list of gene-set members (\code{collection}), (2) a data frame
#' of per-pathway direction statistics (\code{results_set}), and (3) the
#' unique set of DE genes (\code{de_genes}). These three components are the
#' required inputs for \code{\link{generate_overlap_network}}.
#'
#' @param ego An \code{enrichResult} object (or coercible via
#'   \code{as.data.frame}) produced by \code{clusterProfiler::enrichGO} or a
#'   similar ORA function. Must contain columns \code{Description},
#'   \code{p.adjust}, and \code{geneID} (slash-delimited gene symbols).
#' @param de_df A data frame of differentially expressed genes with at least
#'   columns \code{gene} and \code{avg_log2FC}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{collection}{Named list of character vectors; each element contains
#'     the gene symbols associated with one GO term.}
#'   \item{results_set}{Data frame with one row per GO term (row names are
#'     the term descriptions), computed from the DE genes in each term:
#'     \describe{
#'       \item{avg_log2FC}{Mean log2 fold change.}
#'       \item{median_log2FC}{Median log2 fold change (less sensitive to
#'         single genes with extreme fold changes).}
#'       \item{n_up, n_down}{Number of the term's DE genes with positive and
#'         negative fold change.}
#'       \item{go_id}{GO ID of the term (\code{NA} if \code{ego} has no
#'         \code{ID} column).}
#'       \item{zscore}{Direction z-score,
#'         \eqn{(n_{up} - n_{down}) / \sqrt{n_{up} + n_{down}}}, as used by
#'         GOplot. Large positive values mean the term's genes are
#'         consistently up; values near 0 mean mixed.}
#'     }
#'     All values are 0 for a term with no DE genes.}
#'   \item{de_genes}{Character vector of unique DE gene symbols.}
#' }
#' Returns \code{NULL} if no significant terms pass filtering.
#'
#' @seealso \code{\link{generate_overlap_network}}
#'
#' @importFrom stats setNames median
#'
#' @export
build_network_inputs <- function(ego, de_df) {

  .extract <- function(ego_obj) {
    if (is.null(ego_obj) || nrow(as.data.frame(ego_obj)) == 0) return(NULL)
    res <- as.data.frame(ego_obj)
    res[!is.na(res$Description) & res$p.adjust < 0.05, ]
  }

  res <- .extract(ego)

  if (is.null(res) || nrow(res) == 0) return(NULL)

  res <- res[!duplicated(res$Description), ]

  collection <- stats::setNames(
    lapply(res$geneID, function(gids) unlist(strsplit(gids, "/"))),
    res$Description
  )

  de_genes <- unique(de_df$gene)
  fc_lookup <- stats::setNames(de_df$avg_log2FC, de_df$gene)

  results_set <- do.call(rbind, lapply(collection, function(genes) {
    fcs <- fc_lookup[intersect(genes, names(fc_lookup))]
    n_up <- sum(fcs > 0)
    n_down <- sum(fcs < 0)
    n <- n_up + n_down
    data.frame(
      avg_log2FC    = if (length(fcs) > 0) mean(fcs) else 0,
      median_log2FC = if (length(fcs) > 0) stats::median(fcs) else 0,
      n_up          = n_up,
      n_down        = n_down,
      zscore        = if (n > 0) (n_up - n_down) / sqrt(n) else 0
    )
  }))
  rownames(results_set) <- names(collection)
  # GO IDs, used to group a cluster's terms under a shared GO ancestor
  results_set$go_id <- if (!is.null(res$ID)) as.character(res$ID) else NA_character_

  list(collection = collection, results_set = results_set, de_genes = de_genes)
}

#' Generate a Jaccard-Overlap Network of GO Terms with Leiden Clustering
#'
#' Constructs an undirected, weighted network where nodes are enriched GO terms
#' and edge weights are Jaccard overlap coefficients between their DE gene
#' sets. Edges below \code{overlap_threshold} are removed. Leiden community
#' detection is applied to identify clusters of functionally related terms.
#' The function produces two PNG plots (a plain up/down-regulation view and a
#' Leiden-clustered view) and exports CSV summaries of node-level and
#' cluster-level statistics.
#'
#' @param results_set Data frame with row names equal to GO term descriptions
#'   and at least the column \code{avg_log2FC} (mean fold-change across member
#'   genes), as returned by \code{\link{build_network_inputs}}. If it also has
#'   \code{zscore}, \code{median_log2FC}, \code{n_up}, and \code{n_down}
#'   columns, node direction is called from the z-score (see
#'   \code{zscore_threshold}); otherwise from the sign of \code{avg_log2FC}.
#' @param collection Named list of character vectors. Each element is a set of
#'   gene symbols belonging to one GO term. Names must match
#'   \code{rownames(results_set)}.
#' @param de_genes Character vector of all differentially expressed gene
#'   symbols. Gene sets in \code{collection} are intersected with this vector
#'   before computing overlaps.
#' @param collection_name Character label used in plot titles and console
#'   output. Typically formatted as
#'   \code{"<comparison>_<celltype>"} (e.g.,
#'   \code{"CD206_Car_vs_Bleo_Fibroblast_1"}).
#' @param overlap_threshold Numeric \code{[0, 1]}. Jaccard coefficient below
#'   which edges are set to zero. Default \code{0.1}.
#' @param leiden_resolution Numeric. Resolution for Leiden clustering. With
#'   the default modularity objective, values around 1 are standard; higher
#'   values give more, smaller clusters. Use \code{\link{plot_leiden_sweep}} to
#'   choose a value. Default \code{1}.
#' @param objective_function Character. Leiden objective passed to
#'   \code{igraph::cluster_leiden}: \code{"modularity"} (default; resolution
#'   is comparable across networks) or \code{"CPM"} (constant Potts model;
#'   resolution acts as an edge-density cutoff and needs much smaller values,
#'   e.g. 0.03).
#' @param seed Integer seed for Leiden clustering and the plot layout, so
#'   results are reproducible. It does not change the random state of your
#'   session. \code{NULL} disables seeding. Default \code{1234}.
#' @param outdir Character path to the output directory. If \code{NULL},
#'   no files are written. The directory is created if it does not exist.
#' @param plot_to_screen Logical. If \code{TRUE} (default), the Leiden-
#'   clustered plot (\code{plots$leiden}) is also printed to the active
#'   graphics device (e.g., the RStudio Plots pane).
#' @param file_stem Character string used as the filename prefix for all
#'   exported files. If \code{NULL} or empty, \code{collection_name} is used.
#' @param zscore_threshold Numeric. A node is called \code{"up"} when its
#'   direction z-score is at least this value, \code{"down"} when it is at
#'   most minus this value, and \code{"mixed"} otherwise. Default \code{1}.
#' @param theme_label Character. Text used to label each cluster in the
#'   themes plot:
#'   \describe{
#'     \item{\code{"ancestor"} (default)}{Up to \code{theme_max_ancestors}
#'       higher-level GO terms that the cluster's terms fall under in the GO
#'       hierarchy (via GO.db), with how many of the cluster's terms they
#'       cover, e.g. \code{"cytokine production / response to cytokine
#'       (22/39)"}. Each is chosen to cover many of the remaining terms while
#'       staying specific. Falls back to the hub term when no
#'       ancestor qualifies or the results have no GO IDs.}
#'     \item{\code{"keywords"}}{The most distinctive words in the cluster's
#'       term names.}
#'     \item{\code{"phrase"}}{The top term name, as used by the wordclouds.}
#'     \item{\code{"hub"}}{The cluster's most connected term.}
#'   }
#' @param theme_ancestor_min_share Numeric \code{(0, 1]}. Minimum share of a
#'   cluster's terms an ancestor must cover to be used as a label. Default
#'   \code{0.15}.
#' @param theme_max_ancestors Integer. Maximum number of ancestors in one
#'   cluster label. Default \code{2}.
#' @param theme_ancestor_max_size Integer. GO terms with more descendants than
#'   this are considered too generic to be a label (for reference: "cell
#'   cycle" has about 800, "response to stress" about 1,300, "cellular
#'   process" about 16,000). Default \code{1000}.
#' @param theme_min_terms Integer. Clusters with fewer terms get a hull but no
#'   callout label in the themes plot, to limit clutter. They are still listed
#'   colour-coded in the cluster legend. Default \code{3}.
#' @param theme_n_words Integer. Number of keywords per cluster label when
#'   \code{theme_label = "keywords"}. Default \code{3}.
#' @param ident1,ident2 Character. The numerator and denominator identities of
#'   the comparison, used in the legend as \code{"Upregulated in <ident1>"}
#'   and \code{"Upregulated in <ident2>"}. If \code{NULL} (default), they are
#'   parsed from \code{collection_name} around \code{"_vs_"}.
#'
#' @return \code{NULL} (with a message) if no two terms overlap at
#'   \code{overlap_threshold}. Otherwise a list with the following elements:
#' \describe{
#'   \item{network}{An \code{igraph} graph object with vertex attributes
#'     \code{fc}, \code{regulation}, \code{size}, \code{cluster},
#'     \code{color}, and \code{shape}.}
#'   \item{clustering}{The \code{communities} object returned by
#'     \code{igraph::cluster_leiden}.}
#'   \item{overlap_matrix}{Square numeric matrix of pairwise Jaccard
#'     coefficients (before thresholding).}
#'   \item{go_term_df}{Data frame with columns \code{go_term},
#'     \code{cluster}, \code{degree}, \code{regulation} (\code{"up"},
#'     \code{"down"}, or \code{"mixed"}), \code{avg_log2FC} (mean),
#'     \code{median_log2FC}, \code{n_up}, \code{n_down}, \code{zscore},
#'     \code{go_id}, and \code{leiden_resolution}. The direction columns are
#'     \code{NA} when \code{results_set} has no z-score columns, and
#'     \code{go_id} is \code{NA} when it has no GO IDs.}
#'   \item{cluster_stats}{Data frame with columns \code{cluster},
#'     \code{avg_degree}, and \code{n_terms}.}
#'   \item{leiden_resolution, objective_function, overlap_threshold, seed}{The
#'     settings used.}
#'   \item{term_genes}{Named list: for each node (GO term), the DE genes it
#'     contains. Used by \code{\link{plot_driver_genes}}.}
#'   \item{cluster_themes}{Data frame with one row per Leiden cluster:
#'     \code{cluster}, \code{n_terms}; the shared GO ancestors
#'     (\code{ancestor}, names joined by \code{" / "}; \code{ancestor_id},
#'     IDs joined by \code{";"}; \code{ancestor_n_covered}, how many of the
#'     cluster's terms fall under them; \code{NA} if none qualifies); and the
#'     other candidate labels \code{keywords}, \code{phrase}, and \code{hub}.}
#'   \item{plots}{A list of \code{ggplot} objects built with \pkg{ggraph}:
#'     \code{plain} (nodes filled by direction: up, down, or mixed),
#'     \code{leiden} (nodes filled by Leiden cluster, with the cluster legend
#'     on the left), and \code{zscore} (nodes filled by a continuous
#'     blue-white-red gradient of the direction z-score). In \code{plain} and
#'     \code{leiden}, node shape shows direction (circle = up in
#'     \code{ident1}, square = up in \code{ident2}, diamond = mixed). Edge
#'     width shows Jaccard overlap. Cluster colours come from the 40-colour
#'     dittoSeq palette. \code{themes} is the Leiden plot annotated with
#'     cluster themes: a shaded hull around each cluster with a callout naming
#'     its theme, and a layout that keeps each cluster together (so
#'     node positions differ from the other plots). \code{zscore} is
#'     \code{NULL} when \code{results_set} has no \code{zscore} column.}
#' }
#'
#' @section Side effects:
#' When \code{outdir} is non-\code{NULL}, the following files are written:
#' \itemize{
#'   \item \code{<file_stem>_network_summary.csv} -- node/edge/cluster counts.
#'   \item \code{<file_stem>-cluster-stats.csv} -- per-cluster average degree.
#'   \item \code{<file_stem>-cluster-membership.csv} -- every GO term in the
#'     network with its Leiden cluster: the columns of \code{go_term_df} plus
#'     \code{genes}, the term's DE genes separated by \code{"/"}. Sorted by
#'     cluster, then by degree.
#'   \item \code{<file_stem>-cluster-themes.csv} -- one row per cluster with
#'     its size and candidate labels (\code{cluster_themes}).
#'   \item \code{<file_stem>-plain.png} -- \code{plots$plain}: network
#'     coloured by direction (red = up in \code{ident1}, blue = up in
#'     \code{ident2}, grey = mixed).
#'   \item \code{<file_stem>-clustered-leiden.png} -- \code{plots$leiden}:
#'     network coloured by Leiden cluster membership.
#'   \item \code{<file_stem>-zscore.png} -- \code{plots$zscore}: network
#'     coloured by direction z-score.
#'   \item \code{<file_stem>-themes.png} -- \code{plots$themes}: network
#'     with labelled cluster-theme hulls.
#' }
#'
#' @seealso \code{\link{build_network_inputs}},
#'   \code{\link{write_cluster_results}}
#'
#' @importFrom igraph graph_from_adjacency_matrix simplify V E cluster_leiden
#'   membership degree delete_vertices vcount ecount layout_with_fr
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @importFrom ggplot2 .data
#'
#' @export
generate_overlap_network <- function(
    results_set,
    collection,
    de_genes,
    collection_name = "network",
    overlap_threshold = 0.1,
    leiden_resolution = 1,
    objective_function = c("modularity", "CPM"),
    seed = 1234,
    outdir = NULL,
    plot_to_screen = TRUE,
    file_stem = NULL,
    ident1 = NULL,
    ident2 = NULL,
    zscore_threshold = 1,
    theme_label = c("ancestor", "keywords", "phrase", "hub"),
    theme_min_terms = 3,
    theme_n_words = 3,
    theme_ancestor_min_share = 0.15,
    theme_ancestor_max_size = 1000,
    theme_max_ancestors = 2
) {

  objective_function <- match.arg(objective_function)
  theme_label <- match.arg(theme_label)

  parts <- strsplit(collection_name, "_vs_")[[1]]
  if (is.null(ident1)) {
    ident1 <- if (length(parts) >= 1) gsub("_", " ", parts[1]) else "ident1"
  }
  if (is.null(ident2)) {
    ident2 <- if (length(parts) >= 2) gsub("_", " ", parts[2]) else "ident2"
  }
  reg_labels <- c(paste0("Upregulated in ", ident1), paste0("Upregulated in ", ident2), "Mixed")

  rownames(results_set) <- gsub("-", "_", rownames(results_set))
  de_genes <- unique(gsub("-", "_", de_genes))
  names(collection) <- gsub("-", "_", names(collection))

  collection_sub <- collection[names(collection) %in% rownames(results_set)]

  collection_sub <- lapply(collection_sub, function(g) {
    intersect(gsub("-", "_", g), de_genes)
  })
  collection_sub <- collection_sub[sapply(collection_sub, length) > 0]

  if (length(collection_sub) < 2) {
    warning(paste("Not enough pathways for:", collection_name))
  }

  nms <- names(collection_sub)
  mat <- matrix(0, nrow = length(nms), ncol = length(nms),
                dimnames = list(nms, nms))

  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      gi <- unique(collection_sub[[i]])
      gj <- unique(collection_sub[[j]])
      overlap    <- length(intersect(gi, gj))
      union_size <- length(unique(c(gi, gj)))
      mat[i, j]  <- ifelse(union_size > 0, overlap / union_size, 0)
    }
  }

  mat_thres <- mat
  mat_thres[mat_thres < overlap_threshold] <- 0

  network <- igraph::graph_from_adjacency_matrix(
    mat_thres, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  network <- igraph::simplify(network, remove.loops = TRUE, remove.multiple = TRUE)

  idx <- match(igraph::V(network)$name, rownames(results_set))
  .node_stat <- function(col) {
    if (col %in% colnames(results_set)) results_set[[col]][idx] else rep(NA_real_, length(idx))
  }
  fc <- results_set$avg_log2FC[idx]
  z  <- .node_stat("zscore")
  igraph::V(network)$fc        <- fc
  igraph::V(network)$median_fc <- .node_stat("median_log2FC")
  igraph::V(network)$n_up      <- .node_stat("n_up")
  igraph::V(network)$n_down    <- .node_stat("n_down")
  igraph::V(network)$zscore    <- z
  igraph::V(network)$go_id     <- if ("go_id" %in% colnames(results_set)) {
    as.character(results_set$go_id[idx])
  } else {
    rep(NA_character_, length(idx))
  }
  igraph::V(network)$regulation <- if (all(is.na(z))) {
    ifelse(fc > 0, "up", "down")
  } else {
    ifelse(z >= zscore_threshold, "up", ifelse(z <= -zscore_threshold, "down", "mixed"))
  }
  igraph::V(network)$size <- log(sapply(collection_sub[igraph::V(network)$name], length)) * 2

  isolates <- igraph::degree(network) == 0
  if (any(isolates)) {
    cli::cli_alert_info(
      "Dropped {sum(isolates)} of {length(isolates)} term{?s} with no overlap >= {overlap_threshold} (isolates)."
    )
  }
  network <- igraph::delete_vertices(network, isolates)

  # Isolates are removed, so any remaining node has an edge: fewer than 2 nodes
  # means no term overlaps another at this threshold.
  if (igraph::vcount(network) < 2) {
    cli::cli_alert_warning(
      "No network for {.val {collection_name}}: no terms overlap at >= {overlap_threshold}."
    )
    return(NULL)
  }

  igraph::E(network)$width <- igraph::E(network)$weight

  ceb <- .run_leiden(network, leiden_resolution, objective_function, seed)
  mem <- igraph::membership(ceb)
  igraph::V(network)$cluster <- as.integer(mem)

  n_clusters <- length(unique(mem))
  cluster_palette <- rep_len(.ditto_colors, n_clusters)
  igraph::V(network)$color <- cluster_palette[mem]
  igraph::V(network)$shape <- c(up = "circle", down = "square", mixed = "rectangle")[
    igraph::V(network)$regulation
  ]

  node_degrees <- igraph::degree(network)
  go_term_df <- data.frame(
    go_term    = igraph::V(network)$name,
    cluster    = igraph::V(network)$cluster,
    degree     = node_degrees,
    regulation = igraph::V(network)$regulation,
    avg_log2FC = igraph::V(network)$fc,
    median_log2FC = igraph::V(network)$median_fc,
    n_up       = igraph::V(network)$n_up,
    n_down     = igraph::V(network)$n_down,
    zscore     = igraph::V(network)$zscore,
    go_id      = igraph::V(network)$go_id,
    leiden_resolution = leiden_resolution,
    stringsAsFactors = FALSE
  )

  cluster_stats <- stats::aggregate(degree ~ cluster, data = go_term_df, FUN = mean)
  colnames(cluster_stats) <- c("cluster", "avg_degree")
  cluster_stats$n_terms <- as.integer(table(go_term_df$cluster))

  cli::cli_h2("Network Summary: {collection_name}")
  cli::cli_dl(c(
    "Nodes"             = igraph::vcount(network),
    "Edges"             = igraph::ecount(network),
    "Clusters"          = n_clusters,
    "Leiden resolution" = leiden_resolution,
    "Leiden objective"  = objective_function,
    "Overlap threshold" = overlap_threshold
  ))
  cli::cli_h3("Average degree per cluster")
  for (i in seq_len(nrow(cluster_stats))) {
    cli::cli_bullets(c(
      "*" = "Cluster {cluster_stats$cluster[i]}: avg_degree = {round(cluster_stats$avg_degree[i], 2)} ({cluster_stats$n_terms[i]} terms)"
    ))
  }
  cli::cli_rule()

  cluster_themes <- .cluster_themes(
    go_term_df, n_words = theme_n_words,
    ancestor_min_share = theme_ancestor_min_share,
    ancestor_max_size = theme_ancestor_max_size,
    max_ancestors = theme_max_ancestors
  )
  term_genes <- collection_sub[igraph::V(network)$name]

  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    if (is.null(file_stem) || !nzchar(file_stem)) {
      file_stem <- collection_name
    }

    summary_df <- data.frame(
      metric = c("nodes", "edges", "clusters"),
      count  = c(igraph::vcount(network), igraph::ecount(network), n_clusters),
      stringsAsFactors = FALSE
    )
    summary_file <- paste0(outdir, file_stem, "_network_summary.csv")
    utils::write.csv(summary_df, file = summary_file, row.names = FALSE)
    cli::cli_alert_success(
      "Saved network summary CSV: {.file {normalizePath(summary_file, winslash = '/', mustWork = FALSE)}} (exists={file.exists(summary_file)})"
    )

    cluster_stats_out <- cluster_stats[, c("cluster", "avg_degree", "n_terms")]
    cluster_stats_file <- paste0(outdir, file_stem, "-cluster-stats.csv")
    utils::write.csv(cluster_stats_out, file = cluster_stats_file, row.names = FALSE)
    cli::cli_alert_success(
      "Saved cluster stats CSV: {.file {normalizePath(cluster_stats_file, winslash = '/', mustWork = FALSE)}} (exists={file.exists(cluster_stats_file)})"
    )

    # Which GO terms belong to which cluster, with each term's DE genes
    membership <- go_term_df
    membership$genes <- vapply(term_genes[membership$go_term],
                               function(g) paste(g, collapse = "/"), character(1))
    membership <- membership[order(membership$cluster, -membership$degree), ]
    membership_file <- paste0(outdir, file_stem, "-cluster-membership.csv")
    utils::write.csv(membership, file = membership_file, row.names = FALSE)
    cli::cli_alert_success(
      "Saved cluster membership CSV: {.file {normalizePath(membership_file, winslash = '/', mustWork = FALSE)}}"
    )

    themes_file <- paste0(outdir, file_stem, "-cluster-themes.csv")
    utils::write.csv(cluster_themes, file = themes_file, row.names = FALSE)
    cli::cli_alert_success(
      "Saved cluster themes CSV: {.file {normalizePath(themes_file, winslash = '/', mustWork = FALSE)}}"
    )
  }

  lo <- .with_optional_seed(seed, igraph::layout_with_fr(network))

  plots <- list(plain = NULL, leiden = NULL, zscore = NULL, themes = NULL)
  if (igraph::vcount(network) > 0) {
    plots$plain <- .network_plot(
      network, lo, reg_labels, fill_by = "regulation",
      title = paste0(collection_name, " (pathway overlap)")
    )
    plots$leiden <- .network_plot(
      network, lo, reg_labels, fill_by = "cluster",
      title = paste0(collection_name, " (Leiden)"),
      cluster_palette = cluster_palette
    )
    if (!all(is.na(igraph::V(network)$zscore))) {
      plots$zscore <- .network_plot(
        network, lo, reg_labels, fill_by = "zscore",
        title = paste0(collection_name, " (direction z-score)")
      )
    }
    plots$themes <- .themes_plot(
      network, .grouped_layout(network, seed), reg_labels,
      title = paste0(collection_name, " (cluster themes)"),
      cluster_palette = cluster_palette, themes = cluster_themes,
      label_by = theme_label, min_terms = theme_min_terms
    )
  }

  if (!is.null(outdir) && !is.null(plots$plain)) {
    .save_network_plot(plots$plain, paste0(outdir, file_stem, "-plain.png"),
                       width = 10, height = 10)
    .save_network_plot(plots$leiden, paste0(outdir, file_stem, "-clustered-leiden.png"),
                       width = 12.5, height = 10)
    if (!is.null(plots$zscore)) {
      .save_network_plot(plots$zscore, paste0(outdir, file_stem, "-zscore.png"),
                         width = 11, height = 10)
    }
    .save_network_plot(plots$themes, paste0(outdir, file_stem, "-themes.png"),
                       width = 15, height = 11)
  }

  if (isTRUE(plot_to_screen) && !is.null(plots$leiden)) {
    print(plots$leiden)
  }

  list(
    network          = network,
    clustering       = ceb,
    overlap_matrix   = mat,
    go_term_df       = go_term_df,
    cluster_stats    = cluster_stats,
    leiden_resolution = leiden_resolution,
    objective_function = objective_function,
    overlap_threshold = overlap_threshold,
    seed             = seed,
    term_genes       = term_genes,
    cluster_themes   = cluster_themes,
    plots            = plots
  )
}

# Leiden clustering with the package's settings: edge weights, iteration until
# the partition is stable (n_iterations = -1), and an optional seed.
.run_leiden <- function(network, resolution, objective_function, seed) {
  .with_optional_seed(seed, igraph::cluster_leiden(
    network,
    objective_function = objective_function,
    weights = igraph::E(network)$weight,
    resolution = resolution,
    n_iterations = -1
  ))
}

# Build a ggraph plot of a GO-term overlap network. Edge width shows Jaccard
# overlap. `fill_by` chooses the node fill:
#   "regulation": up / down / mixed colours; shape also shows direction.
#   "cluster":    Leiden cluster (legend on the left); shape shows direction.
#   "zscore":     continuous blue-white-red gradient of the direction z-score.
.network_plot <- function(network, layout, reg_labels, title,
                          fill_by = c("regulation", "cluster", "zscore"),
                          cluster_palette = NULL) {
  fill_by <- match.arg(fill_by)
  reg_breaks <- c("up", "down", "mixed")
  reg_shapes <- c(up = 21, down = 22, mixed = 23)
  reg_names  <- stats::setNames(reg_labels, reg_breaks)

  p <- ggraph::ggraph(network, layout = "manual", x = layout[, 1], y = layout[, 2]) +
    ggraph::geom_edge_link(ggplot2::aes(edge_width = .data$weight),
                           colour = "grey70", alpha = 0.6) +
    ggraph::scale_edge_width(range = c(0.2, 1.5), guide = "none")

  if (fill_by == "regulation") {
    p <- p +
      ggraph::geom_node_point(
        ggplot2::aes(fill = .data$regulation, shape = .data$regulation,
                     size = .data$size),
        colour = "grey25", stroke = 0.4
      ) +
      ggplot2::scale_fill_manual(
        name = NULL, breaks = reg_breaks, labels = reg_names,
        values = c(up = "#E67A6E", down = "#278EEA", mixed = "grey75")
      ) +
      ggplot2::scale_shape_manual(
        name = NULL, values = reg_shapes, breaks = reg_breaks, labels = reg_names,
        guide = ggplot2::guide_legend(override.aes = list(size = 4))
      )
  } else if (fill_by == "cluster") {
    p <- p +
      ggraph::geom_node_point(
        ggplot2::aes(fill = factor(.data$cluster), shape = .data$regulation,
                     size = .data$size),
        colour = "grey25", stroke = 0.4
      ) +
      ggplot2::scale_fill_manual(
        name = "Cluster",
        values = stats::setNames(cluster_palette, seq_along(cluster_palette)),
        labels = function(x) paste("Cluster", x),
        guide = ggplot2::guide_legend(
          position = "left",
          override.aes = list(shape = 21, size = 4)
        )
      ) +
      # Fill shows cluster here, so the direction key gets a neutral fill
      ggplot2::scale_shape_manual(
        name = NULL, values = reg_shapes, breaks = reg_breaks, labels = reg_names,
        guide = ggplot2::guide_legend(override.aes = list(size = 4, fill = "grey60"))
      )
  } else {
    z_max <- max(abs(igraph::V(network)$zscore), na.rm = TRUE)
    if (!is.finite(z_max) || z_max == 0) z_max <- 1
    p <- p +
      ggraph::geom_node_point(
        ggplot2::aes(fill = .data$zscore, size = .data$size),
        shape = 21, colour = "grey25", stroke = 0.4
      ) +
      ggplot2::scale_fill_gradient2(
        name = paste0("Direction z-score\n(+ = ", reg_labels[1],
                      ",\n- = ", reg_labels[2], ")"),
        low = "#278EEA", mid = "white", high = "#E67A6E", midpoint = 0,
        limits = c(-z_max, z_max),
        guide = ggplot2::guide_colourbar(
          theme = ggplot2::theme(legend.key.height = ggplot2::unit(12, "lines"))
        )
      )
  }

  p +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

# Save a network plot as a PNG and report the result with cli.
.save_network_plot <- function(plot, file, width, height) {
  tryCatch({
    ggplot2::ggsave(file, plot, width = width, height = height, dpi = 300, bg = "white")
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
      try(system2("chflags", c("nohidden", file)), silent = TRUE)
    }
    cli::cli_alert_success(
      "Saved PNG: {.file {normalizePath(file, winslash = '/', mustWork = FALSE)}} (exists={file.exists(file)})"
    )
  }, error = function(e) {
    cli::cli_alert_danger("Failed to save {.file {file}}: {conditionMessage(e)}")
  })
}

#' Extract and Export GO Terms from a Specific Leiden Cluster
#'
#' Given a network result object and a single GO term of interest, identifies
#' which Leiden cluster that term belongs to and returns (and optionally writes)
#' a data frame of all GO terms sharing that cluster. This is useful for
#' inspecting the biological theme captured by a particular network community.
#'
#' @param net_result A list returned by \code{\link{generate_overlap_network}},
#'   which must contain at least the element \code{go_term_df}.
#' @param go_term Character string. The GO term description to look up (e.g.,
#'   \code{"extracellular matrix organization"}). Hyphens are internally
#'   replaced with underscores to match network node naming conventions.
#' @param outdir Character path to the output directory. If non-\code{NULL},
#'   a CSV file is written. The directory is created if it does not exist.
#' @param label Optional character label used in the output filename. If
#'   \code{NULL}, defaults to \code{"cluster-<id>"}.
#'
#' @return A data frame containing all GO terms in the same Leiden cluster as
#'   \code{go_term}, with the same columns as \code{go_term_df} in
#'   \code{\link{generate_overlap_network}} (including \code{regulation},
#'   \code{avg_log2FC}, \code{median_log2FC}, and \code{zscore}). Returns \code{NULL} with a warning if the
#'   requested term is not found in the network.
#'
#' @seealso \code{\link{generate_overlap_network}}
#'
#' @importFrom utils write.csv
#'
#' @export
write_cluster_results <- function(net_result, go_term, outdir, label = NULL) {

  go_term_df <- net_result$go_term_df
  target_cluster <- go_term_df$cluster[go_term_df$go_term == gsub("-", "_", go_term)]

  if (length(target_cluster) == 0) {
    warning("GO term '", go_term, "' not found in network")
    return(NULL)
  }

  cluster_terms <- go_term_df[go_term_df$cluster == target_cluster[1], ]

  if (!is.null(outdir)) {
    if (is.null(label)) label <- paste0("cluster-", target_cluster[1])
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    utils::write.csv(cluster_terms, paste0(outdir, "network-cluster-", label, ".csv"),
                     row.names = FALSE)
  }

  cluster_terms
}
