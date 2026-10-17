#' Generate Patchwork Word Clouds for All Leiden Clusters
#'
#' For each cell type in a network result list, this function creates a
#' multi-panel word cloud figure where each panel represents one Leiden
#' cluster. GO term labels are sized by their network degree (optionally
#' weighted by absolute log2 fold-change) and semantically similar terms are
#' merged via hierarchical clustering on Jaccard distance of constituent
#' words. The combined figure is saved as a PNG file.
#'
#' @param net_list A named list of network result lists, one per cell type.
#'   Each element must be a list containing at least \code{go_term_df} (a
#'   data frame with columns \code{go_term}, \code{cluster}, \code{degree},
#'   and \code{avg_log2FC}), as returned by
#'   \code{\link{generate_overlap_network}}.
#' @param outputDir Character path to the root output directory. Word cloud
#'   PNGs are saved to \code{<outputDir>/<comparison>/networks/}.
#' @param comparison Character string identifying the comparison (e.g.,
#'   \code{"CD206_Car_vs_Bleo"}). Used to construct the output sub-directory.
#' @param min_terms Integer. Minimum number of (post-merge) GO term phrases
#'   required to generate a word cloud for a cluster. Clusters with fewer
#'   terms are skipped. Default \code{5}.
#' @param max_terms Integer. Maximum number of GO term phrases to include per
#'   cluster word cloud (after sorting by weight). Default \code{50}.
#' @param use_logfc_weight Logical. If \code{TRUE}, term weights are scaled by
#'   \code{1 + abs(avg_log2FC)}, emphasizing terms backed by larger fold-
#'   changes. Default \code{FALSE}.
#' @param ncol Integer. Number of columns in the patchwork grid layout.
#'   Default \code{3}.
#' @param merge_similarity Numeric \code{[0, 1]} or \code{NULL}. Jaccard-
#'   distance threshold for merging semantically similar GO term labels. A
#'   value of \code{0.6} means terms sharing >= 40\% of their words are merged
#'   (the most frequent label is kept and frequencies are summed). Set to
#'   \code{NULL} to disable merging. Default \code{0.6}.
#'
#' @return Called for its side effect (writing PNG files). Returns
#'   \code{invisible(NULL)}.
#'
#' @section Output files:
#' One PNG per cell type:
#' \code{<outputDir>/<comparison>/networks/<celltype>_all_clusters_wordclouds.png}
#'
#' @seealso \code{\link{generate_overlap_network}}
#'
#' @importFrom dplyr filter mutate group_by summarise arrange desc slice_head
#'   select
#' @importFrom stringr str_replace_all str_trim str_wrap
#' @importFrom ggplot2 ggplot aes scale_size_area scale_color_manual
#'   theme_minimal theme ggtitle ggsave
#' @importFrom ggwordcloud geom_text_wordcloud
#' @importFrom scales hue_pal
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom stats as.dist cutree hclust
#'
#' @export
write_cluster_wordclouds_patchwork <- function(
    net_list,
    outputDir,
    comparison = "CD206_Car_vs_Bleo",
    min_terms = 5,
    max_terms = 50,
    use_logfc_weight = FALSE,
    ncol = 3,
    merge_similarity = 0.6
) {

  for (celltype in names(net_list)) {

    df <- net_list[[celltype]]$go_term_df

    if (is.null(df) || nrow(df) == 0) {
      cli::cli_alert_warning("Skipping {.val {celltype}} (no GO terms)")
      next
    }

    celltype_safe <- .safe_filename(celltype)

    plot_list <- list()

    for (cluster_id in sort(unique(df$cluster))) {

      sub_df <- df |>
        dplyr::filter(.data$cluster == cluster_id)

      if (nrow(sub_df) == 0) next

      sub_df <- sub_df |>
        dplyr::mutate(
          go_term_clean = stringr::str_replace_all(.data$go_term, "_", " "),
          go_term_clean = stringr::str_trim(.data$go_term_clean)
        )

      sub_df <- sub_df |>
        dplyr::mutate(
          term_weight = pmax(.data$degree, 1),
          term_weight = if (use_logfc_weight) {
            .data$term_weight * (1 + abs(.data$avg_log2FC))
          } else {
            .data$term_weight
          }
        )

      phrase_freq <- sub_df |>
        dplyr::group_by(.data$go_term_clean) |>
        dplyr::summarise(freq = sum(.data$term_weight), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(.data$freq)) |>
        dplyr::slice_head(n = max_terms)

      if (!is.null(merge_similarity) && nrow(phrase_freq) > 1) {
        tokens <- strsplit(tolower(phrase_freq$go_term_clean), "\\s+")
        n <- length(tokens)
        jdist <- matrix(0, n, n)
        for (i in seq_len(n - 1)) {
          for (j in (i + 1):n) {
            inter <- length(intersect(tokens[[i]], tokens[[j]]))
            union_val <- length(union(tokens[[i]], tokens[[j]]))
            jdist[i, j] <- jdist[j, i] <- if (union_val == 0) 1 else 1 - inter / union_val
          }
        }
        hc <- stats::hclust(stats::as.dist(jdist), method = "average")
        hc$height <- cummax(hc$height)
        phrase_freq$group <- stats::cutree(hc, h = merge_similarity)
        phrase_freq <- phrase_freq |>
          dplyr::group_by(.data$group) |>
          dplyr::summarise(
            go_term_clean = .data$go_term_clean[which.max(.data$freq)],
            freq = sum(.data$freq),
            .groups = "drop"
          ) |>
          dplyr::select(-"group") |>
          dplyr::arrange(dplyr::desc(.data$freq))
      }

      if (nrow(phrase_freq) < min_terms) {
        cli::cli_alert_warning("Skipping {.val {celltype}} cluster {cluster_id} (too few phrases)")
        next
      }

      p <- ggplot2::ggplot(
        phrase_freq,
        ggplot2::aes(
          label = stringr::str_wrap(.data$go_term_clean, 30),
          size = sqrt(.data$freq),
          color = .data$go_term_clean
        )
      ) +
        ggwordcloud::geom_text_wordcloud(area_corr = TRUE) +
        ggplot2::scale_size_area(max_size = 10) +
        ggplot2::scale_color_manual(
          values = scales::hue_pal()(nrow(phrase_freq))
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(paste0("Cluster ", cluster_id))

      plot_list[[paste0("cluster_", cluster_id)]] <- p
    }

    if (length(plot_list) == 0) next

    combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) +
      patchwork::plot_annotation(title = celltype)

    out_file <- paste0(
      outputDir, comparison, "/networks/",
      celltype_safe, "_all_clusters_wordclouds.png"
    )

    ggplot2::ggsave(out_file,
                    plot = combined_plot,
                    width = 4 * ncol,
                    height = 4 * ceiling(length(plot_list) / ncol))
  }

  invisible(NULL)
}
