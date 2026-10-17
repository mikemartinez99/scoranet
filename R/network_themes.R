# Cluster themes: short text labels describing each Leiden cluster, shared by
# the network themes plot and the cluster wordclouds.

# Words that carry little meaning in GO term names, dropped from keywords.
.go_stopwords <- c(
  "the", "and", "for", "via", "from", "with", "into", "onto", "involved",
  "regulation", "regulated", "regulating", "positive", "negative", "process",
  "processes", "cellular", "cell", "cells", "response", "biological",
  "activity", "pathway", "mediated", "dependent", "establishment",
  "maintenance", "mechanism", "system", "type", "level", "levels", "other",
  "organism", "multicellular", "specific", "based", "involved", "part"
)

# Informative lower-case words in a GO term name (network names use "_" for "-").
.term_words <- function(term) {
  words <- strsplit(tolower(gsub("_", " ", term)), "[^a-z0-9]+")[[1]]
  words[nchar(words) >= 3 & !words %in% .go_stopwords & !grepl("^[0-9]+$", words)]
}

# Rank one cluster's GO term names for a wordcloud: weight each term by its
# network degree (optionally scaled by 1 + |avg_log2FC|), keep the top
# `max_terms`, and merge terms whose names share most of their words.
# Returns a data frame with `go_term_clean` and `freq`, highest first.
.rank_cluster_phrases <- function(sub_df, max_terms = 50, use_logfc_weight = FALSE,
                                  merge_similarity = 0.6) {
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

  phrase_freq
}

# --- GO hierarchy ------------------------------------------------------------

# Root terms of the three GO ontologies; never useful as a label.
.go_roots <- c("GO:0008150", "GO:0003674", "GO:0005575", "all")

# Look up GO ancestors or offspring for a set of GO IDs, using the GO.db table
# for each ID's ontology. Returns a named list (unknown IDs are dropped).
.go_relatives <- function(ids, which = c("ANCESTOR", "OFFSPRING")) {
  which <- match.arg(which)
  ids <- unique(ids[!is.na(ids)])
  ids <- ids[ids %in% AnnotationDbi::keys(GO.db::GOTERM)]
  if (length(ids) == 0) return(list())
  ont <- AnnotationDbi::Ontology(GO.db::GOTERM[ids])
  tables <- list(
    ANCESTOR  = list(BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR),
    OFFSPRING = list(BP = GO.db::GOBPOFFSPRING, MF = GO.db::GOMFOFFSPRING, CC = GO.db::GOCCOFFSPRING)
  )[[which]]
  out <- list()
  for (o in c("BP", "MF", "CC")) {
    sel <- names(ont)[ont == o]
    if (length(sel) > 0) {
      out <- c(out, AnnotationDbi::mget(sel, tables[[o]], ifnotfound = NA))
    }
  }
  lapply(out, function(x) x[!is.na(x)])
}

# Pick up to `max_ancestors` GO terms that together describe a cluster's terms.
# Each pick must newly cover at least `min_share` of the cluster's terms, must
# not be an ontology root, and must have at most `max_size` descendants. Picks are greedy: each maximises
#   (terms it covers that earlier picks did not) x information content,
# where information content = log(ontology size / (descendants + 1)), so a
# specific term covering fewer terms can beat a generic one covering more.
# Returns the picks' IDs and names, the terms each newly covers, and how many
# of the cluster's terms are covered in total. `id` is empty if nothing
# qualifies.
.cluster_ancestor <- function(go_ids, min_share = 0.15, max_size = 1000,
                              max_ancestors = 2) {
  none <- list(id = character(0), term = character(0), n_new = integer(0),
               n_covered = 0L)
  go_ids <- unique(go_ids[!is.na(go_ids)])
  if (length(go_ids) == 0) return(none)

  ancestors <- .go_relatives(go_ids, "ANCESTOR")
  covered_by <- lapply(go_ids, function(id) unique(c(id, ancestors[[id]])))
  candidates <- setdiff(unique(unlist(covered_by)), .go_roots)
  if (length(candidates) == 0) return(none)

  # covers[i, j]: candidate j is term i or one of its ancestors
  covers <- vapply(candidates, function(cand) {
    vapply(covered_by, function(set) cand %in% set, logical(1))
  }, logical(length(go_ids)))
  covers <- matrix(covers, nrow = length(go_ids), dimnames = list(go_ids, candidates))

  offspring <- .go_relatives(candidates, "OFFSPRING")
  size <- vapply(candidates, function(cand) length(offspring[[cand]]), numeric(1))
  ont <- AnnotationDbi::Ontology(GO.db::GOTERM[candidates])
  root_size <- c(BP = length(GO.db::GOBPOFFSPRING[["GO:0008150"]]),
                 MF = length(GO.db::GOMFOFFSPRING[["GO:0003674"]]),
                 CC = length(GO.db::GOCCOFFSPRING[["GO:0005575"]]))
  info <- log(root_size[ont[candidates]] / (size + 1))

  eligible <- size <= max_size
  uncovered <- rep(TRUE, length(go_ids))
  picks <- character(0)
  n_new <- integer(0)
  for (k in seq_len(max_ancestors)) {
    gain <- colSums(covers[uncovered, , drop = FALSE])
    ok <- eligible & gain / length(go_ids) >= min_share & !(candidates %in% picks)
    score <- ifelse(ok, gain * info, -Inf)
    if (all(!is.finite(score))) break
    best <- candidates[which.max(score)]
    picks <- c(picks, best)
    n_new <- c(n_new, as.integer(gain[[best]]))
    uncovered <- uncovered & !covers[, best]
    if (!any(uncovered)) break
  }
  if (length(picks) == 0) return(none)

  list(
    id = picks,
    term = unname(AnnotationDbi::Term(GO.db::GOTERM[picks])),
    n_new = n_new,
    n_covered = as.integer(sum(!uncovered))
  )
}

# One row per Leiden cluster with three candidate labels:
#   ancestor: up to `max_ancestors` shared GO ancestors (see
#             .cluster_ancestor), when go_term_df has GO IDs;
#   keywords: the `n_words` most distinctive words across the cluster's term
#             names (degree-weighted counts, down-weighting words that appear
#             in many clusters, like TF-IDF);
#   phrase:   the top merged term name, as used by the wordclouds;
#   hub:      the cluster's most connected term.
.cluster_themes <- function(go_term_df, n_words = 3, ancestor_min_share = 0.15,
                            ancestor_max_size = 1000, max_ancestors = 2) {
  has_ids <- !is.null(go_term_df$go_id) && any(!is.na(go_term_df$go_id))
  clusters <- sort(unique(go_term_df$cluster))
  k <- length(clusters)

  word_rows <- do.call(rbind, lapply(seq_len(nrow(go_term_df)), function(i) {
    words <- unique(.term_words(go_term_df$go_term[i]))
    if (length(words) == 0) return(NULL)
    data.frame(cluster = go_term_df$cluster[i], word = words,
               weight = max(go_term_df$degree[i], 1), stringsAsFactors = FALSE)
  }))

  keywords <- stats::setNames(rep("", k), clusters)
  if (!is.null(word_rows)) {
    tf <- stats::aggregate(weight ~ cluster + word, data = word_rows, FUN = sum)
    n_clusters_with_word <- table(tf$word)
    tf$score <- tf$weight * log(1 + k / as.numeric(n_clusters_with_word[tf$word]))
    for (cl in clusters) {
      sub <- tf[tf$cluster == cl, ]
      sub <- sub[order(-sub$score, sub$word), ]
      keywords[as.character(cl)] <- paste(utils::head(sub$word, n_words), collapse = ", ")
    }
  }

  do.call(rbind, lapply(clusters, function(cl) {
    sub <- go_term_df[go_term_df$cluster == cl, ]
    phrase <- .rank_cluster_phrases(sub, max_terms = 50)$go_term_clean[1]
    hub <- sub$go_term[which.max(sub$degree)]
    anc <- if (has_ids) {
      .cluster_ancestor(sub$go_id, ancestor_min_share, ancestor_max_size, max_ancestors)
    } else {
      list(id = character(0), term = character(0), n_covered = 0L)
    }
    data.frame(
      cluster  = cl,
      n_terms  = nrow(sub),
      ancestor = if (length(anc$term) > 0) paste(anc$term, collapse = " / ") else NA_character_,
      ancestor_id = if (length(anc$id) > 0) paste(anc$id, collapse = ";") else NA_character_,
      ancestor_n_covered = anc$n_covered,
      keywords = unname(keywords[as.character(cl)]),
      phrase   = phrase,
      hub      = trimws(gsub("_", " ", hub)),
      stringsAsFactors = FALSE
    )
  }))
}

# Layout that keeps each Leiden cluster together: edges within a cluster pull
# much harder than edges between clusters.
.grouped_layout <- function(network, seed, between_weight = 0.01) {
  ends <- igraph::ends(network, igraph::E(network), names = FALSE)
  cl <- igraph::V(network)$cluster
  same <- cl[ends[, 1]] == cl[ends[, 2]]
  w <- igraph::E(network)$weight * ifelse(same, 1, between_weight)
  .with_optional_seed(seed, igraph::layout_with_fr(network, weights = w))
}

# Network plot annotated with cluster themes: a shaded hull around each Leiden
# cluster, labelled with its theme in a callout.
.themes_plot <- function(network, layout, reg_labels, title, cluster_palette,
                         themes, label_by = "keywords", min_terms = 3) {
  reg_breaks <- c("up", "down", "mixed")
  reg_shapes <- c(up = 21, down = 22, mixed = 23)
  reg_names  <- stats::setNames(reg_labels, reg_breaks)

  theme_text <- .theme_text(themes, label_by)

  idx <- match(igraph::V(network)$cluster, themes$cluster)
  labelled <- themes$n_terms[idx] >= min_terms
  igraph::V(network)$cluster_f <- factor(igraph::V(network)$cluster, levels = themes$cluster)
  # Callouts put each ancestor on its own line to keep them compact
  callout_text <- gsub(" / ", "\n", theme_text, fixed = TRUE)
  igraph::V(network)$hull_label <- ifelse(
    labelled, paste0(themes$cluster[idx], ": ", callout_text[idx]), NA_character_
  )

  fill_values <- stats::setNames(cluster_palette[seq_len(nrow(themes))], themes$cluster)
  legend_labels <- stats::setNames(paste("Cluster", themes$cluster), themes$cluster)

  # Fade edges between clusters so the clusters stand out
  ends <- igraph::ends(network, igraph::E(network), names = FALSE)
  cl <- igraph::V(network)$cluster
  igraph::E(network)$within <- cl[ends[, 1]] == cl[ends[, 2]]

  hull_args <- list(
    concavity = 4, expand = ggplot2::unit(3, "mm"), radius = ggplot2::unit(3, "mm")
  )

  ggraph::ggraph(network, layout = "manual", x = layout[, 1], y = layout[, 2]) +
    # Shaded hulls underneath everything
    do.call(ggforce::geom_mark_hull, c(list(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, group = .data$cluster_f,
                             fill = .data$cluster_f),
      alpha = 0.15, colour = "grey55", linewidth = 0.3, show.legend = FALSE
    ), hull_args)) +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_width = .data$weight, edge_alpha = .data$within),
      colour = "grey60"
    ) +
    ggraph::scale_edge_width(range = c(0.2, 1.2), guide = "none") +
    ggraph::scale_edge_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.06), guide = "none") +
    ggraph::geom_node_point(
      ggplot2::aes(fill = .data$cluster_f, shape = .data$regulation, size = .data$size),
      colour = "grey25", stroke = 0.4
    ) +
    # Invisible hulls on top, only to draw the callout labels above the edges
    do.call(ggforce::geom_mark_hull, c(list(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, group = .data$cluster_f,
                             label = .data$hull_label),
      fill = NA, colour = NA,
      label.fontsize = 8, label.fontface = "plain", label.fill = "white",
      label.buffer = ggplot2::unit(6, "mm"), label.minwidth = ggplot2::unit(25, "mm"),
      con.colour = "grey40", con.cap = 0, show.legend = FALSE
    ), hull_args)) +
    ggplot2::scale_fill_manual(
      name = "Cluster", values = fill_values, labels = legend_labels,
      guide = ggplot2::guide_legend(position = "left",
                                    override.aes = list(shape = 21, size = 4, alpha = 1))
    ) +
    ggplot2::scale_shape_manual(
      name = NULL, values = reg_shapes, breaks = reg_breaks, labels = reg_names,
      guide = ggplot2::guide_legend(override.aes = list(size = 4, fill = "grey60"))
    ) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right",
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    )
}

# Label text for each cluster. For "ancestor", the shared GO ancestors and how
# many of the cluster's terms they cover, falling back to the hub term when no
# ancestor qualifies (often a sign the cluster mixes unrelated terms).
.theme_text <- function(themes, label_by) {
  if (label_by == "ancestor") {
    has_anc <- !is.na(themes$ancestor)
    text <- ifelse(
      has_anc,
      paste0(themes$ancestor, " (", themes$ancestor_n_covered, "/", themes$n_terms, ")"),
      paste0(themes$hub, " (hub; no shared GO parent)")
    )
  } else {
    text <- themes[[label_by]]
  }
  text[is.na(text) | text == ""] <- "(no label)"
  text
}
