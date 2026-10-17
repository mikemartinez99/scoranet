test_that(".term_words drops filler words and short tokens", {
  words <- .term_words("positive_regulation_of_interferon_beta_production")
  expect_equal(words, c("interferon", "beta", "production"))
})

test_that(".cluster_themes returns keywords, phrase, and hub per cluster", {
  go_term_df <- data.frame(
    go_term = c("collagen fibril organization", "collagen metabolic process",
                "extracellular matrix organization",
                "type I interferon production", "response to interferon beta"),
    cluster = c(1, 1, 1, 2, 2),
    degree  = c(5, 3, 4, 2, 6),
    avg_log2FC = 1,
    stringsAsFactors = FALSE
  )

  themes <- .cluster_themes(go_term_df, n_words = 2)

  expect_named(themes, c("cluster", "n_terms", "ancestor", "ancestor_id",
                         "ancestor_n_covered", "keywords", "phrase", "hub"))
  expect_true(all(is.na(themes$ancestor)))   # no GO IDs in this input
  expect_equal(themes$n_terms, c(3, 2))
  expect_match(themes$keywords[1], "collagen")
  expect_match(themes$keywords[2], "interferon")
  expect_equal(themes$hub, c("collagen fibril organization", "response to interferon beta"))
})

test_that(".rank_cluster_phrases merges near-duplicate term names", {
  sub_df <- data.frame(
    go_term = c("collagen fibril organization", "collagen fibril assembly", "cell migration"),
    degree = c(4, 2, 1), avg_log2FC = 1, stringsAsFactors = FALSE
  )

  merged <- .rank_cluster_phrases(sub_df, merge_similarity = 0.6)
  unmerged <- .rank_cluster_phrases(sub_df, merge_similarity = NULL)

  expect_equal(nrow(unmerged), 3)
  expect_equal(nrow(merged), 2)
  expect_equal(merged$go_term_clean[1], "collagen fibril organization")
  expect_equal(merged$freq[1], 6)
})

test_that("generate_overlap_network returns cluster themes and a themes plot", {
  collection <- list(
    collagen_fibril_organization = paste0("C", 1:6),
    collagen_metabolic_process   = paste0("C", 2:7),
    interferon_beta_production   = paste0("I", 1:6),
    response_to_interferon_beta  = paste0("I", 2:7)
  )
  results_set <- data.frame(avg_log2FC = c(1, 1, -1, -1), row.names = names(collection))

  res <- suppressMessages(generate_overlap_network(
    results_set, collection, unlist(collection),
    outdir = NULL
  ))

  expect_equal(nrow(res$cluster_themes), length(unique(res$go_term_df$cluster)))
  expect_s3_class(plot_network(res, "themes", label_by = "hub"), "ggplot")
  expect_true(any(grepl("collagen", res$cluster_themes$keywords)))
})

test_that(".cluster_ancestor finds a shared GO parent", {
  skip_if_not_installed("GO.db")
  ids <- c(
    "GO:0032481",  # positive regulation of type I interferon production
    "GO:0032479",  # regulation of type I interferon production
    "GO:0032606"   # type I interferon production
  )

  anc <- .cluster_ancestor(ids, min_share = 0.15, max_size = 1000)

  expect_equal(anc$n_covered, 3L)
  expect_match(anc$term[1], "interferon|cytokine")
  expect_length(anc$term, 1)   # one parent covers all three, so no second pick
})

test_that(".cluster_ancestor returns nothing when only generic parents are shared", {
  skip_if_not_installed("GO.db")
  ids <- c("GO:0030198", "GO:0032606")  # ECM organization, type I IFN production

  anc <- .cluster_ancestor(ids, min_share = 0.9, max_size = 100)

  expect_length(anc$id, 0)
  expect_equal(anc$n_covered, 0L)
})

test_that(".theme_text falls back to the hub term without an ancestor", {
  themes <- data.frame(
    cluster = 1:2, n_terms = c(5, 4),
    ancestor = c("cell cycle", NA), ancestor_n_covered = c(4, 0),
    hub = c("mitotic spindle", "response to UV"),
    stringsAsFactors = FALSE
  )
  expect_equal(.theme_text(themes, "ancestor"),
               c("cell cycle (4/5)", "response to UV (hub; no shared GO parent)"))
})
