#' Write Differential Expression Results to CSV Files
#'
#' Saves the output of \code{\link{run_cluster_de}} into cluster-specific CSV
#' files: the filtered DEGs and the full table of tested genes for each
#' cluster.
#'
#' @param de_list The output of \code{\link{run_cluster_de}}: a named list
#'   whose elements each contain \code{$degs} and \code{$tested} data frames.
#' @param outSubDir Character path to the output directory, with or without a
#'   trailing \code{/}. Created recursively if it does not exist.
#'
#' @return \code{invisible(NULL)}. Called for its side effect of writing CSV
#'   files.
#'
#' @section Output files:
#' \itemize{
#'   \item \code{<outSubDir>/<cluster>_DE_Results.csv}: DEGs (skipped when a
#'     cluster has none).
#'   \item \code{<outSubDir>/<cluster>_tested_genes.csv}: all tested genes.
#' }
#'
#' @seealso \code{\link{write_go}} for GO enrichment results.
#'
#' @export
write_de <- function(de_list, outSubDir) {

  outSubDir <- .as_dir(outSubDir)

  for (i in names(de_list)) {

    celltypeName <- .safe_filename(i)
    res <- as.data.frame(de_list[[i]]$degs)
    tested <- as.data.frame(de_list[[i]]$tested)

    if (nrow(tested) > 0) {
      outfile_tested <- paste0(outSubDir, celltypeName, "_tested_genes.csv")
      utils::write.csv(tested, file = outfile_tested, row.names = FALSE)
      cli::cli_alert_success("Wrote {.file {outfile_tested}}")
    }

    if (nrow(res) == 0) {
      cli::cli_alert_warning("{.val {i}}: no DEGs obtained.")
    } else {
      outfile_degs <- paste0(outSubDir, celltypeName, "_DE_Results.csv")
      utils::write.csv(res, file = outfile_degs, row.names = FALSE)
      cli::cli_alert_success("Wrote {.file {outfile_degs}}")
    }

    cli::cli_text("")
  }

  invisible(NULL)
}

#' Write GO Over-Representation Results to CSV Files
#'
#' Saves a named list of GO enrichment results into cluster-specific CSV
#' files. The \code{GeneRatio} column is parsed into numeric
#' numerator/denominator columns, and results can optionally be split by
#' regulation direction.
#'
#' @param go_list A named list of \code{enrichResult} objects, typically the
#'   output of \code{\link{run_go_ora}}. When \code{split_direction = TRUE},
#'   each element must instead be a list with \code{$up} and \code{$down}
#'   \code{enrichResult} objects.
#' @param outSubDir Character path to the output directory, with or without a
#'   trailing \code{/}. Created recursively if it does not exist.
#' @param split_direction Logical. If \code{TRUE}, write separate files for
#'   the \code{$up} and \code{$down} results of each cluster. Must match the
#'   \code{split_direction} used in \code{\link{run_go_ora}}; a mismatch is an
#'   error. Default \code{FALSE}.
#'
#' @return \code{invisible(NULL)}. Called for its side effect of writing CSV
#'   files.
#'
#' @section Output files:
#' \itemize{
#'   \item No split: \code{<outSubDir>/<cluster>_GO_ORA_Results.csv}
#'   \item Split: \code{<outSubDir>/<cluster>_GO_ORA_up_Results.csv} and
#'     \code{..._down_Results.csv}
#' }
#'
#' @seealso \code{\link{write_de}} for differential expression results.
#'
#' @export
write_go <- function(go_list, outSubDir, split_direction = FALSE) {

  outSubDir <- .as_dir(outSubDir)

  .write_go <- function(ego, celltypeName, direction_label = NULL) {
    if (is.null(ego)) return(invisible(NULL))
    res <- as.data.frame(ego)
    if (nrow(res) == 0) {
      label <- if (is.null(direction_label)) "" else paste0(" ", direction_label)
      cli::cli_alert_warning("Skipping {.val {celltypeName}}{label} (no results)")
      return(invisible(NULL))
    }
    tmp <- strsplit(as.character(res$GeneRatio), "/")
    res$geneRatio_num <- as.numeric(sapply(tmp, `[`, 1))
    res$geneRatio_den <- as.numeric(sapply(tmp, `[`, 2))

    suffix <- if (!is.null(direction_label)) {
      paste0("_GO_ORA_", direction_label, "_Results.csv")
    } else {
      "_GO_ORA_Results.csv"
    }

    outfile <- paste0(outSubDir, celltypeName, suffix)
    utils::write.csv(res, file = outfile, row.names = FALSE)
    cli::cli_alert_success("Wrote {.file {outfile}}")
  }

  for (i in names(go_list)) {

    celltypeName <- .safe_filename(i)
    x <- go_list[[i]]

    if (!split_direction) {
      if (is.list(x) && !is.data.frame(x) && any(c("up", "down") %in% names(x))) {
        cli::cli_abort(c(
          "{.code go_list[[\"{i}\"]]} contains {.field up}/{.field down} results from {.code run_go_ora(split_direction = TRUE)}.",
          "i" = "Write them with {.code write_go(..., split_direction = TRUE)}."
        ))
      }
      .write_go(x, celltypeName)
      next
    }

    if (is.null(x)) {
      cli::cli_alert_warning("Skipping {.val {i}} (no results)")
      next
    }

    if (is.data.frame(x) || !is.list(x) || !all(c("up", "down") %in% names(x))) {
      cli::cli_abort(c(
        "{.arg split_direction = TRUE} needs {.code go_list[[\"{i}\"]]} to be a list with {.field up} and {.field down} results.",
        "i" = "{.fn run_go_ora} returns one combined result per cluster; use {.code split_direction = FALSE} for its output."
      ))
    }

    .write_go(x$up,   celltypeName, "up")
    .write_go(x$down, celltypeName, "down")
  }

  invisible(NULL)
}
