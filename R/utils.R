# Replace spaces and slashes so cluster names are safe to use in file names.
.safe_filename <- function(x) {
  gsub("[ /]", "_", x)
}

# Report the Seurat assay version of one assay in a Seurat object: "v5" for
# Seurat 5 assays (class Assay5), otherwise "v3" (the Assay class used by
# Seurat 3 and 4).
.seurat_assay_version <- function(obj, assay = "RNA") {
  if (inherits(obj[[assay]], "Assay5")) "v5" else "v3"
}

# The 40 distinct colours of dittoSeq::dittoColors() (dittoSeq, MIT licence;
# Bunis et al. 2021, Bioinformatics), copied so the package does not depend on
# dittoSeq. Starts with the colour-blind-friendly Okabe-Ito palette, followed by
# darker and lighter variants of it.
.ditto_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700",
  "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF",
  "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
  "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45",
  "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D"
)

# Evaluate `expr` with a fixed random seed without changing the session's random
# state. `seed = NULL` evaluates `expr` without seeding.
.with_optional_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  withr::with_seed(seed, expr)
}

# `x` unless it is NULL, otherwise `y` (base R's `%||%` needs R >= 4.4).
.or_default <- function(x, y) {
  if (is.null(x)) y else x
}

# Output directory built from one or more path parts, created if missing, and
# returned with exactly one trailing "/" so file names can be appended with
# paste0(). Parts may be given with or without trailing slashes.
.as_dir <- function(...) {
  parts <- c(...)
  parts[-1] <- sub("^/+", "", parts[-1])
  parts <- sub("/+$", "", parts)
  path <- do.call(file.path, as.list(parts[nzchar(parts) | seq_along(parts) == 1]))
  if (!nzchar(path)) path <- "/"
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  if (grepl("/$", path)) path else paste0(path, "/")
}
