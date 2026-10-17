# scoranet <img src="man/figures/logo.png" align="right" height="120" alt="scoranet logo" />

**Single-cell differential expression, GO over-representation analysis, and pathway overlap networks.** Developed by Mike Martinez.

 ![Version](https://img.shields.io/badge/version-0.1.0-blue)
 [![docs: pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mikemartinez99.github.io/scoranet/)
 ![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
[![Codecov test coverage](https://codecov.io/gh/mikemartinez99/scoranet/graph/badge.svg)](https://app.codecov.io/gh/mikemartinez99/scoranet)


📖 **Documentation and vignette:** [mikemartinez99.github.io/scoranet](https://mikemartinez99.github.io/scoranet/)

> [!WARNING]
> **🚧 In development:** functions and defaults may still change between versions.

## Workflow

| Step | Functions |
|---|---|
| Differential expression per cell-type cluster | `run_cluster_de()`, `write_de()` |
| GO over-representation analysis | `run_go_ora()`, `write_go()` |
| Pathway overlap networks with Leiden clustering and cluster themes | `run_comparison_networks()` |
| Diagnostics for choosing network settings | `run_network_diagnostics()` |
| Module interpretation | `write_cluster_wordclouds_patchwork()`, `write_cluster_results()` |

See the [interferon-beta PBMC vignette](https://mikemartinez99.github.io/scoranet/articles/ifnb-example.html) for a full walkthrough on real data, and [DIAGNOSTICS.md](DIAGNOSTICS.md) for how to read the diagnostic figures.

## Installation

Install the Bioconductor dependencies first if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "GO.db", "AnnotationDbi", "org.Mm.eg.db"))
```

Then install scoranet from GitHub:

```r
devtools::install_github("https://github.com/mikemartinez99/scoranet/")
```

## Dependencies

**CRAN:** Seurat, cli, dplyr, ggforce, ggplot2, ggraph, ggwordcloud, igraph (>= 2.1.0), patchwork, scales, stringr, withr

**Bioconductor:** clusterProfiler, GO.db, AnnotationDbi, and an organism annotation package (org.Mm.eg.db for mouse or org.Hs.eg.db for human)

## License

MIT
