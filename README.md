# scoranet

**Single-cell differential expression, GO over-representation analysis, and pathway overlap networks.** Developed by the Dartmouth Genomic Data Science Core.

 ![Version](https://img.shields.io/badge/version-0.1.0-blue)
 [![pkgdown documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://dartmouth-data-analytics-core.github.io/GDSCtools/index.html)
 ![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
 [![codecov](https://codecov.io/gh/dartmouth-data-analytics-core/GDSCtools/branch/main/graph/badge.svg)](https://codecov.io/gh/dartmouth-data-analytics-core/GDSCtools)

> [!WARNING]
> **🚧 In development:**  
> This package currently supports basic single-cell quality control, differential expression analysis, ORA analysis, and pathway network analysis. More functionalities and utilities will be added, stay tuned!

## Installation

```r
# Install from github
devtools::install_github("https://github.com/Dartmouth-Data-Analytics-Core/GDSCtools")
```

## Dependencies

**CRAN:** igraph, dplyr, ggplot2, ggwordcloud, patchwork, scales, stringr, purrr, ggExtra

**Bioconductor:** Seurat, clusterProfiler, org.Mm.eg.db (or org.Hs.eg.db)

Install Bioconductor dependencies first if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))
```

## Vignettes

- [Pathway Network Analysis](https://dartmouth-data-analytics-core.github.io/GDSCtools/articles/network-analysis.html)

## License

MIT
