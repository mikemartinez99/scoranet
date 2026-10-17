# scoranet <img src="img/logo.png" align="right" height="120" alt="scoranet logo" />

**Single-cell differential expression, GO over-representation analysis, and pathway overlap networks.** Developed by the Mike Martinez

 ![Version](https://img.shields.io/badge/version-0.1.0-blue)
 [![pkgdown documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://dartmouth-data-analytics-core.github.io/GDSCtools/index.html)
 ![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
 [![codecov](https://codecov.io/gh/dartmouth-data-analytics-core/GDSCtools/branch/main/graph/badge.svg)](https://codecov.io/gh/dartmouth-data-analytics-core/GDSCtools)

> [!WARNING]
> **🚧 In development:**


## Installation

```r
# Install from github
devtools::install_github("https://github.com/mikemartinez99/scoranet/")
```

## Dependencies

Install Bioconductor dependencies first if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))
```



## License

MIT
