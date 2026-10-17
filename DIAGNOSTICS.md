# Interpreting the scoranet diagnostic figures

This guide explains what each diagnostic figure shows, what a healthy result
looks like, and which settings to change when it doesn't look healthy. All
figures come from `run_network_diagnostics()`, or from the individual
functions it calls.

```r
dx <- run_network_diagnostics(net_list, de_list, go_list, outdir = "results/diagnostics/")
```

## Contents

1. [Recommended order](#recommended-order)
2. [Comparison overview](#1-comparison-overview)
3. [Overlap threshold diagnostic](#2-overlap-threshold-diagnostic)
4. [Leiden resolution sweep](#3-leiden-resolution-sweep)
5. [Driver genes per cluster](#4-driver-genes-per-cluster)
6. [The network plots themselves](#5-the-network-plots)
7. [Quick reference](#quick-reference)

---

## Recommended order

The settings depend on each other, so check them in this order:

1. **Overview.** Which cell types have enough signal to bother with?
2. **Overlap threshold.** It decides which terms are in the network at all.
3. **Leiden resolution.** It decides how the remaining terms are grouped.
   Choose it *after* the threshold, because changing the threshold changes the
   network the resolution acts on.
4. **Driver genes.** Read these last, to interpret the clusters you settled on.

---

## 1. Comparison overview

`plot_comparison_overview(de_list, go_list, net_list)`

**What it shows.** For each cell type, three panels of bar counts:

| Panel | Counts |
|---|---|
| DEGs | Genes passing `padj_thresh` and `logfc_thresh`, split into up in the numerator and up in the denominator |
| Enriched GO terms | Significant terms from `run_go_ora()`; split by direction if ORA was run with `split_direction = TRUE`, otherwise one grey "Combined" bar |
| Network clusters | Leiden clusters in each cell type's network, split by direction for split networks |

**How to read it.**

- **Where is the biology?** Cell types with many DEGs *and* many enriched terms
  are where the comparison has the clearest effect.
- **Direction imbalance.** A cell type whose DEGs are almost all one colour
  responds mostly in one direction. Consider running ORA with
  `split_direction = TRUE` for it, so the minority direction isn't drowned out.
- **Many DEGs but few terms.** The DEGs may be biologically scattered, may be
  poorly annotated in GO, or the universe may be too narrow. Check the
  `run_go_ora()` QC table.
- **Few DEGs but many terms.** A small set of genes falls into many related GO
  terms (often one large family, such as ribosomal or histone genes). Check the
  driver-gene plot before trusting those terms.
- **Terms but no network clusters.** The terms don't share enough genes to
  connect at the current overlap threshold (see section 2).

**Caveats.**

- DEG counts depend strongly on **cells per group**. A cell type with 2,000
  cells will find more DEGs than one with 50, even for the same effect size.
  Don't compare raw counts across cell types as if they measured effect
  strength.
- The number of network clusters depends on the resolution. More clusters is
  not "more biology".

---

## 2. Overlap threshold diagnostic

`plot_overlap_threshold(net_result)`

**What it shows.** Two panels.

**Top: histogram of pairwise Jaccard overlaps.** For every pair of enriched
terms that share at least one DE gene, the Jaccard overlap is
`shared genes / total distinct genes`. The dashed line is the
`overlap_threshold` used to build the network: pairs at or above it become
edges. The subtitle reports how many of all term pairs share any gene.

**Bottom: what each threshold would keep.** Across a range of thresholds:

- **Kept** (green): terms with at least one edge.
- **Dropped (isolates)** (orange): terms with no edge. These are **removed from
  the network entirely**: they don't appear in the network plots,
  `go_term_df`, driver genes, or wordclouds.
- **Edges** (right panel): number of connections.

**How to read it.**

- **The threshold should sit before the orange line rises.** If many terms are
  dropped at your threshold, you are discarding enriched biology. Lower the
  threshold, or accept that those terms are unrelated to the rest.
- **Watch the edge count.** Too low a threshold connects almost everything into
  one hairball, and Leiden then has little structure to work with. Too high a
  threshold leaves a few disconnected fragments.
- **A good threshold** keeps nearly all terms while cutting the edges well below
  their maximum. In the histogram it usually falls just above the large mass of
  very small overlaps (pairs sharing one gene out of many).
- **Steps in the curves are normal.** With small gene sets, Jaccard values take
  a limited set of values (1/5, 1/4, 1/3, ...), so the histogram has gaps and
  the curves move in jumps.

**Typical values.** 0.1 to 0.3 for most data sets. The default is 0.1.

**What to change.** `overlap_threshold` in `run_comparison_networks()`.

---

## 3. Leiden resolution sweep

`plot_leiden_sweep(net_result)`

**What it shows.** The network is re-clustered at each resolution, with 10
different seeds at each. Three panels share the x axis (resolution). The dashed
line marks the resolution the network was built with.

| Panel | Meaning | Line / band |
|---|---|---|
| **Clusters** | Number of Leiden clusters | Mean / range across seeds |
| **Modularity** | How well separated the clusters are: the share of edge weight inside clusters, compared with what random wiring would give. Always measured at resolution 1, so values are comparable across the x axis. | Mean / range |
| **Stability (ARI)** | Agreement between the partitions from different seeds (adjusted Rand index; 1 = identical) | Mean / range across all seed pairs |

**How to read it.** Look for a **plateau**: a stretch of resolutions where

- the **cluster count is flat**,
- **modularity is high** and roughly flat, and
- **stability is close to 1** with a narrow band.

A resolution inside such a plateau gives clusters that reflect real structure
in the network rather than the particular seed or the exact resolution you
chose. If the dashed line sits on a steep part of the cluster-count curve, small
changes in resolution would change your clusters: move it onto a plateau.

**Interpreting the panels together.**

- **Rising cluster count, falling modularity.** Clusters are being split past
  the network's natural structure. Go lower.
- **Stability dips with a wide band.** At that resolution, different seeds give
  different clusters: the partition is ambiguous. Avoid those resolutions.
- **Stability of 1 at the extremes is not meaningful on its own.** At very low
  resolution everything merges into one or a few clusters, and at very high
  resolution clusters approach single terms. Both are perfectly reproducible
  but uninformative. Always read stability together with the cluster count and
  modularity.
- **Modularity peak.** The resolution with the highest modularity is a
  reasonable default if there is no clear plateau.

**Modularity values.** Above about 0.3 suggests meaningful community
structure; GO overlap networks are often 0.5 to 0.8, because related terms form
tight groups.

**What to change.** `leiden_resolution` in `run_comparison_networks()`.

**Notes.**

- The default objective is **modularity**, where resolutions around 1 are
  standard (the default is 1). With `objective_function = "CPM"`, resolution
  means something different (an edge-density cutoff) and needs much smaller
  values, such as 0.01 to 0.1. Pass `resolutions = seq(0.005, 0.1, by = 0.005)`
  to the sweep in that case.
- Clustering is seeded (`seed = 1234`), so rerunning gives the same clusters.
  The sweep varies the seed deliberately, to test how much the result would
  change without that fixed seed.
- Networks with very few nodes can't be swept meaningfully; the wrapper skips
  them with a message.

---

## 4. Driver genes per cluster

`plot_driver_genes(net_result, de_df)`

**What it shows.** For each Leiden cluster (columns), the `top_n` DE genes
(rows) that appear in the most of that cluster's GO terms.

- **Dot size** is the share of the cluster's terms that contain the gene. 100%
  means the gene is in every term of the cluster.
- **Dot colour** is the gene's log2 fold change (red = up in the numerator,
  blue = up in the denominator). Grey if no `de_df` was supplied.

**How to read it.**

- **Name the module from its genes, not only its terms.** GO term names are
  generic ("regulation of cellular process"); the driver genes say what is
  actually changing.
- **Large dots for one or two genes** mean the cluster exists because a few
  genes are annotated to many terms. The cluster may represent one gene's
  annotation breadth rather than a coordinated program. Treat it with care.
- **Many medium dots** mean the cluster is supported by a broad set of genes
  shared across its terms: a more convincing module.
- **Mixed colours within a cluster** mean its genes move in both directions.
  This matches a grey (mixed) node colour in the network plot and is common for
  regulatory programs.
- **The same gene driving several clusters** (it appears in more than one
  column) links those modules. They may be one program that Leiden split; try a
  lower resolution.
- **Watch for housekeeping families.** Ribosomal (Rpl/Rps), mitochondrial (mt-),
  and histone genes can dominate clusters through sheer annotation volume, and
  can reflect technical differences between samples.

**Rows** are ordered by the first cluster each gene drives, so the plot reads
roughly along a diagonal.

**What to change.** `top_n` for more or fewer genes per cluster.

---

## 5. The network plots

Each network gets three plots with the same layout (`$plots` in each network
result, also saved as PNGs).

| Plot | Node fill | Use it to |
|---|---|---|
| `plain` | Direction: red = up in numerator, blue = up in denominator, grey = mixed | See which regions of the network are up or down |
| `leiden` | Leiden cluster (legend on the left) | See the modules; match them to the driver-gene plot and wordclouds |
| `zscore` | Continuous direction z-score, blue to white to red | See *how strongly* each term leans one way |

In `plain` and `leiden`, node shape also shows direction (circle = up,
square = down, diamond = mixed). In all three, edge width shows Jaccard overlap
and node size shows the number of DE genes in the term.

**The direction z-score.** For each term,
`z = (n_up - n_down) / sqrt(n_up + n_down)`, using the term's DE genes. It
combines *how consistently* the genes move with *how many* genes support it:

- A term with 9 up and 1 down genes: z = 2.5 (clearly up).
- A term with 3 up and 2 down: z = 0.45 (mixed).
- A term with a single up gene: z = 1.0 (up, but weakly supported).

Terms are called up or down when |z| is at least `zscore_threshold` (default 1)
and mixed otherwise. Raise it (for example to 2) for stricter direction calls.

The output tables also keep the **mean** and **median** log2 fold change for
each term. The mean can be pulled in the opposite direction by one gene with an
extreme fold change; the median and z-score are more robust.

**With `split_direction = TRUE`** every term in the `_up` network is up and every
term in the `_down` network is down, so the `plain` and `zscore` plots are a
single colour. Use the `leiden` plot and driver genes for those networks.

---

## Quick reference

| Symptom | Likely cause | Try |
|---|---|---|
| Many terms dropped as isolates | Overlap threshold too high | Lower `overlap_threshold` |
| One giant hairball cluster | Threshold too low, or resolution too low | Raise `overlap_threshold`, or raise `leiden_resolution` |
| Dozens of tiny clusters | Resolution too high | Lower `leiden_resolution` (check the sweep) |
| Clusters change between runs with different seeds | Resolution in an unstable region | Move the resolution onto a stability plateau |
| A cluster driven by one gene | Annotation breadth of that gene | Interpret with care; check its terms |
| Driver genes are ribosomal or mitochondrial | Housekeeping or technical signal | Check QC; consider whether the module is biological |
| Many DEGs, few terms | Scattered biology, or a narrow universe | Check the ORA QC table; try `use_tested_universe` both ways |
| Most terms are mixed | Up and down genes share annotations | Run ORA with `split_direction = TRUE` |
| Terms enriched, but no network | No two terms overlap at the threshold | Lower `overlap_threshold` |
