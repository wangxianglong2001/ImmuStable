
# ImmuStable

**ImmuStable** is an R package designed for the analysis of immune-related single-cell transcriptomic data.  
It provides a **Z-score normalization** framework to evaluate whether gene expression in test samples falls within the normal reference range. Combined with enrichment analysis and visualization functions, ImmuStable helps researchers better understand immune cell expression patterns and functional states.

------------------------------------------------------------------------

##  Features Overview

- **Z-score Calculation**  
  Normalize gene expression in test samples against a reference database (`WT_database_data`).  
  Automatically classify gene expression as *Within WT range*, *Above WT range*, or *Below WT range*.

- **Reference Database**  
  Includes the built-in `WT_database_data`, which contains mean and standard deviation of gene expression across multiple immune cell types.  
  Supported cell types include Helper T cells, CTLs, Regulatory T cells, B cells, Monocytes, Macrophages, Dendritic Cells (DCs), and NK cells.

- **Enrichment Analysis**  
  Perform KEGG and GO (BP/CC/MF) pathway enrichment analysis based on up- and down-regulated genes.  
  Automatically generate statistical results and optionally save them to Excel files.

- **Result Visualization**  
  Provides **dot plots** and **bar plots** to visualize enrichment results.  
  Supports customization of the number of pathways displayed (top N), plot titles, and color schemes, making it easy to interpret pathway significance and gene counts.

------------------------------------------------------------------------

##  Data Source

The `WT_database_data` dataset originates from quality-controlled and normalized single-cell RNA-seq data. It contains:

- **Gene**: Gene symbol or ENSEMBL ID  
- **CellType**: Immune cell type  
- **ref_mean**: Reference mean expression  
- **ref_sd**: Reference standard deviation  

This database serves as the baseline for Z-score calculation, enabling evaluation of whether test sample expression levels deviate from the normal range.

------------------------------------------------------------------------

##  Usage Example

```r
library(ImmuStable)

# Load the reference database
data(WT_database_data)

# Compute Z-scores
zmat <- compute_zscore(seurat_obj, celltype_col = "celltype", wt_ref = WT_database_data)

# Perform enrichment analysis
res <- enrich_by_celltype(zmat, out_prefix = "example")

# Plot enrichment results
plots <- enrich_results_plot(res, celltype = "B cells", direction = "up", which = c("KEGG","BP"))
plots$KEGG$dot
plots$BP$bar
```

``` r
# Enrichment analysis
res <- enrich_by_celltype(zmat, out_prefix = "example")
```

``` r
# Plot enrichment results
plots <- enrich_results_plot(res, celltype = "B cells", direction = "up", which = c("KEGG","BP"))
plots$KEGG$dot
plots$BP$bar
```


