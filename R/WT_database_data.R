#' Reference database for Z-score computation
#'
#' A reference dataset containing mean and standard deviation of gene expression
#' across multiple immune cell types. This dataset is used by
#' \code{compute_zscore()} to compare observed expression values against
#' reference distributions.
#'
#' @format A data frame with 144,387 rows and 4 columns:
#' \describe{
#'   \item{Gene}{Character. Gene symbol or ENSEMBL ID.}
#'   \item{CellType}{Character. Immune cell type. One of:
#'   "B cells", "CTLs", "DC", "Helper T cells", "Mono/Macro",
#'   "NK cells", "Regulatory T cells".}
#'   \item{ref_mean}{Numeric. Reference mean expression value for the gene in the given cell type.}
#'   \item{ref_sd}{Numeric. Reference standard deviation of expression for the gene in the given cell type.}
#' }
#'
#' @details
#' This dataset provides baseline statistics for immune-related gene expression
#' across major cell types. It is intended for use in enrichment and Z-score
#' analysis workflows within the \pkg{ImmuStable} package.
#'
#' @source Internal reference data curated for ImmuStable package.
#'
#' @examples
#' data(WT_database_data)
#' head(WT_database_data)
#' summary(WT_database_data)
"WT_database_data"
