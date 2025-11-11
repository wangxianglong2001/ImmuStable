# R/immustable-internal-globals.R
# Internal: declare global variable names used in NSE
# This file prevents "no visible binding" notes for known data-frame column names.
# It is internal only and does not produce user documentation.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # mapping / annotation columns
    "SYMBOL", "ENTREZID", "ENSEMBL", "ENTREZID_bm", "ENTREZID.db",
    "ENSEMBL.bm", "ENSEMBL.db", "ENTREZID_final", "ENSEMBL_final",
    # enrichment / plotting columns
    "EnrichmentVal", "QVal", "Term", "ListHits", "LogQ",
    # zscore / reference columns
    "Gene", "CellType", "ref_mean", "ref_sd", "Zscore", "Status",
    # dplyr .data pronoun
    ".data"
  ))
}
