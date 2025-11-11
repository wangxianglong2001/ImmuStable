#' Compute Z-score for Seurat object against reference database
#'
#' This function calculates average gene expression per cell type in a Seurat object,
#' compares it with a reference database, and computes Z-scores. Optionally, results
#' can be saved to Excel or CSV files.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param celltype_col Character string. Column name in `seurat_obj@meta.data` that
#'   specifies cell type annotations.
#' @param assay Character string. Assay to use (default: "RNA").
#' @param use_layer Optional. Layer name in the assay to extract expression data from.
#'   If NULL, the "data" slot is used.
#' @param wt_ref A data.frame reference database with required columns:
#'   \code{Gene}, \code{CellType}, \code{ref_mean}, \code{ref_sd}.
#'   If NULL, the function will look for `WT_database_data` in the global environment.
#' @param save_path Optional. File path to save results. If provided, results will be
#'   written to Excel (if \pkg{openxlsx} is available) or CSV files.
#' @importFrom magrittr %>%
#' @return A numeric matrix of Z-scores (genes x cell types).
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts expression matrix from the Seurat object.
#'   \item Computes mean expression per gene per cell type.
#'   \item Joins with reference database and calculates Z-scores.
#'   \item Classifies results into "Above", "Below", or "Within" database range.
#'   \item Optionally saves results to file.
#' }
#'
#' @examples
#' \dontrun{
#' zmat <- compute_zscore(seurat_obj, celltype_col = "celltype",
#'                        wt_ref = WT_database_data, save_path = "output")
#' }
#'
#' @export
compute_zscore <- function(seurat_obj,
                           celltype_col,
                           assay = "RNA",
                           use_layer = NULL,
                           wt_ref = NULL,
                           save_path = NULL) {

  if (missing(seurat_obj) || is.null(seurat_obj)) stop("Seurat object 'seurat_obj' must be provided")
  if (missing(celltype_col) || !nzchar(celltype_col)) stop("Cell type column name 'celltype_col' must be provided")
  if (!requireNamespace("Seurat", quietly = TRUE) && !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Please load Seurat or SeuratObject package first")
  }
  if (!assay %in% Seurat::Assays(seurat_obj)) stop(paste0("Assay not found in Seurat object: ", assay))

  if (is.null(wt_ref)) {
    if (exists("WT_database_data", envir = .GlobalEnv)) wt_ref <- get("WT_database_data", envir = .GlobalEnv)
  }
  if (is.null(wt_ref) || !is.data.frame(wt_ref)) {
    stop("Missing reference database 'wt_ref' (data.frame). Please provide 'wt_ref' or ensure 'WT_database_data' exists in the environment")
  }
  required_ref_cols <- c("Gene", "CellType", "ref_mean", "ref_sd")
  if (!all(required_ref_cols %in% colnames(wt_ref))) {
    stop(paste0("Reference database must contain columns: ", paste(required_ref_cols, collapse = ", ")))
  }

  Seurat::DefaultAssay(seurat_obj) <- assay
  assay_obj <- seurat_obj[[assay]]
  # Try using layer if available; fall back to slot if layer fails.
  # Prefer explicit calls inside tryCatch so we adapt to S3/S4 wrappers.
  if (!is.null(use_layer)) {
    # first try layer
    expr_mat <- tryCatch({
      as.matrix(SeuratObject::GetAssayData(assay_obj, layer = use_layer))
    }, error = function(e_layer) {
      # if layer failed, try slot with same name (older versions)
      tryCatch({
        as.matrix(SeuratObject::GetAssayData(assay_obj, slot = use_layer))
      }, error = function(e_slot) {
        stop("GetAssayData failed with both layer and slot for provided 'use_layer' (",
             use_layer, "). Original errors:\n layer: ", conditionMessage(e_layer),
             "\n slot: ", conditionMessage(e_slot))
      })
    })
  } else {
    # no use_layer specified: try to read 'data' as a layer first, then as a slot
    expr_mat <- tryCatch({
      as.matrix(SeuratObject::GetAssayData(assay_obj, layer = "data"))
    }, error = function(e_layer) {
      tryCatch({
        as.matrix(SeuratObject::GetAssayData(assay_obj, slot = "data"))
      }, error = function(e_slot) {
        stop("GetAssayData failed for default 'data' with both layer and slot. Original errors:\n layer: ",
             conditionMessage(e_layer), "\n slot: ", conditionMessage(e_slot))
      })
    })
  }

  if (is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) stop("Expression matrix must have row names (Gene) and column names (Cells)")

  meta <- seurat_obj@meta.data
  if (!celltype_col %in% colnames(meta)) stop(paste0("Column not found in meta.data: ", celltype_col))
  celltypes <- unique(meta[[celltype_col]])
  genes <- rownames(expr_mat)

  mean_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(celltypes),
                     dimnames = list(genes, celltypes))
  for (i in seq_along(celltypes)) {
    ct <- celltypes[i]
    cells_in_ct <- rownames(meta)[meta[[celltype_col]] == ct]
    if (length(cells_in_ct) > 0) {
      subm <- expr_mat[, cells_in_ct, drop = FALSE]
      mean_mat[, i] <- rowMeans(subm, na.rm = TRUE)
    }
  }

  mean_df_long <- as.data.frame(mean_mat)
  mean_df_long$Gene <- rownames(mean_df_long)
  mean_df_long <- tidyr::pivot_longer(mean_df_long, cols = -Gene, names_to = "CellType", values_to = "mean")

  final_df <- dplyr::left_join(mean_df_long, wt_ref, by = c("Gene", "CellType")) %>%
    dplyr::mutate(
      Zscore = dplyr::if_else(is.na(ref_sd) | is.na(mean) | is.na(ref_mean),
                              NA_real_,
                              (mean - ref_mean) / ref_sd,
                              missing = NA_real_),
      Status = dplyr::case_when(
        is.na(Zscore) ~ "NA",
        Zscore >  2   ~ "Above database range",
        Zscore < -2   ~ "Below database range",
        TRUE          ~ "Within database range"
      )
    ) %>%
    dplyr::select(Gene, CellType, rf_mean = ref_mean, rf_sd = ref_sd, mean, Zscore, Status) %>%
    dplyr::arrange(CellType, Gene)

  final_df <- dplyr::filter(final_df, Status != "NA")

  zscore_wide <- tidyr::pivot_wider(final_df,
                                    id_cols = Gene,
                                    names_from = CellType,
                                    values_from = Zscore)
  zscore_mat <- as.matrix(zscore_wide[, -1, drop = FALSE])
  rownames(zscore_mat) <- zscore_wide$Gene

  if (!is.null(save_path) && nzchar(save_path)) {
    out_base <- save_path
    if (!grepl("\\.xlsx$", out_base, ignore.case = TRUE)) {
      out_base <- paste0(out_base, "_Zscore.xlsx")
    }
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      utils::write.csv(final_df, file = paste0(tools::file_path_sans_ext(out_base), "_Status.csv"), row.names = FALSE)
      utils::write.csv(data.frame(Gene = rownames(mean_mat), mean_mat, check.names = FALSE),
                       file = paste0(tools::file_path_sans_ext(out_base), "_Mean_matrix.csv"), row.names = FALSE)
      utils::write.csv(data.frame(Gene = rownames(zscore_mat), zscore_mat, check.names = FALSE),
                       file = paste0(tools::file_path_sans_ext(out_base), "_Zscore_matrix.csv"), row.names = FALSE)
    } else {
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Status")
      openxlsx::writeData(wb, "Status", final_df)
      openxlsx::addWorksheet(wb, "Mean_matrix")
      openxlsx::writeData(wb, "Mean_matrix", data.frame(Gene = rownames(mean_mat), mean_mat, check.names = FALSE))
      openxlsx::addWorksheet(wb, "Zscore_matrix")
      openxlsx::writeData(wb, "Zscore_matrix", data.frame(Gene = rownames(zscore_mat), zscore_mat, check.names = FALSE))
      openxlsx::saveWorkbook(wb, out_base, overwrite = TRUE)
    }
    message("File saved: ", normalizePath(out_base))
  }

  return(zscore_mat)
}
