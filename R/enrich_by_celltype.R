#' Pathway enrichment analysis by cell type using Z-score matrix
#'
#' This function performs KEGG and GO enrichment analysis for each cell type
#' based on up- and down-regulated genes identified from a Z-score matrix.
#' Results include statistics, mapping tables, and enrichment outputs, with
#' optional saving to Excel.
#'
#' @param zscore_mat A numeric matrix of Z-scores (genes x cell types).
#' @param out_prefix Optional character string. Prefix for output Excel file.
#'   If NULL, results are returned in R only.
#' @param min_genes_for_enrich Minimum number of genes required to run enrichment
#'   (default: 3).
#' @param kegg_pcutoff Numeric. P-value cutoff for KEGG enrichment (default: 0.05).
#' @param go_pcutoff Numeric. P-value cutoff for GO enrichment (default: 0.05).
#' @param analysis_type Character vector specifying enrichment types. Options:
#'   \code{"all"}, \code{"KEGG"}, \code{"GO"}, \code{"BP"}, \code{"CC"}, \code{"MF"}.
#'   Default: \code{c("all","KEGG","GO","BP","CC","MF")}.
#' @importFrom magrittr %>%
#' @return A list containing enrichment results for each cell type. If
#'   \code{out_prefix} is provided, an Excel file is also saved.
#'
#' @details
#' For each cell type:
#' \enumerate{
#'   \item Identify up- and down-regulated genes (Z-score > 2 or < -2).
#'   \item Map ENSEMBL IDs to SYMBOLs, then to ENTREZ IDs.
#'   \item Perform KEGG enrichment (via \pkg{clusterProfiler}).
#'   \item Perform GO enrichment (BP, CC, MF) using \pkg{clusterProfiler}.
#'   \item Collect statistics and optionally save results to Excel.
#' }
#'
#' @examples
#' \dontrun{
#' res <- enrich_by_celltype(zscore_mat,
#'                           out_prefix = "output",
#'                           analysis_type = c("KEGG","GO"))
#' }
#'
#' @export
enrich_by_celltype <- function(zscore_mat,
                               out_prefix = NULL,
                               min_genes_for_enrich = 3,
                               kegg_pcutoff = 0.05,
                               go_pcutoff = 0.05,
                               analysis_type = c("all","KEGG","GO","BP","CC","MF")) {

  valid_types <- c("all","KEGG","GO","BP","CC","MF")
  if (is.null(analysis_type)) stop("analysis_type cannot be NULL")
  if (!all(analysis_type %in% valid_types)) stop("analysis_type must be one of: ", paste(valid_types, collapse = ", "))
  if ("all" %in% analysis_type) {
    analysis_sel <- c("KEGG","BP","CC","MF")
  } else if ("GO" %in% analysis_type) {
    analysis_sel <- unique(c(setdiff(analysis_type,"GO"), "BP","CC","MF"))
  } else {
    analysis_sel <- unique(analysis_type)
  }

  safe_char_cols <- function(df, cols) {
    for (c in intersect(cols, colnames(df))) {
      df[[c]] <- as.character(df[[c]])
      df[[c]][trimws(df[[c]]) == ""] <- NA_character_
    }
    df
  }

  ensembl_to_symbol <- function(ensembl_ids) {
    if (length(ensembl_ids) == 0) return(data.frame(ENSEMBL=character(0), SYMBOL=character(0), stringsAsFactors = FALSE))
    mart <- tryCatch(biomaRt::useMart("ensembl", dataset = "sscrofa_gene_ensembl"), error = function(e) NULL)
    if (is.null(mart)) return(data.frame(ENSEMBL=character(0), SYMBOL=character(0), stringsAsFactors = FALSE))
    bm <- tryCatch(biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                                  filters = "ensembl_gene_id", values = ensembl_ids, mart = mart),
                   error = function(e) data.frame())
    if (nrow(bm) == 0) return(data.frame(ENSEMBL=character(0), SYMBOL=character(0), stringsAsFactors = FALSE))
    colnames(bm) <- c("ENSEMBL","SYMBOL")
    bm <- safe_char_cols(as.data.frame(bm), c("ENSEMBL","SYMBOL")) %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::distinct()
    return(bm)
  }

  symbol_to_entrez_and_ensembl <- function(symbols) {
    symbols <- unique(stats::na.omit(symbols))
    if (length(symbols) == 0) return(list(collapsed = data.frame(), long = data.frame()))

    # mapped_db: try AnnotationDbi using org.Ss.eg.db if available
    mapped_db <- data.frame()
    if (requireNamespace("org.Ss.eg.db", quietly = TRUE)) {
      OrgDb_obj <- tryCatch(get("org.Ss.eg.db", envir = asNamespace("org.Ss.eg.db")), error = function(e) NULL)
      if (!is.null(OrgDb_obj)) {
        mapped_db <- tryCatch({
          AnnotationDbi::select(OrgDb_obj, keys = symbols, keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL"))
        }, error = function(e) data.frame())
      }
    }

    if (nrow(mapped_db) > 0) {
      mapped_db <- safe_char_cols(as.data.frame(mapped_db), c("ENTREZID","SYMBOL")) %>%
        dplyr::filter(!is.na(.data[["ENTREZID"]])) %>%
        dplyr::distinct() %>%
        dplyr::select(ENTREZID = .data[["ENTREZID"]], SYMBOL = .data[["SYMBOL"]]) %>%
        dplyr::mutate(ENSEMBL = NA_character_)
    }

    # biomaRt mapping (ENSEMBL + entrez from ensembl)
    mart <- tryCatch(biomaRt::useMart("ensembl", dataset = "sscrofa_gene_ensembl"), error = function(e) NULL)
    mapped_bm <- data.frame()
    if (!is.null(mart)) {
      bm <- tryCatch(biomaRt::getBM(attributes = c("external_gene_name","entrezgene_id","ensembl_gene_id"),
                                    filters = "external_gene_name", values = symbols, mart = mart),
                     error = function(e) data.frame())
      if (nrow(bm) > 0) {
        colnames(bm) <- c("SYMBOL","ENTREZID_bm","ENSEMBL")
        bm <- safe_char_cols(as.data.frame(bm), c("SYMBOL","ENTREZID_bm","ENSEMBL"))
        bm <- bm %>%
          dplyr::mutate(ENTREZID = as.character(.data[["ENTREZID_bm"]])) %>%
          dplyr::select(ENTREZID = .data[["ENTREZID"]], SYMBOL = .data[["SYMBOL"]], ENSEMBL = .data[["ENSEMBL"]]) %>%
          dplyr::distinct()
        mapped_bm <- bm
      }
    }

    # if nothing mapped, return empty
    if (nrow(mapped_db) == 0 && nrow(mapped_bm) == 0) {
      return(list(collapsed = data.frame(), long = data.frame()))
    }

    # merge db and bm results when both available
    if (nrow(mapped_db) == 0) {
      merged <- mapped_bm
    } else if (nrow(mapped_bm) == 0) {
      merged <- mapped_db
    } else {
      merged_tmp <- dplyr::full_join(mapped_db, mapped_bm, by = "SYMBOL", suffix = c(".db", ".bm"))
      merged_tmp <- safe_char_cols(as.data.frame(merged_tmp),
                                   c("ENTREZID.db","ENTREZID.bm","ENSEMBL.db","ENSEMBL.bm","SYMBOL"))
      merged <- merged_tmp %>%
        dplyr::mutate(
          ENTREZID_final = ifelse(!is.na(.data[["ENTREZID.db"]]) & .data[["ENTREZID.db"]] != "", .data[["ENTREZID.db"]], .data[["ENTREZID.bm"]]),
          ENSEMBL_final  = ifelse(!is.na(.data[["ENSEMBL.bm"]]) & .data[["ENSEMBL.bm"]] != "", .data[["ENSEMBL.bm"]], .data[["ENSEMBL.db"]])
        ) %>%
        dplyr::transmute(
          ENTREZID = as.character(.data[["ENTREZID_final"]]),
          SYMBOL  = as.character(.data[["SYMBOL"]]),
          ENSEMBL = as.character(.data[["ENSEMBL_final"]])
        ) %>%
        dplyr::distinct()
    }

    # normalize empty strings to NA
    merged <- merged %>%
      safe_char_cols(c("ENTREZID","SYMBOL","ENSEMBL")) %>%
      dplyr::mutate(
        ENTREZID = ifelse(.data[["ENTREZID"]] == "" | is.na(.data[["ENTREZID"]]), NA_character_, .data[["ENTREZID"]]),
        ENSEMBL  = ifelse(.data[["ENSEMBL"]] == ""  | is.na(.data[["ENSEMBL"]]), NA_character_, .data[["ENSEMBL"]])
      )

    # collapsed (one row per SYMBOL) and long (one row per ENTREZID)
    collapsed <- merged %>%
      dplyr::group_by(.data[["SYMBOL"]]) %>%
      dplyr::summarise(
        ENTREZID = paste(unique(stats::na.omit(.data[["ENTREZID"]])), collapse = ";"),
        ENSEMBL  = paste(unique(stats::na.omit(.data[["ENSEMBL"]])), collapse = ";"),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        ENTREZID = ifelse(.data[["ENTREZID"]] == "", NA_character_, .data[["ENTREZID"]]),
        ENSEMBL  = ifelse(.data[["ENSEMBL"]] == "", NA_character_, .data[["ENSEMBL"]])
      ) %>%
      dplyr::select(SYMBOL = .data[["SYMBOL"]], ENTREZID = .data[["ENTREZID"]], ENSEMBL = .data[["ENSEMBL"]])

    long <- merged %>%
      dplyr::filter(!is.na(.data[["ENTREZID"]])) %>%
      dplyr::distinct() %>%
      dplyr::select(ENTREZID = .data[["ENTREZID"]], SYMBOL = .data[["SYMBOL"]], ENSEMBL = .data[["ENSEMBL"]])

    return(list(collapsed = collapsed, long = long))
  }


  convert_geneid_col_to_symbol <- function(df, long_map) {
    if (!("geneID" %in% colnames(df))) return(df)
    if (is.null(long_map) || nrow(long_map) == 0) return(df)
    dict <- long_map %>%
      dplyr::distinct(ENTREZID, SYMBOL) %>%
      dplyr::group_by(ENTREZID) %>%
      dplyr::summarise(SYMBOL = dplyr::first(SYMBOL), .groups = "drop")
    entrez2sym <- stats::setNames(dict$SYMBOL, dict$ENTREZID)
    df$geneID <- sapply(df$geneID, function(g) {
      ids <- unlist(strsplit(as.character(g), "/"))
      syms <- unname(entrez2sym[ids])
      syms[is.na(syms)] <- ids[is.na(syms)]
      paste(unique(syms), collapse = "/")
    }, USE.NAMES = FALSE)
    df
  }

  celltypes <- colnames(zscore_mat)
  results_list <- list()
  wb <- if (!is.null(out_prefix)) openxlsx::createWorkbook() else NULL

  for (ct in celltypes) {
    zvec <- zscore_mat[, ct]
    valid_idx <- which(!is.na(zvec))
    if (length(valid_idx) == 0) {
      message("Cell type ", ct, " has all NA values, skipped")
      next
    }
    zvec <- zvec[valid_idx]
    genes_all <- rownames(zscore_mat)[valid_idx]

    total_genes <- length(genes_all)
    genes_up <- names(zvec)[zvec > 2]
    genes_dn <- names(zvec)[zvec < -2]
    n_up <- length(genes_up)
    n_dn <- length(genes_dn)

    ensembl_up <- genes_up[grepl("^ENSSSCG", genes_up)]
    ensembl_dn <- genes_dn[grepl("^ENSSSCG", genes_dn)]
    symbol_up <- setdiff(genes_up, ensembl_up)
    symbol_dn <- setdiff(genes_dn, ensembl_dn)
    n_ensembl_up <- length(ensembl_up)
    n_ensembl_dn <- length(ensembl_dn)
    n_symbol_up <- length(symbol_up)
    n_symbol_dn <- length(symbol_dn)

    map_up_ensembl <- ensembl_to_symbol(ensembl_up)
    map_dn_ensembl <- ensembl_to_symbol(ensembl_dn)

    all_symbols_up <- unique(c(symbol_up, if (nrow(map_up_ensembl) > 0) map_up_ensembl$SYMBOL else character(0)))
    all_symbols_dn <- unique(c(symbol_dn, if (nrow(map_dn_ensembl) > 0) map_dn_ensembl$SYMBOL else character(0)))

    up_maps <- symbol_to_entrez_and_ensembl(all_symbols_up)
    dn_maps <- symbol_to_entrez_and_ensembl(all_symbols_dn)
    up_collapsed <- up_maps$collapsed
    up_long <- up_maps$long
    dn_collapsed <- dn_maps$collapsed
    dn_long <- dn_maps$long

    n_mapped_entrez_up <- if (nrow(up_long) > 0) length(unique(up_long$ENTREZID)) else 0
    n_mapped_entrez_dn <- if (nrow(dn_long) > 0) length(unique(dn_long$ENTREZID)) else 0

    message(sprintf("Cell type %s - total genes %d; up %d; down %d",
                    ct, total_genes, length(genes_up), length(genes_dn)))
    message(sprintf("Mapped to ENTREZ IDs: up = %d; down = %d",
                    length(unique(up_maps$long$ENTREZID)),
                    length(unique(dn_maps$long$ENTREZID))))

    entrez_up <- if (nrow(up_long) > 0) unique(up_long$ENTREZID) else character(0)
    entrez_dn <- if (nrow(dn_long) > 0) unique(dn_long$ENTREZID) else character(0)

    kegg_up_df <- data.frame(); kegg_dn_df <- data.frame()
    go_up_list <- list(); go_dn_list <- list()

    if ("KEGG" %in% analysis_sel && length(entrez_up) >= min_genes_for_enrich) {
      kegg_up <- tryCatch(clusterProfiler::enrichKEGG(gene = entrez_up, organism = "ssc", pvalueCutoff = kegg_pcutoff), error = function(e) NULL)
      if (!is.null(kegg_up)) kegg_up_df <- convert_geneid_col_to_symbol(as.data.frame(kegg_up), up_long)
    }
    if ("KEGG" %in% analysis_sel && length(entrez_dn) >= min_genes_for_enrich) {
      kegg_dn <- tryCatch(clusterProfiler::enrichKEGG(gene = entrez_dn, organism = "ssc", pvalueCutoff = kegg_pcutoff), error = function(e) NULL)
      if (!is.null(kegg_dn)) kegg_dn_df <- convert_geneid_col_to_symbol(as.data.frame(kegg_dn), dn_long)
    }

    go_onts <- intersect(c("BP", "CC", "MF"), analysis_sel)
    if ("GO" %in% analysis_sel) go_onts <- c("BP", "CC", "MF")
    if ("all" %in% analysis_sel) go_onts <- c("BP", "CC", "MF")


    if ((length(entrez_up) < min_genes_for_enrich) && (length(entrez_dn) < min_genes_for_enrich)) {

    } else {

      if (!requireNamespace("org.Ss.eg.db", quietly = TRUE)) {
        warning(sprintf("Package 'org.Ss.eg.db' not installed; GO enrichment will be skipped for cell type %s", ct))
      } else {
        OrgDb <- get("org.Ss.eg.db", envir = asNamespace("org.Ss.eg.db"))
        for (ont in unique(go_onts)) {
          if (length(entrez_up) >= min_genes_for_enrich) {
            go_up <- tryCatch(
              clusterProfiler::enrichGO(
                gene = entrez_up,
                OrgDb = OrgDb,
                ont = ont,
                keyType = "ENTREZID",
                pvalueCutoff = go_pcutoff
              ),
              error = function(e) NULL
            )
            if (!is.null(go_up)) {
              go_up_list[[ont]] <- convert_geneid_col_to_symbol(as.data.frame(go_up), up_long)
            }
          }

          if (length(entrez_dn) >= min_genes_for_enrich) {
            go_dn <- tryCatch(
              clusterProfiler::enrichGO(
                gene = entrez_dn,
                OrgDb = OrgDb,
                ont = ont,
                keyType = "ENTREZID",
                pvalueCutoff = go_pcutoff
              ),
              error = function(e) NULL
            )
            if (!is.null(go_dn)) {
              go_dn_list[[ont]] <- convert_geneid_col_to_symbol(as.data.frame(go_dn), dn_long)
            }
          }
        }
      }
    }

    stats <- list(
      total_genes = total_genes,
      n_up = n_up,
      n_dn = n_dn,
      n_ensembl_up = n_ensembl_up,
      n_ensembl_dn = n_ensembl_dn,
      n_symbol_up = n_symbol_up,
      n_symbol_dn = n_symbol_dn,
      n_mapped_entrez_up = n_mapped_entrez_up,
      n_mapped_entrez_dn = n_mapped_entrez_dn
    )

    results_list[[ct]] <- list(
      stats = stats,
      up = list(genes_raw = genes_up, symbol2entrez_collapsed = up_collapsed, symbol2entrez_long = up_long, kegg = kegg_up_df, go = go_up_list),
      down = list(genes_raw = genes_dn, symbol2entrez_collapsed = dn_collapsed, symbol2entrez_long = dn_long, kegg = kegg_dn_df, go = go_dn_list)
    )

    if (!is.null(wb)) {
      safe <- function(x) substr(gsub("[^A-Za-z0-9_]","_", x), 1, 31)
      write_safe <- function(sheet, df) {
        sn <- safe(sheet); openxlsx::addWorksheet(wb, sn)
        if (!is.data.frame(df) || nrow(df) == 0) {
          openxlsx::writeData(wb, sn, data.frame(note = "no results"))
        } else {
          openxlsx::writeData(wb, sn, df)
        }
      }
      if (nrow(kegg_up_df) > 0) write_safe(paste0(ct, "_up_KEGG"), kegg_up_df)
      if (nrow(kegg_dn_df) > 0) write_safe(paste0(ct, "_dn_KEGG"), kegg_dn_df)
      for (ont in names(go_up_list)) write_safe(paste0(ct, "_up_GO_", ont), go_up_list[[ont]])
      for (ont in names(go_dn_list)) write_safe(paste0(ct, "_dn_GO_", ont), go_dn_list[[ont]])
    }
  }

  if (!is.null(out_prefix) && !is.null(wb)) {
    out_file <- paste0(out_prefix, "_byCelltype_pathway_enrichment.xlsx")
    openxlsx::saveWorkbook(wb, file = out_file, overwrite = TRUE)
    message("Results saved: ", normalizePath(out_file))
    return(list(results = results_list, out_file = out_file))
  } else {
    return(list(results = results_list))
  }
}
