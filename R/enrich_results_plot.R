#' Plot enrichment results for a given cell type
#'
#' This function generates dot plots and bar plots for enrichment results
#' (KEGG or GO categories) from the output of \code{enrich_by_celltype}.
#' It supports flexible column naming and automatically detects enrichment
#' score, q-value, and gene count columns.
#'
#' @param res A result list from \code{enrich_by_celltype}, or a compatible
#'   list structure containing enrichment results.
#' @param celltype Character string. The cell type to plot (must exist in \code{res}).
#' @param direction Character. Either \code{"up"} or \code{"down"} (default: both options).
#' @param which Character vector. Which enrichment types to plot. Options:
#'   \code{"KEGG"}, \code{"BP"}, \code{"CC"}, \code{"MF"}.
#' @param top_n Integer. Number of top terms to display (default: 20).
#' @param title Optional character string. Plot title. If NULL, no title is shown.
#' @param enrichment_col Candidate column names for enrichment score
#'   (default: \code{c("FoldEnrichment","RichFactor")}).
#' @param qvalue_col Candidate column names for q-value or adjusted p-value
#'   (default: \code{c("qvalue","p.adjust","pvalue")}).
#' @param geneid_col Candidate column names for term description or ID
#'   (default: \code{c("Description","ID")}).
#' @param list_count_col Candidate column names for gene counts
#'   (default: \code{c("Count","GeneRatio","BgRatio")}).
#'
#' @importFrom magrittr %>%
#' @return A list of plots for each enrichment type requested. Each element
#'   contains:
#'   \itemize{
#'     \item \code{dot}: ggplot2 dot plot
#'     \item \code{bar}: ggplot2 bar plot
#'     \item \code{df}: processed data frame used for plotting
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts enrichment results for the specified cell type and direction.
#'   \item Detects enrichment score, q-value, and gene count columns.
#'   \item Prepares a ranked data frame of top terms.
#'   \item Generates dot and bar plots with ggplot2.
#' }
#'
#' @examples
#' \dontrun{
#' plots <- enrich_results_plot(res,
#'                              celltype = "B cells",
#'                              direction = "up",
#'                              which = c("KEGG","BP"),
#'                              top_n = 15)
#' plots$KEGG$dot
#' plots$BP$bar
#' }
#'
#' @export
enrich_results_plot <- function(res,
                                celltype = NULL,
                                direction = c("up","down"),
                                which = c("KEGG","BP","CC","MF"),
                                top_n = 20,
                                title = NULL,
                                enrichment_col = c("FoldEnrichment","RichFactor"),
                                qvalue_col = c("qvalue","p.adjust","pvalue"),
                                geneid_col = c("Description","ID"),
                                list_count_col = c("Count","GeneRatio","BgRatio")) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required")

  results_list <- if ("results" %in% names(res)) res$results else res
  if (is.null(celltype)) stop("Please specify 'celltype', e.g., 'B cells'")
  if (!(celltype %in% names(results_list))) stop("Specified 'celltype' not found in results")
  direction <- match.arg(direction)
  which <- unique(which)

  pick_first_existing <- function(cands, df) {
    cands <- as.character(cands)
    for (c in cands) if (c %in% colnames(df)) return(c)
    return(NA_character_)
  }

  get_df_for <- function(ct, dir, typ) {
    obj <- results_list[[ct]][[dir]]
    if (is.null(obj)) return(NULL)

    if (typ == "KEGG") {
      df <- obj$kegg
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
      df[] <- lapply(df, function(x) if (is.factor(x)) as.character(x) else x)
      return(as.data.frame(df))
    }

    go_obj <- obj$go
    if (is.null(go_obj)) return(NULL)

    if (is.list(go_obj) && !is.data.frame(go_obj)) {
      if (typ %in% names(go_obj)) {
        df <- go_obj[[typ]]
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
          df[] <- lapply(df, function(x) if (is.factor(x)) as.character(x) else x)
          return(as.data.frame(df))
        } else return(NULL)
      } else {
        nm <- names(go_obj)
        match_idx <- which(toupper(nm) == toupper(typ))
        if (length(match_idx) == 1) {
          df <- go_obj[[match_idx]]
          df[] <- lapply(df, function(x) if (is.factor(x)) as.character(x) else x)
          return(as.data.frame(df))
        }
        return(NULL)
      }
    }

    if (is.data.frame(go_obj)) {
      df_all <- as.data.frame(go_obj, stringsAsFactors = FALSE)
      df_all[] <- lapply(df_all, function(x) if (is.factor(x)) as.character(x) else x)

      ont_cols <- c("Ontology","ONTOLOGY","ONT","ont","Category","subcategory")
      ont_col <- intersect(ont_cols, colnames(df_all))[1]
      if (!is.na(ont_col)) {
        df_sub <- dplyr::filter(df_all, toupper(.data[[ont_col]]) == toupper(typ))
        if (nrow(df_sub) > 0) return(df_sub)
      }

      if ("Description" %in% colnames(df_all)) {
        # use explicit column access to avoid NSE notes
        desc_vec <- as.character(df_all[["Description"]])
        hit_idx <- which(stringr::str_detect(toupper(desc_vec), toupper(typ)))
        if (length(hit_idx) > 0) {
          df_sub <- df_all[hit_idx, , drop = FALSE]
          return(df_sub)
        }
      }
      return(NULL)
    }

    return(NULL)
  }

  prepare_plot_df <- function(df, top_n) {
    ec <- pick_first_existing(enrichment_col, df)
    qc <- pick_first_existing(qvalue_col, df)
    termc <- pick_first_existing(geneid_col, df)
    cc <- pick_first_existing(list_count_col, df)

    if (is.na(ec)) stop("Enrichment column not found. Please check FoldEnrichment or RichFactor")
    if (is.na(qc)) stop("q-value column not found. Please check qvalue/p.adjust/pvalue")

    if (is.na(termc)) {
      termc <- pick_first_existing(c("Description","Term","ID"), df)
      if (is.na(termc)) termc <- colnames(df)[1]
    }

    dfw <- dplyr::mutate(df,
                         EnrichmentVal = as.numeric(as.character(.data[[ec]])),
                         QVal = as.numeric(as.character(.data[[qc]])),
                         Term = as.character(.data[[termc]])
    )

    if (!is.na(cc)) {
      dfw$ListHits <- suppressWarnings(as.numeric(as.character(dfw[[cc]])))
      if (all(is.na(dfw$ListHits))) {
        dfw$ListHits <- sapply(dfw[[cc]], function(x) {
          if (is.na(x) || x == "") return(NA_integer_)
          nums <- stringr::str_extract_all(as.character(x), "\\d+")[[1]]
          if (length(nums) >= 1) return(as.integer(nums[1])) else return(NA_integer_)
        })
      }
    } else {
      geneid_cand <- pick_first_existing(c("geneID","GeneRatio","gene_id","gene"), df)
      if (!is.na(geneid_cand) && geneid_cand %in% colnames(df)) {
        dfw$ListHits <- sapply(dfw[[geneid_cand]], function(g) {
          if (is.na(g) || g == "") return(NA_integer_)
          length(unlist(strsplit(as.character(g), "/")))
        })
      } else {
        dfw$ListHits <- NA_integer_
      }
    }

    dfw <- dplyr::filter(dfw, !is.na(EnrichmentVal) & !is.na(QVal))
    if (nrow(dfw) == 0) return(NULL)

    dfw <- dplyr::arrange(dfw, dplyr::desc(EnrichmentVal))
    if (!is.null(top_n) && is.numeric(top_n) && nrow(dfw) > top_n) dfw <- dfw[1:top_n, , drop = FALSE]

    dfw <- dplyr::arrange(dfw, EnrichmentVal)
    dfw$Term <- factor(dfw$Term, levels = unique(dfw$Term))
    dfw$LogQ <- -log10(dfw$QVal)

    return(dfw)
  }


  out <- list()
  for (typ in which) {
    df <- get_df_for(celltype, direction, typ)
    if (is.null(df)) {
      out[[typ]] <- NULL
      next
    }

    df_plot <- prepare_plot_df(df, top_n = top_n)
    if (is.null(df_plot)) {
      out[[typ]] <- NULL
      next
    }


    if (is.null(title)) {
      title_args_dot <- list(x = quote(EnrichmentVal), y = quote(Term))
      label_args_dot <- ggplot2::labs(x = "Enrichment Score", y = NULL, color = expression(-log[10]("Q value")), size = "Gene Count")
      label_args_bar <- ggplot2::labs(x = "Enrichment Score", y = NULL, fill = expression(-log[10]("Q value")))
    } else {
      label_args_dot <- ggplot2::labs(title = title, x = "Enrichment Score", y = NULL, color = expression(-log[10]("Q value")), size = "Gene Count")
      label_args_bar <- ggplot2::labs(title = title, x = "Enrichment Score", y = NULL, fill = expression(-log[10]("Q value")))
    }


    p_dot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = EnrichmentVal, y = Term)) +
      ggplot2::geom_point(ggplot2::aes(size = ListHits, color = LogQ), alpha = 0.85) +
      ggplot2::scale_color_gradient(low = "#5DADE2", high = "#E74C3C") +
      ggplot2::theme_minimal(base_size = 14) +
      label_args_dot +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.y = ggplot2::element_text(size = 11, face = "bold", color = "black"),
        axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "gray90")
      )


    p_bar <- ggplot2::ggplot(df_plot, ggplot2::aes(x = EnrichmentVal, y = Term, fill = LogQ)) +
      ggplot2::geom_col(alpha = 0.9, width = 0.7) +
      ggplot2::scale_fill_gradient(low = "#5DADE2", high = "#E74C3C") +
      ggplot2::theme_minimal(base_size = 14) +
      label_args_bar +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.y = ggplot2::element_text(size = 11, face = "bold", color = "black"),
        axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "gray90")
      )

    out[[typ]] <- list(dot = p_dot, bar = p_bar, df = df_plot)
  }

  return(out)
}

