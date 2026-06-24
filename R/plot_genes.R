# ---------------------------------------------------------------------------
# Per-gene and multi-gene ggplot helpers.
# Requires Suggests: ggplot2 + (for grids) patchwork
# ---------------------------------------------------------------------------

# Internal: ensure a Suggests package is loaded, else error politely.
.require_suggest <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this function. ",
         "Install with: install.packages('", pkg, "').", call. = FALSE)
  }
}

#' Plot one gene over time, faceted by group
#'
#' Looks up the gene by external symbol in \code{ncounts_annot}, then plots
#' \code{log2(counts + 1)} against the time-of-day column, faceted by the
#' grouping column. Subject lines connect repeated samples from the same ID.
#'
#' @param external_gene_name Character symbol.
#' @param ncounts Numeric gene x sample matrix.
#' @param ncounts_annot Data.frame with row-names matching \code{ncounts}
#'   and a column \code{external_gene_name}.
#' @param metadata_filt Data.frame with one row per sample, containing the
#'   sample-name, time-of-day, group, and per-subject columns.
#' @param tp_col Time-of-day column in \code{metadata_filt}; default
#'   \code{"TP.Hour"}.
#' @param sample_col Sample-identifier column; default \code{"Sample.ID"}.
#' @param colour_col Per-subject column to colour by; default \code{"ppid"}.
#' @param group_col Group column to facet on; default \code{"gender"}.
#' @param scale_within_ppid Z-score within subject before plotting.
#' @return A ggplot object.
#' @export
plot_gene_over_TP <- function(external_gene_name, ncounts, ncounts_annot,
                               metadata_filt, tp_col = "TP.Hour",
                               sample_col = "Sample.ID",
                               colour_col = "ppid",
                               group_col  = "gender",
                               scale_within_ppid = FALSE) {
  .require_suggest("ggplot2")
  stopifnot("external_gene_name" %in% colnames(ncounts_annot))
  ncounts  <- as.matrix(ncounts)
  cand_ids <- rownames(ncounts_annot)[ncounts_annot$external_gene_name == external_gene_name]
  cand_ids <- intersect(cand_ids, rownames(ncounts))
  if (length(cand_ids) == 0) stop("Gene not found: ", external_gene_name)

  if (length(cand_ids) > 1) {
    means   <- rowMeans(ncounts[cand_ids, , drop = FALSE], na.rm = TRUE)
    gene_id <- cand_ids[which.max(means)]
  } else gene_id <- cand_ids[1]

  samp <- metadata_filt[[sample_col]]
  keep <- samp %in% colnames(ncounts)
  if (!any(keep)) stop("No matching samples between metadata_filt and ncounts.")

  df <- data.frame(
    TP     = metadata_filt[[tp_col]][keep],
    ppid   = factor(metadata_filt[[colour_col]][keep]),
    group  = factor(metadata_filt[[group_col]][keep]),
    y      = log2(as.numeric(ncounts[gene_id, samp[keep]]) + 1)
  )

  if (scale_within_ppid) {
    df$y_plot <- stats::ave(df$y, df$ppid, FUN = function(z) {
      if (stats::sd(z, na.rm = TRUE) == 0) return(rep(0, length(z)))
      as.numeric(scale(z))
    })
    ylab_txt     <- "Scaled within PPID (z-score of log2 counts + 1)"
    title_suffix <- "within-PPID scaled"
  } else {
    df$y_plot    <- df$y
    ylab_txt     <- "log2(normalised counts + 1)"
    title_suffix <- "raw"
  }

  # Aliased to avoid CMD check 'no visible binding for global variable':
  TP <- y_plot <- ppid <- NULL
  ggplot2::ggplot(df, ggplot2::aes(x = TP, y = y_plot, colour = ppid)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15),
                        alpha = 0.7, size = 2) +
    ggplot2::stat_summary(ggplot2::aes(group = ppid), fun = mean,
                          geom = "line", linewidth = 0.8) +
    ggplot2::facet_grid(stats::reformulate("group")) +
    ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1) +
    ggplot2::labs(x = tp_col, y = ylab_txt, colour = colour_col,
                  title = paste0(external_gene_name, " (", gene_id, "; ",
                                  title_suffix, ")"))
}


#' Single-gene panel for a labelled gene grid
#'
#' Wraps \code{\link{plot_gene_over_TP}} into a square panel sized for
#' \code{plot_gene_set_grid()}. Returns a placeholder if the gene is
#' absent from the dataset.
#'
#' @param gene Gene symbol.
#' @param ncounts,ncounts_annot,metadata_filt See \code{\link{plot_gene_over_TP}}.
#' @param scale_within_ppid Z-score within subject.
#' @param group_col Group / facet column.
#' @param ... Extra arguments to \code{\link{plot_gene_over_TP}}.
#' @return A ggplot object.
#' @export
plot_gene_panel <- function(gene, ncounts, ncounts_annot, metadata_filt,
                             scale_within_ppid = FALSE,
                             group_col = "gender", ...) {
  .require_suggest("ggplot2")
  p <- tryCatch(
    plot_gene_over_TP(gene, ncounts, ncounts_annot, metadata_filt,
                       scale_within_ppid = scale_within_ppid,
                       group_col = group_col, ...),
    error = function(e) NULL)

  if (is.null(p)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = sprintf("%s\n(not in dataset)", gene),
                          size = 4) +
        ggplot2::theme_void() +
        ggplot2::theme(aspect.ratio = 0.7))
  }

  y_label <- if (scale_within_ppid) "Within-subject z-score"
             else                   "log2 normalised counts"

  p +
    ggplot2::facet_grid(stats::reformulate(group_col)) +
    ggplot2::labs(title = gene, x = "Time of day (h)", y = y_label) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(aspect.ratio    = 1,
                   legend.position = "none",
                   plot.title      = ggplot2::element_text(face = "bold",
                                                            hjust = 0.5,
                                                            size = 13),
                   strip.text      = ggplot2::element_text(face = "bold",
                                                            size = 11),
                   axis.title      = ggplot2::element_text(size = 11),
                   axis.text       = ggplot2::element_text(size  = 9))
}


# Choose a sensible grid width
.ncol_for <- function(n) {
  if (n <= 2) 2L
  else if (n <= 4) 2L
  else if (n <= 6) 3L
  else if (n <= 9) 3L
  else 4L
}


#' Build a multi-gene PDF grid (one panel per gene, group as facet)
#'
#' @param set_name Title of the grid (used in the filename, sanitised).
#' @param gene_list Character vector of gene symbols.
#' @param ncounts,ncounts_annot,metadata_filt See \code{\link{plot_gene_over_TP}}.
#' @param scale_within_ppid Z-score within subject.
#' @param group_col Group / facet column.
#' @param out_dir Directory to write the PDF.
#' @return The combined patchwork object (invisibly).
#' @export
plot_gene_set_grid <- function(set_name, gene_list,
                                ncounts, ncounts_annot, metadata_filt,
                                scale_within_ppid = FALSE,
                                group_col = "gender",
                                out_dir = "./output/figures/gene_grids") {
  .require_suggest("ggplot2")
  .require_suggest("patchwork")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  panels <- lapply(gene_list, plot_gene_panel,
                   ncounts        = ncounts,
                   ncounts_annot  = ncounts_annot,
                   metadata_filt  = metadata_filt,
                   scale_within_ppid = scale_within_ppid,
                   group_col      = group_col)

  ncol_g <- .ncol_for(length(gene_list))
  nrow_g <- ceiling(length(gene_list) / ncol_g)

  body <- patchwork::wrap_plots(panels, ncol = ncol_g)

  scale_tag  <- if (scale_within_ppid) "within-PPID scaled" else "raw log-normalised"
  full_title <- sprintf("%s  -  %s", set_name, scale_tag)

  out <- body + patchwork::plot_annotation(
    title = full_title,
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                size = 16,
                                                                hjust = 0.5))
  )

  fname <- sprintf("%s/%s_%s.pdf",
                    out_dir,
                    gsub("[^A-Za-z0-9]+", "_", set_name),
                    if (scale_within_ppid) "scaled" else "raw")
  ggplot2::ggsave(fname, out,
                  width  = 4.5 * ncol_g + 1.5,
                  height = 3.0 * nrow_g + 1.0,
                  limitsize = FALSE)
  message(sprintf("[gene-grid] %-55s -> %s", set_name, fname))
  invisible(out)
}
