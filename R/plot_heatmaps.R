# ---------------------------------------------------------------------------
# Per-partition expression heatmaps (Female | Male facets).
# Requires Suggests: ggplot2, dplyr, tidyr, patchwork, viridis, scales, grid
# ---------------------------------------------------------------------------

# Row-z-score then clip to a window so a few outliers don't crush the palette.
.row_zscore <- function(M, clip = 2.5) {
  Z <- t(scale(t(M)))
  Z[!is.finite(Z)] <- 0
  Z[Z >  clip] <-  clip
  Z[Z < -clip] <- -clip
  Z
}

# One panel: row = gene (sorted by phase), col = sample (sorted by time within sex).
.partition_panel <- function(z_mat, sample_meta, phase_h,
                              title, ylab = NULL,
                              clip = 2.5,
                              show_y_axis = FALSE) {

  ord    <- order(phase_h, na.last = TRUE)
  z_mat  <- z_mat[ord, , drop = FALSE]

  df       <- as.data.frame(z_mat)
  df$gene  <- factor(rownames(z_mat), levels = rownames(z_mat))
  long     <- tidyr::pivot_longer(df, cols = -gene,
                                   names_to = "sample", values_to = "z")
  long     <- dplyr::left_join(long, sample_meta, by = "sample")

  ord_samples <- with(sample_meta,
                      order(as.character(sex),
                            as.numeric(as.character(time)),
                            as.character(ppid)))
  sample_levels <- sample_meta$sample[ord_samples]
  long$sample   <- factor(long$sample, levels = sample_levels)

  brks_df     <- sample_meta[ord_samples, , drop = FALSE]
  brks_df$pos <- seq_len(nrow(brks_df))
  pos <- sex <- time <- z <- gene <- sample <- NULL
  brks <- brks_df %>%
    dplyr::group_by(sex, time) %>%
    dplyr::summarise(pos = stats::median(pos), .groups = "drop")
  brks$label <- sprintf("%02d", as.integer(as.character(brks$time)))

  ggplot2::ggplot(long, ggplot2::aes(x = sample, y = gene, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(. ~ sex, scales = "free_x", space = "free_x",
                        switch = "x") +
    viridis::scale_fill_viridis(option = "magma", limits = c(-clip, clip),
                                 oob = scales::squish, name = "Row z-score") +
    ggplot2::scale_x_discrete(breaks = sample_levels[brks$pos],
                              labels = brks$label, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.spacing.x = grid::unit(0.4, "lines"),
                   axis.text.x     = ggplot2::element_text(angle = 0,
                                                            vjust = 0.5,
                                                            size = 7),
                   axis.text.y     = if (show_y_axis)
                                       ggplot2::element_text(size = 5)
                                     else ggplot2::element_blank(),
                   axis.ticks.y    = if (show_y_axis) ggplot2::element_line()
                                     else ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text.x    = ggplot2::element_text(face = "bold",
                                                            size = 11),
                   panel.grid      = ggplot2::element_blank(),
                   legend.position = "right") +
    ggplot2::labs(x = "Time of day (h)", y = ylab, title = title)
}


#' Per-partition expression heatmaps (Female | Male facets)
#'
#' Builds one heatmap per named gene set. Rows are sorted by acrophase
#' (via \code{phase_lookup}), columns by (group, time, subject). Combined
#' into a vertical patchwork PDF.
#'
#' @param results_df Data.frame with at least \code{ensembl_id} (or matching
#'   row names) and a per-gene \code{f_ph1} column (used by the default
#'   \code{phase_lookup}).
#' @param gene_sets Named list of character vectors, one per partition.
#' @param ncounts_filt_log Numeric gene x sample matrix (log-normalised).
#' @param metadata_filt Sample metadata with rownames matching
#'   \code{ncounts_filt_log} columns and columns \code{gender}, \code{TP.Hour},
#'   \code{ppid}.
#' @param out_pdf Output PDF path.
#' @param out_png_dir Optional directory for per-panel PNGs.
#' @param max_genes_per_panel Cap on rows per panel; keeps the most variable.
#' @param phase_lookup Optional \code{function(ids)} returning per-gene phase.
#'   Default uses \code{results_df$f_ph1}.
#' @param clip Z-score clip for the colour scale.
#' @return The combined patchwork object (invisibly).
#' @export
make_partition_heatmaps <- function(results_df,
                                     gene_sets,
                                     ncounts_filt_log,
                                     metadata_filt,
                                     out_pdf,
                                     out_png_dir = NULL,
                                     max_genes_per_panel = 400,
                                     phase_lookup = NULL,
                                     clip = 2.5) {

  for (pkg in c("ggplot2", "dplyr", "tidyr", "patchwork",
                "viridis", "scales", "grid"))
    .require_suggest(pkg)

  stopifnot(is.list(gene_sets))
  if (is.null(rownames(metadata_filt)))
    stop("metadata_filt needs rownames matching ncounts_filt_log columns.")

  sample_meta <- data.frame(
    sample = colnames(ncounts_filt_log),
    sex    = metadata_filt$gender[match(colnames(ncounts_filt_log),
                                          rownames(metadata_filt))],
    time   = metadata_filt$TP.Hour[match(colnames(ncounts_filt_log),
                                           rownames(metadata_filt))],
    ppid   = metadata_filt$ppid[match(colnames(ncounts_filt_log),
                                        rownames(metadata_filt))],
    stringsAsFactors = FALSE
  )
  sample_meta <- sample_meta[!is.na(sample_meta$sex) & !is.na(sample_meta$time), ]

  if (is.null(phase_lookup)) {
    if (!"f_ph1" %in% colnames(results_df))
      stop("phase_lookup not given and results_df has no f_ph1.")
    phase_lookup <- function(ids) {
      results_df$f_ph1[match(ids, results_df$ensembl_id)]
    }
  }

  panels <- list()
  for (nm in names(gene_sets)) {
    ids   <- intersect(gene_sets[[nm]], rownames(ncounts_filt_log))
    n_keep <- length(ids)
    if (n_keep < 2L) {
      message(sprintf("[heatmap] %s: only %d genes - skipping panel", nm, n_keep))
      panels[[nm]] <- NULL
      next
    }
    if (n_keep > max_genes_per_panel) {
      sds <- apply(ncounts_filt_log[ids, , drop = FALSE], 1, stats::sd, na.rm = TRUE)
      ids <- ids[order(-sds)][seq_len(max_genes_per_panel)]
    }
    mat   <- ncounts_filt_log[ids, sample_meta$sample, drop = FALSE]
    z_mat <- .row_zscore(mat, clip = clip)
    phs   <- phase_lookup(ids)
    panels[[nm]] <- .partition_panel(z_mat, sample_meta, phs,
                                     title = sprintf("%s   (n = %d)", nm, n_keep),
                                     clip = clip)
  }

  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0L) {
    message("[heatmap] no panels to plot - all gene sets empty.")
    return(invisible(NULL))
  }

  combined <- patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = "Per-partition expression heatmaps (group facets)",
      subtitle = "Rows: genes (by acrophase). Columns: samples (group -> time -> subject)."
    ) &
    ggplot2::theme(legend.position = "right")

  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(out_pdf, combined,
                  width = 9, height = 3 + 3 * length(panels),
                  limitsize = FALSE)
  message(sprintf("[heatmap] wrote %s", out_pdf))

  if (!is.null(out_png_dir)) {
    dir.create(out_png_dir, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_along(panels)) {
      f_png <- file.path(out_png_dir,
                         sprintf("partition_heatmap_%02d_%s.png",
                                  i, gsub("[^A-Za-z0-9]+", "_", names(panels)[i])))
      ggplot2::ggsave(f_png, panels[[i]], width = 9, height = 4.2, dpi = 150)
    }
  }
  invisible(combined)
}
