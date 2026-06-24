# ---------------------------------------------------------------------------
# dry_lmm_full()
# Mixed-model extension of dryR with model selection over rhythm and mean
# partitions and a 3-step fallback ladder on the random structure.
# ---------------------------------------------------------------------------

#' Fit dryR-style mixed-model partition selection
#'
#' Per gene, fits the 5 rhythm models and 2 mean models from dryR using
#' \code{glmmTMB::glmmTMB} with random intercept and random slopes on
#' \eqn{\sin(2\pi t/T)} and \eqn{\cos(2\pi t/T)}. The random structure
#' uses a 3-step fallback ladder (full -> diagonal -> intercept-only) to
#' protect against singular fits and non-convergence (see
#' \code{\link{compute_glmmTMB_full}}).
#'
#' Returns BIC weights, the chosen partition per gene, fixed-effect cosinor
#' parameters per group, per-subject BLUPs, and a variance summary.
#'
#' @param data Gene x sample numeric matrix (rownames = gene IDs).
#' @param group Sample-level group label (e.g. sex). Length = ncol(data).
#' @param time Sample-level time of day (numeric, same units as
#'   \code{period}).
#' @param ID Sample-level subject identifier.
#' @param period Rhythm period; default 24.
#' @param sample_name Column names for the response matrix.
#' @param batch Optional batch covariate (length = ncol(data)).
#' @param weights Optional gene x sample weight matrix.
#' @param extract_per_subject If \code{TRUE} (default), build per-subject
#'   BLUP table and variance summary.
#' @param parallel_type One of \code{"auto"}, \code{"fork"}, \code{"psock"},
#'   \code{"serial"}.
#' @param n.cores Number of cores; default 60\% of available cores.
#' @return List with elements:
#'   \describe{
#'     \item{time, group, ID}{Re-ordered sample-level vectors.}
#'     \item{results}{Wide data.frame: response values + parameters per gene.}
#'     \item{BICW_rhythm, BICW_mean}{Per-gene BIC weights for the 5 rhythm
#'       and 2 mean partitions.}
#'     \item{values}{The input \code{data} (re-ordered to match \code{time}).}
#'     \item{parameters}{Per-gene cosinor parameters and chosen partition.}
#'     \item{per_subject}{Optional per-subject BLUP table (\code{NULL} if
#'       \code{extract_per_subject = FALSE}).}
#'     \item{variance_summary}{Optional per-gene variance summary.}
#'     \item{rstruct}{Which random structure won per gene (rhythm + mean).}
#'   }
#' @export
dry_lmm_full <- function(data, group, time, ID, period = 24,
                          sample_name = colnames(data),
                          batch = rep("A", length(sample_name)),
                          weights = NULL,
                          extract_per_subject = TRUE,
                          parallel_type = c("auto", "fork", "psock", "serial"),
                          n.cores = round(parallel::detectCores() * 0.6, 0)) {

  parallel_type <- .resolve_parallel_type(parallel_type)
  message("parallel_type = '", parallel_type, "'")

  # ---- set up parallel backend ----------------------------------------
  cl <- NULL
  if (parallel_type == "psock" && n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(glmmTMB))
      suppressPackageStartupMessages(library(dryRlmm))
    })
  }
  papply <- .make_papply(parallel_type, n.cores, cl)

  vec <- FALSE
  if (is.vector(data)) {
    data <- rbind(data, data)
    rownames(data) <- c("X1", "X2")
    vec <- TRUE
  }

  if (!is.null(weights)) {
    if (!identical(dim(weights), dim(data)))
      stop("`weights` must have the same dimensions as `data`.")
    if (is.null(rownames(weights))) rownames(weights) <- rownames(data)
  }

  # ---- align sample order ----------------------------------------------
  sel         <- order(group, time)
  time        <- time[sel]
  group       <- group[sel]
  data        <- data[, sel]
  batch       <- batch[sel]
  ID          <- ID[sel]
  sample_name <- as.character(sample_name[sel])
  if (!is.null(weights)) weights <- weights[, sel, drop = FALSE]

  s1 <- sin(2 * pi * time / period)
  c1 <- cos(2 * pi * time / period)
  conds <- cbind(group, s1, c1, batch)
  colnames(conds) <- c("group", "s1", "c1", "batch")
  colData <- data.frame(row.names = colnames(data), conds)

  N <- length(unique(group))

  # ---- build rhythm design matrices ------------------------------------
  message("fitting rhythmic models (random slopes on s1, c1)")
  models <- create_matrix_list(time, group, N, period)
  models <- lapply(models, function(l)
    l[, c(grep("u", colnames(l)), grep("a|b", colnames(l)))])
  for (i in seq_along(models)) rownames(models[[i]]) <- rownames(colData)

  if (length(unique(batch)) > 1) {
    model_b <- as.matrix(stats::model.matrix(~batch),
                         contrasts.arg = NULL)[, 2:length(unique(batch)), drop = FALSE]
    colnames(model_b) <- paste0("BATCH_", unique(batch)[-1])
    models <- lapply(models, function(l) cbind(model_b, l))
    models <- lapply(models, function(l)
      l[, c(grep("^u",     colnames(l)),
            grep("^BATCH", colnames(l)),
            grep("^a|^b",  colnames(l)))])
  }

  # ---- rhythm-side fits -------------------------------------------------
  fit <- papply(rownames(data), function(g) {
    x <- as.numeric(data[g, ])
    w <- if (!is.null(weights)) weights[g, ] else NULL
    dryRlmm:::.do_all_glmmTMB_full(x, my_mat = models, ID = ID,
                                   s1 = s1, c1 = c1, w = w)
  })
  names(fit) <- rownames(data)

  BIC <- t(sapply(fit, function(gf) vapply(gf, function(f) f$BIC, numeric(1))))
  rownames(BIC) <- names(fit)
  BIC  <- BIC[rownames(data), , drop = FALSE]
  BICW <- t(apply(BIC, 1, compute_BICW))
  chosen_model      <- apply(BIC,  1, which.min)
  chosen_model_BICW <- apply(BICW, 1, max)

  rstruct_rhythm <- vapply(seq_along(fit), function(g) {
    cm <- chosen_model[g]
    if (length(cm) == 0L || is.na(cm)) return(NA_character_)
    fit[[g]][[cm]]$rstruct %||% NA_character_
  }, character(1))
  names(rstruct_rhythm) <- names(fit)
  rstruct_rhythm <- rstruct_rhythm[rownames(data)]

  # ---- mean-side fits ---------------------------------------------------
  message("fitting mean models")
  model_mean_cond <- create_matrix_list_mean(N, group)
  model_mean_cond <- lapply(model_mean_cond, annotate_matrix, group)
  for (i in seq_along(model_mean_cond))
    rownames(model_mean_cond[[i]]) <- rownames(colData)

  gene.list <- as.list(rownames(data))
  fit_m <- papply(gene.list, function(g) {
    dryRlmm:::.do_all_glmmTMB_mr_full(
      x            = g,
      countData    = data,
      my_mat_r     = models,
      my_mat_m     = model_mean_cond,
      ID           = ID,
      s1           = s1,
      c1           = c1,
      chosen_model = chosen_model,
      weights_mat  = weights
    )
  })
  names(fit_m) <- rownames(data)

  BIC_mean <- t(sapply(fit_m, function(gf) vapply(gf, function(f) f$BIC, numeric(1))))
  rownames(BIC_mean) <- rownames(data)
  BICW_mean <- t(apply(BIC_mean, 1, compute_BICW))
  chosen_model_mean      <- apply(BIC_mean,  1, which.min)
  chosen_model_mean_BICW <- apply(BICW_mean, 1, max)

  rstruct_mean <- vapply(seq_along(fit_m), function(g) {
    cm <- chosen_model_mean[g]
    if (length(cm) == 0L || is.na(cm)) return(NA_character_)
    fit_m[[g]][[cm]]$rstruct %||% NA_character_
  }, character(1))
  names(rstruct_mean) <- rownames(data)

  # ---- group-level cosinor parameters ----------------------------------
  message("extracting rhythmic parameters")
  parameters <- lapply(seq_len(nrow(data)), function(i) {
    cm <- chosen_model_mean[i]
    if (length(cm) == 0L || is.na(cm)) return(rep(NA_real_, N * 6))
    dds <- fit_m[[i]][[cm]]$param
    compute_param_l_mixed(dds, period, N)
  })
  parameters <- data.frame(t(do.call(cbind, parameters)))
  colnames(parameters) <- paste(c("mean", "a", "b", "amp", "relamp", "phase"),
                                 rep(unique(group), each = 6), sep = "_")
  rownames(parameters) <- rownames(data)

  ncounts_RF <- data
  complete_parameters <- cbind(parameters,
                               chosen_model, chosen_model_BICW,
                               chosen_model_mean, chosen_model_mean_BICW,
                               rstruct_rhythm, rstruct_mean)
  global_table_df <- merge(
    data.frame(Row = rownames(ncounts_RF),       ncounts_RF,          check.names = FALSE),
    data.frame(Row = rownames(complete_parameters), complete_parameters, check.names = FALSE),
    by = "Row"
  )
  rownames(global_table_df) <- global_table_df$Row
  global_table_df <- global_table_df[, -1]

  # ---- per-subject parameters and variance summary ---------------------
  per_subject      <- NULL
  variance_summary <- NULL
  if (isTRUE(extract_per_subject)) {
    message("extracting per-subject parameters and variance summaries")
    uniq_groups   <- unique(as.character(group))
    per_subj_list <- lapply(seq_len(nrow(data)), function(i) {
      cm <- chosen_model_mean[i]
      if (length(cm) == 0L || is.na(cm)) return(NULL)
      .per_subject_from_fit(fit_m[[i]][[cm]], ID, group, period,
                            gene_name = rownames(data)[i])
    })
    var_list <- lapply(seq_len(nrow(data)), function(i) {
      cm <- chosen_model_mean[i]
      vc <- if (length(cm) && !is.na(cm)) fit_m[[i]][[cm]]$varcorr else NULL
      .variance_summary_one_gene(per_subj_list[[i]], vc, period,
                                 gene_name   = rownames(data)[i],
                                 uniq_groups = uniq_groups)
    })

    nonnull     <- !vapply(per_subj_list, is.null, logical(1))
    per_subject <- if (any(nonnull))
                     do.call(rbind, per_subj_list[nonnull])
                   else
                     data.frame()
    variance_summary <- do.call(rbind, var_list)
    rownames(variance_summary) <- variance_summary$gene
  }

  message("finished")
  list(
    time             = time,
    group            = group,
    ID               = ID,
    results          = global_table_df,
    BICW_rhythm      = BICW,
    BICW_mean        = BICW_mean,
    values           = ncounts_RF,
    parameters       = complete_parameters,
    per_subject      = per_subject,
    variance_summary = variance_summary,
    rstruct          = data.frame(rhythm = rstruct_rhythm, mean = rstruct_mean,
                                  row.names = rownames(data))
  )
}
