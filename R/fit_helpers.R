# ---------------------------------------------------------------------------
# Per-gene, per-model glmmTMB fitting with a 3-step fallback ladder.
# ---------------------------------------------------------------------------
# The fallback ladder protects the per-gene pipeline against singular fits
# and non-convergence:
#
#   1) (1 + s1 + c1 | ID)                # full correlated random effects
#   2) (1|ID) + (0+s1|ID) + (0+c1|ID)    # diagonal (independent) random effects
#   3) (1 | ID)                          # intercept-only (original dryRlmm)
#
# Group-level rhythmicity differences remain identifiable: fixed effects
# capture group means of mesor/amplitude/phase, random effects capture
# individual deviations (mean zero across IDs).
# ---------------------------------------------------------------------------

# Extract per-ID random effects as a data.frame (rows = IDs, cols = terms).
.extract_ranef_df <- function(m) {
  rf <- tryCatch(glmmTMB::ranef(m)$cond$ID, error = function(e) NULL)
  if (is.null(rf)) return(NULL)
  as.data.frame(rf)
}

# Extract per-term random-effect SDs (population-level variability between
# IDs within a group). Returns a named numeric.
.extract_varcorr <- function(m) {
  vc <- tryCatch(glmmTMB::VarCorr(m)$cond$ID, error = function(e) NULL)
  out <- c(sigma_intercept = NA_real_, sigma_s1 = NA_real_,
           sigma_c1        = NA_real_, sigma_resid = NA_real_)
  if (is.null(vc)) {
    out["sigma_resid"] <- tryCatch(stats::sigma(m), error = function(e) NA_real_)
    return(out)
  }
  d      <- diag(as.matrix(vc))
  sigmas <- sqrt(pmax(d, 0))
  if ("(Intercept)" %in% names(sigmas))
    out["sigma_intercept"] <- unname(sigmas["(Intercept)"])
  if ("s1" %in% names(sigmas))
    out["sigma_s1"]        <- unname(sigmas["s1"])
  if ("c1" %in% names(sigmas))
    out["sigma_c1"]        <- unname(sigmas["c1"])
  out["sigma_resid"] <-
    tryCatch(stats::sigma(m), error = function(e) NA_real_)
  out
}

# TRUE if the fit converged and no variance component is on the boundary.
.is_glmmTMB_ok <- function(m) {
  if (inherits(m, "try-error") || is.null(m))            return(FALSE)
  conv <- tryCatch(m$sdr$pdHess, error = function(e) FALSE)
  if (is.null(conv) || !isTRUE(conv))                    return(FALSE)
  vc <- tryCatch(glmmTMB::VarCorr(m)$cond, error = function(e) NULL)
  if (is.null(vc))                                       return(FALSE)
  for (g in vc) {
    if (length(g) && any(diag(as.matrix(g)) < 1e-8))     return(FALSE)
  }
  TRUE
}

#' Fit one fixed-effect matrix with the random-effect fallback ladder
#'
#' Internal worker used by [dry_lmm_full()]. Tries the full
#' \code{(1 + s1 + c1 | ID)} random structure first; falls back to diagonal,
#' then to intercept-only on singular / non-converging fits. Returns BIC,
#' fixed-effect coefficients, per-subject BLUPs, variance components, and a
#' string tagging which random structure won.
#'
#' @param x Numeric response vector (length = nsamples).
#' @param matX Fixed-effect design matrix.
#' @param ID Subject identifier (length = nsamples).
#' @param s1,c1 Cosinor regressors (\code{sin / cos(2*pi*t/T)}).
#' @param w Optional weights vector (default \code{NULL} = unit weights).
#' @return List with BIC, param (fixed-effect coefficients), ranef_df,
#'   varcorr (named numeric of sigmas), and rstruct ("full" / "diag" /
#'   "intercept" / "failed").
#' @export
compute_glmmTMB_full <- function(x, matX, ID, s1, c1, w = NULL) {

  if (is.null(w)) w <- rep(1, length(x))

  # 1) full correlated random effects
  m <- try(suppressWarnings(
    glmmTMB::glmmTMB(x ~ 0 + matX + (1 + s1 + c1 | ID),
                     weights = w, REML = FALSE,
                     family  = stats::gaussian())
  ), silent = TRUE)
  rstruct <- "full"

  # 2) diagonal random effects
  if (!.is_glmmTMB_ok(m)) {
    m <- try(suppressWarnings(
      glmmTMB::glmmTMB(x ~ 0 + matX + (1 | ID) + (0 + s1 | ID) + (0 + c1 | ID),
                       weights = w, REML = FALSE,
                       family  = stats::gaussian())
    ), silent = TRUE)
    rstruct <- "diag"
  }

  # 3) intercept-only (original dryRlmm behaviour)
  if (!.is_glmmTMB_ok(m)) {
    m <- try(suppressWarnings(
      glmmTMB::glmmTMB(x ~ 0 + matX + (1 | ID),
                       weights = w, REML = FALSE,
                       family  = stats::gaussian())
    ), silent = TRUE)
    rstruct <- "intercept"
  }

  if (inherits(m, "try-error") || is.null(m)) {
    return(list(BIC      = Inf,
                param    = NA_real_,
                ranef_df = NULL,
                varcorr  = c(sigma_intercept = NA_real_,
                             sigma_s1        = NA_real_,
                             sigma_c1        = NA_real_,
                             sigma_resid     = NA_real_),
                rstruct  = "failed"))
  }

  b <- tryCatch(stats::BIC(m), error = function(e) Inf)
  if (!is.finite(b)) b <- Inf

  list(BIC      = b,
       param    = glmmTMB::fixef(m)$cond,
       ranef_df = .extract_ranef_df(m),
       varcorr  = .extract_varcorr(m),
       rstruct  = rstruct)
}

# Apply compute_glmmTMB_full() over a list of fixed-effect matrices for one
# response. Used by dry_lmm_full() for both rhythm and mean-side fits.
.do_all_glmmTMB_full <- function(x, my_mat, ID, s1, c1, w = NULL) {
  x <- as.numeric(x)
  lapply(my_mat, compute_glmmTMB_full, x = x, ID = ID, s1 = s1, c1 = c1, w = w)
}

# Mean-side variant: nest the rhythm side's chosen fixed-effect block into
# each mean-model matrix, then fit. `chosen_model[i]` indexes which rhythm
# model was selected for gene `i`.
.do_all_glmmTMB_mr_full <- function(x, countData, my_mat_r, my_mat_m,
                                    ID, s1, c1, chosen_model,
                                    weights_mat = NULL) {
  i         <- match(x, rownames(countData))
  gene_name <- x
  x         <- countData[x, ]
  M         <- my_mat_r[[chosen_model[i]]]
  gene_specific_mean_models <- lapply(my_mat_m, function(z)
    cbind(z, M[, -grep("u", colnames(M))]))
  x <- as.numeric(x)
  w <- if (!is.null(weights_mat)) weights_mat[gene_name, ] else NULL
  lapply(gene_specific_mean_models, compute_glmmTMB_full,
         x = x, ID = ID, s1 = s1, c1 = c1, w = w)
}
