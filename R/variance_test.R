# ---------------------------------------------------------------------------
# Boundary-corrected LRTs for differences in between-subject variance.
#
# Each test compares the null (pooled random-effect variances on intercept,
# S and C) against an alternative where one or more variance components
# are split by group. The omnibus is the 3-df joint test; the marginal
# tests target single components.
# ---------------------------------------------------------------------------

#' Boundary-corrected p-value for a likelihood-ratio statistic
#'
#' Self & Liang (1987) mixture: when the alternative hypothesis sits on the
#' parameter-space boundary (variance components fixed at zero under H0),
#' the asymptotic null distribution of the LRT is a 50:50 mixture of
#' chi-squared distributions with 0..q degrees of freedom.
#'
#' @param LR Observed \eqn{-2(\ell_0 - \ell_1)}.
#' @param q Number of variance components on the boundary under H0.
#' @return Numeric scalar p-value in \eqn{[0, 1]}.
#' @export
boundary_lrt_pvalue <- function(LR, q) {
  if (!is.finite(LR) || LR <= 0) return(1)
  weights <- choose(q, 0:q) / 2^q
  pvals <- sapply(0:q, function(k)
    if (k == 0) 0 else stats::pchisq(LR, df = k, lower.tail = FALSE))
  sum(weights * pvals)
}

#' Test whether between-subject random-effect variances differ between groups
#'
#' Boundary-corrected LRTs on a glmmTMB cosinor fit. Returns 5 p-values:
#'
#' \describe{
#'   \item{\code{p_v_omnibus}}{3-df joint test of all three random-effect
#'         variances (intercept, S, C) split by group. Natural omnibus.}
#'   \item{\code{p_v_mesor}}{1-df test of the random-intercept variance
#'         differing by group. Cleanest of the three.}
#'   \item{\code{p_v_rhythmic}}{2-df joint test on the random-slope variances
#'         (\eqn{\sigma_S, \sigma_C}). Right answer to "does rhythm spread
#'         differ by group" — amplitude and phase spread cannot be separately
#'         identified at the random-effect level.}
#'   \item{\code{p_v_amp}}{1-df test on \eqn{\sigma_C} alone (heuristic).}
#'   \item{\code{p_v_phase}}{1-df test on \eqn{\sigma_S} alone (heuristic).}
#' }
#'
#' @param df A long-format data.frame with columns:
#'   \describe{
#'     \item{\code{measure}}{numeric response (e.g. log-expression).}
#'     \item{\code{group}}{2-level factor (e.g. sex).}
#'     \item{\code{time}}{numeric hours within \code{period}.}
#'     \item{\code{id}}{subject identifier (factor or character).}
#'   }
#' @param period Rhythm period; default 24.
#' @param timeout_s Per-fit timeout in seconds (default 60).
#' @param drop_ZT24 If TRUE, drop \code{time == period} (treat as ZT0 dup).
#' @return 1-row data.frame with the five p-values.
#' @export
drylmm_variance_test <- function(df,
                                  period    = 24,
                                  timeout_s = 60,
                                  drop_ZT24 = FALSE) {

  na_out <- data.frame(p_v_omnibus  = NA_real_,
                       p_v_mesor    = NA_real_,
                       p_v_rhythmic = NA_real_,
                       p_v_amp      = NA_real_,
                       p_v_phase    = NA_real_)

  d <- data.frame(
    y   = as.numeric(df$measure),
    sex = factor(df$group),
    ZT  = as.numeric(df$time),
    id  = factor(paste0("p", df$id))
  )
  d <- d[stats::complete.cases(d), ]
  if (drop_ZT24) d <- d[d$ZT < period, ]
  if (nrow(d) < 10 || nlevels(droplevels(d$sex)) < 2) return(na_out)

  omega    <- 2 * pi / period
  d$S      <- sin(omega * d$ZT)
  d$C      <- cos(omega * d$ZT)
  d$id_sex <- factor(paste(d$sex, d$id, sep = ":"))

  ctrl <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500,
                                                 eval.max = 500))

  fit_one <- function(formula) {
    tryCatch(
      R.utils::withTimeout(
        suppressWarnings(
          glmmTMB::glmmTMB(formula, data = d, REML = FALSE, control = ctrl)
        ),
        timeout = timeout_s, onTimeout = "silent"
      ),
      error = function(e) NULL
    )
  }

  out <- na_out

  # NOTE: `diag(...)` inside the formula below is glmmTMB's covariance-
  # structure DSL (a diagonal RE covariance matrix), not the base R
  # diag() function. Don't namespace it.

  # NULL: pooled random-effect variances on (1, S, C)
  fit_null <- fit_one(
    y ~ sex + S + C + sex:S + sex:C + diag(1 + S + C | id)
  )
  if (is.null(fit_null)) return(out)
  ll_null <- as.numeric(stats::logLik(fit_null))

  # p_v_omnibus: all three variances differ by sex (3 df)
  fit_omnibus <- fit_one(
    y ~ sex + S + C + sex:S + sex:C +
      diag(0 + sex | id_sex) + diag(0 + sex:S + sex:C | id_sex)
  )
  if (!is.null(fit_omnibus)) {
    LR <- 2 * (as.numeric(stats::logLik(fit_omnibus)) - ll_null)
    if (is.finite(LR) && LR >= 0)
      out$p_v_omnibus <- boundary_lrt_pvalue(LR, q = 3)
  }

  # p_v_mesor: intercept variance differs by sex (1 df)
  fit_mesor <- fit_one(
    y ~ sex + S + C + sex:S + sex:C +
      diag(0 + sex | id_sex) + diag(0 + S + C | id)
  )
  if (!is.null(fit_mesor)) {
    LR <- 2 * (as.numeric(stats::logLik(fit_mesor)) - ll_null)
    if (is.finite(LR) && LR >= 0)
      out$p_v_mesor <- boundary_lrt_pvalue(LR, q = 1)
  }

  # p_v_rhythmic: joint (sigma_S, sigma_C) by sex (2 df)
  fit_rhythmic <- fit_one(
    y ~ sex + S + C + sex:S + sex:C +
      diag(1 | id) + diag(0 + sex:S + sex:C | id_sex)
  )
  if (!is.null(fit_rhythmic)) {
    LR <- 2 * (as.numeric(stats::logLik(fit_rhythmic)) - ll_null)
    if (is.finite(LR) && LR >= 0)
      out$p_v_rhythmic <- boundary_lrt_pvalue(LR, q = 2)
  }

  # p_v_amp: sigma_C alone (1 df, heuristic)
  fit_amp <- fit_one(
    y ~ sex + S + C + sex:S + sex:C +
      diag(1 + S | id) + diag(0 + sex:C | id_sex)
  )
  if (!is.null(fit_amp)) {
    LR <- 2 * (as.numeric(stats::logLik(fit_amp)) - ll_null)
    if (is.finite(LR) && LR >= 0)
      out$p_v_amp <- boundary_lrt_pvalue(LR, q = 1)
  }

  # p_v_phase: sigma_S alone (1 df, heuristic)
  fit_phase <- fit_one(
    y ~ sex + S + C + sex:S + sex:C +
      diag(1 + C | id) + diag(0 + sex:S | id_sex)
  )
  if (!is.null(fit_phase)) {
    LR <- 2 * (as.numeric(stats::logLik(fit_phase)) - ll_null)
    if (is.finite(LR) && LR >= 0)
      out$p_v_phase <- boundary_lrt_pvalue(LR, q = 1)
  }

  out
}
