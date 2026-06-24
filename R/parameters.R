# ---------------------------------------------------------------------------
# Parameter extraction from fitted dryRlmm models.
#
# - compute_param_l_mixed: parse dryR-style coefficient names ("u1", "a1",
#   "b1", "u2", "a2", "b2", ...) into a per-group (u, a, b, amp, relamp,
#   phase) row.
# - .parse_fixef_by_group:  parse the same naming into per-group named
#   numerics, used when reconstructing per-subject parameters from BLUPs.
# - .per_subject_from_fit:  add BLUPs to group fixed-effect parameters to
#   give one row per subject.
# - .variance_summary_one_gene: gene-level summary combining model-based
#   sigmas with descriptive subject-level SDs.
# ---------------------------------------------------------------------------

#' Convert a dryR-style fixed-effect vector to per-group cosinor parameters
#'
#' dryR encodes fixed-effect coefficients with names like \code{u1}, \code{a1},
#' \code{b1}, \code{u2}, ... where the digit is the group index. This helper
#' parses those names and returns a length-\code{6 * N} vector ordered as
#' \code{c(mean, a, b, amp, relamp, phase)} repeated per group.
#'
#' @param dds Named numeric vector (or 1-column matrix) of fixed-effect
#'   coefficients.
#' @param period Period of the rhythm in the same units as \code{phase}.
#'   Default 24 (hours).
#' @param N Number of groups.
#' @param verbose If \code{TRUE}, print parsing debug info.
#' @return Numeric vector of length \code{6 * N}.
#' @export
compute_param_l_mixed <- function(dds, period = 24, N, verbose = FALSE) {
  log <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  dds <- unlist(dds, use.names = TRUE)

  if (length(dds) == 0L) return(rep(NA_real_, N * 6))

  if (length(dds) == 1L) {
    uval <- as.numeric(dds)
    out  <- rep(NA_real_, N * 6)
    for (i in seq_len(N)) out[(1:6) + 6L * (i - 1L)] <- c(uval, NA, NA, NA, NA, NA)
    return(out)
  }

  if (is.null(dim(dds))) {
    nm  <- names(dds)
    mat <- cbind(value = as.numeric(dds))
    if (!is.null(nm)) rownames(mat) <- nm
    dds <- mat
  } else {
    dds <- as.matrix(dds); storage.mode(dds) <- "double"
  }

  if (is.null(rownames(dds))) return(rep(NA_real_, N * 6))
  rn <- trimws(rownames(dds))
  if (length(rn) == 0L) return(rep(NA_real_, N * 6))
  rownames(dds) <- rn

  pat   <- "(?i)([uab])\\W*(\\d+)"
  map_u <- vector("list", N); map_a <- vector("list", N); map_b <- vector("list", N)

  for (r in seq_along(rn)) {
    toks <- regmatches(rn[r], gregexpr(pat, rn[r], perl = TRUE))[[1]]
    if (!length(toks)) next
    syms <- tolower(sub(pat, "\\1", toks, perl = TRUE))
    idxs <- as.integer(sub(pat, "\\2", toks, perl = TRUE))
    keep <- !is.na(idxs) & idxs >= 1L & idxs <= N
    syms <- syms[keep]; idxs <- idxs[keep]
    for (k in seq_along(idxs)) {
      if      (syms[k] == "u") map_u[[idxs[k]]] <- r
      else if (syms[k] == "a") map_a[[idxs[k]]] <- r
      else if (syms[k] == "b") map_b[[idxs[k]]] <- r
    }
  }

  out <- rep(NA_real_, N * 6)
  for (i in seq_len(N)) {
    iu <- if (!is.null(map_u[[i]])) map_u[[i]] else integer(0)
    ia <- if (!is.null(map_a[[i]])) map_a[[i]] else integer(0)
    ib <- if (!is.null(map_b[[i]])) map_b[[i]] else integer(0)

    u <- if (length(iu)) dds[iu[1], 1] else NA_real_
    a <- if (length(ia)) dds[ia[1], 1] else NA_real_
    b <- if (length(ib)) dds[ib[1], 1] else NA_real_

    amp    <- 2 * sqrt(a^2 + b^2)
    relamp <- if (!is.na(u) && u != 0) 0.5 * amp / u else NA_real_
    phase  <- period / (2 * pi) * atan2(b, a)
    if (!is.na(phase)) phase <- phase %% period

    out[(1:6) + 6L * (i - 1L)] <- c(u, a, b, amp, relamp, phase)
  }
  out
}


# Parse fixed-effect names into per-group (u, a, b) named numerics.
.parse_fixef_by_group <- function(fixef_vec, uniq_groups) {
  N   <- length(uniq_groups)
  pat <- "(?i)([uab])\\W*(\\d+)"
  out <- list(u = stats::setNames(rep(NA_real_, N), uniq_groups),
              a = stats::setNames(rep(NA_real_, N), uniq_groups),
              b = stats::setNames(rep(NA_real_, N), uniq_groups))
  if (is.null(fixef_vec) || all(is.na(fixef_vec))) return(out)
  nm <- names(fixef_vec)
  if (is.null(nm)) return(out)
  for (k in seq_along(nm)) {
    matches <- regmatches(nm[k], gregexpr(pat, nm[k], perl = TRUE))[[1]]
    if (!length(matches)) next
    m   <- matches[1]
    sym <- tolower(sub(pat, "\\1", m, perl = TRUE))
    idx <- as.integer(sub(pat, "\\2", m, perl = TRUE))
    if (is.na(idx) || idx < 1L || idx > N) next
    gname <- uniq_groups[idx]
    val   <- unname(fixef_vec[k])
    if      (sym == "u") out$u[gname] <- val
    else if (sym == "a") out$a[gname] <- val
    else if (sym == "b") out$b[gname] <- val
  }
  out
}


# Build a per-subject data frame for one gene from its chosen fit object.
.per_subject_from_fit <- function(fit_obj, ID_vec, group_vec, period,
                                  gene_name = NA_character_) {
  if (is.null(fit_obj) || is.null(fit_obj$ranef_df)) return(NULL)
  rf <- fit_obj$ranef_df
  id_levels <- rownames(rf)
  if (is.null(id_levels) || !length(id_levels)) return(NULL)

  ID_chr    <- as.character(ID_vec)
  group_chr <- as.character(group_vec)
  id_group  <- vapply(id_levels, function(id) {
    idx <- which(ID_chr == id)[1]
    if (is.na(idx)) NA_character_ else group_chr[idx]
  }, character(1))

  uniq_groups <- unique(group_chr)
  parsed      <- .parse_fixef_by_group(fit_obj$param, uniq_groups)

  u_blup <- if ("(Intercept)" %in% colnames(rf)) rf[, "(Intercept)"] else rep(0, nrow(rf))
  a_blup <- if ("s1"          %in% colnames(rf)) rf[, "s1"]          else rep(0, nrow(rf))
  b_blup <- if ("c1"          %in% colnames(rf)) rf[, "c1"]          else rep(0, nrow(rf))

  u_fix <- parsed$u[id_group]; u_fix[is.na(u_fix)] <- 0
  a_fix <- parsed$a[id_group]; a_fix[is.na(a_fix)] <- 0
  b_fix <- parsed$b[id_group]; b_fix[is.na(b_fix)] <- 0

  u      <- u_fix + u_blup
  a      <- a_fix + a_blup
  b      <- b_fix + b_blup
  amp    <- 2 * sqrt(a^2 + b^2)
  phase  <- (period / (2 * pi) * atan2(b, a)) %% period
  relamp <- ifelse(!is.na(u) & u != 0, 0.5 * amp / u, NA_real_)

  data.frame(
    gene   = gene_name,
    ID     = id_levels,
    group  = id_group,
    mesor  = u,
    a      = a,
    b      = b,
    amp    = amp,
    relamp = relamp,
    phase  = phase,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# One-row variance summary for one gene from its per-subject data frame
# and its model-based variance components.
.variance_summary_one_gene <- function(per_subj_df, varcorr_vec, period,
                                       gene_name = NA_character_,
                                       uniq_groups = NULL) {

  base <- data.frame(gene = gene_name, stringsAsFactors = FALSE)

  vc <- if (is.null(varcorr_vec))
          c(sigma_intercept = NA_real_, sigma_s1 = NA_real_,
            sigma_c1        = NA_real_, sigma_resid = NA_real_)
        else varcorr_vec
  base$sigma_intercept <- unname(vc["sigma_intercept"])
  base$sigma_s1        <- unname(vc["sigma_s1"])
  base$sigma_c1        <- unname(vc["sigma_c1"])
  base$sigma_resid     <- unname(vc["sigma_resid"])

  if (is.null(per_subj_df) || !nrow(per_subj_df)) {
    if (!is.null(uniq_groups)) {
      for (g in uniq_groups) {
        base[[paste0("sd_mesor_",     g)]] <- NA_real_
        base[[paste0("sd_amp_",       g)]] <- NA_real_
        base[[paste0("circ_sd_phase_",g)]] <- NA_real_
        base[[paste0("R_",            g)]] <- NA_real_
      }
    }
    return(base)
  }

  if (is.null(uniq_groups)) uniq_groups <- unique(per_subj_df$group)

  for (g in uniq_groups) {
    sel <- per_subj_df$group == g & !is.na(per_subj_df$group)
    base[[paste0("sd_mesor_",      g)]] <- stats::sd(per_subj_df$mesor[sel], na.rm = TRUE)
    base[[paste0("sd_amp_",        g)]] <- stats::sd(per_subj_df$amp[sel],   na.rm = TRUE)
    base[[paste0("circ_sd_phase_", g)]] <- .circ_sd(per_subj_df$phase[sel], period)
    base[[paste0("R_",             g)]] <- .mean_resultant_length(per_subj_df$phase[sel], period)
  }
  base
}
