# ---------------------------------------------------------------------------
# Internal utilities used across files
# ---------------------------------------------------------------------------

# Null-coalescing operator (so we don't pull in rlang)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Pick a parallel backend. Auto-detects fork vs PSOCK based on OS / RStudio.
.resolve_parallel_type <- function(parallel_type = c("auto", "fork",
                                                    "psock", "serial")) {
  parallel_type <- match.arg(parallel_type)
  if (parallel_type != "auto") return(parallel_type)
  is_rstudio <- isTRUE(Sys.getenv("RSTUDIO") == "1")
  is_mac     <- Sys.info()[["sysname"]] == "Darwin"
  is_windows <- .Platform$OS.type == "windows"
  if (is_windows || (is_mac && is_rstudio)) "psock" else "fork"
}

# Build a papply() closure that dispatches lapply / parLapply / mclapply.
# Caller owns the cluster lifecycle.
.make_papply <- function(parallel_type, n.cores, cl = NULL) {
  function(X, FUN, ...) {
    if (parallel_type == "serial" || n.cores <= 1L) {
      lapply(X, FUN, ...)
    } else if (parallel_type == "psock") {
      parallel::parLapply(cl, X, FUN, ...)
    } else {
      parallel::mclapply(X, FUN, mc.cores = n.cores, ...)
    }
  }
}

# Mean resultant length R: 0 = uniform, 1 = perfect sync. Phases in `period`.
.mean_resultant_length <- function(phases_hours, period) {
  ph <- phases_hours[!is.na(phases_hours)]
  if (length(ph) < 1L) return(NA_real_)
  rad <- 2 * pi * ph / period
  sqrt(sum(cos(rad))^2 + sum(sin(rad))^2) / length(ph)
}

# Circular SD in hours (sqrt(-2 ln R) on radians, converted back to hours).
.circ_sd <- function(phases_hours, period) {
  ph <- phases_hours[!is.na(phases_hours)]
  if (length(ph) < 2L) return(NA_real_)
  R <- .mean_resultant_length(ph, period)
  R <- min(max(R, 1e-12), 1 - 1e-12)
  sqrt(-2 * log(R)) * period / (2 * pi)
}
