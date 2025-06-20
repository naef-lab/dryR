# Test a time series for rhythmicity using a harmonic regression method
# This function has been used in PMID:24344304. Please cite this paper if used.
# 
# @param x Numeric vector to test for rhythmicity
# @param t Numeric vector of time (in hours)
# @param period Period of oscillation to test
# @param offset Phase offset (if needed)
# @return A vector with number of timepoints, mean, amplitude, relative amplitude, phase, and p-value of rhythmicity
# 
f24_R2_cycling <- function(x, t = 2 * (0:(length(x) - 1)), period = 24, offset = 0) {
  # remove NAs
  valid <- !is.na(x)
  x <- x[valid]
  t <- t[valid]
  n <- length(x)
  nb.timepoints <- length(valid)

  # handle insufficient data
  if (n < 4) {
    stats <- c(
      nb.timepoints = nb.timepoints,
      mean          = if (n > 0) mean(x) else NA,
      amp           = NA,
      relamp        = NA,
      phase         = NA,
      pval          = NA,
      c1            = NA,
      s1            = NA
    )
    return(stats)
  }

  # compute sine/cosine terms
  c <- cos(2 * pi * t / period)
  s <- sin(2 * pi * t / period)

  # regression coefficients
  A  <- cov(x, c)
  B  <- cov(x, s)
  C1 <- var(c)
  C2 <- cov(c, s)
  C3 <- var(s)

  b  <- (A * C2 - B * C1) / (C2^2 - C1 * C3)
  a  <- (A - b * C2) / C1
  mu <- mean(x) - a * mean(c) - b * mean(s)

  x_hat <- mu + a * c + b * s
  R2    <- if (var(x) > 0) 1 - var(x - x_hat) / var(x) else NA

  # amplitude and phase
  amp    <- 2 * sqrt(a^2 + b^2)
  phi    <- (period / (2 * pi)) * atan2(b, a)
  phase  <- (phi %% period + offset) %% period

  # p-value from beta distribution
  pval <- pbeta(R2, (3 - 1) / 2, (n - 3) / 2, lower.tail = FALSE)

  c(
    nb.timepoints = nb.timepoints,
    mean          = mean(x),
    amp           = amp,
    relamp        = amp / mu,
    phase         = phase,
    pval          = pval,
    c1            = a,
    s1            = b
  )
}



# This function performs rhythmicity analysis for one condition for a dataframe that 
#' @param data	matrix or vector containing linear data; if a matrix is provided each column represents a sample, each row represents a feature.
#' @param time	vector containing numeric values of the zeitgeber/circadian time for each sample.
#' @param period	numeric value to indicate period length of the oscillation. Default: period = 24 h.
#' @param sample_name	vector containing sample names. Default: colnames are sample names.
#' @return a list that contains the following data.frames/vectors: results (summary of results), parameters (rhythmic parameters), time ((timepoints), period (period length)
#
f_24 = function (data, time, period = 24, sample_name = names(data))
{
  if (is.vector(data)) {
    data = rbind(data, data)
    rownames(data) = c("X1", "X2")
  }
  res_tmp = apply(data, 1, function(x) f24_R2_cycling(x, t = time,
                                                      period = period))
  res = t(res_tmp)
  padj = p.adjust(res[, "pval"], method = "BH")
  res_complete = cbind(res, padj)
  colnames(res_complete)[ncol(res_complete)] = "padj"
  global_df = as.data.frame(cbind(data, res_complete))
  out = list(time = time, period = period, results = global_df,
             parameters = res_complete)
  message("finished!")
  return(out)
}
