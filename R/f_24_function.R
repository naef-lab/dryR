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
  valid_data_indices = which(!is.na(x))
  x = x[valid_data_indices]
  t = t[valid_data_indices]
  n = length(x)
  
  if (n < 4) {
    return(c(Intercept = mean(x), amp = NA, phase = NA, s1 = NA, c1 = NA, pval = NA))
  }
  
  x_mean = mean(x)
  x_var = var(x)
  cos_component = cos(2 * pi * t / period)
  sin_component = sin(2 * pi * t / period)
  A = mean(x * cos_component) - x_mean * mean(cos_component)
  B = mean(x * sin_component) - x_mean * mean(sin_component)
  cos_sq_mean = mean(cos_component^2)
  sin_sq_mean = mean(sin_component^2)
  cos_sin_mean = mean(cos_component * sin_component)
  cos_mean = mean(cos_component)
  sin_mean = mean(sin_component)
  
  c1 = cos_sq_mean - cos_mean^2
  c2 = cos_sin_mean - cos_mean * sin_mean
  c3 = sin_sq_mean - sin_mean^2
  b = (A * c2 - B * c1) / (c2^2 - c1 * c3)
  a = (A - b * c2) / c1
  mu = x_mean - a * cos_mean - b * sin_mean
  
  x_hat = mu + a * cos_component + b * sin_component
  residual_var = var(x - x_hat)
  
  if (is.na(a) || is.na(b)) {
    return(c(Intercept = x_mean, amp = NA, phase = NA, s1 = NA, c1 = NA, pval = NA))
  }
  
  p = 3
  R2 = ifelse(x_var > 0, 1 - residual_var / x_var, 0)
  amp = max(x) - min(x)
  phase = period / (2 * pi) * atan2(b, a)
  phase = (phase + offset) %% period
  phase = ifelse(phase < 0, phase + period, phase)
  phase = ifelse(phase > period, phase - period, phase)
  pval = pbeta(R2, (p - 1) / 2, (n - p) / 2, lower.tail = FALSE, log.p = FALSE)
  
  return(c(Intercept = x_mean, amp = 2 * sqrt(a^2 + b^2), phase = phase, s1 = b, c1 = a, pval = pval))
}



# This function performs rhythmicity analysis for one condition for a dataframe that 
#' @param data	matrix or vector containing linear data; if a matrix is provided each column represents a sample, each row represents a feature.
#' @param time	vector containing numeric values of the zeitgeber/circadian time for each sample.
#' @param period	numeric value to indicate period length of the oscillation. Default: period = 24 h.
#' @param sample_name	vector containing sample names. Default: colnames are sample names.
#' @return a list that contains the following data.frames/vectors: results (summary of results), parameters (rhythmic parameters), time ((timepoints), period (period length)

f_24 <- function(data, time, period = 24, sample_name = names(data)) {
  
  if (is.vector(data)) {
    data = rbind(data, data)
    rownames(data) = c("X1", "X2")
  }
  
  res_tmp = apply(data, 1, function(x) f24_R2_cycling(x, t = time, period = period))
  res = t(res_tmp)
  padj = p.adjust(res[, "pval"], method = "BH")
  res_complete = cbind(res, padj)
  colnames(res_complete)[7] = "padj"
  
  global_df = as.data.frame(cbind(data, res_complete))
  
  out = list(time = time, period = period, results = global_df, parameters = res_complete)
  
  message("finished!")
  
  return(out)
}
