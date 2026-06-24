# A tiny rhythmic / non-rhythmic / sex-different simulated dataset
# Two sexes, 6 subjects per sex, 4 timepoints per subject, period 24.
make_sim <- function(n_per_sex = 6, period = 24,
                      times = c(0, 6, 12, 18),
                      seed = 42) {
  set.seed(seed)
  subj_F <- paste0("F", seq_len(n_per_sex))
  subj_M <- paste0("M", seq_len(n_per_sex))

  meta <- do.call(rbind, lapply(c(subj_F, subj_M), function(s) {
    data.frame(id = s,
                sex = ifelse(startsWith(s, "F"), "Female", "Male"),
                time = times,
                stringsAsFactors = FALSE)
  }))
  meta$sample <- sprintf("%s_t%d", meta$id, meta$time)

  S <- sin(2 * pi * meta$time / period)
  C <- cos(2 * pi * meta$time / period)
  is_F <- meta$sex == "Female"

  # Six genes with known partitions
  gene_y <- function(label) {
    re_int <- rnorm(length(unique(meta$id)), 0, 0.3)
    names(re_int) <- unique(meta$id)
    e <- rnorm(nrow(meta), 0, 0.5)
    if (label == "flat") {
      y <- 5 + re_int[meta$id] + e
    } else if (label == "shared") {
      y <- 5 + 1.0 * C + 0.5 * S + re_int[meta$id] + e
    } else if (label == "female_only") {
      y <- 5 + ifelse(is_F, 1.0 * C + 0.5 * S, 0) + re_int[meta$id] + e
    } else if (label == "male_only") {
      y <- 5 + ifelse(!is_F, 1.2 * C + 0.4 * S, 0) + re_int[meta$id] + e
    } else if (label == "different_amp") {
      y <- 5 + ifelse(is_F, 1.5, 0.6) * C + 0.4 * S + re_int[meta$id] + e
    } else if (label == "mesor_diff") {
      y <- 5 + ifelse(is_F, 0, 1) + 0.8 * C + 0.3 * S + re_int[meta$id] + e
    }
    y
  }

  genes  <- c("flat", "shared", "female_only", "male_only",
              "different_amp", "mesor_diff")
  data   <- t(vapply(genes, gene_y, numeric(nrow(meta))))
  rownames(data) <- genes
  colnames(data) <- meta$sample

  list(data = data, meta = meta, period = period)
}
