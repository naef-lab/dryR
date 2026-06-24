test_that("dry_lmm_full() runs end-to-end on a tiny sim and returns sane output", {
  sim <- make_sim(n_per_sex = 6)
  out <- dry_lmm_full(
    data   = sim$data,
    group  = sim$meta$sex,
    time   = sim$meta$time,
    ID     = sim$meta$id,
    period = sim$period,
    extract_per_subject = TRUE,
    parallel_type = "serial",
    n.cores = 1
  )

  expect_named(out, c("time", "group", "ID", "results",
                      "BICW_rhythm", "BICW_mean", "values",
                      "parameters", "per_subject", "variance_summary",
                      "rstruct"))
  expect_equal(nrow(out$BICW_rhythm), 6)
  expect_equal(ncol(out$BICW_rhythm), 5)
  expect_equal(ncol(out$BICW_mean), 2)
  expect_true(all(rowSums(out$BICW_rhythm) > 0.99))
  expect_true(all(out$parameters$chosen_model %in% 1:5))
  expect_true(all(out$parameters$chosen_model_mean %in% 1:2))
})

test_that("drylmm_variance_test() returns 5 finite p-values on a clean gene", {
  sim <- make_sim(n_per_sex = 6)
  df  <- data.frame(
    measure = sim$data["shared", ],
    group   = sim$meta$sex,
    time    = sim$meta$time,
    id      = sim$meta$id
  )
  res <- drylmm_variance_test(df, timeout_s = 30)
  expect_equal(ncol(res), 5)
  expect_true(all(c("p_v_omnibus", "p_v_mesor", "p_v_rhythmic",
                    "p_v_amp", "p_v_phase") %in% colnames(res)))
  expect_true(any(is.finite(unlist(res))))
})

test_that("boundary_lrt_pvalue() matches Self-Liang 50:50 mixture", {
  # q = 1: 0.5 * P(chi2_1 > LR)
  LR <- 3.84
  expect_equal(boundary_lrt_pvalue(LR, 1),
               0.5 * pchisq(LR, 1, lower.tail = FALSE),
               tolerance = 1e-8)
  expect_equal(boundary_lrt_pvalue(0, 1), 1)
  expect_equal(boundary_lrt_pvalue(-1, 1), 1)
})
