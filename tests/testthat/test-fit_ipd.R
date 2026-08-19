test_that("short fit_ipd run returns expected structure", {
  fit <- fit_ipd(burnin = 5, mainrun = 10, verbose = FALSE, seed = 1001L)
  expect_s3_class(fit, "bayesmetaipd_fit")
  expect_equal(dim(fit$posterior_mu), c(10, 4))
  expect_equal(dim(fit$posterior_Sigma_diag), c(10, 4))
  expect_true(is.finite(sum(fit$posterior_mu)))
})

test_that("sim2_rep1 data is available", {
  data(sim2_rep1, package = "bayesmetaipd")
  expect_equal(dim(sim2_rep1$X_cube), c(40, 200, 4))
  expect_equal(dim(sim2_rep1$Y_mat), c(40, 200))
})

test_that("short fit_ipd_ad run returns expected structure", {
  fit <- fit_ipd_ad(burnin = 2, mainrun = 3, verbose = FALSE, seed = 1001L)
  expect_s3_class(fit, "bayesmetaipd_fit")
  expect_equal(dim(fit$posterior_mu), c(3, 4))
  expect_equal(dim(fit$posterior_Sigma_diag), c(3, 4))
  expect_true(fit$settings$J >= 1)
  expect_true(fit$settings$K >= 1)
  expect_true(is.finite(sum(fit$posterior_mu)))
})

test_that("sim2_ipdad_rep1 data is available", {
  data(sim2_ipdad_rep1, package = "bayesmetaipd")
  expect_equal(dim(sim2_ipdad_rep1$X_cube), c(40, 200, 4))
  expect_equal(length(sim2_ipdad_rep1$delta_biased_access), 40)
  expect_true(all(sim2_ipdad_rep1$delta_biased_access %in% c(0, 1)))
})

test_that("example_application_data and IPD-only Application fit", {
  toy <- example_application_data(seed = 2L)
  expect_equal(length(toy$ipd), 8)
  expect_equal(ncol(toy$ipd[[1]]$X), 7)
  fit <- fit_ipd_gaussian(toy$ipd, burnin = 3, mainrun = 5, verbose = FALSE, seed = 117L)
  expect_s3_class(fit, "bayesmetaipd_fit")
  expect_equal(dim(fit$posterior_mu), c(5, 7))
  expect_equal(dim(fit$posterior_sig2), c(5, 8))
  expect_equal(fit$settings$outcome, "gaussian_application")
  expect_equal(fit$settings$J_type3, 0)
  expect_true(all(fit$posterior_sig2 > 0))
})

test_that("Application IPD+AD Type1/2/3 short run", {
  toy <- example_application_data(seed = 3L)
  fit <- fit_ipd_ad_gaussian(
    ipd = toy$ipd,
    ad_type1 = toy$ad_type1,
    ad_type2 = toy$ad_type2,
    ad_type3 = toy$ad_type3,
    burnin = 1, mainrun = 2, verbose = FALSE, seed = 117L,
    tau_update_after = 0L
  )
  expect_equal(dim(fit$posterior_mu), c(2, 7))
  expect_equal(fit$settings$J, 8)
  expect_equal(fit$settings$J_type1, 1)
  expect_equal(fit$settings$J_type2, 1)
  expect_equal(fit$settings$J_type3, 1)
  expect_equal(ncol(fit$posterior_sig2), 9) # 8 IPD + 1 Type3
  expect_true(is.finite(sum(fit$posterior_mu)))
})

test_that("sim1_ipdad_rep1 data is available", {
  data(sim1_ipdad_rep1, package = "bayesmetaipd")
  expect_equal(dim(sim1_ipdad_rep1$X_cube), c(40, 200, 4))
  expect_equal(dim(sim1_ipdad_rep1$Y_mat), c(40, 200))
  expect_equal(sim1_ipdad_rep1$type_vec, c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10)))
})

test_that("short fit_ipd_ad_sim1 run returns expected structure", {
  fit <- fit_ipd_ad_sim1(burnin = 1, mainrun = 2, verbose = FALSE, seed = 1001L)
  expect_s3_class(fit, "bayesmetaipd_fit")
  expect_equal(dim(fit$posterior_mu), c(2, 4))
  expect_equal(dim(fit$posterior_Sigma_diag), c(2, 4))
  expect_equal(length(fit$posterior_sig2), 2)
  expect_equal(fit$settings$outcome, "gaussian_sim1")
  expect_equal(fit$settings$J, 10)
  expect_equal(fit$settings$J_type1, 10)
  expect_equal(fit$settings$J_type2, 10)
  expect_equal(fit$settings$J_type3, 10)
  expect_true(all(fit$posterior_sig2 > 0))
  expect_true(is.finite(sum(fit$posterior_mu)))
})

test_that("short fit_ipd_ad_sim1 C++ engine returns expected structure", {
  fit <- fit_ipd_ad_sim1(
    burnin = 1, mainrun = 2, verbose = FALSE, seed = 1001L, engine = "cpp"
  )
  expect_s3_class(fit, "bayesmetaipd_fit")
  expect_equal(fit$settings$engine, "cpp")
  expect_equal(dim(fit$posterior_mu), c(2, 4))
  expect_equal(dim(fit$posterior_Sigma_diag), c(2, 4))
  expect_equal(length(fit$posterior_sig2), 2)
  expect_true(all(fit$posterior_sig2 > 0))
  expect_true(is.finite(sum(fit$posterior_mu)))
})

test_that("sim1_as_formula_data matches cube dimensions", {
  d <- sim1_as_formula_data()
  expect_equal(nrow(d$ipd), 10 * 200)
  expect_equal(nrow(d$ad_nested), 10)
  expect_equal(nrow(d$ad_subgroup), 10)
  expect_equal(nrow(d$ad_partial), 10)
  expect_true(all(c("X1", "X2", "drm_mean", "drm_var") %in% names(d$ad_nested)))
  expect_true(all(c("ind.1", "ind.2", "ind.3", "ind.4") %in% names(d$ad_subgroup)))
})

test_that("fit_ipd_ad_lm formula interface short run", {
  d <- sim1_as_formula_data()
  fit <- fit_ipd_ad_lm(
    formula = d$formula,
    ipd = d$ipd,
    study = d$study,
    nested_formula = d$nested_formula,
    ad_nested = d$ad_nested,
    nested_reported = d$nested_reported,
    subgroup = d$subgroup,
    ad_subgroup = d$ad_subgroup,
    partial_terms = d$partial_terms,
    ad_partial = d$ad_partial,
    drm_formula = d$drm_formula,
    theta_init_ipd = d$theta_init_ipd,
    theta_init_nested = d$theta_init_nested,
    theta_init_subgroup = d$theta_init_subgroup,
    theta_init_partial = d$theta_init_partial,
    mu_init = d$mu_init,
    sig2_init = d$sig2_init,
    burnin = 1, mainrun = 2, verbose = FALSE, seed = 2L
  )
  expect_s3_class(fit, "bayesmetaipd_fit")
  expect_equal(fit$settings$outcome, "gaussian_lm")
  expect_equal(dim(fit$posterior_mu), c(2, 4))
  expect_equal(colnames(fit$posterior_mu), c("(Intercept)", "X1", "X2", "X1:X2"))
  expect_equal(fit$settings$J, 10)
  expect_equal(fit$settings$J_type1, 10)
  expect_equal(fit$settings$J_type2, 10)
  expect_equal(fit$settings$J_type3, 10)
  expect_true(all(fit$posterior_sig2 > 0))
  expect_true(is.finite(sum(fit$posterior_mu)))
})

test_that("fit_ipd_ad_lm C++ engine short run", {
  d <- sim1_as_formula_data()
  fit <- fit_ipd_ad_lm(
    formula = d$formula,
    ipd = d$ipd,
    study = d$study,
    nested_formula = d$nested_formula,
    ad_nested = d$ad_nested,
    nested_reported = d$nested_reported,
    subgroup = d$subgroup,
    ad_subgroup = d$ad_subgroup,
    partial_terms = d$partial_terms,
    ad_partial = d$ad_partial,
    drm_formula = d$drm_formula,
    theta_init_ipd = d$theta_init_ipd,
    theta_init_nested = d$theta_init_nested,
    theta_init_subgroup = d$theta_init_subgroup,
    theta_init_partial = d$theta_init_partial,
    mu_init = d$mu_init,
    sig2_init = d$sig2_init,
    burnin = 1, mainrun = 2, verbose = FALSE, seed = 2L, engine = "cpp"
  )
  expect_equal(fit$settings$engine, "cpp")
  expect_equal(dim(fit$posterior_mu), c(2, 4))
  expect_equal(colnames(fit$posterior_mu), c("(Intercept)", "X1", "X2", "X1:X2"))
  expect_true(all(fit$posterior_sig2 > 0))
  expect_true(is.finite(sum(fit$posterior_mu)))
})