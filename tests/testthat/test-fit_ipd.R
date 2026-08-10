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