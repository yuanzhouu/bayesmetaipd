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
