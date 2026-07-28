# bayesmetaipd
#
# Bayesian random-effects meta-analysis for logistic models using
# individual participant data (IPD) and aggregate data (AD).
# Defaults reproduce Simulation Study 2 results from
# https://github.com/hang-kim-stat/Bayesian-Meta
#
# Install
# -------
#
#   install.packages("remotes")
#   remotes::install_github("yuanzhouu/bayesmeta")
#
# Quick start: IPD-only Benchmark
# -------------------------------
#
#   library(bayesmetaipd)
#   fit <- fit_ipd()
#   fit <- reproduce_sim2_benchmark(compare_official = TRUE)
#
# Quick start: IPD + AD
# ---------------------
#
#   fit_ad <- fit_ipd_ad()
#   fit_ad <- reproduce_sim2_ipdad(compare_official = TRUE)
#   attr(fit_ad, "comparison")
#
# Short demo
# ----------
#
#   fit_short <- fit_ipd(burnin = 50, mainrun = 100, verbose = FALSE)
#   fit_ad_short <- fit_ipd_ad(burnin = 5, mainrun = 10, verbose = FALSE)
#
# Custom data
# -----------
#
#   # IPD only: X (L x n x p), Y (L x n)
#   fit <- fit_ipd(X = my_X, Y = my_Y, seed = 1)
#
#   # IPD + AD: also provide beta_mat, V_beta_cube, is_ipd
#   fit <- fit_ipd_ad(
#     X = my_X, Y = my_Y,
#     beta_mat = my_beta, V_beta_cube = my_V,
#     is_ipd = my_is_ipd, seed = 1
#   )
#
# Why default seed matches the official script
# --------------------------------------------
#
# Official code does set.seed(1001) then load(SimulationData_2.RData).
# That file was saved with save.image(), so load() restores .Random.seed.
# This package restores that state when seed = 1001 and bundled data are used.
#
# License: MIT
