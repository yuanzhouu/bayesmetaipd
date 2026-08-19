# Compare fit_ipd_ad_lm R vs C++ on Simulation 1 formula data.
# From package root:
#   Rscript tools/compare_r_cpp_lm.R
# Also kept at paper_repo/tmp_code/compare_r_cpp_lm.R

suppressPackageStartupMessages(library(bayesmetaipd))

rel <- function(a, b) abs(a - b) / pmax(abs(a), 1e-12)
col_sd <- function(m) apply(m, 2, sd)
qfun <- function(x) quantile(x, c(0.025, 0.975), names = FALSE)

lm_args <- function() {
  d <- sim1_as_formula_data()
  list(
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
    sig2_init = d$sig2_init
  )
}

run_fit <- function(engine, burnin, mainrun, seed = 1001L) {
  args <- c(lm_args(), list(
    burnin = burnin, mainrun = mainrun, verbose = FALSE, seed = seed, engine = engine
  ))
  do.call(fit_ipd_ad_lm, args)
}

cat("package:", as.character(packageVersion("bayesmetaipd")), "\n")

fit_r_s <- run_fit("r", 2L, 3L)
fit_c_s <- run_fit("cpp", 2L, 3L)
cat("short identical mu?", isTRUE(all.equal(fit_r_s$posterior_mu, fit_c_s$posterior_mu)), "\n")
cat("short max |d mu|:", max(abs(fit_r_s$posterior_mu - fit_c_s$posterior_mu)), "\n")
cat("short max |d sig2|:", max(abs(fit_r_s$posterior_sig2 - fit_c_s$posterior_sig2)), "\n")
cat("short C++ finite?", all(is.finite(fit_c_s$posterior_mu)), all(fit_c_s$posterior_sig2 > 0), "\n")

n_burn <- 100L
n_main <- 300L
t0 <- proc.time()[[3]]
fit_r <- run_fit("r", n_burn, n_main)
t_r <- proc.time()[[3]] - t0
t0 <- proc.time()[[3]]
fit_c <- run_fit("cpp", n_burn, n_main)
t_c <- proc.time()[[3]] - t0

mu_r <- colMeans(fit_r$posterior_mu)
mu_c <- colMeans(fit_c$posterior_mu)
s2_r <- mean(fit_r$posterior_sig2)
s2_c <- mean(fit_c$posterior_sig2)
sdm_r <- col_sd(fit_r$posterior_mu)
sdm_c <- col_sd(fit_c$posterior_mu)
sds_r <- col_sd(fit_r$posterior_Sigma_diag)
sds_c <- col_sd(fit_c$posterior_Sigma_diag)

cat("\n=== posterior means (burnin=", n_burn, " main=", n_main, ") ===\n", sep = "")
print(rbind(R = mu_r, Cpp = mu_c, abs_diff = abs(mu_r - mu_c), rel_diff = rel(mu_r, mu_c)))
cat("sig2 mean R=", s2_r, " Cpp=", s2_c, " abs_diff=", abs(s2_r - s2_c),
    " rel=", rel(s2_r, s2_c), "\n")

cat("\n=== posterior SD of mu ===\n")
print(rbind(R = sdm_r, Cpp = sdm_c, abs_diff = abs(sdm_r - sdm_c), rel_diff = rel(sdm_r, sdm_c)))
cat("sig2 SD R=", sd(fit_r$posterior_sig2), " Cpp=", sd(fit_c$posterior_sig2),
    " rel=", rel(sd(fit_r$posterior_sig2), sd(fit_c$posterior_sig2)), "\n")

cat("\n=== posterior SD of Sigma_diag ===\n")
print(rbind(R = sds_r, Cpp = sds_c, rel_diff = rel(sds_r, sds_c)))

cat("\n=== 95% intervals for mu ===\n")
nm <- colnames(fit_r$posterior_mu)
for (j in seq_along(nm)) {
  qr <- qfun(fit_r$posterior_mu[, j])
  qc <- qfun(fit_c$posterior_mu[, j])
  cat(sprintf("  %s  R [%.3f, %.3f]  Cpp [%.3f, %.3f]  width_rel=%.3f\n",
              nm[j], qr[1], qr[2], qc[1], qc[2], rel(diff(qr), diff(qc))))
}

cat("\n=== timing (same ", n_burn + n_main, " total iters) ===\n", sep = "")
cat("R  :", round(t_r, 2), "s\n")
cat("C++:", round(t_c, 2), "s\n")
cat("speedup R/C++:", round(t_r / t_c, 2), "x\n")
cat("iters/sec R :", round((n_burn + n_main) / t_r, 2), "\n")
cat("iters/sec C++:", round((n_burn + n_main) / t_c, 2), "\n")
