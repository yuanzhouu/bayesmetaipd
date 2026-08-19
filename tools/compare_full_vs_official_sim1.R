# Full-length compare: official Simulation_1/2_IPD-AD rep_1.RData
# (10000 posterior draws; official burnin=10000, mainrun=10000)
# vs package fit_ipd_ad_sim1 and fit_ipd_ad_lm, engines r and cpp.
#
# Re-run:
#   Rscript tmp_code/compare_full_vs_official_sim1.R

opt_lib <- "C:/Users/Yuan Zhou/AppData/Local/Temp/bayesmetaipd_optlib"
if (dir.exists(opt_lib)) .libPaths(c(opt_lib, .libPaths()))

suppressPackageStartupMessages(library(bayesmetaipd))

official_path <- "C:/Users/Yuan Zhou/OneDrive - University of Cincinnati/UCsemester/survey_sampling/paper_repo/Bayesian-Meta/Output/Simulation_1/2_IPD-AD/RData/rep_1.RData"
out_dir <- "C:/Users/Yuan Zhou/OneDrive - University of Cincinnati/UCsemester/survey_sampling/paper_repo/tmp_code"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rel <- function(a, b) abs(a - b) / pmax(abs(a), 1e-12)
col_sd <- function(m) apply(m, 2, sd)

off_env <- new.env(parent = emptyenv())
load(official_path, envir = off_env)
n_main <- nrow(off_env$posterior_mu)
n_burn <- n_main
stopifnot(n_main >= 1L, length(off_env$posterior_sig2) == n_main)

off <- list(
  mu_mean = colMeans(off_env$posterior_mu),
  mu_sd = col_sd(off_env$posterior_mu),
  sig2_mean = mean(off_env$posterior_sig2),
  sig2_sd = sd(off_env$posterior_sig2),
  Sigma_mean = colMeans(off_env$posterior_Sigma_theta),
  Sigma_sd = col_sd(off_env$posterior_Sigma_theta)
)

mu_names <- c("(Intercept)", "X1", "X2", "X1X2")
names(off$mu_mean) <- names(off$mu_sd) <- mu_names
names(off$Sigma_mean) <- names(off$Sigma_sd) <- mu_names

summarize_fit <- function(fit) {
  mu <- fit$posterior_mu
  if (!is.null(colnames(mu))) {
    # formula API uses X1:X2; cube uses unlabeled / X1X2
    colnames(mu) <- mu_names
  }
  sig <- fit$posterior_Sigma_diag
  colnames(sig) <- mu_names
  list(
    n = nrow(mu),
    mu_mean = setNames(colMeans(mu), mu_names),
    mu_sd = setNames(col_sd(mu), mu_names),
    sig2_mean = mean(fit$posterior_sig2),
    sig2_sd = sd(fit$posterior_sig2),
    Sigma_mean = setNames(colMeans(sig), mu_names),
    Sigma_sd = setNames(col_sd(sig), mu_names)
  )
}

cmp_one <- function(lab, sm, elapsed_sec) {
  list(
    label = lab,
    elapsed_sec = elapsed_sec,
    n = sm$n,
    mu_mean = sm$mu_mean,
    mu_sd = sm$mu_sd,
    sig2_mean = sm$sig2_mean,
    sig2_sd = sm$sig2_sd,
    Sigma_mean = sm$Sigma_mean,
    Sigma_sd = sm$Sigma_sd,
    d_mu_mean = sm$mu_mean - off$mu_mean,
    rel_mu_mean = rel(sm$mu_mean, off$mu_mean),
    d_mu_sd = sm$mu_sd - off$mu_sd,
    rel_mu_sd = rel(sm$mu_sd, off$mu_sd),
    d_sig2_mean = sm$sig2_mean - off$sig2_mean,
    rel_sig2_mean = rel(sm$sig2_mean, off$sig2_mean),
    d_sig2_sd = sm$sig2_sd - off$sig2_sd,
    rel_sig2_sd = rel(sm$sig2_sd, off$sig2_sd),
    d_Sigma_mean = sm$Sigma_mean - off$Sigma_mean,
    rel_Sigma_mean = rel(sm$Sigma_mean, off$Sigma_mean)
  )
}

print_block <- function(x) {
  cat("\n======== ", x$label, " ========\n", sep = "")
  cat(sprintf("elapsed: %.1f s (%.2f min) | draws: %d\n", x$elapsed_sec, x$elapsed_sec / 60, x$n))
  cat("posterior mean of mu vs official:\n")
  print(rbind(pkg = x$mu_mean, official = off$mu_mean, abs_diff = abs(x$d_mu_mean), rel_diff = x$rel_mu_mean))
  cat("posterior SD of mu vs official:\n")
  print(rbind(pkg = x$mu_sd, official = off$mu_sd, abs_diff = abs(x$d_mu_sd), rel_diff = x$rel_mu_sd))
  cat(sprintf(
    "sig2 mean pkg=%.6f official=%.6f abs=%.6f rel=%.4f\n",
    x$sig2_mean, off$sig2_mean, abs(x$d_sig2_mean), x$rel_sig2_mean
  ))
  cat(sprintf(
    "sig2 SD   pkg=%.6f official=%.6f abs=%.6f rel=%.4f\n",
    x$sig2_sd, off$sig2_sd, abs(x$d_sig2_sd), x$rel_sig2_sd
  ))
  cat("posterior mean of Sigma_diag vs official:\n")
  print(rbind(pkg = x$Sigma_mean, official = off$Sigma_mean, abs_diff = abs(x$d_Sigma_mean), rel_diff = x$rel_Sigma_mean))
}

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

run_one <- function(kind, engine) {
  cat(sprintf("\n[%s] starting %s engine=%s  burnin=%d mainrun=%d\n",
              format(Sys.time()), kind, engine, n_burn, n_main))
  flush.console()
  t0 <- proc.time()[[3]]
  if (identical(kind, "sim1")) {
    fit <- fit_ipd_ad_sim1(
      burnin = n_burn, mainrun = n_main, verbose = TRUE, seed = 1001L, engine = engine
    )
  } else {
    args <- c(lm_args(), list(
      burnin = n_burn, mainrun = n_main, verbose = TRUE, seed = 1001L, engine = engine
    ))
    fit <- do.call(fit_ipd_ad_lm, args)
  }
  elapsed <- proc.time()[[3]] - t0
  cat(sprintf("[%s] finished %s engine=%s in %.1f s\n", format(Sys.time()), kind, engine, elapsed))
  flush.console()
  list(fit = fit, elapsed = elapsed)
}

jobs <- list(
  c("sim1", "cpp"),
  c("lm", "cpp"),
  c("sim1", "r"),
  c("lm", "r")
)

results <- list()
partial_path <- file.path(out_dir, "full_engine_vs_official_partial.rds")
report_path <- file.path(out_dir, "full_engine_vs_official_report.txt")

zz <- file(report_path, open = "wt")
sink(zz, split = TRUE)
on.exit({ sink(); close(zz) }, add = TRUE)

cat("package:", as.character(packageVersion("bayesmetaipd")), "\n")
cat("has sim1_mcmc_cpp:", exists("sim1_mcmc_cpp", asNamespace("bayesmetaipd"), inherits = FALSE), "\n")
cat("has lm_mcmc_cpp:", exists("lm_mcmc_cpp", asNamespace("bayesmetaipd"), inherits = FALSE), "\n")
cat("official:", official_path, "\n")
cat("official draws:", n_main, " => package burnin=mainrun=", n_burn, "\n")
cat("official mu means:\n")
print(off$mu_mean)

for (job in jobs) {
  kind <- job[[1]]
  engine <- job[[2]]
  lab <- paste(kind, engine, sep = "_")
  run <- run_one(kind, engine)
  sm <- summarize_fit(run$fit)
  results[[lab]] <- cmp_one(lab, sm, run$elapsed)
  results[[lab]]$settings <- run$fit$settings
  saveRDS(results, partial_path)
  print_block(results[[lab]])
}

cat("\n======== timing summary (seconds) ========\n")
times <- vapply(results, `[[`, numeric(1), "elapsed_sec")
print(times)
if ("sim1_r" %in% names(results) && "sim1_cpp" %in% names(results)) {
  cat(sprintf("sim1 speedup R/C++: %.2fx\n", times[["sim1_r"]] / times[["sim1_cpp"]]))
}
if ("lm_r" %in% names(results) && "lm_cpp" %in% names(results)) {
  cat(sprintf("lm   speedup R/C++: %.2fx\n", times[["lm_r"]] / times[["lm_cpp"]]))
}

saveRDS(list(official = off, official_path = official_path, n_main = n_main, results = results),
        file.path(out_dir, "full_engine_vs_official.rds"))
cat("\nDONE_FULL_COMPARE\n")
cat("saved:", file.path(out_dir, "full_engine_vs_official.rds"), "\n")
