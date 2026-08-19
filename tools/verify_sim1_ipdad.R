# Compare package fit_ipd_ad_sim1 vs official Code/Simulation_1/2_IPD-AD.R
# Short MCMC (burnin=2, mainrun=3), run official as Rscript so load() hits .GlobalEnv.

root <- "C:/Users/Yuan Zhou/OneDrive - University of Cincinnati/UCsemester/survey_sampling/paper_repo"
official <- file.path(root, "Bayesian-Meta/Code/Simulation_1/2_IPD-AD.R")
txt <- readLines(official, warn = FALSE)
txt <- gsub("DrawDiagnostics = TRUE", "DrawDiagnostics = FALSE", txt, fixed = TRUE)
txt <- gsub("burnin = 10000 ; mainrun = 10000", "burnin = 2 ; mainrun = 3", txt, fixed = TRUE)
txt <- gsub(
  "library(ModelMetrics) ; library(mvtnorm) ; library(invgamma) ; library(MCMCpack)",
  "library(mvtnorm); library(MCMCpack)",
  txt,
  fixed = TRUE
)
out_dir <- file.path(root, "Bayesian-Meta/Output/Simulation_1/2_IPD-AD_verify")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
txt <- gsub(
  'OutputFolder = "../../Output/Simulation_1/2_IPD-AD"',
  paste0("OutputFolder = ", deparse(out_dir)),
  txt,
  fixed = TRUE
)

tmp_script <- file.path(out_dir, "2_IPD-AD_short.R")
writeLines(txt, tmp_script)

oldwd <- setwd(file.path(root, "Bayesian-Meta/Code/Simulation_1"))
on.exit(setwd(oldwd), add = TRUE)
rscript <- "C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
status <- system2(rscript, shQuote(tmp_script), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")

official_path <- file.path(out_dir, "RData", "rep_1.RData")
if (!file.exists(official_path)) {
  stop("Official short run did not write ", official_path)
}
env <- new.env(parent = emptyenv())
load(official_path, envir = env)

devtools::load_all(file.path(root, "bayesmetaipd"), quiet = TRUE)
pkg <- fit_ipd_ad_sim1(burnin = 2, mainrun = 3, verbose = FALSE, seed = 1001L)

cmp <- list(
  mu_all_equal = isTRUE(all.equal(pkg$posterior_mu, env$posterior_mu)),
  sigma_all_equal = isTRUE(all.equal(pkg$posterior_Sigma_diag, env$posterior_Sigma_theta)),
  sig2_all_equal = isTRUE(all.equal(pkg$posterior_sig2, env$posterior_sig2)),
  mu_max_abs_diff = max(abs(pkg$posterior_mu - env$posterior_mu)),
  sigma_max_abs_diff = max(abs(pkg$posterior_Sigma_diag - env$posterior_Sigma_theta)),
  sig2_max_abs_diff = max(abs(pkg$posterior_sig2 - env$posterior_sig2))
)
print(cmp)
if (!all(c(cmp$mu_all_equal, cmp$sigma_all_equal, cmp$sig2_all_equal))) {
  stop("Package draws do not match official Simulation_1/2_IPD-AD.R")
}
cat("OK: package matches official short run\n")
