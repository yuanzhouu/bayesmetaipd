devtools::load_all(
  "C:/Users/Yuan Zhou/OneDrive - University of Cincinnati/UCsemester/survey_sampling/paper_repo/bayesmetaipd",
  quiet = TRUE
)
d <- sim1_as_formula_data()
dat <- load_sim1_ipdad_rep1()

s <- 31
Xcube <- dat$X_cube[s, , ]
mf <- data.frame(Y = dat$Y_mat[s, ], X1 = dat$X_cube[s, , "X1"], X2 = dat$X_cube[s, , "X2"])
Xmm <- model.matrix(Y ~ X1 * X2, mf)
cat("IPD X max|cube - model.matrix| =", max(abs(Xcube - Xmm)), "\n")

# Same RNG stream as cube API (bundled .Random.seed)
assign(".Random.seed", d$random_seed, envir = .GlobalEnv)
fit_cube <- fit_ipd_ad_sim1(burnin = 5, mainrun = 10, verbose = FALSE, seed = 1001L)

assign(".Random.seed", d$random_seed, envir = .GlobalEnv)
fit_lm <- fit_ipd_ad_lm(
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
  burnin = 5,
  mainrun = 10,
  verbose = FALSE,
  seed = NULL
)

mu_diff <- max(abs(unname(fit_lm$posterior_mu) - fit_cube$posterior_mu))
sig_diff <- max(abs(unname(fit_lm$posterior_Sigma_diag) - fit_cube$posterior_Sigma_diag))
s2_diff <- max(abs(fit_lm$posterior_sig2 - fit_cube$posterior_sig2))
cat("short-run max|mu| =", mu_diff, "\n")
cat("short-run max|Sigma| =", sig_diff, "\n")
cat("short-run max|sig2| =", s2_diff, "\n")
cat("mu all.equal =", isTRUE(all.equal(unname(fit_lm$posterior_mu), fit_cube$posterior_mu)), "\n")
cat("cube mean mu =", paste(round(colMeans(fit_cube$posterior_mu), 4), collapse = ", "), "\n")
cat("lm   mean mu =", paste(round(colMeans(fit_lm$posterior_mu), 4), collapse = ", "), "\n")
cat("mean mu abs diff =", paste(round(abs(colMeans(unname(fit_lm$posterior_mu)) - colMeans(fit_cube$posterior_mu)), 4), collapse = ", "), "\n")
