\name{fit_ipd_ad_lm}
\alias{fit_ipd_ad_lm}
\title{Fit linear IPD + Type 1/2/3 AD from formulas}
\description{
General continuous-outcome IPD+AD sampler. Supply a full IPD formula, then
describe the three AD report types and provide their published estimates.
}
\usage{
fit_ipd_ad_lm(
  formula,
  ipd,
  study,
  nested_formula = NULL,
  ad_nested = NULL,
  nested_reported = NULL,
  subgroup = NULL,
  ad_subgroup = NULL,
  partial_terms = NULL,
  ad_partial = NULL,
  drm_formula = NULL,
  use_drm = TRUE,
  diagonal_V = TRUE,
  burnin = 10000L,
  mainrun = 10000L,
  step_theta = 0.2,
  step_alpha = 0.01,
  step_tau = 0.02,
  lambda = 1e4,
  nu0 = 0.1,
  phi0 = 0.1,
  theta_init_ipd = NULL,
  theta_init_nested = NULL,
  theta_init_subgroup = NULL,
  theta_init_partial = NULL,
  mu_init = NULL,
  Sigma_init = NULL,
  sig2_init = 1,
  seed = 1001L,
  verbose = TRUE,
  engine = c("r", "cpp")
)
}
\arguments{
\item{formula}{Full IPD model, e.g. \code{y ~ x1 * x2}.}

\item{ipd}{IPD data frame.}

\item{study}{Study id column in \code{ipd}.}

\item{nested_formula}{Type 1 working model, e.g. \code{~ x1 + x2}.}

\item{ad_nested}{Type 1 AD table: reported nested coefficients, \code{se_*} or
list-column \code{V}, plus \code{drm_mean}, \code{drm_var}.}

\item{nested_reported}{Names/indices of nested coefficients that were published.
Default: all except intercept.}

\item{subgroup}{Named list of Type 2 subgroup formulas.}

\item{ad_subgroup}{Type 2 AD table; columns match \code{names(subgroup)}.}

\item{partial_terms}{Full-model terms reported by Type 3 studies.}

\item{ad_partial}{Type 3 AD table for those terms.}

\item{drm_formula}{Density-ratio formula (default first covariate).}

\item{use_drm}{If \code{FALSE}, skip density-ratio MH updates.}

\item{diagonal_V}{If \code{TRUE}, keep only the diagonal of each AD covariance.}

\item{burnin, mainrun}{MCMC lengths.}

\item{step_theta, step_alpha, step_tau}{MH scales.}

\item{lambda, nu0, phi0}{Prior hyperparameters.}

\item{theta_init_ipd, theta_init_nested, theta_init_subgroup, theta_init_partial}{Optional starting \code{theta} matrices.}

\item{mu_init, Sigma_init, sig2_init}{Starting values.}

\item{seed}{RNG seed.}

\item{verbose}{Print progress.}

\item{engine}{\code{"r"} (default) R sampler; \code{"cpp"} faster Rcpp loop
(same model, not bit-identical).}
}
\value{
A \code{bayesmetaipd_fit} object.
}
\examples{
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
  burnin = 1, mainrun = 2, verbose = FALSE, seed = 1
)
colMeans(fit$posterior_mu)
}
