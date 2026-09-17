\name{fit_ipd_ad_lm}
\alias{fit_ipd_ad_lm}
\title{Fit gaussian linear model to IPD + Type 1/2/3 AD data}
\description{
Fit gaussian linear model to IPD + Type 1/2/3 AD data within a Bayesian
random-effects framework via estimating equations, multiplier bootstrap,
and density ratio models.
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
\item{formula}{Two-sided formula specifying the full target model for Individual
Participant Data (IPD), e.g., \code{Y ~ X1 + X2 + X1:X2} (or equivalently
\code{Y ~ X1 * X2}). This model represents the primary inferential target of
interest across all studies: for each study \eqn{i}, the individual-level
continuous outcome is modeled as \eqn{Y_{ij} = X_{ij}^\top \beta_i + \epsilon_{ij}}
with a study-specific coefficient vector \eqn{\beta_i} (including a random
intercept by default) and residual error \eqn{\epsilon_{ij} \sim \mathcal{N}(0, \sigma^2)},
where \eqn{\beta_i} follows a hierarchical random-effects distribution.}

\item{ipd}{Data frame of Individual Patient Data (Individual Participant Data),
containing individual-level records from all available IPD studies. It must
contain the continuous response variable (\code{Y}), all predictors/covariates
specified in \code{formula} (e.g., \code{X1}, \code{X2}), and a study
identifier column identifying which study each participant belongs to.}

\item{study}{Character string specifying the name of the study identifier column
in the \code{ipd} data frame (e.g., \code{"study"}).}

\item{nested_formula}{One-sided formula defining the Type 1 working model
(nested/reduced model), e.g., \code{~ X1 + X2}. Type 1 Aggregate Data (AD)
studies fit a reduced working model that omits certain terms from the full
target model (e.g., omitting the interaction term \code{X1:X2}). For
simplicity, the response variable \code{Y} on the left-hand side of the tilde
(\code{~}) is omitted. Ignored if \code{ad_nested} is \code{NULL}.}

\item{ad_nested}{Data frame of summary statistics from Type 1 AD studies (one
row per study). Columns must include: (1) an optional study identifier column
(e.g., \code{"study"}); (2) columns for reported point estimates of the nested
model coefficients matching \code{nested_reported} or terms in
\code{nested_formula} (e.g., columns named \code{"coef_X1"}, \code{"coef_X2"});
(3) standard errors in columns named \code{se_<term>} (e.g., \code{"se_X1"},
\code{"se_X2"} or \code{"se_coef_X1"}, \code{"se_coef_X2"}); and (4)
density-ratio summary columns \code{drm_mean} and \code{drm_var} containing the
sample mean and variance-of-the-mean of the baseline covariate specified in
\code{drm_formula}.}

\item{nested_reported}{Character vector (or integer indices) indicating which
coefficients from the Type 1 nested working model were published by the AD
studies. By default, all coefficients in \code{nested_formula} except the
intercept are assumed to be reported.}

\item{subgroup}{Named list of one-sided formulas defining the subgroup
indicators for Type 2 AD studies. Type 2 AD studies report sample subgroup
outcome means and standard errors across partitions of the covariate space
rather than regression coefficients. For example:
\code{list(g1 = ~ X1 > 0 & X2 == 0, g2 = ~ X1 <= 0 & X2 == 0, g3 = ~ X1 > 0 & X2 == 1, g4 = ~ X1 <= 0 & X2 == 1)}.
Each formula evaluates to a binary indicator on the participant-level
covariate space.}

\item{ad_subgroup}{Data frame of summary statistics from Type 2 AD studies (one
row per study). Columns must include: (1) an optional study identifier column
(e.g., \code{"study"}); (2) sample subgroup mean estimates, where column names
match \code{names(subgroup)} (e.g., columns named \code{"mean_g1"},
\code{"mean_g2"}); (3) corresponding standard errors in columns named
\code{se_<subgroup>} (e.g., \code{"se_g1"}, \code{"se_g2"} or
\code{"se_mean_g1"}, \code{"se_mean_g2"}); and (4) density-ratio summary
columns \code{drm_mean} and \code{drm_var} for baseline covariate shift
adjustment.}

\item{partial_terms}{Character vector (or integer indices) specifying which
full-model terms were published by Type 3 AD studies, e.g.,
\code{c("X2", "X1:X2")}. Type 3 AD studies fit the complete target model but
report only a subset of the estimated coefficients.}

\item{ad_partial}{Data frame of summary statistics from Type 3 AD studies (one
row per study). Columns must include: (1) an optional study identifier column
(e.g., \code{"study"}); (2) reported coefficient estimates for the terms in
\code{partial_terms} (e.g., columns named \code{"coef_X2"},
\code{"coef_X1:X2"}); (3) corresponding standard errors in columns named
\code{se_<term>} (e.g., \code{"se_X2"}, \code{"se_X1:X2"} or
\code{"se_coef_X2"}, \code{"se_coef_X1:X2"}); and (4) density-ratio summary
columns \code{drm_mean} and \code{drm_var} for baseline covariate shift
adjustment.}

\item{drm_formula}{One-sided formula specifying the baseline covariate used in
the semi-parametric Density-Ratio Model (DRM), e.g., \code{~ X1}. The DRM uses
exponential tilting to account for covariate shift (population heterogeneity)
between AD and IPD study populations based on published aggregate moments
(\code{drm_mean} and \code{drm_var}). If \code{NULL}, defaults to the first
non-intercept covariate in the full target model.}

\item{use_drm}{Logical; if \code{TRUE} (default), applies semi-parametric
density-ratio modeling to adjust for covariate shift between IPD and AD
populations. If \code{FALSE}, assumes homogeneous covariate distributions
across studies, fixing density-ratio weights to 1 and skipping
Metropolis-Hastings updates for DRM tilt parameters.}

\item{diagonal_V}{Logical; if \code{TRUE} (default), keeps only the diagonal
variances of each AD covariance matrix (assuming zero off-diagonal sampling
covariances).}

\item{burnin, mainrun}{Positive integers specifying MCMC sampling lengths:
\code{burnin} is the number of initial warm-up/burn-in iterations to discard,
and \code{mainrun} is the number of post-burn-in iterations retained for
posterior inference.}

\item{step_theta, step_alpha, step_tau}{Numeric proposal standard deviations for
random-walk Metropolis-Hastings steps: \code{step_theta} for study-specific
effect vectors (\eqn{\theta_l}), \code{step_alpha} for DRM tilt parameters
(\eqn{\alpha_l}), and \code{step_tau} for baseline covariate moments
(\eqn{\tau_l}).}

\item{lambda, nu0, phi0}{Prior hyperparameters for the Bayesian hierarchical
model: \code{lambda} is the variance multiplier for the population mean prior
\eqn{\mu \sim \mathcal{N}(0, \lambda I)}; \code{nu0} and \code{phi0} are the
degrees of freedom and scale multiplier for the Inverse-Wishart hyperprior on
between-study covariance \eqn{\Sigma \sim \text{Inv-Wishart}(\nu_0, \phi_0 I)}.}

\item{theta_init_ipd, theta_init_nested, theta_init_subgroup, theta_init_partial}{Optional
numeric matrices of starting values for study-specific parameters
(\eqn{\theta_l}) across IPD studies, Type 1 nested AD studies, Type 2 subgroup
AD studies, and Type 3 partial AD studies, respectively. Each matrix must have
rows equal to the number of studies in that category and columns equal to the
number of full-model terms.}

\item{mu_init, Sigma_init, sig2_init}{Optional starting values for the Markov
chain: \code{mu_init} is a numeric vector for the population mean effect
\eqn{\mu}; \code{Sigma_init} is a positive-definite matrix for between-study
covariance \eqn{\Sigma}; and \code{sig2_init} is a positive scalar for residual
error variance \eqn{\sigma^2}.}

\item{seed}{An integer random seed passed to \code{set.seed()} for exact
reproducibility of MCMC chains. If \code{NULL}, the current RNG state is left
unchanged.}

\item{verbose}{Logical; if \code{TRUE}, displays MCMC iteration progress and
diagnostic messages to the console.}

\item{engine}{Character string specifying the computation engine: \code{"r"}
(default) uses the pure R MCMC sampler; \code{"cpp"} uses an optimized C++
(Rcpp/RcppArmadillo) inner loop providing significant speedup (same
mathematical model, but not bit-identical due to floating-point differences).}
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
