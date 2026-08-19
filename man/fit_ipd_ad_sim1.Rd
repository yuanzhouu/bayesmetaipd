\name{fit_ipd_ad_sim1}
\alias{fit_ipd_ad_sim1}
\title{Fit Simulation Study 1 IPD + AD (continuous / Gaussian)}
\description{
Implements \code{Code/Simulation_1/2_IPD-AD.R}: Gaussian linear IPD, three AD
report types (nested model, subgroup means, full model with partial
coefficients), and density-ratio adjustment for covariate shift.

Defaults use bundled \code{\link{sim1_ipdad_rep1}} (10 Type1 + 10 Type2 +
10 Type3 AD, 10 IPD). \code{seed = 1001} restores the generator
\code{.Random.seed} stored with that replicate.
}
\usage{
fit_ipd_ad_sim1(
  X = NULL,
  Y = NULL,
  beta_mat = NULL,
  V_beta_cube = NULL,
  type_vec = NULL,
  theta_init = NULL,
  mu_init = NULL,
  Sigma_init = NULL,
  sig2_init = NULL,
  burnin = 10000L,
  mainrun = 10000L,
  step_theta = 0.2,
  step_alpha = 0.01,
  step_tau = 0.02,
  lambda = 1e4,
  nu0 = 0.1,
  phi0 = 0.1,
  seed = 1001L,
  verbose = TRUE
)
}
\arguments{
\item{X}{Optional \code{L x n x p} covariate array. If \code{NULL}, uses
\code{\link{sim1_ipdad_rep1}}. Columns should be named \code{Int}, \code{X1},
\code{X2}, \code{X1X2}.}

\item{Y}{Optional \code{L x n} continuous response matrix.}

\item{beta_mat}{Optional \code{L x p} AD working-model estimates.}

\item{V_beta_cube}{Optional \code{L x p x p} AD working covariances.}

\item{type_vec}{Optional length-\code{L} integer vector: \code{1} Type1 AD,
\code{2} Type2 AD, \code{3} Type3 AD, \code{4} IPD.}

\item{theta_init}{Optional \code{L x p} starting coefficients.}

\item{mu_init}{Optional length-\code{p} starting \code{mu}.}

\item{Sigma_init}{Optional \code{p x p} starting covariance.}

\item{sig2_init}{Starting shared residual variance.}

\item{burnin}{Integer burn-in iterations. Default \code{10000}.}

\item{mainrun}{Integer post-burn-in iterations retained. Default \code{10000}.}

\item{step_theta}{MH SD for \code{theta} (IPD and AD). Default \code{0.2}.}

\item{step_alpha}{MH SD for density-ratio \code{alpha}. Default \code{0.01}.}

\item{step_tau}{MH SD for \code{tau}. Default \code{0.02}.}

\item{lambda}{Prior variance for \code{mu ~ N(0, lambda * I)}.}

\item{nu0}{Inverse-Wishart degrees of freedom.}

\item{phi0}{Inverse-Wishart scale multiplier.}

\item{seed}{Integer RNG seed. Default \code{1001}.}

\item{verbose}{Logical; print progress every 1000 iterations.}
}
\value{
A list of class \code{bayesmetaipd_fit} with \code{posterior_mu},
\code{posterior_Sigma_diag}, \code{posterior_sig2}, \code{call}, and
\code{settings}.
}
\examples{
fit_short <- fit_ipd_ad_sim1(burnin = 2, mainrun = 3, verbose = FALSE)
colMeans(fit_short$posterior_mu)
}
