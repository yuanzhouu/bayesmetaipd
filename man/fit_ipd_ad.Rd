\name{fit_ipd_ad}
\alias{fit_ipd_ad}
\title{Fit IPD + AD Bayesian random-effects meta-analysis (logistic)}
\description{
Implements the Simulation Study 2 proposed method
(\code{Code/Simulation_2/2_IPD-AD.R}): integrate IPD and AD under a logistic
random-effects model with density-ratio adjustment for covariate shift.

Defaults reproduce the official Simulation 2 IPD+AD \code{rep_1} result when
called as \code{fit_ipd_ad()} (bundled \code{sim2_ipdad_rep1} data,
\code{seed = 1001}, burn-in/main = 10000).
}
\usage{
fit_ipd_ad(
  X = NULL,
  Y = NULL,
  beta_mat = NULL,
  V_beta_cube = NULL,
  is_ipd = NULL,
  theta_init = NULL,
  mu_init = NULL,
  Sigma_init = NULL,
  burnin = 10000L,
  mainrun = 10000L,
  step_theta_ipd = 0.15,
  step_theta_ad = c(0.2, 0.1, 0.1, 0.2) * 1.5,
  step_alpha = 0.005,
  step_tau = 0.1,
  lambda = 1e4,
  nu0 = 0.1,
  phi0 = 0.1,
  seed = 1001L,
  verbose = TRUE
)
}
\arguments{
\item{X}{Optional \code{L x n x p} covariate array. If \code{NULL}, uses
\code{\link{sim2_ipdad_rep1}}.}

\item{Y}{Optional \code{L x n} binary response matrix.}

\item{beta_mat}{Optional \code{L x p} matrix of AD working-model estimates.}

\item{V_beta_cube}{Optional \code{L x p x p} array of AD working variances.}

\item{is_ipd}{Optional length-\code{L} logical/0-1 vector; \code{TRUE}/\code{1}
= IPD study.}

\item{theta_init}{Optional \code{L x p} starting values for study-specific
coefficients.}

\item{mu_init}{Optional length-\code{p} starting value for \code{mu}.}

\item{Sigma_init}{Optional \code{p x p} starting covariance.}

\item{burnin}{Integer burn-in iterations. Default \code{10000}.}

\item{mainrun}{Integer post-burn-in iterations retained. Default \code{10000}.}

\item{step_theta_ipd}{MH proposal SD for IPD \code{theta}. Default \code{0.15}.}

\item{step_theta_ad}{Length-\code{p} MH proposal SDs for AD \code{theta}.}

\item{step_alpha}{MH proposal SD for density-ratio \code{alpha}.}

\item{step_tau}{MH proposal SD for \code{tau}.}

\item{lambda}{Prior variance for \code{mu ~ N(0, lambda * I)}.}

\item{nu0}{Inverse-Wishart degrees of freedom.}

\item{phi0}{Inverse-Wishart scale multiplier.}

\item{seed}{Integer RNG seed. Default \code{1001}. With default data,
restores the official \code{.Random.seed}.}

\item{verbose}{Logical; print progress every 1000 iterations.}
}
\value{
A list of class \code{bayesmetaipd_fit} with \code{posterior_mu},
\code{posterior_Sigma_diag}, \code{call}, and \code{settings}.
}
\examples{
fit_short <- fit_ipd_ad(burnin = 5, mainrun = 10, verbose = FALSE)
colMeans(fit_short$posterior_mu)
}
