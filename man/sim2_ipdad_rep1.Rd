\name{sim2_ipdad_rep1}
\alias{sim2_ipdad_rep1}
\docType{data}
\title{Simulation Study 2, replicate 1 (IPD + AD data)}
\description{
Bundled first replicate from Simulation Study 2, including biased-access
indicators and AD working-model summaries. Used by default in
\code{\link{fit_ipd_ad}} and \code{\link{reproduce_sim2_ipdad}}.
}
\format{
A list with components:
\describe{
  \item{X_cube}{Numeric array \code{L x n x p} of covariates.}
  \item{Y_mat}{Numeric matrix \code{L x n} of binary responses.}
  \item{theta_l_mat}{Numeric matrix \code{L x p} of true study-specific coefficients.}
  \item{beta_mat}{Numeric matrix \code{L x p} of AD working-model estimates.}
  \item{V_beta_cube}{Numeric array \code{L x p x p} of AD working variances.}
  \item{delta_biased_access}{Length-\code{L} indicator: \code{1}=IPD, \code{0}=AD.}
  \item{true_mu}{Numeric vector of length \code{p}.}
  \item{random_seed}{\code{.Random.seed} from \code{SimulationData_2.RData}.}
}
}
\source{
Extracted from \code{SimulationData_2.RData} in
\url{https://github.com/hang-kim-stat/Bayesian-Meta}.
}
\examples{
data(sim2_ipdad_rep1)
str(sim2_ipdad_rep1, max.level = 1)
}
\keyword{datasets}
