\name{reproduce_sim1_ipdad}
\alias{reproduce_sim1_ipdad}
\title{Reproduce Simulation Study 1 IPD+AD (rep_1) draws}
\description{
Wrapper around \code{\link{fit_ipd_ad_sim1}} with the official default
settings from \code{Code/Simulation_1/2_IPD-AD.R}.
}
\usage{
reproduce_sim1_ipdad(...)
}
\arguments{
\item{...}{Passed to \code{\link{fit_ipd_ad_sim1}} (e.g. \code{verbose = FALSE}).}
}
\value{
A \code{bayesmetaipd_fit} object.
}
\examples{
\dontrun{
fit <- reproduce_sim1_ipdad()
colMeans(fit$posterior_mu)
}
}
