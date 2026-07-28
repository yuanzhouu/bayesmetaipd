\name{reproduce_sim2_ipdad}
\alias{reproduce_sim2_ipdad}
\title{Reproduce official Simulation 2 IPD+AD (rep_1) draws}
\description{
Convenience wrapper around \code{\link{fit_ipd_ad}} with the exact default
settings used in \code{Code/Simulation_2/2_IPD-AD.R}.
}
\usage{
reproduce_sim2_ipdad(..., compare_official = FALSE)
}
\arguments{
\item{...}{Passed to \code{\link{fit_ipd_ad}}.}

\item{compare_official}{If \code{TRUE}, compare to bundled official posterior draws.}
}
\value{
A \code{bayesmetaipd_fit} fit object. If \code{compare_official = TRUE}, an
attribute \code{comparison} is attached.
}
\examples{
\dontrun{
fit <- reproduce_sim2_ipdad(compare_official = TRUE)
attr(fit, "comparison")
}
}
