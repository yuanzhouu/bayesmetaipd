\name{sim1_as_formula_data}
\alias{sim1_as_formula_data}
\title{Convert Simulation 1 replicate 1 into formula-style inputs}
\description{
Builds an IPD data frame plus Type 1/2/3 AD tables from
\code{\link{sim1_ipdad_rep1}}.
}
\usage{
sim1_as_formula_data()
}
\value{
A list with \code{ipd}, formulas, AD tables, and MCMC starting values.
}
\examples{
d <- sim1_as_formula_data()
names(d$ad_nested)
}
