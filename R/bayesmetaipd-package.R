# Bayesian IPD / IPD+AD meta-analysis
#
# Install from GitHub:
#   install.packages("remotes")
#   remotes::install_github("yuanzhouu/bayesmetaipd")
#
# Quick start:
#   library(bayesmetaipd)
#   fit <- fit_ipd()                       # logistic IPD-only (Sim2)
#   fit_ad <- fit_ipd_ad()                 # logistic IPD+AD (Sim2)
#   toy <- example_application_data()
#   fit_g <- fit_ipd_gaussian(toy$ipd)     # Application continuous IPD
#   fit_gad <- fit_ipd_ad_gaussian(        # Application continuous IPD+AD
#     toy$ipd, toy$ad_type1, toy$ad_type2, toy$ad_type3
#   )

#' @keywords internal
"_PACKAGE"
