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
#   fit_s1 <- fit_ipd_ad_sim1()            # Gaussian IPD+AD Types 1/2/3 (Sim1)
#   d <- sim1_as_formula_data()
#   fit_lm <- fit_ipd_ad_lm(d$formula, d$ipd, d$study,
#     nested_formula = d$nested_formula, ad_nested = d$ad_nested,
#     subgroup = d$subgroup, ad_subgroup = d$ad_subgroup,
#     partial_terms = d$partial_terms, ad_partial = d$ad_partial)
#   toy <- example_application_data()
#   fit_g <- fit_ipd_gaussian(toy$ipd)     # Application continuous IPD
#   fit_gad <- fit_ipd_ad_gaussian(        # Application continuous IPD+AD
#     toy$ipd, toy$ad_type1, toy$ad_type2, toy$ad_type3
#   )

#' @keywords internal
"_PACKAGE"
