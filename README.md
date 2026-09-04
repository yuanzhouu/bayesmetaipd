# bayesmetaipd

<p align="center">
  <img src="man/figures/logo.png" alt="bayesmetaipd logo" width="180"/>
</p>

Bayesian hierarchical random-effects meta-analysis combining **Individual Participant Data (IPD)** and **Aggregate Data (AD)** for continuous outcomes via `fit_ipd_ad_lm()`.

The package accommodates three primary aggregate data reporting paradigms:
1. **Type 1 AD (Nested Working Model)**: AD studies reporting estimates from a reduced or misspecified model (e.g., omitting interaction terms).
2. **Type 2 AD (Subgroup Means)**: AD studies reporting sample means and standard errors across partitions of the covariate space.
3. **Type 3 AD (Partial Full Model)**: AD studies fitting the full model but publishing only a subset of estimated coefficients.

To address potential population heterogeneity and covariate distribution shifts between AD and IPD studies, the framework incorporates a semi-parametric **Density-Ratio Model (DRM)** via exponential tilting, along with high-performance **C++ (Rcpp)** computational acceleration (~28x speedup).

---

## Installation

```r
install.packages("remotes")
remotes::install_github("yuanzhouu/bayesmetaipd")
```

---

## Quick Start

```r
library(bayesmetaipd)

# Fit Bayesian IPD + AD linear model via formula interface
fit <- fit_ipd_ad_lm(
  formula = Y ~ X1 * X2,            # Full IPD regression specification
  ipd = ipd_df,                     # Individual participant dataset
  study = "study",                  # Study identifier column
  nested_formula = ~ X1 + X2,       # Type 1 AD: reduced working model
  ad_nested = ad_nested_df,         # Type 1 AD estimates and variances
  subgroup = subgroup_list,         # Type 2 AD: subgroup definitions
  ad_subgroup = ad_subgroup_df,     # Type 2 AD estimates and variances
  partial_terms = c("X2", "X1:X2"), # Type 3 AD: reported subset of coefficients
  ad_partial = ad_partial_df,       # Type 3 AD estimates and variances
  drm_formula = ~ X1,               # Density-ratio model covariates
  burnin = 5000, 
  mainrun = 10000, 
  engine = "cpp"                    # "cpp" (fast C++ sampler) or "r" (native R)
)

print(fit)
summary(fit$posterior_mu)
```

---

## Key Features

- **Intuitive R Formula Interface**: Specify the full study regression model and AD reporting structures using standard R formulas.
- **Covariate Shift Adjustment**: Integrated Density-Ratio Modeling (DRM) matches moments between reference IPD and aggregate studies.
- **Dual-Engine Architecture**: Native R implementation for transparent verification and C++ (Rcpp / Armadillo) for fast MCMC sampling.
- **Detailed Validation Reports**: See [`docs/fit_ipd_ad_lm_results.md`](docs/fit_ipd_ad_lm_results.md) for complete mathematical derivation, benchmark comparisons across random seeds, and simulation tutorials with ridge plots.

---

## License

MIT © Yuan Zhou
