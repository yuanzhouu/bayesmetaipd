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

The following self-contained example illustrates how to prepare data for each argument and fit the model:

```r
library(bayesmetaipd)

# -------------------------------------------------------------------------
# 1. IPD Data: Individual patient data across multiple studies (J = 5)
# -------------------------------------------------------------------------
set.seed(42)
ipd <- data.frame(
  study = rep(1:5, each = 100),
  X1    = rnorm(500),
  X2    = rbinom(500, size = 1, prob = 0.5)
)
ipd$Y <- 1.0 + 0.5 * ipd$X1 - 0.8 * ipd$X2 + 0.6 * (ipd$X1 * ipd$X2) + rnorm(500)

# -------------------------------------------------------------------------
# 2. Type 1 AD: Studies reporting reduced model coefficients (omitting X1:X2)
#    - Columns match terms in nested_formula and their SEs ('se_<term>')
#    - drm_mean & drm_var: published baseline mean and variance of the DRM covariate
# -------------------------------------------------------------------------
ad_nested <- data.frame(
  study    = 6:10,
  X1       = rnorm(5, mean =  0.5, sd = 0.05),
  X2       = rnorm(5, mean = -0.8, sd = 0.05),
  se_X1    = rep(0.08, 5),
  se_X2    = rep(0.08, 5),
  drm_mean = rnorm(5, mean = 0.0, sd = 0.1),
  drm_var  = rep(1.0 / 100, 5)
)

# -------------------------------------------------------------------------
# 3. Type 2 AD: Studies reporting subgroup outcome means and standard errors
#    - 'subgroup' defines covariate partition formulas
#    - 'ad_subgroup' columns match group names and standard errors ('se_<group>')
# -------------------------------------------------------------------------
subgroup_def <- list(
  g1 = ~ X1 > 0 & X2 == 0,
  g2 = ~ X1 <= 0 & X2 == 0
)

ad_subgroup <- data.frame(
  study    = 11:15,
  g1       = rnorm(5, mean = 1.2, sd = 0.1),
  g2       = rnorm(5, mean = 0.4, sd = 0.1),
  se_g1    = rep(0.12, 5),
  se_g2    = rep(0.12, 5),
  drm_mean = rnorm(5, mean = 0.0, sd = 0.1),
  drm_var  = rep(1.0 / 100, 5)
)

# -------------------------------------------------------------------------
# 4. Type 3 AD: Studies fitting full model but reporting only a subset of terms
#    - 'partial_terms' specifies which terms are published
#    - 'ad_partial' contains those estimates and their SEs ('se_<term>')
# -------------------------------------------------------------------------
ad_partial <- data.frame(
  study      = 16:20,
  X2         = rnorm(5, mean = -0.8, sd = 0.06),
  `X1:X2`    = rnorm(5, mean =  0.6, sd = 0.06),
  se_X2      = rep(0.09, 5),
  `se_X1:X2` = rep(0.09, 5),
  drm_mean   = rnorm(5, mean = 0.0, sd = 0.1),
  drm_var    = rep(1.0 / 100, 5),
  check.names = FALSE
)

# -------------------------------------------------------------------------
# 5. Fit the Bayesian Hierarchical Random-Effects Model
# -------------------------------------------------------------------------
fit <- fit_ipd_ad_lm(
  formula        = Y ~ X1 * X2,            # Full IPD regression specification
  ipd            = ipd,                    # Individual participant dataset
  study          = "study",                # Study identifier column in ipd
  nested_formula = ~ X1 + X2,              # Type 1 AD: reduced working formula
  ad_nested      = ad_nested,              # Type 1 AD data frame
  subgroup       = subgroup_def,           # Type 2 AD: partition formulas
  ad_subgroup    = ad_subgroup,            # Type 2 AD data frame
  partial_terms  = c("X2", "X1:X2"),       # Type 3 AD: reported subset of terms
  ad_partial     = ad_partial,             # Type 3 AD data frame
  drm_formula    = ~ X1,                   # Density-ratio model covariate
  burnin         = 1000, 
  mainrun        = 2000, 
  engine         = "cpp"                   # "cpp" (fast C++ sampler) or "r" (native R)
)

# Inspect results
print(fit)
colMeans(fit$posterior_mu)
```

For detailed mathematical specifications and simulation comparisons, see [`docs/fit_ipd_ad_lm_results.md`](docs/fit_ipd_ad_lm_results.md).

---

## License

MIT © Yuan Zhou
