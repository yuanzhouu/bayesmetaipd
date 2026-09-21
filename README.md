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

The package includes the official simulation dataset from the paper (Simulation Study 1 replicate 1), ready to use via `load_example()` or `data("example_data")`.

### 1. Load Simulation Data

```r
library(bayesmetaipd)

# Load official Simulation Study 1 formula dataset
d <- load_example()
```

### 2. Inspect Data Structures

Print the structure and first few rows of the IPD and the three AD tables (including standard errors `se_coef_*` / `se_mean_*`):

```r
# --- 1. IPD: Individual participant data (one row per subject) ---
head(d$ipd, 3)
#>   study          Y         X1 X2
#> 1    31 -6.4342892 -0.9431838  0
#> 2    31 -0.8290351 -0.1559129  0
#> 3    31 -5.9490264 -1.1335461  0

# --- 2. Type 1 AD: Nested working model (estimates, SEs, and DRM summary stats) ---
head(d$ad_nested, 3)
#>   study   coef_X1    coef_X2 se_coef_X1 se_coef_X2    drm_mean     drm_var
#> 1     1 1.9897200 -0.1295904 0.09809379  0.1430922 -0.04094950 0.002663604
#> 2     2 1.5704165 -2.0140107 0.10785919  0.1480977 -0.22959586 0.002367529
#> 3     3 0.4885982  0.9193593 0.09977003  0.1362024  0.01449613 0.002335441

# --- 3. Type 2 AD: Subgroup means (sample means, SEs, and DRM summary stats) ---
head(d$ad_subgroup, 3)
#>   study mean_ind.1 mean_ind.2 mean_ind.3 mean_ind.4 se_mean_ind.1 se_mean_ind.2 se_mean_ind.3 se_mean_ind.4   drm_mean     drm_var
#> 1    11   1.244086   3.863963  0.6635194  0.7109921     0.1551659     0.1515569     0.3534072     0.3778082  0.7275264 0.002923500
#> 2    12   1.798221   2.594988  0.3522915 -0.1966817     0.1526580     0.1574297     0.2054263     0.2175203  0.2764594 0.002345003
#> 3    13   2.690740   3.975390  1.7192349  2.8468055     0.1491089     0.1439359     0.1242996     0.1275856 -0.1694563 0.002561947

# --- 4. Type 3 AD: Partial full model (subset of terms, SEs, and DRM summary stats) ---
head(d$ad_partial, 3)
#>   study   coef_X2  coef_X1:X2 se_coef_X2 se_coef_X1:X2   drm_mean     drm_var
#> 1    21 0.6994860  0.08862467  0.1514120     0.2066515  0.2079521 0.002469874
#> 2    22 2.3032548  0.32059343  0.1913503     0.1810330 -0.7834332 0.002529736
#> 3    23 0.4278566 -2.01219431  0.1872224     0.2184944  0.5393566 0.002233042
```

### 3. Fit the Model

Call `fit_ipd_ad_lm()` directly on the loaded dataset:

```r
# Fit the Bayesian hierarchical random-effects model (L = 40 studies: 10 IPD + 30 AD)
fit <- fit_ipd_ad_lm(
  formula         = d$formula,              # Full model: Y ~ X1 * X2
  ipd             = d$ipd,                  # Individual participant dataset
  study           = d$study,                # Study identifier column
  nested_formula  = d$nested_formula,       # Type 1 AD: nested working formula (~ X1 + X2)
  ad_nested       = d$ad_nested,            # Type 1 AD table (uses default non-intercept terms)
  subgroup        = d$subgroup,             # Type 2 AD: 4 subgroup partition formulas
  ad_subgroup     = d$ad_subgroup,          # Type 2 AD table
  partial_terms   = d$partial_terms,        # Type 3 AD: reported subset c("X2", "X1:X2")
  ad_partial      = d$ad_partial,           # Type 3 AD table
  drm_formula     = d$drm_formula,          # Density-ratio covariate (~ X1)
  burnin          = 1000, 
  mainrun         = 2000, 
  engine          = "cpp"                   # Accelerated C++ MCMC sampler
)

# Inspect posterior summary
print(fit)
colMeans(fit$posterior_mu)
```

---

## License

MIT © Yuan Zhou
