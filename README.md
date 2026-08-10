# bayesmetaipd

<p align="center">
  <img src="man/figures/logo.png" alt="bayesmetaipd logo" width="180"/>
</p>

Bayesian random-effects meta-analysis for **logistic** (Simulation Study 2) and
**continuous Application** outcomes using **individual participant data (IPD)**
and **aggregate data (AD)**.


| Function | What it does |
|----------|----------------|
| `fit_ipd()` | Logistic IPD-only Benchmark (Simulation 2) |
| `fit_ipd_ad()` | Logistic IPD + AD with density-ratio adjustment (Simulation 2) |
| `fit_ipd_gaussian()` | Application continuous IPD-only (study-specific \(\sigma_l^2\)) |
| `fit_ipd_ad_gaussian()` | Application continuous IPD + Type1/2/3 AD |
| `example_application_data()` | Tiny synthetic Application-style demo data |

> Real I-WIP Application CSVs are **not** redistributed (not public in the
> upstream repo). Pass your own prepared IPD / AD objects.

---

## Installation

```r
install.packages("remotes")
remotes::install_github("yuanzhouu/bayesmetaipd")
```


---

## Quick start

```r
library(bayesmetaipd)

# Logistic — Simulation 2 (bundled data defaults)
fit <- fit_ipd()
fit_ad <- fit_ipd_ad()

# Continuous — Application engine (supply your data)
toy <- example_application_data()
fit_g <- fit_ipd_gaussian(toy$ipd, burnin = 100, mainrun = 100, verbose = FALSE)
fit_gad <- fit_ipd_ad_gaussian(
  ipd = toy$ipd,
  ad_type1 = toy$ad_type1,
  ad_type2 = toy$ad_type2,
  ad_type3 = toy$ad_type3,
  burnin = 50, mainrun = 50, verbose = FALSE,
  tau_update_after = 10L
)
colMeans(fit_gad$posterior_mu)
```

---

## Application data shape

**IPD** — list of studies `list(y, X)` with Application coding
`(bmi_cat1, bmi_cat2, bmi_cat3, b_wt, bmi_cat1.trt, bmi_cat2.trt, bmi_cat3.trt)`,
or a data frame via `as_application_ipd()`.

**Type1 AD** — `data.frame(beta_hat, sqrt_hat_V)` bridged to coefficient 7.

**Type2 AD** — proportions `p_01…p_13` plus `beta_hat`, `sqrt_hat_V`.

**Type3 AD** — list with `betahat`/`se` (`J×6`), `use_drm`, `hat_tau`,
`hat_gamma_tau`, and lists `Design_X`, `nu_X`, `threshold`, `DRM_X`
(columns `int`, `DRM_X_type3`, `a_n`).

---

## Model (brief)

**Logistic (Sim2)** — Bernoulli–logit IPD, hierarchical \(\theta_l \sim N(\mu,\Sigma)\).

**Application continuous** — Gaussian IPD with study-specific \(\sigma_l^2\)
(truncated-normal prior), Gibbs \(\theta_l\) for IPD, MH bridges for Type1/2/3 AD
(Type3 uses normal-CDF excess-weight thresholds + optional density ratio).

See [hang-kim-stat/Bayesian-Meta](https://github.com/hang-kim-stat/Bayesian-Meta).

---

## License

MIT © Yuan Zhou
