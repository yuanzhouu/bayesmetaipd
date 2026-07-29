# bayesmetaipd

Bayesian random-effects meta-analysis for **logistic** outcomes using
**individual participant data (IPD)** and **aggregate data (AD)**.

Default settings reproduce Simulation Study 2 from the supplementary materials of
[hang-kim-stat/Bayesian-Meta](https://github.com/hang-kim-stat/Bayesian-Meta).

| Function | What it does |
|----------|----------------|
| `fit_ipd()` | IPD-only Benchmark (all studies as IPD) |
| `fit_ipd_ad()` | Proposed IPD + AD analysis with density-ratio adjustment |
| `reproduce_sim2_benchmark()` | Run IPD Benchmark and optionally compare to official draws |
| `reproduce_sim2_ipdad()` | Run IPD+AD and optionally compare to official draws |

---

## Installation

```r
install.packages("remotes")
remotes::install_github("yuanzhouu/bayesmeta")
```

Requires R ≥ 4.1 and packages `mvtnorm`, `MCMCpack`.

---

## Quick start

```r
library(bayesmetaipd)

# IPD-only Benchmark (Simulation 2, rep_1 defaults)
fit <- fit_ipd()
colMeans(fit$posterior_mu)

# IPD + AD (Simulation 2 proposed method)
fit_ad <- fit_ipd_ad()
colMeans(fit_ad$posterior_mu)

# Optional: compare to bundled official posterior draws
fit <- reproduce_sim2_benchmark(compare_official = TRUE)
fit_ad <- reproduce_sim2_ipdad(compare_official = TRUE)
attr(fit_ad, "comparison")
```

Short demo (not the official long chain):

```r
fit_short <- fit_ipd(burnin = 50, mainrun = 100, verbose = FALSE)
fit_ad_short <- fit_ipd_ad(burnin = 5, mainrun = 10, verbose = FALSE)
print(fit_short)
```

> **Note.** Full default runs use 10,000 burn-in + 10,000 main iterations and can take several minutes (IPD) to about an hour (IPD+AD).

---

## Custom data

**IPD only** — covariate array `X` (`L × n × p`) and binary responses `Y` (`L × n`):

```r
fit <- fit_ipd(X = my_X, Y = my_Y, seed = 1)
```

**IPD + AD** — also supply AD summaries and an IPD indicator:

```r
fit <- fit_ipd_ad(
  X = my_X,
  Y = my_Y,
  beta_mat = my_beta,       # L x p working-model estimates
  V_beta_cube = my_V,       # L x p x p working variances
  is_ipd = my_is_ipd,       # length L; 1 = IPD, 0 = AD
  seed = 1
)
```

Bundled example datasets:

```r
data(sim2_rep1)        # IPD Benchmark data
data(sim2_ipdad_rep1)  # IPD + AD data (biased access)
```

---

## Model (brief)

Hierarchical logistic random-effects model:

$$
\begin{aligned}
y_{li} \mid X_{li}, \theta_l &\sim \mathrm{Bernoulli}\!\bigl(\mathrm{logit}^{-1}(X_{li}^\top \theta_l)\bigr) \\
\theta_l \mid \mu, \Sigma &\sim \mathcal{N}(\mu, \Sigma)
\end{aligned}
$$

- **IPD-only:** Metropolis–Hastings for each $\theta_l$, then Gibbs for $(\mu, \Sigma)$.
- **IPD+AD:** same random-effects layer, plus a moment bridge from AD working-model estimates and a density-ratio adjustment for covariate shift.

See also the original code and data in [hang-kim-stat/Bayesian-Meta](https://github.com/hang-kim-stat/Bayesian-Meta).

---

## Reproducibility

Official scripts do:

```r
set.seed(1001)
load("SimulationData_2.RData")  # written with save.image()
```

`save.image()` stores `.Random.seed`, so `load()` overwrites `set.seed(1001)`.  
When defaults are used with `seed = 1001`, this package restores that saved RNG state so results match the official `rep_1` posterior draws (up to floating-point noise).

---

## License

MIT © Yuan Zhou
