# Simulation 1 official IPD+AD posterior vs `fit_ipd_ad_lm` (R / C++)

This note compares the **official Simulation Study 1 IPD+AD posterior** with the two engines of `fit_ipd_ad_lm()` (`engine = "r"` and `engine = "cpp"`). Official results are compared first, then the two engines are compared with each other. Every table reports **posterior means** and **posterior SDs** (MCMC sample SDs after burn-in).

Numbers come from the full-length run on 2026-08-19, stored in `tmp_code/full_engine_vs_official.rds`. Reproduction script: `tmp_code/compare_full_vs_official_sim1.R`. Extra C++ seeds are in `tmp_code/cpp_seeds_vs_official.rds`.

---

## 1. Setup

| Item | Detail |
|------|--------|
| Official draws | `Bayesian-Meta/Output/Simulation_1/2_IPD-AD/RData/rep_1.RData` (2025-06-08) |
| Official script | `Code/Simulation_1/2_IPD-AD.R`, `rep_no = 1` |
| Package API | `fit_ipd_ad_lm(..., engine = "r"` or `"cpp")` with `sim1_as_formula_data()` |
| Data | Replicate 1 of official `SimulationData_1.RData`. Bundled `sim1_ipdad_rep1` matches that replicate on $X$, $Y$, $\theta$, $\beta$, $V$ (difference at most $10^{-15}$) |
| Chain length | Official and package: **10,000 burn-in + 10,000 main**. Comparisons use the 10,000 post-burn-in draws only |
| RNG | Official: `set.seed(1001)` then `load` the data. Package with `seed = 1001`: restore `.Random.seed` saved at the end of data generation. **Starting states differ; draws are not required to match pointwise** |
| Runtime | R engine **1241.5 s** (20.7 min); C++ engine **44.0 s**; about **28 times** faster |

The official RData stores only population-level parameters:

$$
\mu = (\mu_0,\ \mu_{X_1},\ \mu_{X_2},\ \mu_{X_1 X_2}),
\qquad
\sigma^2,
\qquad
\mathrm{diag}(\Sigma)=(\Sigma_{11},\ \Sigma_{22},\ \Sigma_{33},\ \Sigma_{44}).
$$

Here $\mu$ is (intercept, $X_1$, $X_2$, interaction), $\sigma^2$ is the shared residual variance, and $\mathrm{diag}(\Sigma)$ is the four diagonal entries of the random-effects covariance (official object `posterior_Sigma_theta`). These $\Sigma_{jj}$ values are the posterior of those diagonal entries, not a posterior SD matrix for $\Sigma$.

The official posterior mean of the intercept is near zero, so **relative error is inflated**. Tables report absolute differences. Differences should be read against each parameter’s own posterior SD.

---

## 2. Official vs both `fit_ipd_ad_lm` engines

### 2.1 Posterior means

| Parameter | Official | lm R | lm C++ | abs(R - official) | abs(C++ - official) |
|-----------|----------|------|--------|-------------------|---------------------|
| $\mu_0$ (intercept) | 0.00015 | 0.08108 | 0.05275 | 0.08093 | 0.05260 |
| $\mu_{X_1}$ | 1.17992 | 1.16068 | 1.18460 | 0.01924 | 0.00468 |
| $\mu_{X_2}$ | 0.47457 | 0.46707 | 0.46431 | 0.00750 | 0.01026 |
| $\mu_{X_1 X_2}$ | 0.35727 | 0.34832 | 0.32482 | 0.00896 | 0.03246 |
| $\sigma^2$ | 0.99468 | 0.99524 | 0.99505 | 0.00056 | 0.00037 |
| $\Sigma_{11}$ | 1.26030 | 1.28754 | 1.53343 | 0.02724 | 0.27313 |
| $\Sigma_{22}$ | 0.77485 | 0.78087 | 0.79131 | 0.00602 | 0.01646 |
| $\Sigma_{33}$ | 0.91755 | 0.92258 | 0.94071 | 0.00502 | 0.02316 |
| $\Sigma_{44}$ | 1.24065 | 1.20445 | 1.26683 | 0.03620 | 0.02619 |

Slope components of $\mu$ and $\sigma^2$ are close to the official posterior (absolute gaps much smaller than the posterior SDs). The intercept differs by about 0.05–0.08, still smaller than its posterior SD (about 0.22–0.24). At `seed = 1001`, the largest gap versus official is the C++ posterior mean of $\Sigma_{11}$ (about 0.27 too high); Section 4 shows that this gap depends on the C++ seed.

### 2.2 Posterior SDs

| Parameter | Official | lm R | lm C++ | abs(R - official) | abs(C++ - official) |
|-----------|----------|------|--------|-------------------|---------------------|
| $\mu_0$ | 0.22278 | 0.24341 | 0.23317 | 0.02063 | 0.01039 |
| $\mu_{X_1}$ | 0.17170 | 0.17240 | 0.17820 | 0.00070 | 0.00651 |
| $\mu_{X_2}$ | 0.15343 | 0.15411 | 0.15670 | 0.00068 | 0.00327 |
| $\mu_{X_1 X_2}$ | 0.19474 | 0.19049 | 0.19655 | 0.00425 | 0.00181 |
| $\sigma^2$ | 0.03162 | 0.03195 | 0.03197 | 0.00032 | 0.00035 |
| $\Sigma_{11}$ | 0.46643 | 0.54234 | 0.62089 | 0.07591 | 0.15446 |
| $\Sigma_{22}$ | 0.23761 | 0.25567 | 0.25578 | 0.01806 | 0.01817 |
| $\Sigma_{33}$ | 0.23583 | 0.24372 | 0.24789 | 0.00788 | 0.01206 |
| $\Sigma_{44}$ | 0.40704 | 0.36270 | 0.39590 | 0.04434 | 0.01114 |

Posterior SDs for $\mu$ and $\sigma^2$ match the official chain closely. $\mathrm{diag}(\Sigma)$ is heavier-tailed, so SD gaps are larger. At `seed = 1001`, the C++ posterior SD of $\Sigma_{11}$ (0.62) is above official (0.47) and the R engine (0.54). Other C++ seeds bring this SD back near 0.47–0.55 (Section 4).

---

## 3. R engine vs C++ engine

Same formula API, same data, same `seed = 1001`, same 10,000+10,000. C++ draws $\mu$ with a Cholesky `rmvnorm`; R uses the `mvtnorm` default eigen decomposition. Trajectories are not pointwise identical.

### 3.1 Posterior means

| Parameter | lm R | lm C++ | abs(C++ - R) |
|-----------|------|--------|--------------|
| $\mu_0$ | 0.08108 | 0.05275 | 0.02832 |
| $\mu_{X_1}$ | 1.16068 | 1.18460 | 0.02392 |
| $\mu_{X_2}$ | 0.46707 | 0.46431 | 0.00276 |
| $\mu_{X_1 X_2}$ | 0.34832 | 0.32482 | 0.02350 |
| $\sigma^2$ | 0.99524 | 0.99505 | 0.00019 |
| $\Sigma_{11}$ | 1.28754 | 1.53343 | 0.24589 |
| $\Sigma_{22}$ | 0.78087 | 0.79131 | 0.01044 |
| $\Sigma_{33}$ | 0.92258 | 0.94071 | 0.01813 |
| $\Sigma_{44}$ | 1.20445 | 1.26683 | 0.06238 |

$\sigma^2$ is essentially the same. Components of $\mu$ differ by about 0.003–0.028 (roughly 0.02–0.12 posterior SDs). At `seed = 1001`, the mean of $\Sigma_{11}$ differs by 0.25, the largest discrepancy between engines (see Section 4).

### 3.2 Posterior SDs

| Parameter | lm R | lm C++ | abs(C++ - R) |
|-----------|------|--------|--------------|
| $\mu_0$ | 0.24341 | 0.23317 | 0.01024 |
| $\mu_{X_1}$ | 0.17240 | 0.17820 | 0.00580 |
| $\mu_{X_2}$ | 0.15411 | 0.15670 | 0.00260 |
| $\mu_{X_1 X_2}$ | 0.19049 | 0.19655 | 0.00606 |
| $\sigma^2$ | 0.03195 | 0.03197 | 0.00003 |
| $\Sigma_{11}$ | 0.54234 | 0.62089 | 0.07855 |
| $\Sigma_{22}$ | 0.25567 | 0.25578 | 0.00011 |
| $\Sigma_{33}$ | 0.24372 | 0.24789 | 0.00418 |
| $\Sigma_{44}$ | 0.36270 | 0.39590 | 0.03320 |

Posterior SDs for $\mu$ and $\sigma^2$ agree closely. The posterior SD of $\Sigma_{11}$ differs by about 0.08; the other diagonal entries are close.

---

## 4. C++ engine at other seeds vs official

Sections 2–3 use `seed = 1001` for both package engines. The C++ $\Sigma_{11}$ mean at that seed (1.53 vs official 1.26) looks large. The same C++ engine was rerun at `seed = 42`, `2024`, and `7`, still with 10,000 burn-in + 10,000 main, same data and formula inputs. Draws: `tmp_code/cpp_seeds_vs_official.rds`.

`fit_ipd_ad_lm` uses `set.seed(seed)` (not the official `set.seed(1001)` then `load` protocol), so none of these C++ chains is meant to match the official trajectory pointwise.

### 4.1 Posterior means vs official

| Parameter | Official | C++ 1001 | C++ 42 | C++ 2024 | C++ 7 |
|-----------|----------|----------|--------|----------|-------|
| $\mu_0$ | 0.00015 | 0.05275 | 0.02046 | 0.12717 | −0.02095 |
| $\mu_{X_1}$ | 1.17992 | 1.18460 | 1.15744 | 1.14058 | 1.18415 |
| $\mu_{X_2}$ | 0.47457 | 0.46431 | 0.47212 | 0.46515 | 0.46660 |
| $\mu_{X_1 X_2}$ | 0.35727 | 0.32482 | 0.34232 | 0.35786 | 0.30663 |
| $\sigma^2$ | 0.99468 | 0.99505 | 0.99449 | 0.99461 | 0.99440 |
| $\Sigma_{11}$ | 1.26030 | 1.53343 | 1.32923 | 1.28625 | 1.30506 |
| $\Sigma_{22}$ | 0.77485 | 0.79131 | 0.80922 | 0.79309 | 0.78800 |
| $\Sigma_{33}$ | 0.91755 | 0.94071 | 0.92675 | 0.92457 | 0.94281 |
| $\Sigma_{44}$ | 1.24065 | 1.26683 | 1.15687 | 1.23168 | 1.16818 |

Absolute difference from official:

| Parameter | C++ 1001 | C++ 42 | C++ 2024 | C++ 7 |
|-----------|----------|--------|----------|-------|
| $\mu_0$ | 0.05260 | 0.02031 | 0.12701 | 0.02110 |
| $\mu_{X_1}$ | 0.00468 | 0.02248 | 0.03934 | 0.00423 |
| $\mu_{X_2}$ | 0.01026 | 0.00245 | 0.00942 | 0.00797 |
| $\mu_{X_1 X_2}$ | 0.03246 | 0.01496 | 0.00059 | 0.05064 |
| $\sigma^2$ | 0.00037 | 0.00020 | 0.00007 | 0.00029 |
| $\Sigma_{11}$ | **0.27314** | 0.06893 | **0.02595** | 0.04476 |
| $\Sigma_{22}$ | 0.01646 | 0.03437 | 0.01824 | 0.01314 |
| $\Sigma_{33}$ | 0.02316 | 0.00919 | 0.00702 | 0.02526 |
| $\Sigma_{44}$ | 0.02619 | 0.08377 | 0.00897 | 0.07247 |

At `seed = 2024`, $\Sigma_{11}$ is as close to official as the R engine was (abs 0.026). `seed = 1001` is the outlier among these four C++ runs.

### 4.2 Posterior SDs vs official

| Parameter | Official | C++ 1001 | C++ 42 | C++ 2024 | C++ 7 |
|-----------|----------|----------|--------|----------|-------|
| $\mu_0$ | 0.22278 | 0.23317 | 0.23461 | 0.25121 | 0.22445 |
| $\mu_{X_1}$ | 0.17170 | 0.17820 | 0.17650 | 0.18139 | 0.16687 |
| $\mu_{X_2}$ | 0.15343 | 0.15670 | 0.15359 | 0.15397 | 0.15493 |
| $\mu_{X_1 X_2}$ | 0.19474 | 0.19655 | 0.18891 | 0.19255 | 0.18483 |
| $\sigma^2$ | 0.03162 | 0.03197 | 0.03199 | 0.03164 | 0.03196 |
| $\Sigma_{11}$ | 0.46643 | 0.62089 | 0.53374 | 0.46896 | 0.54946 |
| $\Sigma_{22}$ | 0.23761 | 0.25578 | 0.28277 | 0.26104 | 0.25660 |
| $\Sigma_{33}$ | 0.23583 | 0.24789 | 0.24278 | 0.24309 | 0.24729 |
| $\Sigma_{44}$ | 0.40704 | 0.39590 | 0.33367 | 0.37707 | 0.35578 |

Absolute difference from official (SDs):

| Parameter | C++ 1001 | C++ 42 | C++ 2024 | C++ 7 |
|-----------|----------|--------|----------|-------|
| $\mu_0$ | 0.01039 | 0.01183 | 0.02843 | 0.00166 |
| $\mu_{X_1}$ | 0.00651 | 0.00481 | 0.00969 | 0.00483 |
| $\mu_{X_2}$ | 0.00327 | 0.00016 | 0.00054 | 0.00150 |
| $\mu_{X_1 X_2}$ | 0.00181 | 0.00582 | 0.00219 | 0.00991 |
| $\sigma^2$ | 0.00035 | 0.00036 | 0.00001 | 0.00033 |
| $\Sigma_{11}$ | 0.15446 | 0.06731 | 0.00253 | 0.08303 |
| $\Sigma_{22}$ | 0.01817 | 0.04515 | 0.02342 | 0.01899 |
| $\Sigma_{33}$ | 0.01206 | 0.00694 | 0.00726 | 0.01146 |
| $\Sigma_{44}$ | 0.01114 | 0.07336 | 0.02996 | 0.05125 |

$\sigma^2$ is stable across seeds. $\Sigma_{11}$ is not: the `1001` chain has a heavier right tail (mean 1.53, SD 0.62); the other three seeds sit at means 1.29–1.33, in line with official 1.26 and the R engine 1.29. The official posterior SD of $\Sigma_{11}$ is already 0.47, so a 0.27 shift in the mean at a single seed is compatible with MCMC noise on this functional. A single C++ run at `seed = 1001` should not be read as a systematic engine bias.

---

## 5. Summary

1. **Same data.** The official `rep_1` posterior and `fit_ipd_ad_lm` both use replicate 1 of `SimulationData_1.RData`.
2. **Slope $\mu$ and $\sigma^2$.** Both engines match official posterior means and SDs closely.
3. **Intercept $\mu_0$.** Official mean is near 0; package means are about 0.05–0.08. The gap is smaller than the posterior SD, and the official vs package MCMC seeds are not the same protocol.
4. **$\mathrm{diag}(\Sigma)$.** At `seed = 1001`, C++ $\Sigma_{11}$ is high in both mean and SD. Other C++ seeds (`42`, `2024`, `7`) match official as closely as the R engine; the `1001` gap is seed-sensitive MCMC variation, not a consistent C++ bias (Section 4).
5. **Between engines.** $\sigma^2$ is nearly identical; $\mu$ gaps are at Monte Carlo / implementation scale; the largest split at `seed = 1001` is $\Sigma_{11}$.
6. **Speed.** For the same 20,000 iterations, C++ is about 28 times faster than R.

---

## 6. Usage

```r
library(bayesmetaipd)
d <- sim1_as_formula_data()
common <- list(
  formula = d$formula, ipd = d$ipd, study = d$study,
  nested_formula = d$nested_formula, ad_nested = d$ad_nested,
  nested_reported = d$nested_reported,
  subgroup = d$subgroup, ad_subgroup = d$ad_subgroup,
  partial_terms = d$partial_terms, ad_partial = d$ad_partial,
  drm_formula = d$drm_formula,
  theta_init_ipd = d$theta_init_ipd,
  theta_init_nested = d$theta_init_nested,
  theta_init_subgroup = d$theta_init_subgroup,
  theta_init_partial = d$theta_init_partial,
  mu_init = d$mu_init, sig2_init = d$sig2_init,
  burnin = 10000L, mainrun = 10000L, seed = 1001L, verbose = FALSE
)
fit_r   <- do.call(fit_ipd_ad_lm, c(common, list(engine = "r")))
fit_cpp <- do.call(fit_ipd_ad_lm, c(common, list(engine = "cpp")))
```
