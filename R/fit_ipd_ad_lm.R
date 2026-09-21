# General linear IPD + 3-type AD (formula interface)

rhs_formula <- function(formula) {
  formula <- stats::as.formula(formula)
  if (length(formula) == 3L) {
    stats::reformulate(deparse(formula[[3]]), intercept = TRUE)
  } else {
    formula
  }
}

eval_subgroup <- function(expr, data) {
  if (inherits(expr, "formula")) {
    val <- eval(expr[[length(expr)]], data, parent.frame())
  } else if (is.character(expr) && length(expr) == 1L) {
    val <- eval(parse(text = expr), data, parent.frame())
  } else if (is.function(expr)) {
    val <- expr(data)
  } else {
    stop("Each subgroup must be a formula, character expression, or function.", call. = FALSE)
  }
  as.numeric(val)
}

MB.est.linear <- function(theta, X, Z, w, drm_psi, alpha) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  theta <- as.vector(theta)
  alpha <- as.vector(alpha)
  w <- as.numeric(w)
  drm <- as.vector(exp(as.matrix(drm_psi) %*% alpha))
  ww <- w * drm
  Wz <- Z * ww
  Mn1 <- t(Z) %*% Wz
  Mn2 <- t(Z) %*% (X * ww) %*% theta
  solve(Mn1, Mn2)
}

parse_ad_estimates <- function(ad, reported, label) {
  if (is.null(ad)) return(NULL)
  if (is.list(ad) && !is.data.frame(ad)) {
    if (is.null(ad$beta)) {
      stop(label, ": list input must contain `beta`.", call. = FALSE)
    }
    beta <- as.matrix(ad$beta)
    if (!is.null(colnames(beta))) {
      clean_names <- gsub("^coef_|^mean_", "", colnames(beta))
      if (all(reported %in% clean_names)) {
        colnames(beta) <- clean_names
      }
    }
    if (is.null(colnames(beta))) {
      if (ncol(beta) != length(reported)) {
        stop(label, ": `beta` needs colnames or ncol = length(reported).", call. = FALSE)
      }
      colnames(beta) <- reported
    }
    miss <- setdiff(reported, colnames(beta))
    if (length(miss)) {
      stop(label, ": `beta` missing columns ", paste(miss, collapse = ", "), call. = FALSE)
    }
    beta <- beta[, reported, drop = FALSE]
    K <- nrow(beta)
    if (!is.null(ad$V)) {
      if (length(ad$V) != K) stop(label, ": `V` must have one matrix per study.", call. = FALSE)
      V <- ad$V
    } else if (!is.null(ad$se)) {
      se <- as.matrix(ad$se)
      if (!is.null(colnames(se))) {
        clean_se_names <- gsub("^se_|^coef_|^mean_", "", colnames(se))
        if (all(reported %in% clean_se_names)) {
          colnames(se) <- clean_se_names
        }
      }
      if (is.null(colnames(se))) colnames(se) <- reported
      se <- se[, reported, drop = FALSE]
      V <- lapply(seq_len(K), function(k) diag(se[k, ]^2, length(reported)))
    } else {
      stop(label, ": provide `V` (list of matrices) or `se`.", call. = FALSE)
    }
    drm_mean <- if (is.null(ad$drm_mean)) rep(NA_real_, K) else as.numeric(ad$drm_mean)
    drm_var <- if (is.null(ad$drm_var)) rep(NA_real_, K) else as.numeric(ad$drm_var)
    if (length(drm_mean) != K || length(drm_var) != K) {
      stop(label, ": `drm_mean` / `drm_var` length must equal number of studies.", call. = FALSE)
    }
    return(list(beta = beta, V = V, drm_mean = drm_mean, drm_var = drm_var, K = K))
  }

  ad <- as.data.frame(ad)
  nms <- names(ad)

  # Match coefficient / subgroup mean columns:
  # Accepts: "coef_<term>", "mean_<term>", or "<term>"
  est_cols <- character(length(reported))
  for (j in seq_along(reported)) {
    r <- reported[j]
    cands <- c(paste0("coef_", r), paste0("mean_", r), r)
    hit <- cands[cands %in% nms]
    if (length(hit) > 0L) {
      est_cols[j] <- hit[1L]
    } else {
      est_cols[j] <- NA_character_
    }
  }

  if (anyNA(est_cols)) {
    miss <- reported[is.na(est_cols)]
    stop(label, ": missing coefficient/mean columns for: ", paste(miss, collapse = ", "),
         " (expected e.g. `coef_", miss[1L], "`, `mean_", miss[1L], "`, or `", miss[1L], "`).",
         call. = FALSE)
  }

  beta <- as.matrix(ad[, est_cols, drop = FALSE])
  colnames(beta) <- reported
  storage.mode(beta) <- "double"
  K <- nrow(beta)

  if (is.list(ad$V)) {
    if (length(ad$V) != K) stop(label, ": `V` must have one matrix per row.", call. = FALSE)
    V <- ad$V
  } else {
    # Match standard error columns:
    # Accepts: "se_<est_col>", "se_<term>", "se_coef_<term>", "se_mean_<term>"
    se_cols <- character(length(reported))
    for (j in seq_along(reported)) {
      r <- reported[j]
      matched_est <- est_cols[j]
      se_cands <- unique(c(
        paste0("se_", matched_est),
        paste0("se_", r),
        paste0("se_coef_", r),
        paste0("se_mean_", r)
      ))
      hit_se <- se_cands[se_cands %in% nms]
      if (length(hit_se) > 0L) {
        se_cols[j] <- hit_se[1L]
      } else {
        se_cols[j] <- NA_character_
      }
    }
    if (anyNA(se_cols)) {
      miss_se <- reported[is.na(se_cols)]
      stop(label, ": missing standard error columns for: ", paste(miss_se, collapse = ", "),
           " (expected e.g. `se_", miss_se[1L], "` or `se_", est_cols[which(reported == miss_se[1L])], "`).",
           call. = FALSE)
    }
    se <- as.matrix(ad[, se_cols, drop = FALSE])
    V <- lapply(seq_len(K), function(k) diag(as.numeric(se[k, ])^2, length(reported)))
  }

  if (!("drm_mean" %in% names(ad)) || !("drm_var" %in% names(ad))) {
    stop(label, ": need columns `drm_mean` and `drm_var` for the density-ratio covariate.", call. = FALSE)
  }
  list(
    beta = beta,
    V = V,
    drm_mean = as.numeric(ad$drm_mean),
    drm_var = as.numeric(ad$drm_var),
    K = K
  )
}

prepare_ipd_lm <- function(formula, ipd, study) {
  if (!is.data.frame(ipd)) stop("`ipd` must be a data.frame.", call. = FALSE)
  if (!study %in% names(ipd)) stop("`study` column not found in `ipd`.", call. = FALSE)
  formula <- stats::as.formula(formula)
  if (length(formula) != 3L) stop("`formula` must be like y ~ x1 * x2.", call. = FALSE)
  ids <- unique(ipd[[study]])
  out <- vector("list", length(ids))
  p <- NULL
  xnames <- NULL
  for (i in seq_along(ids)) {
    rows <- ipd[ipd[[study]] == ids[[i]], , drop = FALSE]
    mf <- stats::model.frame(formula, data = rows, na.action = stats::na.pass)
    if (any(is.na(mf))) {
      stop("IPD study ", ids[[i]], " has missing values in the model frame.", call. = FALSE)
    }
    y <- as.numeric(stats::model.response(mf))
    X <- stats::model.matrix(formula, data = mf)
    storage.mode(X) <- "double"
    if (is.null(p)) {
      p <- ncol(X)
      xnames <- colnames(X)
    } else if (!identical(colnames(X), xnames)) {
      stop("All IPD studies must share the same design matrix columns.", call. = FALSE)
    }
    out[[i]] <- list(
      y = y,
      X = X,
      data = rows,
      n = length(y),
      study = as.character(ids[[i]])
    )
  }
  list(studies = out, p = p, xnames = xnames, formula = formula)
}

match_reported <- function(requested, available, label) {
  if (is.numeric(requested)) {
    idx <- as.integer(requested)
    if (any(idx < 1L | idx > length(available))) {
      stop(label, ": coefficient index out of range.", call. = FALSE)
    }
    return(list(names = available[idx], index = idx))
  }
  requested <- as.character(requested)
  idx <- match(requested, available)
  if (anyNA(idx)) {
    alt_avail <- gsub(":", "", available, fixed = TRUE)
    idx2 <- match(gsub(":", "", requested, fixed = TRUE), alt_avail)
    idx[is.na(idx)] <- idx2[is.na(idx)]
  }
  if (anyNA(idx)) {
    stop(
      label, ": cannot match ", paste(requested[is.na(idx)], collapse = ", "),
      ". Available: ", paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  list(names = available[idx], index = as.integer(idx))
}


#' Load example formula-style dataset from Simulation Study 1
#'
#' Builds an IPD data frame plus Type 1/2/3 AD tables with explicit coefficients
#' (`coef_*` / `mean_*`), standard errors (`se_coef_*` / `se_mean_*`), and density-ratio
#' moments (`drm_mean`, `drm_var`) from [sim1_ipdad_rep1].
#' Uses the official full model `Y ~ X1 * X2`, nested model `~ X1 + X2`,
#' four `(X1>0) x X2` subgroups, and partial terms `X2` and `X1:X2`.
#'
#' @return A list with `ipd`, `study`, `formula`, `nested_formula`, `ad_nested`,
#'   `nested_reported`, `subgroup`, `ad_subgroup`, `partial_terms`, `ad_partial`, `drm_formula`.
#' @export
#' @examples
#' d <- load_example()
#' head(d$ipd)
#' head(d$ad_nested)
load_example <- function() {
  dat <- load_sim1_ipdad_rep1()
  L <- dim(dat$X_cube)[1]
  n <- dim(dat$X_cube)[2]
  type_vec <- as.integer(dat$type_vec)
  ipd_ids <- which(type_vec == 4L)
  rows <- vector("list", length(ipd_ids))
  for (t in seq_along(ipd_ids)) {
    s <- ipd_ids[t]
    rows[[t]] <- data.frame(
      study = s,
      Y = dat$Y_mat[s, ],
      X1 = dat$X_cube[s, , "X1"],
      X2 = dat$X_cube[s, , "X2"],
      stringsAsFactors = FALSE
    )
  }
  ipd <- do.call(rbind, rows)

  drm_stats <- function(idx) {
    x1 <- dat$X_cube[idx, , "X1"]
    c(mean = mean(x1), var = stats::var(x1) / length(x1))
  }

  t1 <- which(type_vec == 1L)
  V_nested <- lapply(t1, function(s) {
    as.matrix(dat$V_beta_cube[s, 2:3, 2:3])
  })
  ad_nested <- data.frame(
    study      = t1,
    coef_X1    = dat$beta_mat[t1, 2],
    coef_X2    = dat$beta_mat[t1, 3],
    se_coef_X1 = vapply(V_nested, function(m) sqrt(m[1, 1]), numeric(1)),
    se_coef_X2 = vapply(V_nested, function(m) sqrt(m[2, 2]), numeric(1)),
    drm_mean   = vapply(t1, function(s) drm_stats(s)[["mean"]], 1),
    drm_var    = vapply(t1, function(s) drm_stats(s)[["var"]], 1)
  )

  t2 <- which(type_vec == 2L)
  V_subgroup <- lapply(t2, function(s) {
    as.matrix(dat$V_beta_cube[s, , ])
  })
  ad_subgroup <- data.frame(
    study         = t2,
    mean_ind.1    = dat$beta_mat[t2, 1],
    mean_ind.2    = dat$beta_mat[t2, 2],
    mean_ind.3    = dat$beta_mat[t2, 3],
    mean_ind.4    = dat$beta_mat[t2, 4],
    se_mean_ind.1 = vapply(V_subgroup, function(m) sqrt(m[1, 1]), numeric(1)),
    se_mean_ind.2 = vapply(V_subgroup, function(m) sqrt(m[2, 2]), numeric(1)),
    se_mean_ind.3 = vapply(V_subgroup, function(m) sqrt(m[3, 3]), numeric(1)),
    se_mean_ind.4 = vapply(V_subgroup, function(m) sqrt(m[4, 4]), numeric(1)),
    drm_mean      = vapply(t2, function(s) drm_stats(s)[["mean"]], 1),
    drm_var       = vapply(t2, function(s) drm_stats(s)[["var"]], 1)
  )

  t3 <- which(type_vec == 3L)
  V_partial <- lapply(t3, function(s) {
    as.matrix(dat$V_beta_cube[s, 3:4, 3:4])
  })
  ad_partial <- data.frame(
    study           = t3,
    coef_X2         = dat$beta_mat[t3, 3],
    `coef_X1:X2`    = dat$beta_mat[t3, 4],
    se_coef_X2      = vapply(V_partial, function(m) sqrt(m[1, 1]), numeric(1)),
    `se_coef_X1:X2` = vapply(V_partial, function(m) sqrt(m[2, 2]), numeric(1)),
    drm_mean        = vapply(t3, function(s) drm_stats(s)[["mean"]], 1),
    drm_var         = vapply(t3, function(s) drm_stats(s)[["var"]], 1),
    check.names     = FALSE
  )

  list(
    formula = Y ~ X1 * X2,
    ipd = ipd,
    study = "study",
    nested_formula = ~ X1 + X2,
    ad_nested = ad_nested,
    nested_reported = c("X1", "X2"),
    subgroup = list(
      ind.1 = ~ X1 > 0 & X2 == 0,
      ind.2 = ~ X1 > 0 & X2 == 1,
      ind.3 = ~ X1 <= 0 & X2 == 0,
      ind.4 = ~ X1 <= 0 & X2 == 1
    ),
    ad_subgroup = ad_subgroup,
    partial_terms = c("X2", "X1:X2"),
    ad_partial = ad_partial,
    drm_formula = ~ X1,
    theta_init_ipd = dat$theta_l_mat[ipd_ids, , drop = FALSE],
    theta_init_nested = dat$theta_l_mat[t1, , drop = FALSE],
    theta_init_subgroup = dat$theta_l_mat[t2, , drop = FALSE],
    theta_init_partial = dat$theta_l_mat[t3, , drop = FALSE],
    mu_init = dat$true_mu,
    sig2_init = dat$true_kappa,
    random_seed = dat$random_seed
  )
}

#' @export
sim1_as_formula_data <- load_example


#' Fit gaussian linear model to IPD + Type 1/2/3 AD data
#'
#' Fit gaussian linear model to IPD + Type 1/2/3 AD data within a Bayesian
#' random-effects framework via estimating equations, multiplier bootstrap,
#' and density ratio models.
#'
#' The user supplies:
#' * a **full IPD formula**
#' * a **nested working formula** (Type 1) and its reported coefficients
#' * **subgroup definitions** (Type 2) and subgroup means
#' * **partial terms** from the full model (Type 3) and those coefficients
#'
#' Each AD table also needs density-ratio summaries `drm_mean` and `drm_var`
#' for the covariate in `drm_formula` (paper: mean and variance-of-the-mean of
#' the baseline covariate).
#'
#' @param formula Two-sided formula specifying the full target model for
#'   Individual Participant Data (IPD), e.g., `Y ~ X1 + X2 + X1:X2` (or
#'   `Y ~ X1 * X2`). This model represents the primary inferential target of
#'   interest across all studies: for each study \eqn{i}, the individual-level
#'   continuous outcome is modeled as \eqn{Y_{ij} = X_{ij}^\top \beta_i + \epsilon_{ij}}
#'   with a study-specific coefficient vector \eqn{\beta_i} (including a random
#'   intercept by default) and residual error
#'   \eqn{\epsilon_{ij} \sim \mathcal{N}(0, \sigma^2)}, where \eqn{\beta_i} follows
#'   a hierarchical random-effects distribution.
#' @param ipd Data frame of Individual Patient Data (Individual Participant Data),
#'   containing individual-level records from all available IPD studies. It must
#'   contain the continuous response variable (`Y`), all predictors/covariates
#'   specified in `formula` (e.g., `X1`, `X2`), and a study identifier column
#'   identifying which study each participant belongs to.
#' @param study Character string specifying the name of the study identifier
#'   column in the `ipd` data frame (e.g., `"study"`).
#' @param nested_formula One-sided formula defining the Type 1 working model
#'   (nested/reduced model), e.g., `~ X1 + X2`. Type 1 Aggregate Data (AD)
#'   studies fit a reduced working model that omits certain terms from the full
#'   target model (e.g., omitting the interaction term `X1:X2`). For simplicity,
#'   the response variable `Y` on the left-hand side of the tilde (`~`) is
#'   omitted. Ignored if `ad_nested` is `NULL`.
#' @param ad_nested Data frame of summary statistics from Type 1 AD studies (one
#'   row per study). Columns must include: (1) an optional study identifier
#'   column (e.g., `"study"`); (2) columns for reported point estimates of the
#'   nested model coefficients matching `nested_reported` or terms in
#'   `nested_formula` (e.g., columns named `"coef_X1"`, `"coef_X2"`); (3)
#'   standard errors in columns named `se_<term>` (e.g., `"se_X1"`, `"se_X2"`
#'   or `"se_coef_X1"`, `"se_coef_X2"`); and (4) density-ratio summary columns
#'   `drm_mean` (sample mean \eqn{\bar{x}}) and `drm_var` (variance of the
#'   sample mean \eqn{\widehat{\mathrm{Var}}(\bar{x}) = s^2 / n = \mathrm{SE}^2},
#'   where \eqn{s} is the sample standard deviation and \eqn{n} is the study
#'   sample size; note that it is *not* the sample variance \eqn{s^2}) of the
#'   baseline covariate specified in `drm_formula`.
#' @param nested_reported Character vector (or integer indices) indicating which
#'   coefficients from the Type 1 nested working model were published by the AD
#'   studies. By default, all coefficients in `nested_formula` except the
#'   intercept are assumed to be reported.
#' @param subgroup Named list of one-sided formulas defining the subgroup
#'   indicators for Type 2 AD studies. Type 2 AD studies report sample subgroup
#'   outcome means and standard errors across partitions of the covariate space
#'   rather than regression coefficients. For example:
#'   `list(g1 = ~ X1 > 0 & X2 == 0, g2 = ~ X1 <= 0 & X2 == 0, g3 = ~ X1 > 0 & X2 == 1, g4 = ~ X1 <= 0 & X2 == 1)`.
#'   Each formula evaluates to a binary indicator on the participant-level
#'   covariate space.
#' @param ad_subgroup Data frame of summary statistics from Type 2 AD studies
#'   (one row per study). Columns must include: (1) an optional study identifier
#'   column (e.g., `"study"`); (2) sample subgroup mean estimates, where column
#'   names match `names(subgroup)` (e.g., columns named `"mean_g1"`,
#'   `"mean_g2"`); (3) corresponding standard errors in columns named
#'   `se_<subgroup>` (e.g., `"se_g1"`, `"se_g2"` or `"se_mean_g1"`,
#'   `"se_mean_g2"`); and (4) density-ratio summary columns `drm_mean` and
#'   `drm_var` for baseline covariate shift adjustment.
#' @param partial_terms Character vector (or integer indices) specifying which
#'   full-model terms were published by Type 3 AD studies, e.g.,
#'   `c("X2", "X1:X2")`. Type 3 AD studies fit the complete target model but
#'   report only a subset of the estimated coefficients.
#' @param ad_partial Data frame of summary statistics from Type 3 AD studies
#'   (one row per study). Columns must include: (1) an optional study identifier
#'   column (e.g., `"study"`); (2) reported coefficient estimates for the terms
#'   in `partial_terms` (e.g., columns named `"coef_X2"`, `"coef_X1:X2"`); (3)
#'   corresponding standard errors in columns named `se_<term>` (e.g.,
#'   `"se_X2"`, `"se_X1:X2"` or `"se_coef_X2"`, `"se_coef_X1:X2"`); and (4)
#'   density-ratio summary columns `drm_mean` and `drm_var` for baseline
#'   covariate shift adjustment.
#' @param drm_formula One-sided formula specifying the baseline covariate used
#'   in the semi-parametric Density-Ratio Model (DRM), e.g., `~ X1`. The DRM uses
#'   exponential tilting to account for covariate shift (population heterogeneity)
#'   between AD and IPD study populations based on published aggregate moments
#'   `drm_mean` and `drm_var`. Here, `drm_mean` is the reported sample mean of
#'   the covariate (\eqn{\bar{x}}), and `drm_var` is the variance of the sample
#'   mean (\eqn{\widehat{\mathrm{Var}}(\bar{x}) = s^2 / n = \mathrm{SE}^2}, computed as
#'   sample variance divided by study sample size, or squared standard error;
#'   it should *not* be the individual-level sample variance \eqn{s^2}). If `NULL`,
#'   defaults to the first non-intercept covariate in the full target model.
#' @param use_drm Logical; if `TRUE` (default), applies semi-parametric
#'   density-ratio modeling to adjust for covariate shift between IPD and AD
#'   populations. If `FALSE`, assumes homogeneous covariate distributions
#'   across studies, fixing density-ratio weights to 1 and skipping
#'   Metropolis-Hastings updates for DRM tilt parameters.
#' @param diagonal_V Logical; if `TRUE` (default), keeps only the diagonal
#'   variances of each AD covariance matrix (assuming zero off-diagonal
#'   sampling covariances).
#' @param burnin,mainrun Positive integers specifying MCMC sampling lengths:
#'   `burnin` is the number of initial warm-up/burn-in iterations to discard
#'   (default: `10000`), and `mainrun` is the number of post-burn-in iterations
#'   retained for posterior inference (default: `10000`).
#' @param step_theta,step_alpha,step_tau Numeric proposal standard deviations
#'   for random-walk Metropolis-Hastings steps: `step_theta` for study-specific
#'   effect vectors (\eqn{\theta_l}, default: `0.2`); `step_alpha` for DRM tilt
#'   parameters (\eqn{\alpha_l}, default: `0.01`); and `step_tau` for baseline
#'   covariate moments (\eqn{\tau_l}, default: `0.02`).
#' @param lambda,nu0,phi0 Prior hyperparameters for the Bayesian hierarchical
#'   model: `lambda` is the variance multiplier for the population mean prior
#'   \eqn{\mu \sim \mathcal{N}(0, \lambda I)} (default: `1e4`); `nu0` and
#'   `phi0` are the degrees of freedom and scale multiplier for the
#'   Inverse-Wishart hyperprior on between-study covariance
#'   \eqn{\Sigma \sim \text{Inv-Wishart}(\nu_0, \phi_0 I)} (defaults:
#'   `nu0 = 0.1`, `phi0 = 0.1`).
#' @param theta_init_ipd,theta_init_nested,theta_init_subgroup,theta_init_partial
#'   Optional numeric matrices of starting values for study-specific parameters
#'   (\eqn{\theta_l}) across IPD studies, Type 1 nested AD studies, Type 2
#'   subgroup AD studies, and Type 3 partial AD studies, respectively. Each
#'   matrix must have rows equal to the number of studies in that category and
#'   columns equal to the number of full-model terms. Defaults to `NULL` for
#'   each, which initializes all starting values to zeros.
#' @param mu_init,Sigma_init,sig2_init Optional starting values for the Markov
#'   chain: `mu_init` is a numeric vector for the population mean effect
#'   \eqn{\mu} (default: `NULL`, initialized to a zero vector); `Sigma_init` is
#'   a positive-definite matrix for between-study covariance \eqn{\Sigma}
#'   (default: `NULL`, initialized to the identity matrix); and `sig2_init` is
#'   a positive scalar for residual error variance \eqn{\sigma^2} (default:
#'   `1`).
#' @param seed An integer random seed passed to `set.seed()` for exact
#'   reproducibility of MCMC chains. If `NULL`, the current RNG state is left
#'   unchanged.
#' @param verbose Logical; if `TRUE`, displays MCMC iteration progress and
#'   diagnostic messages to the console.
#' @param engine Character string specifying the computation engine: `"r"`
#'   (default) uses the pure R MCMC sampler; `"cpp"` uses an optimized C++
#'   (Rcpp/RcppArmadillo) inner loop providing significant speedup (same
#'   mathematical model, but not bit-identical due to floating-point
#'   differences).
#' @return An S3 object of class `bayesmetaipd_fit` (inheriting from `list`)
#'   containing post-burn-in MCMC samples and model metadata:
#'   \describe{
#'     \item{`posterior_mu`}{A numeric matrix of dimension `mainrun x p` containing
#'       posterior draws of the population-level mean coefficients (\eqn{\mu}).
#'       Column names match the full-model terms.}
#'     \item{`posterior_Sigma_diag`}{A numeric matrix of dimension `mainrun x p`
#'       containing posterior draws of the between-study variance components
#'       (\eqn{\text{diag}(\Sigma)}).}
#'     \item{`posterior_sig2`}{A numeric vector of length `mainrun` containing
#'       posterior draws of the residual error variance (\eqn{\sigma^2}).}
#'     \item{`call`}{The matched function call.}
#'     \item{`settings`}{A list recording execution configurations and summary
#'       metadata, including MCMC parameters (`burnin`, `mainrun`, MH step sizes),
#'       prior hyperparameters (`lambda`, `nu0`, `phi0`), random `seed`, study
#'       counts (`L`, `J`, `K`, `J_type1`, `J_type2`, `J_type3`), coefficient names,
#'       and model formulas.}
#'   }
#'   Objects of class `bayesmetaipd_fit` have a dedicated `print` method.
#'
#' @examples
#' \donttest{
#' d <- load_example()
#' fit <- fit_ipd_ad_lm(
#'   formula         = d$formula,              # Full model: Y ~ X1 * X2
#'   ipd             = d$ipd,                  # Individual participant dataset
#'   study           = d$study,                # Study identifier column
#'   nested_formula  = d$nested_formula,       # Type 1 AD: nested working formula (~ X1 + X2)
#'   ad_nested       = d$ad_nested,            # Type 1 AD table (uses default non-intercept terms)
#'   subgroup        = d$subgroup,             # Type 2 AD: 4 subgroup partition formulas
#'   ad_subgroup     = d$ad_subgroup,          # Type 2 AD table
#'   partial_terms   = d$partial_terms,        # Type 3 AD: reported subset c("X2", "X1:X2")
#'   ad_partial      = d$ad_partial,           # Type 3 AD table
#'   drm_formula     = d$drm_formula,          # Density-ratio covariate (~ X1)
#'   burnin          = 1000, 
#'   mainrun         = 2000, 
#'   engine          = "cpp"                   # Accelerated C++ MCMC sampler
#' )
#' print(fit)
#' colMeans(fit$posterior_mu)
#' }
#'
#' @export
fit_ipd_ad_lm <- function(formula,
                          ipd,
                          study,
                          nested_formula = NULL,
                          ad_nested = NULL,
                          nested_reported = NULL,
                          subgroup = NULL,
                          ad_subgroup = NULL,
                          partial_terms = NULL,
                          ad_partial = NULL,
                          drm_formula = NULL,
                          use_drm = TRUE,
                          diagonal_V = TRUE,
                          burnin = 10000L,
                          mainrun = 10000L,
                          step_theta = 0.2,
                          step_alpha = 0.01,
                          step_tau = 0.02,
                          lambda = 1e4,
                          nu0 = 0.1,
                          phi0 = 0.1,
                          theta_init_ipd = NULL,
                          theta_init_nested = NULL,
                          theta_init_subgroup = NULL,
                          theta_init_partial = NULL,
                          mu_init = NULL,
                          Sigma_init = NULL,
                          sig2_init = 1,
                          seed = 1001L,
                          verbose = TRUE,
                          engine = c("r", "cpp")) {
  engine <- match.arg(engine)
  ipd_prep <- prepare_ipd_lm(formula, ipd, study)
  ipd_st <- ipd_prep$studies
  J <- length(ipd_st)
  p_theta <- ipd_prep$p
  xnames <- ipd_prep$xnames
  if (J < 1L) stop("Need at least one IPD study.", call. = FALSE)

  if (isTRUE(use_drm)) {
    if (is.null(drm_formula)) {
      non_int <- xnames[xnames != "(Intercept)"]
      if (!length(non_int)) {
        stop("Cannot infer `drm_formula`; provide it explicitly.", call. = FALSE)
      }
      drm_formula <- stats::reformulate(non_int[[1]])
    }
    drm_formula <- rhs_formula(drm_formula)
    drm_probe <- stats::model.matrix(drm_formula, data = ipd_st[[1]]$data)
    p_alpha <- ncol(drm_probe)
    if (p_alpha < 2L) {
      stop("`drm_formula` should include intercept + at least one covariate (paper: ~ x1).", call. = FALSE)
    }
    drm_cov_col <- 2L
  } else {
    drm_formula <- ~1
    p_alpha <- 2L
    drm_cov_col <- 2L
  }

  build_drm <- function(data) {
    if (!isTRUE(use_drm)) {
      return(cbind(1, rep(0, nrow(data))))
    }
    mm <- stats::model.matrix(drm_formula, data = data)
    if (ncol(mm) != p_alpha) {
      stop("DRM design rank changed across studies.", call. = FALSE)
    }
    mm
  }
  drm_mean_ipd <- vapply(ipd_st, function(st) {
    mean(build_drm(st$data)[, drm_cov_col])
  }, numeric(1))

  pick_reference <- function(hat_tau) {
    ipd_st[[which.min(abs(hat_tau - drm_mean_ipd))]]
  }

  # ---- Type 1 nested ----
  t1 <- list()
  if (!is.null(ad_nested)) {
    if (is.null(nested_formula)) {
      stop("`nested_formula` is required when `ad_nested` is supplied.", call. = FALSE)
    }
    nested_formula <- rhs_formula(nested_formula)
    Z_probe <- stats::model.matrix(nested_formula, data = ipd_st[[1]]$data)
    nested_names <- colnames(Z_probe)
    if (is.null(nested_reported)) {
      nested_reported <- setdiff(nested_names, "(Intercept)")
    }
    hit <- match_reported(nested_reported, nested_names, "Type 1")
    parsed <- parse_ad_estimates(ad_nested, hit$names, "ad_nested")
    for (k in seq_len(parsed$K)) {
      Vk <- as.matrix(parsed$V[[k]])
      if (isTRUE(diagonal_V)) Vk <- diag(diag(Vk), nrow(Vk))
      t1[[k]] <- list(
        beta_tilde = as.numeric(parsed$beta[k, ]),
        V = Vk,
        reported = hit$index,
        hat_tau = parsed$drm_mean[k],
        hat_gamma = parsed$drm_var[k],
        ref = pick_reference(parsed$drm_mean[k])
      )
      t1[[k]]$X <- t1[[k]]$ref$X
      t1[[k]]$Z <- stats::model.matrix(nested_formula, data = t1[[k]]$ref$data)
      t1[[k]]$drm_psi <- build_drm(t1[[k]]$ref$data)
    }
  }

  # ---- Type 2 subgroup ----
  t2 <- list()
  if (!is.null(ad_subgroup)) {
    if (is.null(subgroup) || !length(subgroup)) {
      stop("`subgroup` is required when `ad_subgroup` is supplied.", call. = FALSE)
    }
    if (is.null(names(subgroup)) || any(!nzchar(names(subgroup)))) {
      names(subgroup) <- paste0("g", seq_along(subgroup))
    }
    gnames <- names(subgroup)
    parsed <- parse_ad_estimates(ad_subgroup, gnames, "ad_subgroup")
    for (k in seq_len(parsed$K)) {
      Vk <- as.matrix(parsed$V[[k]])
      if (isTRUE(diagonal_V)) Vk <- diag(diag(Vk), nrow(Vk))
      ref <- pick_reference(parsed$drm_mean[k])
      Z <- sapply(subgroup, eval_subgroup, data = ref$data)
      if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
      colnames(Z) <- gnames
      t2[[k]] <- list(
        beta_tilde = as.numeric(parsed$beta[k, ]),
        V = Vk,
        reported = seq_along(gnames),
        hat_tau = parsed$drm_mean[k],
        hat_gamma = parsed$drm_var[k],
        ref = ref,
        X = ref$X,
        Z = Z,
        drm_psi = build_drm(ref$data)
      )
    }
  }

  # ---- Type 3 partial full-model coefficients ----
  t3 <- list()
  if (!is.null(ad_partial)) {
    if (is.null(partial_terms)) {
      stop("`partial_terms` is required when `ad_partial` is supplied.", call. = FALSE)
    }
    hit <- match_reported(partial_terms, xnames, "Type 3")
    parsed <- parse_ad_estimates(ad_partial, hit$names, "ad_partial")
    for (k in seq_len(parsed$K)) {
      Vk <- as.matrix(parsed$V[[k]])
      if (isTRUE(diagonal_V)) Vk <- diag(diag(Vk), nrow(Vk))
      ref <- pick_reference(parsed$drm_mean[k])
      t3[[k]] <- list(
        beta_tilde = as.numeric(parsed$beta[k, ]),
        V = Vk,
        reported = hit$index,
        hat_tau = parsed$drm_mean[k],
        hat_gamma = parsed$drm_var[k],
        ref = ref,
        X = ref$X,
        Z = ref$X,
        drm_psi = build_drm(ref$data)
      )
    }
  }

  K1 <- length(t1)
  K2 <- length(t2)
  K3 <- length(t3)
  K <- K1 + K2 + K3
  if (K < 1L) {
    warning("No AD data provided (`ad_nested`, `ad_subgroup`, and `ad_partial` are all NULL). Running in IPD-only mode.", call. = FALSE)
  }
  L <- J + K
  if ((nu0 + L) < p_theta) {
    stop("Need nu0 + L >= p for Inverse-Wishart. Got nu0+L=", nu0 + L, ", p=", p_theta, ".", call. = FALSE)
  }

  init_theta <- function(mat, n_st, label) {
    if (is.null(mat)) return(matrix(0, n_st, p_theta))
    mat <- as.matrix(mat)
    if (!all(dim(mat) == c(n_st, p_theta))) {
      stop(label, " must be ", n_st, " x ", p_theta, ".", call. = FALSE)
    }
    mat
  }
  theta_mat_IPD <- init_theta(theta_init_ipd, J, "`theta_init_ipd`")
  theta_t1 <- init_theta(theta_init_nested, K1, "`theta_init_nested`")
  theta_t2 <- init_theta(theta_init_subgroup, K2, "`theta_init_subgroup`")
  theta_t3 <- init_theta(theta_init_partial, K3, "`theta_init_partial`")
  if (is.null(mu_init)) mu_init <- rep(0, p_theta)
  if (is.null(Sigma_init)) Sigma_init <- diag(1, p_theta)
  mu_vec <- as.numeric(mu_init)
  Sigma_theta_mat <- as.matrix(Sigma_init)
  sig2 <- as.numeric(sig2_init)[1]

  alpha0 <- if (isTRUE(use_drm)) 0.1 else 0
  start_beta <- function(studies, theta_mat) {
    beta <- vector("list", length(studies))
    for (k in seq_along(studies)) {
      st <- studies[[k]]
      beta[[k]] <- tryCatch(
        as.numeric(MB.est.linear(
          theta_mat[k, ], st$X, st$Z, rep(1, nrow(st$X)), st$drm_psi, rep(alpha0, p_alpha)
        )),
        error = function(e) {
          if (ncol(st$Z) == p_theta) as.numeric(theta_mat[k, ]) else rep(0, ncol(st$Z))
        }
      )
    }
    beta
  }
  beta_t1 <- start_beta(t1, theta_t1)
  beta_t2 <- start_beta(t2, theta_t2)
  beta_t3 <- start_beta(t3, theta_t3)

  alpha_t1 <- array(alpha0, c(max(K1, 1L), p_alpha))
  alpha_t2 <- array(alpha0, c(max(K2, 1L), p_alpha))
  alpha_t3 <- array(alpha0, c(max(K3, 1L), p_alpha))
  tau_t1 <- if (K1) vapply(t1, `[[`, numeric(1), "hat_gamma") else numeric(0)
  tau_t2 <- if (K2) vapply(t2, `[[`, numeric(1), "hat_gamma") else numeric(0)
  tau_t3 <- if (K3) vapply(t3, `[[`, numeric(1), "hat_gamma") else numeric(0)

  burnin <- as.integer(burnin)
  mainrun <- as.integer(mainrun)
  if (burnin < 0L || mainrun < 1L) stop("`burnin` >= 0 and `mainrun` >= 1 required.", call. = FALSE)
  if (!is.null(seed)) set.seed(as.integer(seed))

  invLambda_theta <- diag(1 / lambda, p_theta)
  Phi0 <- diag(phi0, p_theta)
  n_iter <- burnin + mainrun
  draw_mu <- array(0, c(n_iter, p_theta))
  colnames(draw_mu) <- xnames
  draw_Sigma <- array(0, c(n_iter, p_theta))
  colnames(draw_Sigma) <- xnames
  draw_sig2 <- rep(0, n_iter)

  drm_moments <- function(drm_psi, alpha, tau) {
    exp_a <- exp(as.matrix(drm_psi) %*% as.numeric(alpha))
    q <- cbind(exp_a - 1, drm_psi[, drm_cov_col] * exp_a - tau)
    list(bar = apply(q, 2, mean), S = stats::var(q) / nrow(q))
  }

  if (identical(engine, "cpp")) {
    all_st <- c(t1, t2, t3)
    theta_ad0 <- rbind(theta_t1, theta_t2, theta_t3)
    alpha_rows <- function(a, n) {
      if (n == 0L) {
        return(matrix(numeric(), nrow = 0L, ncol = p_alpha))
      }
      a[seq_len(n), , drop = FALSE]
    }
    alpha_ad0 <- rbind(
      alpha_rows(alpha_t1, K1),
      alpha_rows(alpha_t2, K2),
      alpha_rows(alpha_t3, K3)
    )
    if (verbose) {
      message(sprintf(
        "Starting formula IPD+AD MCMC (C++): %d iters; J=%d IPD, nested=%d, subgroup=%d, partial=%d, p=%d",
        n_iter, J, K1, K2, K3, p_theta
      ))
    }
    cpp_out <- lm_mcmc_cpp(
      X_ipd = lapply(ipd_st, function(st) unname(as.matrix(st$X))),
      y_ipd = lapply(ipd_st, function(st) as.numeric(st$y)),
      X_ad = lapply(all_st, function(st) unname(as.matrix(st$X))),
      Z_ad = lapply(all_st, function(st) unname(as.matrix(st$Z))),
      drm_ad = lapply(all_st, function(st) unname(as.matrix(st$drm_psi))),
      V_ad = lapply(all_st, function(st) unname(as.matrix(st$V))),
      beta_tilde = lapply(all_st, function(st) as.numeric(st$beta_tilde)),
      reported = lapply(all_st, function(st) as.integer(st$reported - 1L)),
      hat_tau = vapply(all_st, `[[`, numeric(1), "hat_tau"),
      hat_gamma = vapply(all_st, `[[`, numeric(1), "hat_gamma"),
      theta_ipd0 = unname(theta_mat_IPD),
      theta_ad0 = unname(theta_ad0),
      beta_ad0 = c(beta_t1, beta_t2, beta_t3),
      alpha0 = unname(alpha_ad0),
      tau0 = as.numeric(c(tau_t1, tau_t2, tau_t3)),
      extra_theta = as.logical(c(rep(FALSE, K1), rep(FALSE, K2), rep(TRUE, K3))),
      mu0 = as.numeric(mu_vec),
      Sigma0 = unname(Sigma_theta_mat),
      sig2 = sig2,
      invLambda = unname(invLambda_theta),
      Phi0 = unname(Phi0),
      n_iter = as.integer(n_iter),
      burnin = as.integer(burnin),
      step_theta = step_theta,
      step_alpha = step_alpha,
      step_tau = step_tau,
      nu0 = nu0,
      use_drm = isTRUE(use_drm),
      drm_cov_col = as.integer(drm_cov_col - 1L),
      verbose = verbose
    )
    mu_keep <- cpp_out$posterior_mu
    colnames(mu_keep) <- xnames
    sig_keep <- cpp_out$posterior_Sigma_diag
    colnames(sig_keep) <- xnames
    out <- list(
      posterior_mu = mu_keep,
      posterior_Sigma_diag = sig_keep,
      posterior_sig2 = as.numeric(cpp_out$posterior_sig2),
      call = match.call(),
      settings = list(
        burnin = burnin,
        mainrun = mainrun,
        step_theta = step_theta,
        step_alpha = step_alpha,
        step_tau = step_tau,
        lambda = lambda,
        nu0 = nu0,
        phi0 = phi0,
        seed = seed,
        L = L,
        J = J,
        K = K,
        J_type1 = K1,
        J_type2 = K2,
        J_type3 = K3,
        n = mean(vapply(ipd_st, `[[`, numeric(1), "n")),
        p = p_theta,
        coef_names = xnames,
        formula = formula,
        nested_formula = nested_formula,
        drm_formula = drm_formula,
        use_drm = isTRUE(use_drm),
        outcome = "gaussian_lm",
        engine = "cpp",
        used_default_data = FALSE
      )
    )
    class(out) <- c("bayesmetaipd_fit", "list")
    return(out)
  }

  mh_ad <- function(studies, theta_mat, beta_list, alpha_mat, tau_vec, extra_theta_in_alpha) {
    for (k in seq_along(studies)) {
      st <- studies[[k]]
      n_ref <- nrow(st$X)
      q_z <- ncol(st$Z)
      rep_idx <- st$reported
      theta_q <- stats::rnorm(p_theta, mean = theta_mat[k, ], sd = step_theta)
      ww <- stats::rnorm(n_ref, 1, 1)
      beta_q <- tryCatch(
        as.numeric(MB.est.linear(theta_q, st$X, st$Z, ww, st$drm_psi, alpha_mat[k, ])),
        error = function(e) NULL
      )
      if (!is.null(beta_q) && length(beta_q) == q_z) {
        logNum <- mvtnorm::dmvnorm(st$beta_tilde, mean = beta_q[rep_idx], sigma = st$V, log = TRUE)
        logNum <- logNum + mvtnorm::dmvnorm(theta_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
        logDen <- mvtnorm::dmvnorm(st$beta_tilde, mean = beta_list[[k]][rep_idx], sigma = st$V, log = TRUE)
        logDen <- logDen + mvtnorm::dmvnorm(theta_mat[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
        if (stats::runif(1) < exp(logNum - logDen)) {
          theta_mat[k, ] <- theta_q
          beta_list[[k]] <- beta_q
        }
      }

      if (!isTRUE(use_drm)) next
      alpha_q <- stats::rnorm(p_alpha, mean = alpha_mat[k, ], sd = step_alpha)
      ww <- stats::rnorm(n_ref, 1, 1)
      beta_q <- tryCatch(
        as.numeric(MB.est.linear(theta_mat[k, ], st$X, st$Z, ww, st$drm_psi, alpha_q)),
        error = function(e) NULL
      )
      if (!is.null(beta_q) && length(beta_q) == q_z) {
        mq <- drm_moments(st$drm_psi, alpha_q, tau_vec[k])
        mc <- drm_moments(st$drm_psi, alpha_mat[k, ], tau_vec[k])
        logNum <- mvtnorm::dmvnorm(st$beta_tilde, mean = beta_q[rep_idx], sigma = st$V, log = TRUE)
        logNum <- logNum + mvtnorm::dmvnorm(mq$bar, mean = rep(0, 2), sigma = mq$S, log = TRUE)
        logDen <- mvtnorm::dmvnorm(st$beta_tilde, mean = beta_list[[k]][rep_idx], sigma = st$V, log = TRUE)
        logDen <- logDen + mvtnorm::dmvnorm(mc$bar, mean = rep(0, 2), sigma = mc$S, log = TRUE)
        if (isTRUE(extra_theta_in_alpha)) {
          logNum <- logNum + mvtnorm::dmvnorm(theta_mat[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
          logDen <- logDen + mvtnorm::dmvnorm(theta_mat[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
        }
        if (stats::runif(1) < exp(logNum - logDen)) {
          alpha_mat[k, ] <- alpha_q
          beta_list[[k]] <- beta_q
        }
      }
    }
    list(theta = theta_mat, beta = beta_list, alpha = alpha_mat)
  }

  if (verbose) {
    message(sprintf(
      "Starting formula IPD+AD MCMC: %d iters; J=%d IPD, nested=%d, subgroup=%d, partial=%d, p=%d",
      n_iter, J, K1, K2, K3, p_theta
    ))
  }
  prev_time <- proc.time()[[3]]

  for (i_iter in seq_len(n_iter)) {
    if (K1) {
      u <- mh_ad(t1, theta_t1, beta_t1, alpha_t1, tau_t1, extra_theta_in_alpha = FALSE)
      theta_t1 <- u$theta
      beta_t1 <- u$beta
      alpha_t1 <- u$alpha
    }
    if (K2) {
      u <- mh_ad(t2, theta_t2, beta_t2, alpha_t2, tau_t2, extra_theta_in_alpha = FALSE)
      theta_t2 <- u$theta
      beta_t2 <- u$beta
      alpha_t2 <- u$alpha
    }
    if (K3) {
      u <- mh_ad(t3, theta_t3, beta_t3, alpha_t3, tau_t3, extra_theta_in_alpha = TRUE)
      theta_t3 <- u$theta
      beta_t3 <- u$beta
      alpha_t3 <- u$alpha
    }

    for (j in seq_len(J)) {
      theta_q <- stats::rnorm(p_theta, mean = theta_mat_IPD[j, ], sd = step_theta)
      x.the_q <- ipd_st[[j]]$X %*% theta_q
      x.the <- ipd_st[[j]]$X %*% theta_mat_IPD[j, ]
      logAcc <- sum(
        stats::dnorm(ipd_st[[j]]$y, x.the_q, sqrt(sig2), log = TRUE) -
          stats::dnorm(ipd_st[[j]]$y, x.the, sqrt(sig2), log = TRUE)
      )
      logAcc <- logAcc +
        mvtnorm::dmvnorm(theta_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE) -
        mvtnorm::dmvnorm(theta_mat_IPD[j, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
      if (stats::runif(1) < exp(logAcc)) theta_mat_IPD[j, ] <- theta_q
    }

    theta_merge <- rbind(theta_mat_IPD, theta_t1, theta_t2, theta_t3)
    inv_Sigma <- solve(Sigma_theta_mat)
    inv_Var <- invLambda_theta + L * inv_Sigma
    Var <- solve(inv_Var)
    Mean <- Var %*% (inv_Sigma %*% apply(theta_merge, 2, sum))
    mu_vec <- as.numeric(mvtnorm::rmvnorm(1, mean = Mean, sigma = Var))

    SS <- array(0, c(p_theta, p_theta))
    for (l in seq_len(L)) {
      dth <- theta_merge[l, ] - mu_vec
      SS <- SS + tcrossprod(dth)
    }
    Sigma_theta_mat <- MCMCpack::riwish(nu0 + L, Phi0 + SS)

    sse <- 0
    ntot <- 0
    for (j in seq_len(J)) {
      r <- ipd_st[[j]]$y - as.numeric(ipd_st[[j]]$X %*% theta_mat_IPD[j, ])
      sse <- sse + sum(r^2)
      ntot <- ntot + ipd_st[[j]]$n
    }
    sig2 <- 1 / stats::rgamma(1, shape = 1 + ntot / 2, rate = 1 + sse / 2)

    update_tau <- function(studies, alpha_mat, tau_vec) {
      for (k in seq_along(studies)) {
        st <- studies[[k]]
        tau_q <- stats::rnorm(1, tau_vec[k], step_tau)
        mq <- drm_moments(st$drm_psi, alpha_mat[k, ], tau_q)
        mc <- drm_moments(st$drm_psi, alpha_mat[k, ], tau_vec[k])
        logNum <- stats::dnorm(st$hat_tau, tau_q, sqrt(st$hat_gamma), log = TRUE) +
          mvtnorm::dmvnorm(mq$bar, mean = rep(0, 2), sigma = mq$S, log = TRUE)
        logDen <- stats::dnorm(st$hat_tau, tau_vec[k], sqrt(st$hat_gamma), log = TRUE) +
          mvtnorm::dmvnorm(mc$bar, mean = rep(0, 2), sigma = mc$S, log = TRUE)
        if (stats::runif(1) < exp(logNum - logDen)) tau_vec[k] <- tau_q
      }
      tau_vec
    }
    if (isTRUE(use_drm)) {
      if (K1) tau_t1 <- update_tau(t1, alpha_t1, tau_t1)
      if (K2) tau_t2 <- update_tau(t2, alpha_t2, tau_t2)
      if (K3) tau_t3 <- update_tau(t3, alpha_t3, tau_t3)
    }

    draw_mu[i_iter, ] <- mu_vec
    draw_Sigma[i_iter, ] <- diag(Sigma_theta_mat)
    draw_sig2[i_iter] <- sig2

    if (verbose && (i_iter %% 1000L == 0L)) {
      cur <- proc.time()[[3]]
      last_batch <- cur - prev_time
      eta <- (n_iter - i_iter) * (last_batch / 1000)
      prev_time <- cur
      message(sprintf(
        "iter %d / %d | last 1000 = %.1f min | ETA = %.1f min",
        i_iter, n_iter, last_batch / 60, eta / 60
      ))
    }
  }

  seq_keep <- if (burnin == 0L) seq_len(mainrun) else (burnin + 1L):(burnin + mainrun)
  out <- list(
    posterior_mu = draw_mu[seq_keep, , drop = FALSE],
    posterior_Sigma_diag = draw_Sigma[seq_keep, , drop = FALSE],
    posterior_sig2 = draw_sig2[seq_keep],
    call = match.call(),
    settings = list(
      burnin = burnin,
      mainrun = mainrun,
      step_theta = step_theta,
      step_alpha = step_alpha,
      step_tau = step_tau,
      lambda = lambda,
      nu0 = nu0,
      phi0 = phi0,
      seed = seed,
      L = L,
      J = J,
      K = K,
      J_type1 = K1,
      J_type2 = K2,
      J_type3 = K3,
      n = mean(vapply(ipd_st, `[[`, numeric(1), "n")),
      p = p_theta,
      coef_names = xnames,
      formula = formula,
      nested_formula = nested_formula,
      drm_formula = drm_formula,
      use_drm = isTRUE(use_drm),
      outcome = "gaussian_lm",
      engine = "r",
      used_default_data = FALSE
    )
  )
  class(out) <- c("bayesmetaipd_fit", "list")
  out
}
