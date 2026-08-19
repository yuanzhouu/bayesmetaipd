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
  miss <- setdiff(reported, names(ad))
  if (length(miss)) {
    stop(label, ": missing coefficient columns ", paste(miss, collapse = ", "), call. = FALSE)
  }
  beta <- as.matrix(ad[, reported, drop = FALSE])
  storage.mode(beta) <- "double"
  K <- nrow(beta)
  se_cols <- paste0("se_", reported)
  if (is.list(ad$V)) {
    if (length(ad$V) != K) stop(label, ": `V` must have one matrix per row.", call. = FALSE)
    V <- ad$V
  } else if (all(se_cols %in% names(ad))) {
    se <- as.matrix(ad[, se_cols, drop = FALSE])
    V <- lapply(seq_len(K), function(k) diag(as.numeric(se[k, ])^2, length(reported)))
  } else {
    stop(
      label, ": provide `se_` columns for each reported coefficient, or a list-column `V`.",
      call. = FALSE
    )
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


#' Convert Simulation 1 replicate 1 into formula-style inputs
#'
#' Builds an IPD data frame plus Type 1/2/3 AD tables from [sim1_ipdad_rep1],
#' using the official full model `Y ~ X1 * X2`, nested model `~ X1 + X2`,
#' four `(X1>0) x X2` subgroups, and partial terms `X2` and `X1:X2`.
#'
#' @return A list with `ipd`, `study`, `formula`, `nested_formula`, `ad_nested`,
#'   `subgroup`, `ad_subgroup`, `partial_terms`, `ad_partial`, `drm_formula`,
#'   and MCMC starting values.
#' @export
#' @examples
#' d <- sim1_as_formula_data()
#' names(d$ad_nested)
sim1_as_formula_data <- function() {
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
  ad_nested <- data.frame(
    study = t1,
    X1 = dat$beta_mat[t1, 2],
    X2 = dat$beta_mat[t1, 3],
    drm_mean = vapply(t1, function(s) drm_stats(s)[["mean"]], 1),
    drm_var = vapply(t1, function(s) drm_stats(s)[["var"]], 1)
  )
  ad_nested$V <- lapply(t1, function(s) {
    as.matrix(dat$V_beta_cube[s, 2:3, 2:3])
  })

  t2 <- which(type_vec == 2L)
  ad_subgroup <- data.frame(
    study = t2,
    ind.1 = dat$beta_mat[t2, 1],
    ind.2 = dat$beta_mat[t2, 2],
    ind.3 = dat$beta_mat[t2, 3],
    ind.4 = dat$beta_mat[t2, 4],
    drm_mean = vapply(t2, function(s) drm_stats(s)[["mean"]], 1),
    drm_var = vapply(t2, function(s) drm_stats(s)[["var"]], 1)
  )
  ad_subgroup$V <- lapply(t2, function(s) {
    as.matrix(dat$V_beta_cube[s, , ])
  })

  t3 <- which(type_vec == 3L)
  ad_partial <- data.frame(
    study = t3,
    X2 = dat$beta_mat[t3, 3],
    `X1:X2` = dat$beta_mat[t3, 4],
    drm_mean = vapply(t3, function(s) drm_stats(s)[["mean"]], 1),
    drm_var = vapply(t3, function(s) drm_stats(s)[["var"]], 1),
    check.names = FALSE
  )
  ad_partial$V <- lapply(t3, function(s) {
    as.matrix(dat$V_beta_cube[s, 3:4, 3:4])
  })

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


#' Fit linear IPD + Type 1/2/3 AD from formulas
#'
#' General continuous-outcome version of Simulation Study 1 IPD+AD. The user
#' supplies:
#' * a **full IPD formula**
#' * a **nested working formula** (Type 1) and its reported coefficients
#' * **subgroup definitions** (Type 2) and subgroup means
#' * **partial terms** from the full model (Type 3) and those coefficients
#'
#' Each AD table also needs density-ratio summaries `drm_mean` and `drm_var`
#' for the covariate in `drm_formula` (paper: mean and variance-of-the-mean of
#' the baseline covariate).
#'
#' @param formula Full model for IPD, e.g. `y ~ x1 * x2`.
#' @param ipd IPD data frame.
#' @param study Study id column in `ipd`.
#' @param nested_formula Type 1 working model, e.g. `~ x1 + x2`. Ignored if
#'   `ad_nested` is `NULL`.
#' @param ad_nested Type 1 AD: data frame with reported nested coefficients,
#'   either `se_<coef>` columns or a list-column `V`, plus `drm_mean`,
#'   `drm_var`. Or a list with `beta`, `V`/`se`, `drm_mean`, `drm_var`.
#' @param nested_reported Names (or indices) of nested coefficients that were
#'   published. Default: all except the intercept.
#' @param subgroup Named list of Type 2 subgroup formulas, e.g.
#'   `list(g1 = ~ x1 > 0 & x2 == 0, ...)`.
#' @param ad_subgroup Type 2 AD table; coefficient columns must match
#'   `names(subgroup)`.
#' @param partial_terms Character vector (or indices) of full-model terms
#'   reported by Type 3 studies, e.g. `c("x2", "x1:x2")`.
#' @param ad_partial Type 3 AD table for those terms.
#' @param drm_formula One-sided formula for the density-ratio covariate,
#'   default `~` the first non-intercept full-model term. Paper uses `~ x1`.
#' @param use_drm If `FALSE`, fix the density ratio at 1 (no `alpha`/`tau` MH).
#' @param diagonal_V If `TRUE` (default, official), replace each AD covariance
#'   with its diagonal.
#' @param burnin,mainrun MCMC lengths.
#' @param step_theta,step_alpha,step_tau MH scales (official Sim1: 0.2, 0.01, 0.02).
#' @param lambda,nu0,phi0 Prior hyperparameters.
#' @param theta_init_ipd,theta_init_nested,theta_init_subgroup,theta_init_partial
#'   Optional starting `theta` matrices.
#' @param mu_init,Sigma_init,sig2_init Starting values.
#' @param seed RNG seed. `NULL` leaves RNG unchanged.
#' @param verbose Progress printing.
#'
#' @return A `bayesmetaipd_fit` with `posterior_mu` (columns = full-model terms),
#'   `posterior_Sigma_diag`, `posterior_sig2`, `call`, `settings`.
#'
#' @examples
#' d <- sim1_as_formula_data()
#' fit <- fit_ipd_ad_lm(
#'   formula = d$formula,
#'   ipd = d$ipd,
#'   study = d$study,
#'   nested_formula = d$nested_formula,
#'   ad_nested = d$ad_nested,
#'   nested_reported = d$nested_reported,
#'   subgroup = d$subgroup,
#'   ad_subgroup = d$ad_subgroup,
#'   partial_terms = d$partial_terms,
#'   ad_partial = d$ad_partial,
#'   drm_formula = d$drm_formula,
#'   burnin = 1, mainrun = 2, verbose = FALSE, seed = 1
#' )
#' colMeans(fit$posterior_mu)
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
                          verbose = TRUE) {
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
  if (K < 1L) stop("Provide at least one of `ad_nested`, `ad_subgroup`, `ad_partial`.", call. = FALSE)
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
      used_default_data = FALSE
    )
  )
  class(out) <- c("bayesmetaipd_fit", "list")
  out
}
