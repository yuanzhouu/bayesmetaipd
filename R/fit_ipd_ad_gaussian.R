# Application Type-3 moment objective (Gaussian threshold probabilities).
# Official script names this MB.logit.est; kept separate from Sim2 logistic helper.
MB.app.type3.est <- function(param,
                             theta_vec,
                             Design_mat,
                             bootstrap_w,
                             nu_x_mat,
                             threshold_type3,
                             sig2_IPD_type3,
                             UseDRM,
                             alpha,
                             REF_X_type3_input) {
  Design_mat <- as.matrix(Design_mat)
  theta_vec <- as.vector(theta_vec)
  param <- as.vector(param)
  nu_x_mat <- as.matrix(nu_x_mat)
  std_dev <- (threshold_type3 - Design_mat %*% theta_vec) / sqrt(sig2_IPD_type3)
  exp_y_above_tau <- 1 - stats::pnorm(std_dev, mean = 0, sd = 1)
  if (isTRUE(UseDRM)) {
    REF_X_type3_mat <- REF_X_type3_input[, c("int", "DRM_X_type3"), drop = FALSE]
    w.DRM <- exp(as.matrix(REF_X_type3_mat) %*% alpha)
    bootstrap_w <- as.vector(w.DRM) * bootstrap_w
  }
  exp_minus_fitted <- bootstrap_w * as.vector(exp_y_above_tau - nu_x_mat %*% param)
  g.i <- cbind(
    nu_x_mat[, 1] * exp_minus_fitted,
    nu_x_mat[, 2] * exp_minus_fitted,
    nu_x_mat[, 3] * exp_minus_fitted,
    nu_x_mat[, 4] * exp_minus_fitted,
    nu_x_mat[, 5] * exp_minus_fitted,
    nu_x_mat[, 6] * exp_minus_fitted
  )
  g.sum <- apply(g.i, 2, sum)
  sum(g.sum^2)
}


validate_ad_type3 <- function(ad_type3, p_theta) {
  req <- c(
    "betahat", "se", "use_drm", "hat_tau", "hat_gamma_tau",
    "Design_X", "nu_X", "threshold", "DRM_X"
  )
  miss <- setdiff(req, names(ad_type3))
  if (length(miss)) {
    stop("`ad_type3` missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  J3 <- nrow(as.matrix(ad_type3$betahat))
  if (J3 < 1L) stop("`ad_type3` must have at least one study.", call. = FALSE)
  if (!all(lengths(list(
    ad_type3$Design_X, ad_type3$nu_X, ad_type3$threshold, ad_type3$DRM_X
  )) == J3)) {
    stop("Type3 list components must all have length J_type3.", call. = FALSE)
  }
  if (ncol(as.matrix(ad_type3$betahat)) != 6L || ncol(as.matrix(ad_type3$se)) != 6L) {
    stop("`ad_type3$betahat` and `$se` must be J x 6.", call. = FALSE)
  }
  if (length(ad_type3$use_drm) != J3) stop("`use_drm` length must be J_type3.", call. = FALSE)
  if (length(ad_type3$hat_tau) != J3 || length(ad_type3$hat_gamma_tau) != J3) {
    stop("`hat_tau` / `hat_gamma_tau` length must be J_type3.", call. = FALSE)
  }
  if (p_theta != 7L) {
    stop("Application Type1–3 currently require p = 7 covariates.", call. = FALSE)
  }
  invisible(NULL)
}


type2_beta_from_theta <- function(theta, p_vec) {
  beta <- (p_vec[["p_11"]] - p_vec[["p_01"]]) * theta[1] + p_vec[["p_11"]] * theta[5]
  beta <- beta + (p_vec[["p_12"]] - p_vec[["p_02"]]) * theta[2] + p_vec[["p_12"]] * theta[6]
  beta <- beta + (p_vec[["p_13"]] - p_vec[["p_03"]]) * theta[3] + p_vec[["p_13"]] * theta[7]
  as.numeric(beta)
}


#' Fit Application-style continuous IPD + AD meta-analysis
#'
#' Implements `Code/Application_CodeOnly.R`: Gaussian IPD with study-specific
#' residual variances, plus Aggregate Data of three report types
#' (treatment-effect Type1, BMI-mixture Type2, threshold-probability Type3 with
#' optional density-ratio adjustment).
#'
#' Real I-WIP Application CSVs are not bundled. Provide prepared objects, or use
#' [example_application_data()] for a small synthetic example.
#'
#' @param ipd IPD via [as_application_ipd()].
#' @param ad_type1 Optional data frame with columns `beta_hat`, `sqrt_hat_V`.
#'   Type1 bridges to coefficient 7 (`bmi_cat3.trt` in Application coding).
#' @param ad_type2 Optional data frame with columns
#'   `p_01,p_02,p_03,p_11,p_12,p_13,beta_hat,sqrt_hat_V`.
#' @param ad_type3 Optional list with Type3 components:
#'   `betahat`, `se` (`J x 6`), `use_drm` (logical), `hat_tau`, `hat_gamma_tau`,
#'   `Design_X`, `nu_X`, `threshold`, `DRM_X` (each length-`J` lists / vectors).
#'   Each `DRM_X[[k]]` needs columns `int`, `DRM_X_type3`, `a_n`.
#' @param burnin,mainrun MCMC lengths. Defaults `10000`.
#' @param mu_sig2,tau_sig2 Truncated-normal prior for residual variances.
#' @param mu0,Sigma0 Prior for `mu`.
#' @param nu0,psi0 Inverse-Wishart hyperparameters.
#' @param step_type1,step_type2,step_type3,step_sig2,step_alpha,step_tau MH scales.
#' @param tau_update_after Start tau MH only after this iteration (official: 1000).
#' @param type1_coef_index Coefficient index bridged by Type1. Default `7`.
#' @param theta_init,mu_init,Sigma_init,sig2_init Starting values.
#' @param seed Default `117`. `NULL` leaves RNG unchanged.
#' @param verbose,progress_every Progress printing controls.
#'
#' @return A `bayesmetaipd_fit` with `posterior_mu`, `posterior_Sigma_diag`,
#'   `posterior_sig2` (matrix over IPD + Type3 studies), `call`, `settings`.
#'
#' @examples
#' toy <- example_application_data()
#' fit <- fit_ipd_ad_gaussian(
#'   ipd = toy$ipd,
#'   ad_type1 = toy$ad_type1,
#'   ad_type2 = toy$ad_type2,
#'   ad_type3 = toy$ad_type3,
#'   burnin = 2, mainrun = 3, verbose = FALSE,
#'   tau_update_after = 0L
#' )
#' colMeans(fit$posterior_mu)
#'
#' @keywords internal
#' @noRd
fit_ipd_ad_gaussian <- function(ipd,
                                ad_type1 = NULL,
                                ad_type2 = NULL,
                                ad_type3 = NULL,
                                burnin = 10000L,
                                mainrun = 10000L,
                                mu_sig2 = 20,
                                tau_sig2 = 20,
                                mu0 = NULL,
                                Sigma0 = NULL,
                                nu0 = 0.1,
                                psi0 = 0.1,
                                step_type1 = 0.2,
                                step_type2 = 0.2,
                                step_type3 = 0.1,
                                step_sig2 = 2,
                                step_alpha = 0.001,
                                step_tau = 0.2,
                                tau_update_after = 1000L,
                                type1_coef_index = 7L,
                                theta_init = NULL,
                                mu_init = NULL,
                                Sigma_init = NULL,
                                sig2_init = NULL,
                                seed = 117L,
                                verbose = TRUE,
                                progress_every = 100L) {
  IPD_listobj <- as_application_ipd(ipd)
  J_IPD <- length(IPD_listobj)
  p_theta <- ncol(IPD_listobj[[1]]$X)
  for (l in seq_len(J_IPD)) {
    if (ncol(IPD_listobj[[l]]$X) != p_theta) {
      stop("All IPD studies must share the same number of covariates.", call. = FALSE)
    }
  }

  if (is.null(ad_type1)) {
    AD_type1 <- NULL
    J_type1 <- 0L
  } else {
    AD_type1 <- as.data.frame(ad_type1)
    need1 <- c("beta_hat", "sqrt_hat_V")
    if (!all(need1 %in% names(AD_type1))) {
      stop("`ad_type1` needs columns beta_hat, sqrt_hat_V.", call. = FALSE)
    }
    J_type1 <- nrow(AD_type1)
  }

  if (is.null(ad_type2)) {
    AD_type2 <- NULL
    J_type2 <- 0L
  } else {
    AD_type2 <- as.data.frame(ad_type2)
    need2 <- c("p_01", "p_02", "p_03", "p_11", "p_12", "p_13", "beta_hat", "sqrt_hat_V")
    if (!all(need2 %in% names(AD_type2))) {
      stop("`ad_type2` needs columns ", paste(need2, collapse = ", "), call. = FALSE)
    }
    J_type2 <- nrow(AD_type2)
  }

  if (is.null(ad_type3)) {
    J_type3 <- 0L
  } else {
    validate_ad_type3(ad_type3, p_theta)
    J_type3 <- nrow(as.matrix(ad_type3$betahat))
  }

  if ((J_type1 + J_type2) > 0L && p_theta < max(7L, type1_coef_index)) {
    stop("Type1/Type2 require p >= 7 (Application covariate coding).", call. = FALSE)
  }

  L <- J_IPD + J_type1 + J_type2 + J_type3
  if ((nu0 + L) < p_theta) {
    stop(
      "Need nu0 + L >= p for Inverse-Wishart (`riwish`). ",
      "Got nu0 + L = ", nu0 + L, ", p = ", p_theta, ".",
      call. = FALSE
    )
  }
  p_alpha <- 2L
  p_beta_type3 <- 6L

  if (is.null(mu0)) mu0 <- rep(0, p_theta)
  if (is.null(Sigma0)) Sigma0 <- diag(100, p_theta)
  mu0 <- as.numeric(mu0)
  Sigma0 <- as.matrix(Sigma0)
  inv_Sigma_0_mat <- solve(Sigma0)
  Psi_mat <- diag(psi0, p_theta)

  burnin <- as.integer(burnin)
  mainrun <- as.integer(mainrun)
  if (burnin < 0L || mainrun < 1L) {
    stop("`burnin` >= 0 and `mainrun` >= 1 required.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  if (is.null(theta_init)) {
    theta_l_mat <- matrix(0.1, nrow = L, ncol = p_theta)
  } else {
    theta_l_mat <- as.matrix(theta_init)
    if (!all(dim(theta_l_mat) == c(L, p_theta))) {
      stop("`theta_init` must be L x p with L = J_IPD + J_type1 + J_type2 + J_type3.", call. = FALSE)
    }
  }

  n_sig2 <- J_IPD + J_type3
  if (is.null(sig2_init)) {
    sig2_IPD_type3 <- rep(10, n_sig2)
  } else {
    sig2_IPD_type3 <- as.numeric(sig2_init)
    if (length(sig2_IPD_type3) == 1L) sig2_IPD_type3 <- rep(sig2_IPD_type3, n_sig2)
    if (length(sig2_IPD_type3) != n_sig2) {
      stop("`sig2_init` must have length J_IPD + J_type3.", call. = FALSE)
    }
  }

  if (is.null(mu_init)) {
    mu_vec <- matrix(0, nrow = p_theta, ncol = 1)
  } else {
    mu_vec <- matrix(as.numeric(mu_init), nrow = p_theta, ncol = 1)
  }
  if (is.null(Sigma_init)) {
    Sigma_mat <- diag(1, p_theta)
  } else {
    Sigma_mat <- as.matrix(Sigma_init)
  }
  inv_Sigma_mat <- solve(Sigma_mat)

  if (J_type3 > 0L) {
    AD_type3_betahat <- as.matrix(ad_type3$betahat)
    AD_type3_SE_betahat <- as.matrix(ad_type3$se)
    use_drm <- as.logical(ad_type3$use_drm)
    hat_tau_type3_vec <- as.numeric(ad_type3$hat_tau)
    hat_Gamma_tau_type3_vec <- as.numeric(ad_type3$hat_gamma_tau)
    Design_X_type3 <- ad_type3$Design_X
    nu_X_type3 <- ad_type3$nu_X
    threshold_type3 <- ad_type3$threshold
    DRM_X_type3 <- ad_type3$DRM_X
    tau_vec <- hat_tau_type3_vec
    alpha_mat_type3 <- array(0.01, c(J_type3, p_alpha))
    beta_l_type3_mat <- array(0.1, c(J_type3, p_beta_type3))

    for (k in seq_len(J_type3)) {
      l <- J_IPD + J_type1 + J_type2 + k
      n_l <- nrow(as.matrix(Design_X_type3[[k]]))
      ww <- stats::rnorm(n_l, mean = 1, sd = 1)
      temp.optim <- stats::optim(
        beta_l_type3_mat[k, ],
        fn = MB.app.type3.est,
        theta_vec = theta_l_mat[l, ],
        Design_mat = Design_X_type3[[k]],
        bootstrap_w = ww,
        nu_x_mat = nu_X_type3[[k]],
        threshold_type3 = threshold_type3[[k]],
        sig2_IPD_type3 = sig2_IPD_type3[J_IPD + k],
        UseDRM = use_drm[k],
        alpha = alpha_mat_type3[k, ],
        REF_X_type3_input = DRM_X_type3[[k]],
        control = list(pgtol = 1e-7, maxit = 1000),
        lower = 0,
        upper = 1,
        method = "L-BFGS-B"
      )
      beta_l_type3_mat[k, ] <- temp.optim$par
    }
  }

  n_iter <- burnin + mainrun
  Draw_mu_vec <- array(0, c(n_iter, p_theta))
  Draw_Sigma_diag <- array(0, c(n_iter, p_theta))
  Draw_sig2 <- array(0, c(n_iter, n_sig2))

  if (verbose) {
    message(sprintf(
      "Starting Application MCMC: %d iters (%d burn-in + %d main); J_IPD=%d, type1=%d, type2=%d, type3=%d",
      n_iter, burnin, mainrun, J_IPD, J_type1, J_type2, J_type3
    ))
  }
  prev_time <- proc.time()[[3]]

  for (i_iter in seq_len(n_iter)) {
    # ---- IPD theta (Gibbs) ----
    for (l in seq_len(J_IPD)) {
      y <- IPD_listobj[[l]]$y
      X <- IPD_listobj[[l]]$X
      inv_Var <- t(X) %*% X / sig2_IPD_type3[l] + inv_Sigma_mat
      Var <- solve(inv_Var)
      Mean_temp_part <- t(X) %*% y / sig2_IPD_type3[l] + inv_Sigma_mat %*% mu_vec
      Mean <- Var %*% Mean_temp_part
      theta_l_mat[l, ] <- mvtnorm::rmvnorm(n = 1, mean = Mean, sigma = Var)
    }

    # ---- Type 1 AD ----
    for (k in seq_len(J_type1)) {
      l <- J_IPD + k
      theta_cur <- theta_l_mat[l, ]
      beta_cur <- theta_cur[type1_coef_index]
      theta_q <- stats::rnorm(n = p_theta, mean = theta_cur, sd = step_type1)
      beta_q <- theta_q[type1_coef_index]
      logNUM <- mvtnorm::dmvnorm(theta_q, mean = as.numeric(mu_vec), sigma = Sigma_mat, log = TRUE) +
        stats::dnorm(AD_type1[k, "beta_hat"], beta_q, AD_type1[k, "sqrt_hat_V"], log = TRUE)
      logDEN <- mvtnorm::dmvnorm(theta_cur, mean = as.numeric(mu_vec), sigma = Sigma_mat, log = TRUE) +
        stats::dnorm(AD_type1[k, "beta_hat"], beta_cur, AD_type1[k, "sqrt_hat_V"], log = TRUE)
      if (stats::runif(1) <= exp(logNUM - logDEN)) {
        theta_l_mat[l, ] <- theta_q
      }
    }

    # ---- Type 2 AD ----
    for (k in seq_len(J_type2)) {
      l <- J_IPD + J_type1 + k
      theta_cur <- theta_l_mat[l, ]
      p_vec <- AD_type2[k, c("p_01", "p_02", "p_03", "p_11", "p_12", "p_13")]
      beta_cur <- type2_beta_from_theta(theta_cur, p_vec)
      theta_q <- stats::rnorm(n = p_theta, mean = theta_cur, sd = step_type2)
      beta_q <- type2_beta_from_theta(theta_q, p_vec)
      logNUM <- mvtnorm::dmvnorm(theta_q, mean = as.numeric(mu_vec), sigma = Sigma_mat, log = TRUE) +
        stats::dnorm(AD_type2[k, "beta_hat"], beta_q, AD_type2[k, "sqrt_hat_V"], log = TRUE)
      logDEN <- mvtnorm::dmvnorm(theta_cur, mean = as.numeric(mu_vec), sigma = Sigma_mat, log = TRUE) +
        stats::dnorm(AD_type2[k, "beta_hat"], beta_cur, AD_type2[k, "sqrt_hat_V"], log = TRUE)
      if (stats::runif(1) <= exp(logNUM - logDEN)) {
        theta_l_mat[l, ] <- theta_q
      }
    }

    # ---- Type 3 AD ----
    for (k in seq_len(J_type3)) {
      n_l <- nrow(as.matrix(Design_X_type3[[k]]))
      l <- J_IPD + J_type1 + J_type2 + k

      theta_type3_vec_q <- stats::rnorm(n = p_theta, mean = theta_l_mat[l, ], sd = step_type3)
      logNUMDEN1 <- mvtnorm::dmvnorm(theta_type3_vec_q, mean = as.numeric(mu_vec), sigma = Sigma_mat, log = TRUE) -
        mvtnorm::dmvnorm(theta_l_mat[l, ], mean = as.numeric(mu_vec), sigma = Sigma_mat, log = TRUE)

      ww <- stats::rnorm(n_l, mean = 1, sd = 1)
      temp.optim <- stats::optim(
        beta_l_type3_mat[k, ],
        fn = MB.app.type3.est,
        theta_vec = theta_type3_vec_q,
        Design_mat = Design_X_type3[[k]],
        bootstrap_w = ww,
        nu_x_mat = nu_X_type3[[k]],
        threshold_type3 = threshold_type3[[k]],
        sig2_IPD_type3 = sig2_IPD_type3[J_IPD + k],
        UseDRM = use_drm[k],
        alpha = alpha_mat_type3[k, ],
        REF_X_type3_input = DRM_X_type3[[k]],
        control = list(pgtol = 1e-7, maxit = 1000),
        lower = 0,
        upper = 1,
        method = "L-BFGS-B"
      )

      if (temp.optim$convergence == 0) {
        beta_type3_vec_q <- temp.optim$par
        hat_beta_temp <- AD_type3_betahat[k, ]
        SEQ <- which(!is.na(hat_beta_temp))
        hat_beta_obs <- as.numeric(hat_beta_temp[SEQ])
        se_hat_beta_temp <- as.numeric(AD_type3_SE_betahat[k, SEQ])
        logNUMDEN2 <- sum(
          stats::dnorm(hat_beta_obs, beta_type3_vec_q[SEQ], se_hat_beta_temp, log = TRUE) -
            stats::dnorm(hat_beta_obs, beta_l_type3_mat[k, SEQ], se_hat_beta_temp, log = TRUE)
        )
        if (stats::runif(n = 1) < exp(logNUMDEN1 + logNUMDEN2)) {
          theta_l_mat[l, ] <- theta_type3_vec_q
          beta_l_type3_mat[k, ] <- beta_type3_vec_q
        }
      }

      if (isTRUE(use_drm[k])) {
        alpha_vec_q <- stats::rnorm(n = p_alpha, mean = alpha_mat_type3[k, ], sd = step_alpha)
        a_n <- DRM_X_type3[[k]][, "a_n"]

        Exp_alpha_psi_q <- exp(as.matrix(DRM_X_type3[[k]][, c("int", "DRM_X_type3")]) %*% alpha_vec_q)
        g_l_mat_q <- cbind(
          a_n * (Exp_alpha_psi_q - 1),
          a_n * (DRM_X_type3[[k]][, "DRM_X_type3"] * Exp_alpha_psi_q - tau_vec[k])
        )
        bar_g_l_q <- apply(g_l_mat_q, 2, mean)
        Sigma_g_l_q <- stats::var(g_l_mat_q) / nrow(g_l_mat_q)

        Exp_alpha_psi <- exp(as.matrix(DRM_X_type3[[k]][, c("int", "DRM_X_type3")]) %*% alpha_mat_type3[k, ])
        g_l_mat <- cbind(
          a_n * (Exp_alpha_psi - 1),
          a_n * (DRM_X_type3[[k]][, "DRM_X_type3"] * Exp_alpha_psi - tau_vec[k])
        )
        bar_g_l <- apply(g_l_mat, 2, mean)
        Sigma_g_l <- stats::var(g_l_mat) / nrow(g_l_mat)

        logNUMDEN1 <- mvtnorm::dmvnorm(bar_g_l_q, mean = rep(0, 2), sigma = Sigma_g_l_q, log = TRUE) -
          mvtnorm::dmvnorm(bar_g_l, mean = rep(0, 2), sigma = Sigma_g_l, log = TRUE)

        ww <- stats::rnorm(n_l, mean = 1, sd = 1)
        temp.optim <- stats::optim(
          beta_l_type3_mat[k, ],
          fn = MB.app.type3.est,
          theta_vec = theta_l_mat[l, ],
          Design_mat = Design_X_type3[[k]],
          bootstrap_w = ww,
          nu_x_mat = nu_X_type3[[k]],
          threshold_type3 = threshold_type3[[k]],
          sig2_IPD_type3 = sig2_IPD_type3[J_IPD + k],
          UseDRM = use_drm[k],
          alpha = alpha_vec_q,
          REF_X_type3_input = DRM_X_type3[[k]],
          control = list(pgtol = 1e-7, maxit = 1000),
          lower = 0,
          upper = 1,
          method = "L-BFGS-B"
        )

        if (temp.optim$convergence == 0) {
          beta_type3_vec_q <- temp.optim$par
          hat_beta_temp <- AD_type3_betahat[k, ]
          SEQ <- which(!is.na(hat_beta_temp))
          hat_beta_obs <- as.numeric(hat_beta_temp[SEQ])
          se_hat_beta_temp <- as.numeric(AD_type3_SE_betahat[k, SEQ])
          logNUMDEN2 <- sum(
            stats::dnorm(hat_beta_obs, beta_type3_vec_q[SEQ], se_hat_beta_temp, log = TRUE) -
              stats::dnorm(hat_beta_obs, beta_l_type3_mat[k, SEQ], se_hat_beta_temp, log = TRUE)
          )
          if (stats::runif(n = 1) < exp(logNUMDEN1 + logNUMDEN2)) {
            alpha_mat_type3[k, ] <- alpha_vec_q
            beta_l_type3_mat[k, ] <- beta_type3_vec_q
          }
        }
      }
    }

    # ---- sig2 IPD ----
    for (l in seq_len(J_IPD)) {
      sig2_q <- stats::rnorm(n = 1, mean = sig2_IPD_type3[l], sd = step_sig2)
      if (sig2_q > 0) {
        y <- IPD_listobj[[l]]$y
        X <- IPD_listobj[[l]]$X
        theta_l_vec <- matrix(theta_l_mat[l, ], ncol = 1)
        logNUMDEN1 <- sum(stats::dnorm(y, X %*% theta_l_vec, sqrt(sig2_q), log = TRUE)) -
          sum(stats::dnorm(y, X %*% theta_l_vec, sqrt(sig2_IPD_type3[l]), log = TRUE))
        logSig2 <- log(truncnorm::dtruncnorm(sig2_q, a = 0, b = Inf, mean = mu_sig2, sd = tau_sig2)) -
          log(truncnorm::dtruncnorm(sig2_IPD_type3[l], a = 0, b = Inf, mean = mu_sig2, sd = tau_sig2))
        if (stats::runif(n = 1) < exp(logNUMDEN1 + logSig2)) {
          sig2_IPD_type3[l] <- sig2_q
        }
      }
    }

    # ---- sig2 Type3 ----
    for (k in seq_len(J_type3)) {
      l <- J_IPD + J_type1 + J_type2 + k
      n_l <- nrow(as.matrix(Design_X_type3[[k]]))
      sig2_q <- stats::rnorm(n = 1, mean = sig2_IPD_type3[J_IPD + k], sd = step_sig2)
      if (sig2_q > 0) {
        ww <- stats::rnorm(n_l, mean = 1, sd = 1)
        temp.optim <- stats::optim(
          beta_l_type3_mat[k, ],
          fn = MB.app.type3.est,
          theta_vec = theta_l_mat[l, ],
          Design_mat = Design_X_type3[[k]],
          bootstrap_w = ww,
          nu_x_mat = nu_X_type3[[k]],
          threshold_type3 = threshold_type3[[k]],
          sig2_IPD_type3 = sig2_q,
          UseDRM = use_drm[k],
          alpha = alpha_mat_type3[k, ],
          REF_X_type3_input = DRM_X_type3[[k]],
          control = list(pgtol = 1e-7, maxit = 1000),
          lower = 0,
          upper = 1,
          method = "L-BFGS-B"
        )
        if (temp.optim$convergence == 0) {
          beta_type3_vec_q <- temp.optim$par
          hat_beta_temp <- AD_type3_betahat[k, ]
          SEQ <- which(!is.na(hat_beta_temp))
          hat_beta_obs <- as.numeric(hat_beta_temp[SEQ])
          se_hat_beta_temp <- as.numeric(AD_type3_SE_betahat[k, SEQ])
          logNUMDEN2 <- sum(
            stats::dnorm(hat_beta_obs, beta_type3_vec_q[SEQ], se_hat_beta_temp, log = TRUE) -
              stats::dnorm(hat_beta_obs, beta_l_type3_mat[k, SEQ], se_hat_beta_temp, log = TRUE)
          )
          logSig2 <- log(truncnorm::dtruncnorm(sig2_q, a = 0, b = Inf, mean = mu_sig2, sd = tau_sig2)) -
            log(truncnorm::dtruncnorm(
              sig2_IPD_type3[J_IPD + k], a = 0, b = Inf, mean = mu_sig2, sd = tau_sig2
            ))
          if (stats::runif(n = 1) < exp(logNUMDEN2 + logSig2)) {
            beta_l_type3_mat[k, ] <- beta_type3_vec_q
            sig2_IPD_type3[J_IPD + k] <- sig2_q
          }
        }
      }
    }

    # ---- mu / Sigma ----
    inv_Var <- inv_Sigma_0_mat + L * inv_Sigma_mat
    Var <- solve(inv_Var)
    theta_sum_vec <- apply(theta_l_mat, 2, sum)
    Mean_temp_part <- inv_Sigma_mat %*% theta_sum_vec
    Mean <- Var %*% Mean_temp_part
    mu_vec <- t(mvtnorm::rmvnorm(n = 1, mean = Mean, sigma = Var))

    df1 <- nu0 + L
    SqSum <- array(0, c(p_theta, p_theta))
    for (l in seq_len(L)) {
      theta_l_vec <- matrix(theta_l_mat[l, ], ncol = 1)
      SqSum <- SqSum + (theta_l_vec - mu_vec) %*% t(theta_l_vec - mu_vec)
    }
    Sigma_mat <- MCMCpack::riwish(v = df1, S = Psi_mat + SqSum)
    inv_Sigma_mat <- solve(Sigma_mat)

    # ---- tau (after burn warmup iterations) ----
    if (J_type3 > 0L && i_iter > tau_update_after) {
      for (k in seq_len(J_type3)) {
        if (isTRUE(use_drm[k])) {
          Exp_alpha_psi <- exp(as.matrix(DRM_X_type3[[k]][, c("int", "DRM_X_type3")]) %*% alpha_mat_type3[k, ])
          a_n <- DRM_X_type3[[k]][, "a_n"]
          tau_q <- stats::rnorm(n = 1, mean = tau_vec[k], sd = step_tau)

          g_l_mat_q <- cbind(
            a_n * (Exp_alpha_psi - 1),
            a_n * (DRM_X_type3[[k]][, "DRM_X_type3"] * Exp_alpha_psi - tau_q)
          )
          bar_g_l_q <- apply(g_l_mat_q, 2, mean)
          Sigma_g_l_q <- stats::var(g_l_mat_q) / nrow(g_l_mat_q)
          logNum <- stats::dnorm(
            hat_tau_type3_vec[k],
            mean = tau_q,
            sd = sqrt(hat_Gamma_tau_type3_vec[k]),
            log = TRUE
          )
          logNum <- logNum + mvtnorm::dmvnorm(bar_g_l_q, mean = rep(0, 2), sigma = Sigma_g_l_q, log = TRUE)

          g_l_mat <- cbind(
            a_n * (Exp_alpha_psi - 1),
            a_n * (DRM_X_type3[[k]][, "DRM_X_type3"] * Exp_alpha_psi - tau_vec[k])
          )
          bar_g_l <- apply(g_l_mat, 2, mean)
          Sigma_g_l <- stats::var(g_l_mat) / nrow(g_l_mat)
          logDen <- stats::dnorm(
            hat_tau_type3_vec[k],
            mean = tau_vec[k],
            sd = sqrt(hat_Gamma_tau_type3_vec[k]),
            log = TRUE
          )
          logDen <- logDen + mvtnorm::dmvnorm(bar_g_l, mean = rep(0, 2), sigma = Sigma_g_l, log = TRUE)

          if (stats::runif(n = 1) < exp(logNum - logDen)) {
            tau_vec[k] <- tau_q
          }
        }
      }
    }

    Draw_mu_vec[i_iter, ] <- as.numeric(mu_vec)
    Draw_Sigma_diag[i_iter, ] <- diag(Sigma_mat)
    Draw_sig2[i_iter, ] <- sig2_IPD_type3

    if (verbose && (i_iter %% as.integer(progress_every) == 0L)) {
      cur <- proc.time()[[3]]
      last_batch <- cur - prev_time
      eta <- (n_iter - i_iter) * (last_batch / progress_every)
      prev_time <- cur
      message(sprintf(
        "iter %d / %d | last %d = %.1f min | ETA = %.1f min",
        i_iter, n_iter, progress_every, last_batch / 60, eta / 60
      ))
    }
  }

  seq_keep <- if (burnin == 0L) seq_len(mainrun) else (burnin + 1L):(burnin + mainrun)
  out <- list(
    posterior_mu = Draw_mu_vec[seq_keep, , drop = FALSE],
    posterior_Sigma_diag = Draw_Sigma_diag[seq_keep, , drop = FALSE],
    posterior_sig2 = Draw_sig2[seq_keep, , drop = FALSE],
    call = match.call(),
    settings = list(
      burnin = burnin,
      mainrun = mainrun,
      seed = seed,
      L = L,
      J = J_IPD,
      K = J_type1 + J_type2 + J_type3,
      J_type1 = J_type1,
      J_type2 = J_type2,
      J_type3 = J_type3,
      n = NA_integer_,
      p = p_theta,
      outcome = "gaussian_application",
      used_default_data = FALSE
    )
  )
  class(out) <- c("bayesmetaipd_fit", "list")
  out
}


#' Small synthetic Application-style dataset
#'
#' Generates a tiny continuous IPD + Type1/2/3 AD example for testing and demos.
#' Not the real I-WIP Application data. Uses enough studies so that
#' `nu0 + L >= p` for the Inverse-Wishart update (`p = 7`).
#'
#' @param seed Integer seed. Default `1`.
#' @param n_ipd Number of IPD studies. Default `8`.
#' @return List with `ipd`, `ad_type1`, `ad_type2`, `ad_type3`.
#' @keywords internal
#' @noRd
example_application_data <- function(seed = 1L, n_ipd = 8L) {
  set.seed(as.integer(seed))
  p <- 7L
  n_ipd <- as.integer(n_ipd)
  if (n_ipd < 7L) {
    stop("`n_ipd` must be >= 7 so that nu0 + L can exceed p = 7.", call. = FALSE)
  }
  mu <- c(12, 11, 10, 0.2, -1.5, -1.2, -1.0)
  make_ipd <- function(n, theta) {
    bmi <- sample(1:3, n, replace = TRUE)
    X <- cbind(
      bmi_cat1 = as.numeric(bmi == 1),
      bmi_cat2 = as.numeric(bmi == 2),
      bmi_cat3 = as.numeric(bmi == 3),
      b_wt = stats::rnorm(n, 70, 8),
      bmi_cat1.trt = 0,
      bmi_cat2.trt = 0,
      bmi_cat3.trt = 0
    )
    trt <- stats::rbinom(n, 1, 0.5)
    X[, "bmi_cat1.trt"] <- X[, "bmi_cat1"] * trt
    X[, "bmi_cat2.trt"] <- X[, "bmi_cat2"] * trt
    X[, "bmi_cat3.trt"] <- X[, "bmi_cat3"] * trt
    y <- as.numeric(X %*% theta + stats::rnorm(n, 0, sqrt(8)))
    list(y = y, X = X, n_l = n)
  }

  ipd <- vector("list", n_ipd)
  for (i in seq_len(n_ipd)) {
    ipd[[i]] <- make_ipd(25L, mu + stats::rnorm(p, 0, 0.25))
  }

  ad_type1 <- data.frame(beta_hat = -0.9, sqrt_hat_V = 0.25)

  ad_type2 <- data.frame(
    p_01 = 0.3, p_02 = 0.4, p_03 = 0.3,
    p_11 = 0.35, p_12 = 0.35, p_13 = 0.3,
    beta_hat = -0.8, sqrt_hat_V = 0.3
  )

  ref <- as.data.frame(ipd[[1]]$X)
  ref$new_trt <- as.numeric(
    ref$bmi_cat1.trt + ref$bmi_cat2.trt + ref$bmi_cat3.trt > 0
  )
  ref$b_bmi <- 22 + 4 * ref$bmi_cat2 + 8 * ref$bmi_cat3 + stats::rnorm(nrow(ref), 0, 1)
  keep <- which(rowSums(ref[, c("bmi_cat1", "bmi_cat2", "bmi_cat3")]) == 1)
  ref <- ref[keep, , drop = FALSE]
  Design_X <- as.matrix(ref[, c(
    "bmi_cat1", "bmi_cat2", "bmi_cat3", "b_wt",
    "bmi_cat1.trt", "bmi_cat2.trt", "bmi_cat3.trt"
  )])
  s1 <- (ref$new_trt == 0) * (ref$bmi_cat1 == 1)
  s2 <- (ref$new_trt == 0) * (ref$bmi_cat2 == 1)
  s3 <- (ref$new_trt == 0) * (ref$bmi_cat3 == 1)
  s4 <- (ref$new_trt == 1) * (ref$bmi_cat1 == 1)
  s5 <- (ref$new_trt == 1) * (ref$bmi_cat2 == 1)
  s6 <- (ref$new_trt == 1) * (ref$bmi_cat3 == 1)
  nu_X <- cbind(s1, s2, s3, s4, s5, s6)
  ok <- which(rowSums(nu_X) == 1)
  Design_X <- Design_X[ok, , drop = FALSE]
  nu_X <- nu_X[ok, , drop = FALSE]
  ref <- ref[ok, , drop = FALSE]
  thr <- ifelse(Design_X[, "bmi_cat1"] == 1, 15.9,
                ifelse(Design_X[, "bmi_cat2"] == 1, 11.3, 9.1))

  n_AD <- c(n_AD = 80, n_AD1 = 25, n_AD2 = 30, n_AD3 = 25)
  prop_ipd <- colMeans(ref[, c("bmi_cat1", "bmi_cat2", "bmi_cat3")])
  prop_ad <- n_AD[c("n_AD1", "n_AD2", "n_AD3")] / n_AD["n_AD"]
  a_n <- rep(0, nrow(ref))
  for (j in 1:3) {
    seqj <- which(ref[, j] == 1)
    if (length(seqj) && prop_ipd[j] > 0) a_n[seqj] <- as.numeric(prop_ad[j] / prop_ipd[j])
  }
  DRM_X <- data.frame(
    bmi_cat1 = ref$bmi_cat1,
    bmi_cat2 = ref$bmi_cat2,
    bmi_cat3 = ref$bmi_cat3,
    DRM_X_type3 = ref$b_bmi,
    int = 1,
    a_n = a_n
  )

  true_beta <- c(0.55, 0.45, 0.4, 0.5, 0.42, 0.35)
  ad_type3 <- list(
    betahat = matrix(true_beta, nrow = 1),
    se = matrix(0.08, nrow = 1, ncol = 6),
    use_drm = TRUE,
    hat_tau = mean(ref$b_bmi),
    hat_gamma_tau = stats::var(ref$b_bmi) / (n_AD["n_AD"] - 1),
    Design_X = list(Design_X),
    nu_X = list(nu_X),
    threshold = list(thr),
    DRM_X = list(DRM_X)
  )

  list(ipd = ipd, ad_type1 = ad_type1, ad_type2 = ad_type2, ad_type3 = ad_type3)
}
