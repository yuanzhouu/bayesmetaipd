# Internal helpers for logistic IPD+AD bridge / density ratio
logit_inv <- function(alpha) {
  exp(alpha) / (exp(alpha) + 1)
}

# Official MB.logit.est from Code/Simulation_2/2_IPD-AD.R
MB.logit.est <- function(param, theta, x, w, alpha) {
  x <- as.matrix(x)
  theta <- as.vector(theta)
  param <- as.vector(param)
  alpha <- as.vector(alpha)
  w.DRM <- exp(x[, 1:2] %*% alpha)
  w <- w.DRM * w
  resid <- w * (logit_inv(x %*% theta) - logit_inv(x[, -4] %*% param))
  g.i <- cbind(x[, 1] * resid, x[, 2] * resid, x[, 3] * resid)
  g.sum <- apply(g.i, 2, sum)
  sum(g.sum^2)
}

load_sim2_ipdad_rep1 <- function() {
  sim2_ipdad_rep1 <- NULL
  if (exists("sim2_ipdad_rep1", envir = asNamespace("bayesmetaipd"), inherits = FALSE)) {
    sim2_ipdad_rep1 <- get("sim2_ipdad_rep1", envir = asNamespace("bayesmetaipd"), inherits = FALSE)
  }
  if (is.null(sim2_ipdad_rep1)) {
    env <- new.env(parent = emptyenv())
    utils::data("sim2_ipdad_rep1", package = "bayesmetaipd", envir = env)
    if (exists("sim2_ipdad_rep1", envir = env, inherits = FALSE)) {
      sim2_ipdad_rep1 <- env$sim2_ipdad_rep1
    }
  }
  if (is.null(sim2_ipdad_rep1)) {
    path <- system.file("extdata", "sim2_ipdad_rep1.rda", package = "bayesmetaipd", mustWork = FALSE)
    if (!nzchar(path) || !file.exists(path)) {
      stop("Could not load bundled dataset sim2_ipdad_rep1.", call. = FALSE)
    }
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    sim2_ipdad_rep1 <- env$sim2_ipdad_rep1
  }
  sim2_ipdad_rep1
}


#' Fit IPD + AD Bayesian random-effects meta-analysis (logistic)
#'
#' Implements the Simulation Study 2 proposed method
#' (`Code/Simulation_2/2_IPD-AD.R`): integrate IPD and AD under a logistic
#' random-effects model with density-ratio adjustment for covariate shift.
#'
#' **Defaults reproduce the official Simulation 2 IPD+AD `rep_1` result** when
#' called as `fit_ipd_ad()` (bundled [sim2_ipdad_rep1] data, `seed = 1001`,
#' burn-in/main = 10000).
#'
#' @param X Optional `L x n x p` covariate array. If `NULL`, uses
#'   [sim2_ipdad_rep1].
#' @param Y Optional `L x n` binary response matrix. If `NULL`, uses
#'   [sim2_ipdad_rep1].
#' @param beta_mat Optional `L x p` matrix of AD working-model estimates.
#' @param V_beta_cube Optional `L x p x p` array of AD working variances.
#' @param is_ipd Optional length-`L` logical/0-1 vector; `TRUE`/`1` = IPD study.
#'   Default uses `sim2_ipdad_rep1$delta_biased_access`.
#' @param theta_init Optional `L x p` starting values for study-specific
#'   coefficients.
#' @param mu_init Optional length-`p` starting value for `mu`.
#' @param Sigma_init Optional `p x p` starting covariance. Default `diag(1, p)`.
#' @param burnin Integer burn-in iterations. Default `10000`.
#' @param mainrun Integer post-burn-in iterations retained. Default `10000`.
#' @param step_theta_ipd MH proposal SD for IPD `theta`. Default `0.15`.
#' @param step_theta_ad Length-`p` MH proposal SDs for AD `theta`.
#'   Default `c(0.2, 0.1, 0.1, 0.2) * 1.5`.
#' @param step_alpha MH proposal SD for density-ratio `alpha`. Default `0.005`.
#' @param step_tau MH proposal SD for `tau`. Default `0.1`.
#' @param lambda Prior variance for `mu ~ N(0, lambda * I)`. Default `1e4`.
#' @param nu0 Inverse-Wishart degrees of freedom. Default `0.1`.
#' @param phi0 Inverse-Wishart scale multiplier. Default `0.1`.
#' @param seed Integer RNG seed. Default `1001`. With default bundled data,
#'   `seed = 1001` restores the official `.Random.seed` from
#'   `SimulationData_2.RData`. Set `NULL` to leave RNG unchanged.
#' @param verbose Logical; print progress every 1000 iterations.
#'
#' @return A list of class `bayesmetaipd_fit` with `posterior_mu`,
#'   `posterior_Sigma_diag`, `call`, and `settings`.
#'
#' @examples
#' \dontrun{
#' fit <- fit_ipd_ad()
#' colMeans(fit$posterior_mu)
#' }
#'
#' fit_short <- fit_ipd_ad(burnin = 5, mainrun = 10, verbose = FALSE)
#' str(fit_short$posterior_mu)
#'
#' @export
fit_ipd_ad <- function(X = NULL,
                       Y = NULL,
                       beta_mat = NULL,
                       V_beta_cube = NULL,
                       is_ipd = NULL,
                       theta_init = NULL,
                       mu_init = NULL,
                       Sigma_init = NULL,
                       burnin = 10000L,
                       mainrun = 10000L,
                       step_theta_ipd = 0.15,
                       step_theta_ad = c(0.2, 0.1, 0.1, 0.2) * 1.5,
                       step_alpha = 0.005,
                       step_tau = 0.1,
                       lambda = 1e4,
                       nu0 = 0.1,
                       phi0 = 0.1,
                       seed = 1001L,
                       verbose = TRUE) {
  using_default_data <- is.null(X) && is.null(Y) &&
    is.null(beta_mat) && is.null(V_beta_cube) && is.null(is_ipd)
  bundled_random_seed <- NULL

  if (using_default_data) {
    dat <- load_sim2_ipdad_rep1()
    X <- dat$X_cube
    Y <- dat$Y_mat
    beta_mat <- dat$beta_mat
    V_beta_cube <- dat$V_beta_cube
    is_ipd <- dat$delta_biased_access
    if (is.null(theta_init)) theta_init <- dat$theta_l_mat
    if (is.null(mu_init)) mu_init <- dat$true_mu
    if (!is.null(dat$random_seed)) bundled_random_seed <- dat$random_seed
  } else {
    if (is.null(X) || is.null(Y) || is.null(beta_mat) ||
        is.null(V_beta_cube) || is.null(is_ipd)) {
      stop(
        "Provide X, Y, beta_mat, V_beta_cube, and is_ipd, ",
        "or leave all NULL to use sim2_ipdad_rep1.",
        call. = FALSE
      )
    }
  }

  if (!is.array(X) || length(dim(X)) != 3L) {
    stop("`X` must be a 3-D array with dimensions L x n x p.", call. = FALSE)
  }
  Y <- as.matrix(Y)
  beta_mat <- as.matrix(beta_mat)

  L <- dim(X)[1]
  n_sample <- dim(X)[2]
  p_theta <- dim(X)[3]
  p_beta <- p_theta - 1L
  p_alpha <- 2L

  if (!all(dim(Y) == c(L, n_sample))) {
    stop("`Y` dimensions must match the first two dimensions of `X`.", call. = FALSE)
  }
  if (!all(dim(beta_mat) == c(L, p_theta))) {
    stop("`beta_mat` must be L x p.", call. = FALSE)
  }
  if (!is.array(V_beta_cube) || !all(dim(V_beta_cube) == c(L, p_theta, p_theta))) {
    stop("`V_beta_cube` must be L x p x p.", call. = FALSE)
  }
  if (length(is_ipd) != L) {
    stop("`is_ipd` must have length L.", call. = FALSE)
  }
  is_ipd <- as.integer(is_ipd)

  if (is.null(theta_init)) theta_init <- matrix(0, nrow = L, ncol = p_theta)
  if (is.null(mu_init)) mu_init <- rep(0, p_theta)
  if (is.null(Sigma_init)) Sigma_init <- diag(1, p_theta)

  theta_init <- as.matrix(theta_init)
  mu_init <- as.numeric(mu_init)
  Sigma_init <- as.matrix(Sigma_init)
  if (!all(dim(theta_init) == c(L, p_theta))) stop("`theta_init` must be L x p.", call. = FALSE)
  if (length(mu_init) != p_theta) stop("`mu_init` must have length p.", call. = FALSE)
  if (!all(dim(Sigma_init) == c(p_theta, p_theta))) stop("`Sigma_init` must be p x p.", call. = FALSE)
  if (length(step_theta_ad) != p_theta) {
    stop("`step_theta_ad` must have length p.", call. = FALSE)
  }

  burnin <- as.integer(burnin)
  mainrun <- as.integer(mainrun)
  if (burnin < 0L || mainrun < 1L) {
    stop("`burnin` must be >= 0 and `mainrun` must be >= 1.", call. = FALSE)
  }

  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (using_default_data && identical(seed, 1001L) && !is.null(bundled_random_seed)) {
      assign(".Random.seed", bundled_random_seed, envir = .GlobalEnv)
    } else {
      set.seed(seed)
    }
  }

  # Split IPD / AD (match official indexing)
  SEQ_IPD <- which(is_ipd == 1L)
  J <- length(SEQ_IPD)
  X_IPD <- X[SEQ_IPD, , , drop = FALSE]
  Y_mat <- Y[SEQ_IPD, , drop = FALSE]
  SEQ_AD <- which(is_ipd == 0L)
  K <- length(SEQ_AD)
  if (J < 1L || K < 1L) {
    stop("Need at least one IPD and one AD study.", call. = FALSE)
  }

  beta_tilde_mat <- beta_mat[SEQ_AD, , drop = FALSE]
  V_tilde_cube <- V_beta_cube[SEQ_AD, , , drop = FALSE]
  for (kk in seq_len(K)) {
    V_tilde_cube[kk, , ] <- diag(diag(V_tilde_cube[kk, , ]))
  }

  # Reference data for density ratio model
  tilde_D_x <- vector("list", K)
  hat_tau_vec <- hat_Gamma_tau_vec <- rep(0, K)
  mean_X1_IPD <- rep(0, J)
  for (j in seq_len(J)) {
    mean_X1_IPD[j] <- mean(X_IPD[j, , "X1"])
  }
  for (k in seq_len(K)) {
    i_study <- SEQ_AD[k]
    X_AD_l <- X[i_study, , ]
    hat_tau_vec[k] <- mean(X_AD_l[, "X1"])
    hat_Gamma_tau_vec[k] <- stats::var(X_AD_l[, "X1"]) / length(X_AD_l[, "X1"])
    WHICH <- which.min(abs(hat_tau_vec[k] - mean_X1_IPD))
    tilde_D_x[[k]] <- X_IPD[WHICH, , ]
  }

  invLambda_theta <- diag(1 / lambda, p_theta)
  Phi0 <- diag(phi0, p_theta)

  theta_mat_IPD <- theta_init[SEQ_IPD, , drop = FALSE]
  mu_vec <- mu_init
  Sigma_theta_mat <- Sigma_init
  tau_vec <- hat_Gamma_tau_vec
  alpha_mat <- array(0.1, c(K, p_alpha))
  theta_mat_AD <- theta_init[SEQ_AD, , drop = FALSE]
  beta_mat_AD <- beta_tilde_mat

  n_iter <- burnin + mainrun
  draw_mu <- array(0, c(n_iter, p_theta))
  draw_Sigma_theta <- array(0, c(n_iter, p_theta))

  if (verbose) {
    message(sprintf(
      "Starting IPD+AD MCMC: %d iterations (%d burn-in + %d main); J=%d IPD, K=%d AD",
      n_iter, burnin, mainrun, J, K
    ))
  }
  prev_time <- proc.time()[[3]]

  for (i_iter in seq_len(n_iter)) {
    # ---- Update AD theta, beta, alpha ----
    for (k in seq_len(K)) {
      # Step 2: update theta and beta
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_AD[k, ], sd = step_theta_ad)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      temp.optim <- stats::optim(
        beta_mat_AD[k, 1:3],
        fn = MB.logit.est,
        theta = theta_vec_q,
        x = tilde_D_x[[k]],
        w = ww,
        control = list(reltol = 1e-7, maxit = 1000),
        alpha = alpha_mat[k, ]
      )

      if (temp.optim$convergence == 0) {
        beta_vec_q <- temp.optim$par
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3],
          mean = beta_vec_q[2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3],
          log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(
          theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3],
          mean = beta_mat_AD[k, 2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3],
          log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(
          theta_mat_AD[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logAcc <- logNum - logDen
        if (stats::runif(n = 1) < exp(logAcc)) {
          theta_mat_AD[k, ] <- theta_vec_q
          beta_mat_AD[k, 1:3] <- beta_vec_q
        }
      }

      # Step 3: update alpha and beta
      alpha_vec_q <- stats::rnorm(n = p_alpha, mean = alpha_mat[k, ], sd = step_alpha)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      temp.optim <- stats::optim(
        beta_mat_AD[k, 1:3],
        fn = MB.logit.est,
        theta = theta_mat_AD[k, ],
        x = tilde_D_x[[k]],
        w = ww,
        control = list(reltol = 1e-7, maxit = 1000),
        alpha = alpha_vec_q
      )

      if (temp.optim$convergence == 0) {
        beta_vec_q <- temp.optim$par
        g_l1_vec_q <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_vec_q)
        g_l_mat_q <- cbind(
          g_l1_vec_q - 1,
          tilde_D_x[[k]][, "X1"] * g_l1_vec_q - tau_vec[k]
        )
        bar_g_l_q <- apply(g_l_mat_q, 2, mean)
        Sigma_g_l_q <- stats::var(g_l_mat_q) / nrow(g_l_mat_q)
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3],
          mean = beta_vec_q[2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3],
          log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(
          bar_g_l_q, mean = rep(0, 2), sigma = Sigma_g_l_q, log = TRUE
        )

        g_l1_vec <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
        g_l_mat <- cbind(
          g_l1_vec - 1,
          tilde_D_x[[k]][, "X1"] * g_l1_vec - tau_vec[k]
        )
        bar_g_l <- apply(g_l_mat, 2, mean)
        Sigma_g_l <- stats::var(g_l_mat) / nrow(g_l_mat)
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3],
          mean = beta_mat_AD[k, 2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3],
          log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(
          bar_g_l, mean = rep(0, 2), sigma = Sigma_g_l, log = TRUE
        )
        logAcc <- logNum - logDen
        if (stats::runif(n = 1) < exp(logAcc)) {
          alpha_mat[k, ] <- alpha_vec_q
          beta_mat_AD[k, 1:3] <- beta_vec_q
        }
      }

      # keep g_l recomputed as in official script (side effect unused later)
      g_l1_vec <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
      g_l_mat <- cbind(
        g_l1_vec - 1,
        tilde_D_x[[k]][, "X1"] * g_l1_vec - tau_vec[k]
      )
      bar_g_l <- apply(g_l_mat, 2, mean)
      Sigma_g_l <- stats::var(g_l_mat) / nrow(g_l_mat)
    }

    # ---- Update IPD theta (Step 1) ----
    for (j in seq_len(J)) {
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_IPD[j, ], sd = step_theta_ipd)
      x.the_q <- X_IPD[j, , ] %*% theta_vec_q
      x.the <- X_IPD[j, , ] %*% theta_mat_IPD[j, ]
      logAcc <- sum(
        x.the_q * Y_mat[j, ] - log(1 + exp(x.the_q)) -
          x.the * Y_mat[j, ] + log(1 + exp(x.the))
      )
      logAcc <- logAcc +
        mvtnorm::dmvnorm(theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE) -
        mvtnorm::dmvnorm(theta_mat_IPD[j, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
      if (stats::runif(n = 1) < exp(logAcc)) {
        theta_mat_IPD[j, ] <- theta_vec_q
      }
    }

    theta_merge <- rbind(theta_mat_IPD, theta_mat_AD)

    # ---- Update mu (Step 5) ----
    inv_Sigma <- solve(Sigma_theta_mat)
    inv_Var <- invLambda_theta + L * inv_Sigma
    Var <- solve(inv_Var)
    Mean_2ndpart <- inv_Sigma %*% apply(theta_merge, 2, sum)
    Mean <- Var %*% Mean_2ndpart
    mu_vec <- mvtnorm::rmvnorm(n = 1, mean = Mean, sigma = Var)

    # ---- Update Sigma (Step 5) ----
    SS <- array(0, c(p_theta, p_theta))
    for (l in seq_len(L)) {
      SS <- SS + t(theta_merge[l, ] - mu_vec) %*% t(t(theta_merge[l, ] - mu_vec))
    }
    Sigma_theta_mat <- MCMCpack::riwish((nu0 + L), (Phi0 + SS))

    # ---- Update tau (Step 4) ----
    for (k in seq_len(K)) {
      g_l1_vec <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
      tau_q <- stats::rnorm(n = 1, mean = tau_vec[k], sd = step_tau)
      g_l_mat_q <- cbind(g_l1_vec - 1, tilde_D_x[[k]][, "X1"] * g_l1_vec - tau_q)
      bar_g_l_q <- apply(g_l_mat_q, 2, mean)
      Sigma_g_l_q <- stats::var(g_l_mat_q) / nrow(g_l_mat_q)
      logNum <- stats::dnorm(hat_tau_vec[k], mean = tau_q, sd = sqrt(hat_Gamma_tau_vec[k]), log = TRUE)
      logNum <- logNum + mvtnorm::dmvnorm(bar_g_l_q, mean = rep(0, 2), sigma = Sigma_g_l_q, log = TRUE)
      g_l_mat <- cbind(g_l1_vec - 1, tilde_D_x[[k]][, "X1"] * g_l1_vec - tau_vec[k])
      bar_g_l <- apply(g_l_mat, 2, mean)
      Sigma_g_l <- stats::var(g_l_mat) / nrow(g_l_mat)
      logDen <- stats::dnorm(hat_tau_vec[k], mean = tau_vec[k], sd = sqrt(hat_Gamma_tau_vec[k]), log = TRUE)
      logDen <- logDen + mvtnorm::dmvnorm(bar_g_l, mean = rep(0, 2), sigma = Sigma_g_l, log = TRUE)
      logAcc <- logNum - logDen
      if (stats::runif(n = 1) < exp(logAcc)) {
        tau_vec[k] <- tau_q
      }
    }

    draw_mu[i_iter, ] <- mu_vec
    draw_Sigma_theta[i_iter, ] <- diag(Sigma_theta_mat)

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

  seq_keep <- if (burnin == 0L) {
    seq_len(mainrun)
  } else {
    (burnin + 1L):(burnin + mainrun)
  }

  out <- list(
    posterior_mu = draw_mu[seq_keep, , drop = FALSE],
    posterior_Sigma_diag = draw_Sigma_theta[seq_keep, , drop = FALSE],
    call = match.call(),
    settings = list(
      burnin = burnin,
      mainrun = mainrun,
      step_theta_ipd = step_theta_ipd,
      step_theta_ad = step_theta_ad,
      step_alpha = step_alpha,
      step_tau = step_tau,
      lambda = lambda,
      nu0 = nu0,
      phi0 = phi0,
      seed = seed,
      L = L,
      J = J,
      K = K,
      n = n_sample,
      p = p_theta,
      used_default_data = using_default_data
    )
  )
  class(out) <- c("bayesmetaipd_fit", "list")
  out
}


#' Reproduce official Simulation 2 IPD+AD (`rep_1`) draws
#'
#' Convenience wrapper around [fit_ipd_ad()] with the exact default settings
#' used in `Code/Simulation_2/2_IPD-AD.R`.
#'
#' @param ... Passed to [fit_ipd_ad()] (e.g. `verbose = FALSE`).
#' @param compare_official If `TRUE`, compare to bundled official posterior draws.
#'
#' @return A `bayesmetaipd_fit` fit object. If `compare_official = TRUE`, an
#'   attribute `comparison` is attached.
#'
#' @examples
#' \dontrun{
#' fit <- reproduce_sim2_ipdad(compare_official = TRUE)
#' attr(fit, "comparison")
#' }
#'
#' @export
reproduce_sim2_ipdad <- function(..., compare_official = FALSE) {
  fit <- fit_ipd_ad(...)

  if (isTRUE(compare_official)) {
    official_path <- system.file(
      "extdata", "official_sim2_ipdad_rep1.RData",
      package = "bayesmetaipd"
    )
    if (!nzchar(official_path) || !file.exists(official_path)) {
      warning("Official IPD+AD reference file not found in the installed package.", call. = FALSE)
    } else {
      env <- new.env(parent = emptyenv())
      load(official_path, envir = env)
      comparison <- list(
        mu_all_equal = isTRUE(all.equal(fit$posterior_mu, env$posterior_mu)),
        sigma_all_equal = isTRUE(all.equal(
          fit$posterior_Sigma_diag, env$posterior_Sigma_theta
        )),
        mu_max_abs_diff = max(abs(fit$posterior_mu - env$posterior_mu)),
        sigma_max_abs_diff = max(abs(
          fit$posterior_Sigma_diag - env$posterior_Sigma_theta
        ))
      )
      attr(fit, "comparison") <- comparison
      message(
        "Comparison vs official IPD+AD rep_1: ",
        "mu all.equal=", comparison$mu_all_equal,
        ", Sigma all.equal=", comparison$sigma_all_equal,
        ", mu max|diff|=", format(comparison$mu_max_abs_diff, digits = 3),
        ", Sigma max|diff|=", format(comparison$sigma_max_abs_diff, digits = 3)
      )
    }
  }

  fit
}
