# Simulation Study 1 moment bridges (Code/Simulation_1/2_IPD-AD.R)

MB.est.ADtype1 <- function(theta, x, w, alpha) {
  Wn <- diag(w)
  Xn <- as.matrix(x)
  subgroup.Xn <- as.matrix(x[, -4])
  theta <- as.vector(theta)
  alpha <- as.vector(alpha)
  w.DRM <- diag(as.vector(exp(as.matrix(x[, 1:2]) %*% alpha)))
  Wn <- Wn * w.DRM
  Mn1 <- t(subgroup.Xn) %*% Wn %*% subgroup.Xn
  Mn2 <- t(subgroup.Xn) %*% Wn %*% Xn %*% theta
  solve(Mn1) %*% Mn2
}

MB.est.ADtype2 <- function(theta, x, w, subgroup.x, alpha) {
  Wn <- diag(w)
  Xn <- as.matrix(x)
  subgroup.Xn <- as.matrix(subgroup.x)
  theta <- as.vector(theta)
  alpha <- as.vector(alpha)
  w.DRM <- diag(as.vector(exp(as.matrix(x[, 1:2]) %*% alpha)))
  Wn <- Wn * w.DRM
  Mn1 <- t(subgroup.Xn) %*% Wn %*% subgroup.Xn
  Mn2 <- t(subgroup.Xn) %*% Wn %*% Xn %*% theta
  solve(Mn1) %*% Mn2
}

MB.est.ADtype3 <- function(theta, x, w, alpha) {
  Wn <- diag(w)
  Xn <- as.matrix(x)
  theta <- as.vector(theta)
  alpha <- as.vector(alpha)
  w.DRM <- diag(as.vector(exp(as.matrix(x[, 1:2]) %*% alpha)))
  Wn <- Wn * w.DRM
  Mn1 <- t(Xn) %*% Wn %*% Xn
  Mn2 <- t(Xn) %*% Wn %*% Xn %*% theta
  solve(Mn1) %*% Mn2
}

sim1_subgroup_indicators <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  phi <- integer(n)
  for (i in seq_len(n)) {
    if (x[i, "X1"] > 0 && x[i, "X2"] == 0) {
      phi[i] <- 1L
    } else if (x[i, "X1"] > 0 && x[i, "X2"] == 1) {
      phi[i] <- 2L
    } else if (x[i, "X1"] <= 0 && x[i, "X2"] == 0) {
      phi[i] <- 3L
    } else {
      phi[i] <- 4L
    }
  }
  as.data.frame(cbind(
    ind.1 = 1 * (phi == 1),
    ind.2 = 1 * (phi == 2),
    ind.3 = 1 * (phi == 3),
    ind.4 = 1 * (phi == 4)
  ))
}

load_sim1_ipdad_rep1 <- function() {
  sim1_ipdad_rep1 <- NULL
  if (exists("sim1_ipdad_rep1", envir = asNamespace("bayesmetaipd"), inherits = FALSE)) {
    sim1_ipdad_rep1 <- get("sim1_ipdad_rep1", envir = asNamespace("bayesmetaipd"), inherits = FALSE)
  }
  if (is.null(sim1_ipdad_rep1)) {
    env <- new.env(parent = emptyenv())
    utils::data("sim1_ipdad_rep1", package = "bayesmetaipd", envir = env)
    if (exists("sim1_ipdad_rep1", envir = env, inherits = FALSE)) {
      sim1_ipdad_rep1 <- env$sim1_ipdad_rep1
    }
  }
  if (is.null(sim1_ipdad_rep1)) {
    path <- system.file("extdata", "sim1_ipdad_rep1.rda", package = "bayesmetaipd", mustWork = FALSE)
    if (!nzchar(path) || !file.exists(path)) {
      stop("Could not load bundled dataset sim1_ipdad_rep1.", call. = FALSE)
    }
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    sim1_ipdad_rep1 <- env$sim1_ipdad_rep1
  }
  sim1_ipdad_rep1
}


#' Fit Simulation Study 1 IPD + AD (continuous / Gaussian)
#'
#' Implements `Code/Simulation_1/2_IPD-AD.R`: Gaussian linear IPD, three AD
#' report types (nested model, subgroup means, full model with partial
#' coefficients), and density-ratio adjustment for covariate shift.
#'
#' **Defaults use bundled [sim1_ipdad_rep1]** (Simulation 1, replicate 1:
#' 10 Type1 + 10 Type2 + 10 Type3 AD, 10 IPD). `seed = 1001` restores the
#' generator `.Random.seed` stored with that replicate.
#'
#' @param X Optional `L x n x p` covariate array. If `NULL`, uses
#'   [sim1_ipdad_rep1]. Columns should be named `Int`, `X1`, `X2`, `X1X2`.
#' @param Y Optional `L x n` continuous response matrix.
#' @param beta_mat Optional `L x p` AD working-model estimates.
#' @param V_beta_cube Optional `L x p x p` AD working covariances.
#' @param type_vec Optional length-`L` integer vector: `1` Type1 AD, `2` Type2
#'   AD, `3` Type3 AD, `4` IPD. Default official split `c(rep(1,10), rep(2,10),
#'   rep(3,10), rep(4,10))`.
#' @param theta_init Optional `L x p` starting coefficients.
#' @param mu_init Optional length-`p` starting `mu`.
#' @param Sigma_init Optional `p x p` starting covariance. Default `diag(1, p)`.
#' @param sig2_init Starting shared residual variance. Default `true_kappa`
#'   from bundled data, otherwise `1`.
#' @param burnin,mainrun MCMC lengths. Defaults `10000`.
#' @param step_theta MH SD for `theta` (IPD and AD). Default `0.2`.
#' @param step_alpha MH SD for density-ratio `alpha`. Default `0.01`.
#' @param step_tau MH SD for `tau`. Default `0.02`.
#' @param lambda Prior variance for `mu ~ N(0, lambda * I)`. Default `1e4`.
#' @param nu0 Inverse-Wishart df. Default `0.1`.
#' @param phi0 Inverse-Wishart scale multiplier. Default `0.1`.
#' @param seed Integer seed. Default `1001`. With bundled data this restores
#'   `sim1_ipdad_rep1$random_seed`. `NULL` leaves RNG unchanged.
#' @param verbose Print progress every 1000 iterations.
#' @param engine `"r"` (default) uses the official-style R sampler (bit-identical
#'   to `Code/Simulation_1/2_IPD-AD.R` when seeds match). `"cpp"` uses the
#'   RcppArmadillo inner loop (same model, faster; multivariate draws use a
#'   Cholesky / Bartlett implementation so chains are not bit-identical to R).
#'
#' @return A `bayesmetaipd_fit` with `posterior_mu`, `posterior_Sigma_diag`,
#'   `posterior_sig2` (vector of shared residual variance), `call`, `settings`.
#'
#' @examples
#' fit_short <- fit_ipd_ad_sim1(burnin = 2, mainrun = 3, verbose = FALSE)
#' colMeans(fit_short$posterior_mu)
#'
#' @keywords internal
#' @noRd
fit_ipd_ad_sim1 <- function(X = NULL,
                            Y = NULL,
                            beta_mat = NULL,
                            V_beta_cube = NULL,
                            type_vec = NULL,
                            theta_init = NULL,
                            mu_init = NULL,
                            Sigma_init = NULL,
                            sig2_init = NULL,
                            burnin = 10000L,
                            mainrun = 10000L,
                            step_theta = 0.2,
                            step_alpha = 0.01,
                            step_tau = 0.02,
    lambda = 1e4,
    nu0 = 0.1,
    phi0 = 0.1,
    seed = 1001L,
    verbose = TRUE,
    engine = c("r", "cpp")) {
  using_default_data <- is.null(X) && is.null(Y) &&
    is.null(beta_mat) && is.null(V_beta_cube) && is.null(type_vec)
  bundled_random_seed <- NULL
  true_kappa <- 1
  engine <- match.arg(engine)

  if (using_default_data) {
    dat <- load_sim1_ipdad_rep1()
    X <- dat$X_cube
    Y <- dat$Y_mat
    beta_mat <- dat$beta_mat
    V_beta_cube <- dat$V_beta_cube
    type_vec <- dat$type_vec
    if (is.null(theta_init)) theta_init <- dat$theta_l_mat
    if (is.null(mu_init)) mu_init <- dat$true_mu
    if (!is.null(dat$true_kappa)) true_kappa <- dat$true_kappa
    if (!is.null(dat$random_seed)) bundled_random_seed <- dat$random_seed
  } else {
    if (is.null(X) || is.null(Y) || is.null(beta_mat) ||
        is.null(V_beta_cube) || is.null(type_vec)) {
      stop(
        "Provide X, Y, beta_mat, V_beta_cube, and type_vec, ",
        "or leave all NULL to use sim1_ipdad_rep1.",
        call. = FALSE
      )
    }
  }

  if (!is.array(X) || length(dim(X)) != 3L) {
    stop("`X` must be a 3-D array with dimensions L x n x p.", call. = FALSE)
  }
  Y <- as.matrix(Y)
  beta_mat <- as.matrix(beta_mat)
  type_vec <- as.integer(type_vec)

  L <- dim(X)[1]
  n_sample <- dim(X)[2]
  p_theta <- dim(X)[3]
  p_alpha <- 2L

  if (is.null(dimnames(X)[[3]]) || !all(c("X1", "X2") %in% dimnames(X)[[3]])) {
    if (p_theta == 4L) {
      dimnames(X)[[3]] <- c("Int", "X1", "X2", "X1X2")
    } else {
      stop("`X` third dimension must be named (need at least X1, X2).", call. = FALSE)
    }
  }

  if (!all(dim(Y) == c(L, n_sample))) {
    stop("`Y` dimensions must match the first two dimensions of `X`.", call. = FALSE)
  }
  if (!all(dim(beta_mat) == c(L, p_theta))) {
    stop("`beta_mat` must be L x p.", call. = FALSE)
  }
  if (!is.array(V_beta_cube) || !all(dim(V_beta_cube) == c(L, p_theta, p_theta))) {
    stop("`V_beta_cube` must be L x p x p.", call. = FALSE)
  }
  if (length(type_vec) != L) stop("`type_vec` must have length L.", call. = FALSE)
  if (!all(type_vec %in% 1:4)) {
    stop("`type_vec` values must be 1, 2, 3 (AD types) or 4 (IPD).", call. = FALSE)
  }

  SEQ_IPD <- which(type_vec == 4L)
  SEQ_AD <- which(type_vec %in% 1:3)
  J <- length(SEQ_IPD)
  K <- length(SEQ_AD)
  if (J < 1L || K < 1L) {
    stop("Need at least one IPD (type 4) and one AD (types 1-3) study.", call. = FALSE)
  }
  ad_type <- type_vec[SEQ_AD]
  idx_t1 <- which(ad_type == 1L)
  idx_t2 <- which(ad_type == 2L)
  idx_t3 <- which(ad_type == 3L)

  if (is.null(theta_init)) theta_init <- matrix(0, nrow = L, ncol = p_theta)
  if (is.null(mu_init)) mu_init <- rep(0, p_theta)
  if (is.null(Sigma_init)) Sigma_init <- diag(1, p_theta)
  if (is.null(sig2_init)) sig2_init <- true_kappa

  theta_init <- as.matrix(theta_init)
  mu_init <- as.numeric(mu_init)
  Sigma_init <- as.matrix(Sigma_init)
  if (!all(dim(theta_init) == c(L, p_theta))) stop("`theta_init` must be L x p.", call. = FALSE)
  if (length(mu_init) != p_theta) stop("`mu_init` must have length p.", call. = FALSE)
  if (!all(dim(Sigma_init) == c(p_theta, p_theta))) stop("`Sigma_init` must be p x p.", call. = FALSE)

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

  X_IPD <- X[SEQ_IPD, , , drop = FALSE]
  Y_mat <- Y[SEQ_IPD, , drop = FALSE]
  beta_tilde_mat <- beta_mat[SEQ_AD, , drop = FALSE]
  V_tilde_cube <- V_beta_cube[SEQ_AD, , , drop = FALSE]
  for (kk in seq_len(K)) {
    V_tilde_cube[kk, , ] <- diag(diag(V_tilde_cube[kk, , ]))
  }

  tilde_D_x <- vector("list", K)
  tilde_D_x_subgroup <- vector("list", K)
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
    tilde_D_x_subgroup[[k]] <- sim1_subgroup_indicators(tilde_D_x[[k]])
  }

  invLambda_theta <- diag(1 / lambda, p_theta)
  Phi0 <- diag(phi0, p_theta)

  theta_mat_IPD <- theta_init[SEQ_IPD, , drop = FALSE]
  mu_vec <- mu_init
  Sigma_theta_mat <- Sigma_init
  sig2 <- as.numeric(sig2_init)[1]
  tau_vec <- hat_Gamma_tau_vec
  alpha_mat <- array(0.1, c(K, p_alpha))
  theta_mat_AD <- theta_init[SEQ_AD, , drop = FALSE]
  beta_mat_AD <- beta_tilde_mat

  for (k in idx_t1) {
    beta_mat_AD[k, 1:3] <- MB.est.ADtype1(
      theta = theta_mat_AD[k, ],
      x = tilde_D_x[[k]],
      w = rep(1, nrow(tilde_D_x[[k]])),
      alpha = alpha_mat[k, ]
    )
  }
  for (k in idx_t2) {
    beta_mat_AD[k, ] <- MB.est.ADtype2(
      theta = theta_mat_AD[k, ],
      x = tilde_D_x[[k]],
      w = rep(1, nrow(tilde_D_x[[k]])),
      subgroup.x = tilde_D_x_subgroup[[k]],
      alpha = alpha_mat[k, ]
    )
  }
  for (k in idx_t3) {
    beta_mat_AD[k, ] <- MB.est.ADtype3(
      theta = theta_mat_AD[k, ],
      x = tilde_D_x[[k]],
      w = rep(1, nrow(tilde_D_x[[k]])),
      alpha = alpha_mat[k, ]
    )
  }

  n_iter <- burnin + mainrun

  if (identical(engine, "cpp")) {
    if (verbose) {
      message(sprintf(
        "Starting Sim1 IPD+AD MCMC (C++): %d iters (%d burn-in + %d main); J=%d IPD, type1=%d, type2=%d, type3=%d",
        n_iter, burnin, mainrun, J, length(idx_t1), length(idx_t2), length(idx_t3)
      ))
    }
    X_ipd_list <- lapply(seq_len(J), function(j) unname(X_IPD[j, , ]))
    X_ref_list <- lapply(tilde_D_x, function(m) unname(as.matrix(m)))
    Z_ref_list <- lapply(tilde_D_x_subgroup, function(m) unname(as.matrix(m)))
    V_list <- lapply(seq_len(K), function(k) unname(V_tilde_cube[k, , ]))
    cpp_out <- sim1_mcmc_cpp(
      X_ipd = X_ipd_list,
      Y_ipd = unname(Y_mat),
      X_ref = X_ref_list,
      Z_ref = Z_ref_list,
      beta_tilde = unname(beta_tilde_mat),
      V_list = V_list,
      idx_t1 = as.integer(idx_t1 - 1L),
      idx_t2 = as.integer(idx_t2 - 1L),
      idx_t3 = as.integer(idx_t3 - 1L),
      theta_ipd0 = unname(theta_mat_IPD),
      theta_ad0 = unname(theta_mat_AD),
      beta_ad0 = unname(beta_mat_AD),
      mu0 = as.numeric(mu_vec),
      Sigma0 = unname(Sigma_theta_mat),
      sig2 = sig2,
      tau0 = as.numeric(tau_vec),
      alpha0 = unname(alpha_mat),
      hat_tau = as.numeric(hat_tau_vec),
      hat_gamma = as.numeric(hat_Gamma_tau_vec),
      invLambda = unname(invLambda_theta),
      Phi0 = unname(Phi0),
      n_iter = as.integer(n_iter),
      burnin = as.integer(burnin),
      step_theta = step_theta,
      step_alpha = step_alpha,
      step_tau = step_tau,
      nu0 = nu0,
      verbose = verbose
    )
    out <- list(
      posterior_mu = cpp_out$posterior_mu,
      posterior_Sigma_diag = cpp_out$posterior_Sigma_diag,
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
        J_type1 = length(idx_t1),
        J_type2 = length(idx_t2),
        J_type3 = length(idx_t3),
        n = n_sample,
        p = p_theta,
        outcome = "gaussian_sim1",
        engine = "cpp",
        used_default_data = using_default_data
      )
    )
    class(out) <- c("bayesmetaipd_fit", "list")
    return(out)
  }
  draw_mu <- array(0, c(n_iter, p_theta))
  draw_Sigma_theta <- array(0, c(n_iter, p_theta))
  draw_sig2 <- rep(0, n_iter)

  if (verbose) {
    message(sprintf(
      "Starting Sim1 IPD+AD MCMC: %d iters (%d burn-in + %d main); J=%d IPD, type1=%d, type2=%d, type3=%d",
      n_iter, burnin, mainrun, J, length(idx_t1), length(idx_t2), length(idx_t3)
    ))
  }
  prev_time <- proc.time()[[3]]

  for (i_iter in seq_len(n_iter)) {
    for (k in idx_t1) {
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_AD[k, ], sd = step_theta)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      beta_vec_q <- tryCatch({
        MB.est.ADtype1(theta = theta_vec_q, x = tilde_D_x[[k]], w = ww, alpha = alpha_mat[k, ])
      }, error = function(e) NULL)
      if (!is.null(beta_vec_q)) {
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3], mean = beta_vec_q[2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3], log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(
          theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3], mean = beta_mat_AD[k, 2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3], log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(
          theta_mat_AD[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        if (stats::runif(n = 1) < exp(logNum - logDen)) {
          theta_mat_AD[k, ] <- theta_vec_q
          beta_mat_AD[k, 1:3] <- beta_vec_q
        }
      }

      alpha_vec_q <- stats::rnorm(n = p_alpha, mean = alpha_mat[k, ], sd = step_alpha)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      beta_vec_q <- tryCatch({
        MB.est.ADtype1(theta = theta_mat_AD[k, ], x = tilde_D_x[[k]], w = ww, alpha = alpha_vec_q)
      }, error = function(e) NULL)
      if (!is.null(beta_vec_q)) {
        Exp_alpha_psi_q <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_vec_q)
        q_l_mat_q <- cbind(
          Exp_alpha_psi_q - 1,
          tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_q - tau_vec[k]
        )
        bar_q_l_q <- apply(q_l_mat_q, 2, mean)
        Sigma_q_l_q <- stats::var(q_l_mat_q) / nrow(q_l_mat_q)
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3], mean = beta_vec_q[2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3], log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(bar_q_l_q, mean = rep(0, 2), sigma = Sigma_q_l_q, log = TRUE)

        Exp_alpha_psi_k <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
        q_l_mat <- cbind(
          Exp_alpha_psi_k - 1,
          tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_k - tau_vec[k]
        )
        bar_q_l <- apply(q_l_mat, 2, mean)
        Sigma_q_l <- stats::var(q_l_mat) / nrow(q_l_mat)
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 2:3], mean = beta_mat_AD[k, 2:3],
          sigma = V_tilde_cube[k, 2:3, 2:3], log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(bar_q_l, mean = rep(0, 2), sigma = Sigma_q_l, log = TRUE)
        if (stats::runif(n = 1) < exp(logNum - logDen)) {
          alpha_mat[k, ] <- alpha_vec_q
          beta_mat_AD[k, 1:3] <- beta_vec_q
        }
      }
    }

    for (k in idx_t2) {
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_AD[k, ], sd = step_theta)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      beta_vec_q <- tryCatch({
        MB.est.ADtype2(
          theta = theta_vec_q, x = tilde_D_x[[k]], w = ww,
          subgroup.x = tilde_D_x_subgroup[[k]], alpha = alpha_mat[k, ]
        )
      }, error = function(e) NULL)
      if (!is.null(beta_vec_q)) {
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, ], mean = as.numeric(beta_vec_q),
          sigma = V_tilde_cube[k, , ], log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(
          theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, ], mean = beta_mat_AD[k, ],
          sigma = V_tilde_cube[k, , ], log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(
          theta_mat_AD[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        if (stats::runif(n = 1) < exp(logNum - logDen)) {
          theta_mat_AD[k, ] <- theta_vec_q
          beta_mat_AD[k, ] <- beta_vec_q
        }
      }

      alpha_vec_q <- stats::rnorm(n = p_alpha, mean = alpha_mat[k, ], sd = step_alpha)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      beta_vec_q <- tryCatch({
        MB.est.ADtype2(
          theta = theta_mat_AD[k, ], x = tilde_D_x[[k]], w = ww,
          subgroup.x = tilde_D_x_subgroup[[k]], alpha = alpha_vec_q
        )
      }, error = function(e) NULL)
      if (!is.null(beta_vec_q)) {
        Exp_alpha_psi_q <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_vec_q)
        q_l_mat_q <- cbind(
          Exp_alpha_psi_q - 1,
          tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_q - tau_vec[k]
        )
        bar_q_l_q <- apply(q_l_mat_q, 2, mean)
        Sigma_q_l_q <- stats::var(q_l_mat_q) / nrow(q_l_mat_q)
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, ], mean = as.numeric(beta_vec_q),
          sigma = V_tilde_cube[k, , ], log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(bar_q_l_q, mean = rep(0, 2), sigma = Sigma_q_l_q, log = TRUE)

        Exp_alpha_psi_k <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
        q_l_mat <- cbind(
          Exp_alpha_psi_k - 1,
          tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_k - tau_vec[k]
        )
        bar_q_l <- apply(q_l_mat, 2, mean)
        Sigma_q_l <- stats::var(q_l_mat) / nrow(q_l_mat)
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, ], mean = beta_mat_AD[k, ],
          sigma = V_tilde_cube[k, , ], log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(bar_q_l, mean = rep(0, 2), sigma = Sigma_q_l, log = TRUE)
        if (stats::runif(n = 1) < exp(logNum - logDen)) {
          alpha_mat[k, ] <- alpha_vec_q
          beta_mat_AD[k, ] <- beta_vec_q
        }
      }
    }

    for (k in idx_t3) {
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_AD[k, ], sd = step_theta)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      beta_vec_q <- tryCatch({
        MB.est.ADtype3(theta = theta_vec_q, x = tilde_D_x[[k]], w = ww, alpha = alpha_mat[k, ])
      }, error = function(e) NULL)
      if (!is.null(beta_vec_q)) {
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 3:4], mean = beta_vec_q[3:4],
          sigma = V_tilde_cube[k, 3:4, 3:4], log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(
          theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 3:4], mean = beta_mat_AD[k, 3:4],
          sigma = V_tilde_cube[k, 3:4, 3:4], log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(
          theta_mat_AD[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        if (stats::runif(n = 1) < exp(logNum - logDen)) {
          theta_mat_AD[k, ] <- theta_vec_q
          beta_mat_AD[k, ] <- beta_vec_q
        }
      }

      alpha_vec_q <- stats::rnorm(n = p_alpha, mean = alpha_mat[k, ], sd = step_alpha)
      ww <- stats::rnorm(nrow(tilde_D_x[[k]]), mean = 1, sd = 1)
      beta_vec_q <- tryCatch({
        MB.est.ADtype3(theta = theta_mat_AD[k, ], x = tilde_D_x[[k]], w = ww, alpha = alpha_vec_q)
      }, error = function(e) NULL)
      if (!is.null(beta_vec_q)) {
        Exp_alpha_psi_q <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_vec_q)
        q_l_mat_q <- cbind(
          Exp_alpha_psi_q - 1,
          tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_q - tau_vec[k]
        )
        bar_q_l_q <- apply(q_l_mat_q, 2, mean)
        Sigma_q_l_q <- stats::var(q_l_mat_q) / nrow(q_l_mat_q)
        logNum <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 3:4], mean = beta_vec_q[3:4],
          sigma = V_tilde_cube[k, 3:4, 3:4], log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(
          theta_mat_AD[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logNum <- logNum + mvtnorm::dmvnorm(bar_q_l_q, mean = rep(0, 2), sigma = Sigma_q_l_q, log = TRUE)

        Exp_alpha_psi_k <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
        q_l_mat <- cbind(
          Exp_alpha_psi_k - 1,
          tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_k - tau_vec[k]
        )
        bar_q_l <- apply(q_l_mat, 2, mean)
        Sigma_q_l <- stats::var(q_l_mat) / nrow(q_l_mat)
        logDen <- mvtnorm::dmvnorm(
          beta_tilde_mat[k, 3:4], mean = beta_mat_AD[k, 3:4],
          sigma = V_tilde_cube[k, 3:4, 3:4], log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(
          theta_mat_AD[k, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE
        )
        logDen <- logDen + mvtnorm::dmvnorm(bar_q_l, mean = rep(0, 2), sigma = Sigma_q_l, log = TRUE)
        if (stats::runif(n = 1) < exp(logNum - logDen)) {
          alpha_mat[k, ] <- alpha_vec_q
          beta_mat_AD[k, ] <- beta_vec_q
        }
      }
    }

    for (j in seq_len(J)) {
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_IPD[j, ], sd = step_theta)
      x.the_q <- X_IPD[j, , 1:4] %*% theta_vec_q
      x.the <- X_IPD[j, , 1:4] %*% theta_mat_IPD[j, ]
      logAcc <- sum(
        stats::dnorm(Y_mat[j, ], x.the_q, sqrt(sig2), log = TRUE) -
          stats::dnorm(Y_mat[j, ], x.the, sqrt(sig2), log = TRUE)
      )
      logAcc <- logAcc +
        mvtnorm::dmvnorm(theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE) -
        mvtnorm::dmvnorm(theta_mat_IPD[j, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)
      if (stats::runif(n = 1) < exp(logAcc)) {
        theta_mat_IPD[j, ] <- theta_vec_q
      }
    }

    theta_merge <- rbind(theta_mat_IPD, theta_mat_AD)

    inv_Sigma <- solve(Sigma_theta_mat)
    inv_Var <- invLambda_theta + (J + K) * inv_Sigma
    Var <- solve(inv_Var)
    Mean_2ndpart <- inv_Sigma %*% apply(theta_merge, 2, sum)
    Mean <- Var %*% Mean_2ndpart
    mu_vec <- mvtnorm::rmvnorm(n = 1, mean = Mean, sigma = Var)

    SS <- array(0, c(p_theta, p_theta))
    for (l in seq_len(J + K)) {
      SS <- SS + t(theta_merge[l, ] - mu_vec) %*% t(t(theta_merge[l, ] - mu_vec))
    }
    Sigma_theta_mat <- MCMCpack::riwish((nu0 + J + K), (Phi0 + SS))

    llik <- rep(0, J)
    for (j in seq_len(J)) {
      x.the <- X_IPD[j, , 1:4] %*% theta_mat_IPD[j, ]
      llik[j] <- sum((Y_mat[j, ] - x.the)^2)
    }
    sig2 <- 1.0 / stats::rgamma(1, shape = (1 + prod(dim(Y_mat)) / 2), rate = (1 + sum(llik) / 2))

    for (k in seq_len(K)) {
      Exp_alpha_psi_k <- exp(tilde_D_x[[k]][, 1:2] %*% alpha_mat[k, ])
      tau_q <- stats::rnorm(n = 1, mean = tau_vec[k], sd = step_tau)
      q_l_mat_q <- cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_k - tau_q)
      bar_q_l_q <- apply(q_l_mat_q, 2, mean)
      Sigma_q_l_q <- stats::var(q_l_mat_q) / nrow(q_l_mat_q)
      logNum <- stats::dnorm(hat_tau_vec[k], mean = tau_q, sd = sqrt(hat_Gamma_tau_vec[k]), log = TRUE)
      logNum <- logNum + mvtnorm::dmvnorm(bar_q_l_q, mean = rep(0, 2), sigma = Sigma_q_l_q, log = TRUE)
      q_l_mat <- cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][, "X1"] * Exp_alpha_psi_k - tau_vec[k])
      bar_q_l <- apply(q_l_mat, 2, mean)
      Sigma_q_l <- stats::var(q_l_mat) / nrow(q_l_mat)
      logDen <- stats::dnorm(hat_tau_vec[k], mean = tau_vec[k], sd = sqrt(hat_Gamma_tau_vec[k]), log = TRUE)
      logDen <- logDen + mvtnorm::dmvnorm(bar_q_l, mean = rep(0, 2), sigma = Sigma_q_l, log = TRUE)
      if (stats::runif(n = 1) < exp(logNum - logDen)) {
        tau_vec[k] <- tau_q
      }
    }

    draw_mu[i_iter, ] <- mu_vec
    draw_Sigma_theta[i_iter, ] <- diag(Sigma_theta_mat)
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
    posterior_Sigma_diag = draw_Sigma_theta[seq_keep, , drop = FALSE],
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
      J_type1 = length(idx_t1),
      J_type2 = length(idx_t2),
      J_type3 = length(idx_t3),
      n = n_sample,
      p = p_theta,
      outcome = "gaussian_sim1",
      engine = "r",
      used_default_data = using_default_data
    )
  )
  class(out) <- c("bayesmetaipd_fit", "list")
  out
}


#' Reproduce Simulation Study 1 IPD+AD (`rep_1`) draws
#'
#' Wrapper around [fit_ipd_ad_sim1()] with the official default settings from
#' `Code/Simulation_1/2_IPD-AD.R`.
#'
#' @param ... Passed to [fit_ipd_ad_sim1()] (e.g. `verbose = FALSE`).
#' @return A `bayesmetaipd_fit` object.
#'
#' @examples
#' \dontrun{
#' fit <- reproduce_sim1_ipdad()
#' colMeans(fit$posterior_mu)
#' }
#'
#' @keywords internal
#' @noRd
reproduce_sim1_ipdad <- function(...) {
  fit_ipd_ad_sim1(...)
}
