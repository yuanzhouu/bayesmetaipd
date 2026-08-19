#' Fit IPD-only Bayesian random-effects meta-analysis (logistic)
#'
#' Implements the Simulation Study 2 Benchmark sampler: all studies are treated
#' as individual participant data (IPD), with a logistic observation model and
#' hierarchical prior `theta_l ~ N(mu, Sigma)`.
#'
#' **Defaults reproduce the official Benchmark `rep_1` result** when called as
#' `fit_ipd()` (bundled `sim2_rep1` data, `seed = 1001`, burn-in/main = 10000).
#'
#' @param X Optional `L x n x p` array of covariates. If `NULL`, uses
#'   [sim2_rep1].
#' @param Y Optional `L x n` matrix of binary responses. If `NULL`, uses
#'   [sim2_rep1].
#' @param theta_init Optional `L x p` starting values for study-specific
#'   coefficients. Default uses `sim2_rep1$theta_l_mat` when data are default,
#'   otherwise `0`.
#' @param mu_init Optional length-`p` starting value for `mu`. Default uses
#'   `sim2_rep1$true_mu` when data are default, otherwise `0`.
#' @param Sigma_init Optional `p x p` starting covariance. Default `diag(1, p)`.
#' @param burnin Integer burn-in iterations. Default `10000`.
#' @param mainrun Integer post-burn-in iterations retained. Default `10000`.
#' @param step_theta MH proposal SD for each component of `theta_l`. Default `0.15`.
#' @param lambda Prior variance for `mu ~ N(0, lambda * I)`. Default `1e4`.
#' @param nu0 Inverse-Wishart degrees of freedom. Default `0.1`.
#' @param phi0 Inverse-Wishart scale multiplier (`Phi0 = phi0 * I`). Default `0.1`.
#' @param seed Integer RNG seed. Default `1001` (= `rep_no + 1000` with `rep_no = 1`
#'   in the official script). With the default bundled data, `seed = 1001` restores
#'   the `.Random.seed` that the official script obtains after
#'   `load(SimulationData_2.RData)` (that file was written with `save.image()`,
#'   so `load()` overwrites `set.seed(1001)`). Set to `NULL` to leave the RNG
#'   state unchanged. Any other integer calls `set.seed(seed)`.
#' @param verbose Logical; print progress every 1000 iterations.
#'
#' @return A list of class `bayesmetaipd_fit` with components:
#' \describe{
#'   \item{posterior_mu}{`mainrun x p` matrix of posterior draws of `mu`.}
#'   \item{posterior_Sigma_diag}{`mainrun x p` matrix of posterior draws of
#'     `diag(Sigma)`.}
#'   \item{call}{Matched call.}
#'   \item{settings}{List of MCMC / prior settings used.}
#' }
#'
#' @examples
#' \dontrun{
#' # Reproduces official Simulation 2 Benchmark rep_1 (takes a few minutes)
#' fit <- fit_ipd()
#' colMeans(fit$posterior_mu)
#' }
#'
#' # Short run for illustration
#' fit_short <- fit_ipd(burnin = 20, mainrun = 30, verbose = FALSE)
#' str(fit_short$posterior_mu)
#'
#' @export
fit_ipd <- function(X = NULL,
                    Y = NULL,
                    theta_init = NULL,
                    mu_init = NULL,
                    Sigma_init = NULL,
                    burnin = 10000L,
                    mainrun = 10000L,
                    step_theta = 0.15,
                    lambda = 1e4,
                    nu0 = 0.1,
                    phi0 = 0.1,
                    seed = 1001L,
                    verbose = TRUE) {
  using_default_data <- is.null(X) && is.null(Y)
  bundled_random_seed <- NULL

  if (using_default_data) {
    # Prefer namespace lazy-data; fall back to bundled extdata file.
    sim2_rep1 <- NULL
    if (exists("sim2_rep1", envir = asNamespace("bayesmetaipd"), inherits = FALSE)) {
      sim2_rep1 <- get("sim2_rep1", envir = asNamespace("bayesmetaipd"), inherits = FALSE)
    }
    if (is.null(sim2_rep1)) {
      env <- new.env(parent = emptyenv())
      utils::data("sim2_rep1", package = "bayesmetaipd", envir = env)
      if (exists("sim2_rep1", envir = env, inherits = FALSE)) {
        sim2_rep1 <- env$sim2_rep1
      }
    }
    if (is.null(sim2_rep1)) {
      path <- system.file("extdata", "sim2_rep1.rda", package = "bayesmetaipd", mustWork = FALSE)
      if (!nzchar(path) || !file.exists(path)) {
        stop("Could not load bundled dataset sim2_rep1.", call. = FALSE)
      }
      env <- new.env(parent = emptyenv())
      load(path, envir = env)
      sim2_rep1 <- env$sim2_rep1
    }
    X <- sim2_rep1$X_cube
    Y <- sim2_rep1$Y_mat
    if (is.null(theta_init)) theta_init <- sim2_rep1$theta_l_mat
    if (is.null(mu_init)) mu_init <- sim2_rep1$true_mu
    if (!is.null(sim2_rep1$random_seed)) {
      bundled_random_seed <- sim2_rep1$random_seed
    }
  } else {
    if (is.null(X) || is.null(Y)) {
      stop("Provide both `X` and `Y`, or leave both NULL to use sim2_rep1.", call. = FALSE)
    }
  }

  if (!is.array(X) || length(dim(X)) != 3L) {
    stop("`X` must be a 3-D array with dimensions L x n x p.", call. = FALSE)
  }
  if (!is.matrix(Y) && !is.array(Y)) {
    stop("`Y` must be a matrix with dimensions L x n.", call. = FALSE)
  }
  Y <- as.matrix(Y)

  L <- dim(X)[1]
  n_sample <- dim(X)[2]
  p_theta <- dim(X)[3]

  if (!all(dim(Y) == c(L, n_sample))) {
    stop("`Y` dimensions must match the first two dimensions of `X`.", call. = FALSE)
  }

  if (is.null(theta_init)) {
    theta_init <- matrix(0, nrow = L, ncol = p_theta)
  }
  if (is.null(mu_init)) {
    mu_init <- rep(0, p_theta)
  }
  if (is.null(Sigma_init)) {
    Sigma_init <- diag(1, p_theta)
  }

  theta_init <- as.matrix(theta_init)
  mu_init <- as.numeric(mu_init)
  Sigma_init <- as.matrix(Sigma_init)

  if (!all(dim(theta_init) == c(L, p_theta))) {
    stop("`theta_init` must be L x p.", call. = FALSE)
  }
  if (length(mu_init) != p_theta) {
    stop("`mu_init` must have length p.", call. = FALSE)
  }
  if (!all(dim(Sigma_init) == c(p_theta, p_theta))) {
    stop("`Sigma_init` must be p x p.", call. = FALSE)
  }

  burnin <- as.integer(burnin)
  mainrun <- as.integer(mainrun)
  if (burnin < 0L || mainrun < 1L) {
    stop("`burnin` must be >= 0 and `mainrun` must be >= 1.", call. = FALSE)
  }

  if (!is.null(seed)) {
    seed <- as.integer(seed)
    # Official Benchmark.R does set.seed(1001) then load(SimulationData_2.RData).
    # That .RData was written with save.image(), so load() restores .Random.seed
    # and overwrites set.seed(1001). Match that RNG state for default data + seed.
    if (using_default_data && identical(seed, 1001L) && !is.null(bundled_random_seed)) {
      assign(".Random.seed", bundled_random_seed, envir = .GlobalEnv)
    } else {
      set.seed(seed)
    }
  }

  invLambda_theta <- diag(1 / lambda, p_theta)
  Phi0 <- diag(phi0, p_theta)

  # Working copies (match official Benchmark indexing: all studies as IPD)
  J <- L
  X_IPD <- X
  Y_mat <- Y
  theta_mat_IPD <- theta_init
  # Keep mu as length-p numeric initially; after the first rmvnorm draw it
  # becomes a 1 x p matrix, matching the official Benchmark script.
  mu_vec <- mu_init
  Sigma_theta_mat <- Sigma_init

  n_iter <- burnin + mainrun
  draw_mu <- array(0, c(n_iter, p_theta))
  draw_Sigma_theta <- array(0, c(n_iter, p_theta))

  if (verbose) {
    message(sprintf("Starting MCMC: %d iterations (%d burn-in + %d main)",
                    n_iter, burnin, mainrun))
  }
  prev_time <- proc.time()[[3]]

  for (i_iter in seq_len(n_iter)) {
    # Step 1: MH update of each theta_l
    for (j in seq_len(J)) {
      theta_vec_q <- stats::rnorm(n = p_theta, mean = theta_mat_IPD[j, ], sd = step_theta)
      # Use columns 1:p explicitly (official script uses 1:4)
      x_the_q <- X_IPD[j, , seq_len(p_theta)] %*% theta_vec_q
      x_the <- X_IPD[j, , seq_len(p_theta)] %*% theta_mat_IPD[j, ]
      logAcc <- sum(
        x_the_q * Y_mat[j, ] - log(1 + exp(x_the_q)) -
          x_the * Y_mat[j, ] + log(1 + exp(x_the))
      )
      logAcc <- logAcc +
        mvtnorm::dmvnorm(theta_vec_q, mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE) -
        mvtnorm::dmvnorm(theta_mat_IPD[j, ], mean = mu_vec, sigma = Sigma_theta_mat, log = TRUE)

      if (stats::runif(n = 1) < exp(logAcc)) {
        theta_mat_IPD[j, ] <- theta_vec_q
      }
    }

    # Step 2: Gibbs update of mu
    inv_Sigma <- solve(Sigma_theta_mat)
    inv_Var <- invLambda_theta + J * inv_Sigma
    Var <- solve(inv_Var)
    Mean_2ndpart <- inv_Sigma %*% apply(theta_mat_IPD, 2, sum)
    Mean <- Var %*% Mean_2ndpart
    # Do not coerce to numeric: official code keeps the 1 x p matrix from rmvnorm
    mu_vec <- mvtnorm::rmvnorm(n = 1, mean = Mean, sigma = Var)

    # Step 3: Gibbs update of Sigma
    # Exact arithmetic form from the official Benchmark script (depends on
    # mu_vec being a 1 x p matrix after the first update).
    SS <- array(0, c(p_theta, p_theta))
    for (l in seq_len(J)) {
      SS <- SS + t(theta_mat_IPD[l, ] - mu_vec) %*% t(t(theta_mat_IPD[l, ] - mu_vec))
    }
    Sigma_theta_mat <- MCMCpack::riwish((nu0 + J), (Phi0 + SS))

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
      step_theta = step_theta,
      lambda = lambda,
      nu0 = nu0,
      phi0 = phi0,
      seed = seed,
      L = L,
      n = n_sample,
      p = p_theta,
      used_default_data = using_default_data
    )
  )
  class(out) <- c("bayesmetaipd_fit", "list")
  out
}


#' Reproduce official Simulation 2 Benchmark (`rep_1`) draws
#'
#' Convenience wrapper around [fit_ipd()] with the exact default settings used
#' in `Code/Simulation_2/1_Benchmark.R` of the Bayesian-Meta repository.
#'
#' @param ... Passed to [fit_ipd()] (e.g. `verbose = FALSE`).
#' @param compare_official If `TRUE`, also load bundled official draws and
#'   report `all.equal` comparisons.
#'
#' @return A `bayesmetaipd_fit` fit object. If `compare_official = TRUE`, an
#'   attribute `comparison` is attached with equality checks.
#'
#' @examples
#' \dontrun{
#' fit <- reproduce_sim2_benchmark(compare_official = TRUE)
#' attr(fit, "comparison")
#' }
#'
#' @export
reproduce_sim2_benchmark <- function(..., compare_official = FALSE) {
  fit <- fit_ipd(...)

  if (isTRUE(compare_official)) {
    official_path <- system.file(
      "extdata", "official_sim2_benchmark_rep1.RData",
      package = "bayesmetaipd"
    )
    if (!nzchar(official_path) || !file.exists(official_path)) {
      warning("Official reference file not found in the installed package.", call. = FALSE)
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
        "Comparison vs official rep_1: ",
        "mu all.equal=", comparison$mu_all_equal,
        ", Sigma all.equal=", comparison$sigma_all_equal,
        ", mu max|diff|=", format(comparison$mu_max_abs_diff, digits = 3),
        ", Sigma max|diff|=", format(comparison$sigma_max_abs_diff, digits = 3)
      )
    }
  }

  fit
}


#' @export
print.bayesmetaipd_fit <- function(x, ...) {
  outcome <- x$settings$outcome
  if (is.null(outcome)) outcome <- "logistic"
  has_ad <- !is.null(x$settings$J) && !is.null(x$settings$K) &&
    isTRUE(x$settings$K > 0)
  if (identical(outcome, "gaussian_application")) {
    lab <- if (has_ad) {
      "Application continuous IPD+AD meta-analysis"
    } else {
      "Application continuous IPD meta-analysis"
    }
    cat(lab, "\n", sep = "")
  } else if (identical(outcome, "gaussian_sim1")) {
    cat("Simulation 1 continuous IPD+AD meta-analysis\n")
  } else if (identical(outcome, "gaussian_lm")) {
    cat("Linear IPD+AD meta-analysis (formula interface)\n")
  } else if (has_ad) {
    cat(sprintf("Bayesian IPD+AD random-effects meta-analysis (%s)\n", outcome))
  } else {
    cat(sprintf("Bayesian IPD random-effects meta-analysis (%s)\n", outcome))
  }
  cat(sprintf(
    "  Draws: %d (after burn-in %d)\n",
    x$settings$mainrun, x$settings$burnin
  ))
  if (identical(outcome, "gaussian_application")) {
    jt1 <- if (is.null(x$settings$J_type1)) 0L else x$settings$J_type1
    jt2 <- if (is.null(x$settings$J_type2)) 0L else x$settings$J_type2
    jt3 <- if (is.null(x$settings$J_type3)) 0L else x$settings$J_type3
    cat(sprintf(
      "  Studies L=%d (IPD=%d, type1=%d, type2=%d, type3=%d), p=%d\n",
      x$settings$L, x$settings$J, jt1, jt2, jt3, x$settings$p
    ))
  } else if (identical(outcome, "gaussian_lm") || identical(outcome, "gaussian_sim1")) {
    jt1 <- if (is.null(x$settings$J_type1)) 0L else x$settings$J_type1
    jt2 <- if (is.null(x$settings$J_type2)) 0L else x$settings$J_type2
    jt3 <- if (is.null(x$settings$J_type3)) 0L else x$settings$J_type3
    cat(sprintf(
      "  Studies L=%d (IPD=%d, nested=%d, subgroup=%d, partial=%d), p=%d\n",
      x$settings$L, x$settings$J, jt1, jt2, jt3, x$settings$p
    ))
  } else if (has_ad) {
    cat(sprintf(
      "  Studies L=%d (J=%d IPD, K=%d AD), n=%d, p=%d\n",
      x$settings$L, x$settings$J, x$settings$K, x$settings$n, x$settings$p
    ))
  } else {
    cat(sprintf("  Studies L=%d, n=%d, p=%d\n", x$settings$L, x$settings$n, x$settings$p))
  }
  cat("  Posterior mean of mu:\n")
  print(colMeans(x$posterior_mu))
  if (!is.null(x$posterior_sig2)) {
    cat(sprintf("  Posterior mean of sig2: %s\n", format(mean(x$posterior_sig2), digits = 4)))
  }
  invisible(x)
}
