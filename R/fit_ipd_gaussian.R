#' Convert IPD input into Application list format
#'
#' @param ipd Either a list of studies (`list(y=, X=)`) or a data frame with a
#'   study id column.
#' @param study_col Study id column name when `ipd` is a data frame.
#' @param y_col Outcome column name when `ipd` is a data frame.
#' @param x_cols Character vector of covariate column names when `ipd` is a
#'   data frame. Default uses Application order.
#' @return List of `list(y, X, n_l)`.
#' @keywords internal
#' @noRd
as_application_ipd <- function(ipd,
                               study_col = "study_name",
                               y_col = "w_gain",
                               x_cols = c(
                                 "bmi_cat1", "bmi_cat2", "bmi_cat3", "b_wt",
                                 "bmi_cat1.trt", "bmi_cat2.trt", "bmi_cat3.trt"
                               )) {
  if (is.data.frame(ipd)) {
    if (!study_col %in% names(ipd)) stop("`study_col` not found in `ipd`.", call. = FALSE)
    if (!y_col %in% names(ipd)) stop("`y_col` not found in `ipd`.", call. = FALSE)
    miss_x <- setdiff(x_cols, names(ipd))
    if (length(miss_x)) {
      stop("Missing covariate columns: ", paste(miss_x, collapse = ", "), call. = FALSE)
    }
    studies <- sort(unique(ipd[[study_col]]))
    out <- vector("list", length(studies))
    for (i in seq_along(studies)) {
      rows <- ipd[ipd[[study_col]] == studies[[i]], , drop = FALSE]
      y <- as.numeric(rows[[y_col]])
      X <- as.matrix(rows[, x_cols, drop = FALSE])
      storage.mode(X) <- "double"
      out[[i]] <- list(y = y, X = X, n_l = length(y), study_name = as.character(studies[[i]]))
    }
    return(out)
  }

  if (!is.list(ipd) || !length(ipd)) {
    stop("`ipd` must be a non-empty list or data frame.", call. = FALSE)
  }
  out <- vector("list", length(ipd))
  for (i in seq_along(ipd)) {
    st <- ipd[[i]]
    if (is.null(st$y) || is.null(st$X)) {
      stop("Each IPD study must provide `y` and `X`.", call. = FALSE)
    }
    y <- as.numeric(st$y)
    X <- as.matrix(st$X)
    storage.mode(X) <- "double"
    if (length(y) != nrow(X)) {
      stop("Study ", i, ": length(y) must equal nrow(X).", call. = FALSE)
    }
    out[[i]] <- list(
      y = y,
      X = X,
      n_l = length(y),
      study_name = if (!is.null(st$study_name)) as.character(st$study_name) else as.character(i)
    )
  }
  out
}


#' Fit Application-style continuous IPD-only meta-analysis
#'
#' Implements the continuous IPD updates from
#' `Code/Application_CodeOnly.R`: Gaussian IPD likelihood, study-specific
#' residual variances with truncated-normal prior, Gibbs updates for study
#' coefficients, and hierarchical random effects for `mu` / `Sigma`.
#'
#' This is the Application model (not Simulation Study 1). Real Application
#' I-WIP data are not redistributed; supply your own IPD.
#'
#' @param ipd IPD studies via [as_application_ipd()] format (list or data frame).
#' @param burnin Burn-in iterations. Default `10000`.
#' @param mainrun Retained iterations. Default `10000`.
#' @param mu_sig2,tau_sig2 Truncated-normal prior mean/sd for each `sigma_l^2`
#'   (support `(0, Inf)`). Defaults `20`.
#' @param mu0 Prior mean for `mu`. Default `0`.
#' @param Sigma0 Prior covariance for `mu`. Default `100 * I`.
#' @param nu0 Inverse-Wishart df. Default `0.1`.
#' @param psi0 Inverse-Wishart scale multiplier. Default `0.1`.
#' @param step_sig2 MH proposal SD for each `sigma_l^2`. Default `2`.
#' @param theta_init Optional `J x p` starting coefficients.
#' @param mu_init Optional length-`p` starting `mu`.
#' @param Sigma_init Optional `p x p` starting covariance.
#' @param sig2_init Optional length-`J` starting residual variances. Default `10`.
#' @param seed RNG seed. Default `117` (official Application script). `NULL`
#'   leaves RNG unchanged.
#' @param verbose Print progress every `progress_every` iterations.
#' @param progress_every Progress print interval. Default `100`.
#'
#' @return A `bayesmetaipd_fit` list with `posterior_mu`,
#'   `posterior_Sigma_diag`, `posterior_sig2` (`mainrun x J`), `call`,
#'   `settings`.
#'
#' @examples
#' toy <- example_application_data()
#' fit <- fit_ipd_gaussian(toy$ipd, burnin = 5, mainrun = 10, verbose = FALSE)
#' colMeans(fit$posterior_mu)
#'
#' @keywords internal
#' @noRd
fit_ipd_gaussian <- function(ipd,
                             burnin = 10000L,
                             mainrun = 10000L,
                             mu_sig2 = 20,
                             tau_sig2 = 20,
                             mu0 = NULL,
                             Sigma0 = NULL,
                             nu0 = 0.1,
                             psi0 = 0.1,
                             step_sig2 = 2,
                             theta_init = NULL,
                             mu_init = NULL,
                             Sigma_init = NULL,
                             sig2_init = NULL,
                             seed = 117L,
                             verbose = TRUE,
                             progress_every = 100L) {
  cl <- match.call()
  out <- fit_ipd_ad_gaussian(
    ipd = ipd,
    ad_type1 = NULL,
    ad_type2 = NULL,
    ad_type3 = NULL,
    burnin = burnin,
    mainrun = mainrun,
    mu_sig2 = mu_sig2,
    tau_sig2 = tau_sig2,
    mu0 = mu0,
    Sigma0 = Sigma0,
    nu0 = nu0,
    psi0 = psi0,
    step_sig2 = step_sig2,
    theta_init = theta_init,
    mu_init = mu_init,
    Sigma_init = Sigma_init,
    sig2_init = sig2_init,
    seed = seed,
    verbose = verbose,
    progress_every = progress_every
  )
  out$call <- cl
  out
}
