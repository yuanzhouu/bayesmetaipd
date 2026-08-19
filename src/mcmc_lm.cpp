// Formula-interface linear IPD+AD MCMC (general X, Z, DRM design)
#include "mcmc_common.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List lm_mcmc_cpp(
    List X_ipd,
    List y_ipd,
    List X_ad,
    List Z_ad,
    List drm_ad,
    List V_ad,
    List beta_tilde,
    List reported,
    NumericVector hat_tau,
    NumericVector hat_gamma,
    NumericMatrix theta_ipd0,
    NumericMatrix theta_ad0,
    List beta_ad0,
    NumericMatrix alpha0,
    NumericVector tau0,
    LogicalVector extra_theta,
    NumericVector mu0,
    NumericMatrix Sigma0,
    double sig2,
    NumericMatrix invLambda,
    NumericMatrix Phi0,
    int n_iter,
    int burnin,
    double step_theta,
    double step_alpha,
    double step_tau,
    double nu0,
    bool use_drm,
    int drm_cov_col,
    bool verbose
) {
  const int J = X_ipd.size();
  const int K = X_ad.size();
  const int p = theta_ipd0.ncol();
  const int Ltot = J + K;
  const uword cov_col = static_cast<uword>(drm_cov_col);

  std::vector<mat> Xj(J), Xa(K), Za(K), Psi(K), V(K);
  std::vector<vec> yj(J), btil(K), beta(K);
  std::vector<uvec> rep(K);
  for (int j = 0; j < J; ++j) {
    Xj[j] = as<mat>(X_ipd[j]);
    yj[j] = as<vec>(y_ipd[j]);
  }
  for (int k = 0; k < K; ++k) {
    Xa[k] = as<mat>(X_ad[k]);
    Za[k] = as<mat>(Z_ad[k]);
    Psi[k] = as<mat>(drm_ad[k]);
    V[k] = as<mat>(V_ad[k]);
    btil[k] = as<vec>(beta_tilde[k]);
    beta[k] = as<vec>(beta_ad0[k]);
    IntegerVector r = reported[k];
    uvec ru(r.size());
    for (int i = 0; i < r.size(); ++i) ru[i] = static_cast<uword>(r[i]);
    rep[k] = ru;
  }

  mat theta_ipd = as<mat>(theta_ipd0);
  mat theta_ad = as<mat>(theta_ad0);
  mat alpha = as<mat>(alpha0);
  vec tau = as<vec>(tau0);
  vec mu = as<vec>(mu0);
  mat Sigma = as<mat>(Sigma0);
  vec htau = as<vec>(hat_tau);
  vec hgam = as<vec>(hat_gamma);
  mat invLam = as<mat>(invLambda);
  mat Phi = as<mat>(Phi0);

  mat draw_mu(n_iter, p, fill::zeros);
  mat draw_Sig(n_iter, p, fill::zeros);
  vec draw_s2(n_iter, fill::zeros);

  GetRNGstate();
  for (int it = 0; it < n_iter; ++it) {
    for (int k = 0; k < K; ++k) {
      const int n_ref = static_cast<int>(Xa[k].n_rows);
      const uword qz = Za[k].n_cols;
      vec th_q = rnorm_vec(vectorise(theta_ad.row(k)), step_theta);
      vec ww(n_ref);
      for (int i = 0; i < n_ref; ++i) ww[i] = R::rnorm(1.0, 1.0);
      vec beta_q;
      bool ok = mb_est_drm(Xa[k], Za[k], Psi[k], ww, th_q, vectorise(alpha.row(k)), beta_q);
      if (ok && beta_q.n_elem == qz) {
        double logNum = dmvnorm_log(btil[k], beta_q.elem(rep[k]), V[k]);
        logNum += dmvnorm_log(th_q, mu, Sigma);
        double logDen = dmvnorm_log(btil[k], beta[k].elem(rep[k]), V[k]);
        logDen += dmvnorm_log(vectorise(theta_ad.row(k)), mu, Sigma);
        if (R::runif(0.0, 1.0) < std::exp(logNum - logDen)) {
          theta_ad.row(k) = th_q.t();
          beta[k] = beta_q;
        }
      }

      if (use_drm) {
        vec al_q = rnorm_vec(vectorise(alpha.row(k)), step_alpha);
        for (int i = 0; i < n_ref; ++i) ww[i] = R::rnorm(1.0, 1.0);
        ok = mb_est_drm(Xa[k], Za[k], Psi[k], ww, vectorise(theta_ad.row(k)), al_q, beta_q);
        if (ok && beta_q.n_elem == qz) {
          vec bar_q, bar_c;
          mat Sq, Sc;
          drm_moments_psi(Psi[k], al_q, tau[k], cov_col, bar_q, Sq);
          drm_moments_psi(Psi[k], vectorise(alpha.row(k)), tau[k], cov_col, bar_c, Sc);
          vec z2(2, fill::zeros);
          double logNum = dmvnorm_log(btil[k], beta_q.elem(rep[k]), V[k]);
          logNum += dmvnorm_log(bar_q, z2, Sq);
          double logDen = dmvnorm_log(btil[k], beta[k].elem(rep[k]), V[k]);
          logDen += dmvnorm_log(bar_c, z2, Sc);
          if (extra_theta[k]) {
            vec th = vectorise(theta_ad.row(k));
            logNum += dmvnorm_log(th, mu, Sigma);
            logDen += dmvnorm_log(th, mu, Sigma);
          }
          if (R::runif(0.0, 1.0) < std::exp(logNum - logDen)) {
            alpha.row(k) = al_q.t();
            beta[k] = beta_q;
          }
        }
      }
    }

    for (int j = 0; j < J; ++j) {
      vec th_q = rnorm_vec(vectorise(theta_ipd.row(j)), step_theta);
      vec th = vectorise(theta_ipd.row(j));
      vec fq = Xj[j] * th_q;
      vec f = Xj[j] * th;
      const double s = std::sqrt(sig2);
      const int n = static_cast<int>(yj[j].n_elem);
      double logAcc = 0.0;
      for (int i = 0; i < n; ++i) {
        logAcc += R::dnorm(yj[j][i], fq[i], s, 1) - R::dnorm(yj[j][i], f[i], s, 1);
      }
      logAcc += dmvnorm_log(th_q, mu, Sigma) - dmvnorm_log(th, mu, Sigma);
      if (R::runif(0.0, 1.0) < std::exp(logAcc)) {
        theta_ipd.row(j) = th_q.t();
      }
    }

    mat theta_merge = join_cols(theta_ipd, theta_ad);
    mat inv_Sigma = inv_sympd(Sigma);
    mat inv_Var = invLam + static_cast<double>(Ltot) * inv_Sigma;
    mat Var = inv_sympd(inv_Var);
    vec Mean = Var * (inv_Sigma * vectorise(sum(theta_merge, 0)));
    mu = rmvnorm_chol(Mean, Var);

    mat SS(p, p, fill::zeros);
    for (int l = 0; l < Ltot; ++l) {
      vec dth = vectorise(theta_merge.row(l)) - mu;
      SS += dth * dth.t();
    }
    Sigma = riwish_bartlett(nu0 + Ltot, Phi + SS);

    double sse = 0.0;
    double ntot = 0.0;
    for (int j = 0; j < J; ++j) {
      vec r = yj[j] - Xj[j] * vectorise(theta_ipd.row(j));
      sse += dot(r, r);
      ntot += static_cast<double>(yj[j].n_elem);
    }
    const double shape = 1.0 + ntot / 2.0;
    const double rate = 1.0 + sse / 2.0;
    sig2 = 1.0 / R::rgamma(shape, 1.0 / rate);

    if (use_drm) {
      for (int k = 0; k < K; ++k) {
        const double tau_q = R::rnorm(tau[k], step_tau);
        vec bar_q, bar_c;
        mat Sq, Sc;
        drm_moments_psi(Psi[k], vectorise(alpha.row(k)), tau_q, cov_col, bar_q, Sq);
        drm_moments_psi(Psi[k], vectorise(alpha.row(k)), tau[k], cov_col, bar_c, Sc);
        vec z2(2, fill::zeros);
        double logNum = R::dnorm(htau[k], tau_q, std::sqrt(hgam[k]), 1);
        logNum += dmvnorm_log(bar_q, z2, Sq);
        double logDen = R::dnorm(htau[k], tau[k], std::sqrt(hgam[k]), 1);
        logDen += dmvnorm_log(bar_c, z2, Sc);
        if (R::runif(0.0, 1.0) < std::exp(logNum - logDen)) tau[k] = tau_q;
      }
    }

    draw_mu.row(it) = mu.t();
    draw_Sig.row(it) = diagvec(Sigma).t();
    draw_s2[it] = sig2;

    if (verbose && ((it + 1) % 1000 == 0)) {
      Rprintf("iter %d / %d\n", it + 1, n_iter);
    }
  }
  PutRNGstate();

  const int keep0 = burnin;
  return List::create(
    Named("posterior_mu") = draw_mu.rows(keep0, n_iter - 1),
    Named("posterior_Sigma_diag") = draw_Sig.rows(keep0, n_iter - 1),
    Named("posterior_sig2") = draw_s2.subvec(keep0, n_iter - 1)
  );
}
