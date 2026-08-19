// Simulation 1 IPD+AD MCMC inner loop (Gaussian + 3 AD types + DRM)
#include "mcmc_common.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List sim1_mcmc_cpp(
    List X_ipd,
    NumericMatrix Y_ipd,
    List X_ref,
    List Z_ref,
    NumericMatrix beta_tilde,
    List V_list,
    IntegerVector idx_t1,
    IntegerVector idx_t2,
    IntegerVector idx_t3,
    NumericMatrix theta_ipd0,
    NumericMatrix theta_ad0,
    NumericMatrix beta_ad0,
    NumericVector mu0,
    NumericMatrix Sigma0,
    double sig2,
    NumericVector tau0,
    NumericMatrix alpha0,
    NumericVector hat_tau,
    NumericVector hat_gamma,
    NumericMatrix invLambda,
    NumericMatrix Phi0,
    int n_iter,
    int burnin,
    double step_theta,
    double step_alpha,
    double step_tau,
    double nu0,
    bool verbose
) {
  const int J = X_ipd.size();
  const int K = X_ref.size();
  const int p = theta_ipd0.ncol();
  const int n = Y_ipd.ncol();
  const int Ltot = J + K;

  std::vector<mat> Xj(J), Xr(K), Zr(K), V(K);
  for (int j = 0; j < J; ++j) Xj[j] = as<mat>(X_ipd[j]);
  for (int k = 0; k < K; ++k) {
    Xr[k] = as<mat>(X_ref[k]);
    Zr[k] = as<mat>(Z_ref[k]);
    V[k] = as<mat>(V_list[k]);
  }

  mat Y = as<mat>(Y_ipd);
  mat theta_ipd = as<mat>(theta_ipd0);
  mat theta_ad = as<mat>(theta_ad0);
  mat beta_ad = as<mat>(beta_ad0);
  vec mu = as<vec>(mu0);
  mat Sigma = as<mat>(Sigma0);
  vec tau = as<vec>(tau0);
  mat alpha = as<mat>(alpha0);
  mat btil = as<mat>(beta_tilde);
  vec htau = as<vec>(hat_tau);
  vec hgam = as<vec>(hat_gamma);
  mat invLam = as<mat>(invLambda);
  mat Phi = as<mat>(Phi0);

  mat draw_mu(n_iter, p, fill::zeros);
  mat draw_Sig(n_iter, p, fill::zeros);
  vec draw_s2(n_iter, fill::zeros);

  auto do_type = [&](int k, int type) {
    const int n_ref = Xr[k].n_rows;
    mat Zuse;
    uvec rep;
    if (type == 1) {
      Zuse = Xr[k].cols(0, 2);
      uvec r = {1, 2};
      rep = r;
    } else if (type == 2) {
      Zuse = Zr[k];
      rep = linspace<uvec>(0, Zuse.n_cols - 1, Zuse.n_cols);
    } else {
      Zuse = Xr[k];
      uvec r = {2, 3};
      rep = r;
    }

    vec th_q = rnorm_vec(vectorise(theta_ad.row(k)), step_theta);
    vec ww(n_ref);
    for (int i = 0; i < n_ref; ++i) ww[i] = R::rnorm(1.0, 1.0);
    vec beta_q;
    bool ok = mb_est(Xr[k], Zuse, ww, th_q, vectorise(alpha.row(k)), beta_q);
    if (ok) {
      vec bt = btil.row(k).t();
      vec bcur = beta_ad.row(k).t();
      mat Vrep = V[k].submat(rep, rep);
      double logNum = dmvnorm_log(bt.elem(rep), beta_q.elem(rep), Vrep);
      logNum += dmvnorm_log(th_q, mu, Sigma);
      double logDen = dmvnorm_log(bt.elem(rep), bcur.elem(rep), Vrep);
      logDen += dmvnorm_log(vectorise(theta_ad.row(k)), mu, Sigma);
      if (R::runif(0.0, 1.0) < std::exp(logNum - logDen)) {
        theta_ad.row(k) = th_q.t();
        if (type == 1) {
          beta_ad(k, 0) = beta_q[0];
          beta_ad(k, 1) = beta_q[1];
          beta_ad(k, 2) = beta_q[2];
        } else {
          beta_ad.row(k) = beta_q.t();
        }
      }
    }

    vec al_q = rnorm_vec(vectorise(alpha.row(k)), step_alpha);
    for (int i = 0; i < n_ref; ++i) ww[i] = R::rnorm(1.0, 1.0);
    ok = mb_est(Xr[k], Zuse, ww, vectorise(theta_ad.row(k)), al_q, beta_q);
    if (ok) {
      vec bt = btil.row(k).t();
      vec bcur = beta_ad.row(k).t();
      mat Vrep = V[k].submat(rep, rep);
      vec bar_q, bar_c;
      mat Sq, Sc;
      drm_moments(Xr[k], al_q, tau[k], bar_q, Sq);
      drm_moments(Xr[k], vectorise(alpha.row(k)), tau[k], bar_c, Sc);
      vec z2(2, fill::zeros);
      double logNum = dmvnorm_log(bt.elem(rep), beta_q.elem(rep), Vrep);
      logNum += dmvnorm_log(bar_q, z2, Sq);
      double logDen = dmvnorm_log(bt.elem(rep), bcur.elem(rep), Vrep);
      logDen += dmvnorm_log(bar_c, z2, Sc);
      if (type == 3) {
        vec th = vectorise(theta_ad.row(k));
        logNum += dmvnorm_log(th, mu, Sigma);
        logDen += dmvnorm_log(th, mu, Sigma);
      }
      if (R::runif(0.0, 1.0) < std::exp(logNum - logDen)) {
        alpha.row(k) = al_q.t();
        if (type == 1) {
          beta_ad(k, 0) = beta_q[0];
          beta_ad(k, 1) = beta_q[1];
          beta_ad(k, 2) = beta_q[2];
        } else {
          beta_ad.row(k) = beta_q.t();
        }
      }
    }
  };

  GetRNGstate();
  for (int it = 0; it < n_iter; ++it) {
    for (int t = 0; t < idx_t1.size(); ++t) do_type(idx_t1[t], 1);
    for (int t = 0; t < idx_t2.size(); ++t) do_type(idx_t2[t], 2);
    for (int t = 0; t < idx_t3.size(); ++t) do_type(idx_t3[t], 3);

    for (int j = 0; j < J; ++j) {
      vec th_q = rnorm_vec(vectorise(theta_ipd.row(j)), step_theta);
      vec th = vectorise(theta_ipd.row(j));
      vec y = Y.row(j).t();
      vec fq = Xj[j] * th_q;
      vec f = Xj[j] * th;
      const double s = std::sqrt(sig2);
      double logAcc = 0.0;
      for (int i = 0; i < n; ++i) {
        logAcc += R::dnorm(y[i], fq[i], s, 1) - R::dnorm(y[i], f[i], s, 1);
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
    for (int j = 0; j < J; ++j) {
      vec r = Y.row(j).t() - Xj[j] * vectorise(theta_ipd.row(j));
      sse += dot(r, r);
    }
    const double ntot = static_cast<double>(J * n);
    const double shape = 1.0 + ntot / 2.0;
    const double rate = 1.0 + sse / 2.0;
    sig2 = 1.0 / R::rgamma(shape, 1.0 / rate);

    for (int k = 0; k < K; ++k) {
      vec ea = exp(Xr[k].cols(0, 1) * vectorise(alpha.row(k)));
      const double tau_q = R::rnorm(tau[k], step_tau);
      mat Qq(Xr[k].n_rows, 2), Qc(Xr[k].n_rows, 2);
      Qq.col(0) = ea - 1.0;
      Qq.col(1) = Xr[k].col(1) % ea - tau_q;
      Qc.col(0) = ea - 1.0;
      Qc.col(1) = Xr[k].col(1) % ea - tau[k];
      vec bar_q = vectorise(mean(Qq, 0));
      vec bar_c = vectorise(mean(Qc, 0));
      mat Sq = sample_cov(Qq) / static_cast<double>(Xr[k].n_rows);
      mat Sc = sample_cov(Qc) / static_cast<double>(Xr[k].n_rows);
      vec z2(2, fill::zeros);
      double logNum = R::dnorm(htau[k], tau_q, std::sqrt(hgam[k]), 1);
      logNum += dmvnorm_log(bar_q, z2, Sq);
      double logDen = R::dnorm(htau[k], tau[k], std::sqrt(hgam[k]), 1);
      logDen += dmvnorm_log(bar_c, z2, Sc);
      if (R::runif(0.0, 1.0) < std::exp(logNum - logDen)) tau[k] = tau_q;
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
  mat mu_keep = draw_mu.rows(keep0, n_iter - 1);
  mat sig_keep = draw_Sig.rows(keep0, n_iter - 1);
  vec s2_keep = draw_s2.subvec(keep0, n_iter - 1);

  return List::create(
    Named("posterior_mu") = mu_keep,
    Named("posterior_Sigma_diag") = sig_keep,
    Named("posterior_sig2") = s2_keep
  );
}
