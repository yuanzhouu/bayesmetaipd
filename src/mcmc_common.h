#ifndef BAYESMETAIPD_MCMC_COMMON_H
#define BAYESMETAIPD_MCMC_COMMON_H

#include <RcppArmadillo.h>

inline double dmvnorm_log(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma) {
  const arma::uword k = x.n_elem;
  arma::vec diff = x - mean;
  arma::mat L;
  if (!arma::chol(L, sigma, "lower")) return R_NegInf;
  arma::vec u = arma::solve(arma::trimatl(L), diff);
  const double quad = arma::dot(u, u);
  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  return -0.5 * (static_cast<double>(k) * std::log(2.0 * M_PI) + logdet + quad);
}

inline bool mb_est_drm(const arma::mat& X, const arma::mat& Z, const arma::mat& drm_psi,
                       const arma::vec& w, const arma::vec& theta, const arma::vec& alpha,
                       arma::vec& beta) {
  const arma::vec drm = arma::exp(drm_psi * alpha);
  const arma::vec ww = w % drm;
  const arma::vec fitted = X * theta;
  arma::mat ZtWZ = Z.t() * (Z.each_col() % ww);
  arma::vec ZtWy = Z.t() * (ww % fitted);
  return arma::solve(beta, ZtWZ, ZtWy, arma::solve_opts::no_approx);
}

// Sim1: DRM design is intercept + X1 = first two columns of X
inline bool mb_est(const arma::mat& X, const arma::mat& Z, const arma::vec& w,
                   const arma::vec& theta, const arma::vec& alpha, arma::vec& beta) {
  return mb_est_drm(X, Z, X.cols(0, 1), w, theta, alpha, beta);
}

inline arma::mat sample_cov(const arma::mat& Q) {
  const arma::uword n = Q.n_rows;
  arma::rowvec mu = arma::mean(Q, 0);
  arma::mat C = Q.each_row() - mu;
  return (C.t() * C) / static_cast<double>(n - 1);
}

inline arma::vec rnorm_vec(const arma::vec& mean, double sd) {
  arma::vec out(mean.n_elem);
  for (arma::uword i = 0; i < mean.n_elem; ++i) out[i] = R::rnorm(mean[i], sd);
  return out;
}

inline arma::mat riwish_bartlett(double v, const arma::mat& S) {
  const arma::uword p = S.n_rows;
  arma::mat Sinv = arma::inv_sympd(S);
  arma::mat CC = arma::chol(Sinv);
  arma::mat Z(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    Z(i, i) = std::sqrt(R::rchisq(v - static_cast<double>(i)));
  }
  if (p > 1) {
    for (arma::uword col = 1; col < p; ++col) {
      for (arma::uword row = 0; row < col; ++row) {
        Z(row, col) = R::rnorm(0.0, 1.0);
      }
    }
  }
  arma::mat W = arma::trans(Z * CC) * (Z * CC);
  return arma::inv_sympd(W);
}

inline arma::vec rmvnorm_chol(const arma::vec& mean, const arma::mat& sigma) {
  arma::mat U = arma::chol(sigma);
  arma::vec z(mean.n_elem);
  for (arma::uword i = 0; i < mean.n_elem; ++i) z[i] = R::rnorm(0.0, 1.0);
  return mean + U.t() * z;
}

inline void drm_moments_psi(const arma::mat& psi, const arma::vec& alpha, double tau,
                            arma::uword cov_col, arma::vec& bar, arma::mat& Sigma_bar) {
  const arma::vec ea = arma::exp(psi * alpha);
  arma::mat Q(psi.n_rows, 2);
  Q.col(0) = ea - 1.0;
  Q.col(1) = psi.col(cov_col) % ea - tau;
  bar = arma::vectorise(arma::mean(Q, 0));
  Sigma_bar = sample_cov(Q) / static_cast<double>(psi.n_rows);
}

inline void drm_moments(const arma::mat& X, const arma::vec& alpha, double tau,
                        arma::vec& bar, arma::mat& Sigma_bar) {
  drm_moments_psi(X.cols(0, 1), alpha, tau, 1, bar, Sigma_bar);
}

#endif
