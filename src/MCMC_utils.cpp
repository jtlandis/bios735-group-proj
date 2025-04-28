#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double dnorm_vec_cpp(NumericVector x, NumericVector mu, NumericVector sigsq, bool log=true) {
  int p = x.size();
  double res;
  if (log == true) {
    res = 0;
    for (int i = 0; i < p; ++i) {
      res += R::dnorm4(x[i], mu[i], pow(sigsq[i], 0.5), 1);
    }
  } else {
    res = 1;
    for (int i = 0; i < p; ++i) {
      res *= R::dnorm4(x[i], mu[i], pow(sigsq[i], 0.5), 0);
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericVector rdirichlet_cpp(NumericVector alpha_m) {
  int K = alpha_m.size();
  NumericVector gamma_draws(K);
  double sum = 0.0;
  
  for (int j = 0; j < K; ++j) {
    gamma_draws[j] = R::rgamma(alpha_m[j], 1.0);
    sum += gamma_draws[j];
  }
  
  for (int j = 0; j < K; ++j) {
    gamma_draws[j] /= sum;
  }
  
  return gamma_draws;
}

// [[Rcpp::export]]
double ddirichlet_cpp(NumericVector x, NumericVector alpha, bool log) {
  int q = x.size();
  double alpha0 = 0.0;
  double log_gam_prod = 0.0;
  double log_x_term = 0.0;
  
  for (int i = 0; i < q; ++i) {
    alpha0 += alpha[i];
    log_gam_prod += std::lgamma(alpha[i]);
    if (x[i] <= 0.0) return R_NegInf;
    log_x_term += (alpha[i] - 1.0) * std::log(x[i]);
  }
  
  double log_beta_term = log_gam_prod - std::lgamma(alpha0);
  double log_density = log_x_term - log_beta_term;
  
  if (log) {
    return log_density;
  } else {
    return std::exp(log_density);
  }
}

// [[Rcpp::export]]
NumericVector dirichlet_to_eta(NumericVector beta_tilde) {
  int q = beta_tilde.size();
  NumericVector eta(q - 1);
  for (int i = 0; i < q - 1; ++i) {
    eta[i] = std::log(beta_tilde[i] / beta_tilde[q - 1]);
  }
  return eta;
}

// [[Rcpp::export]]
NumericVector eta_to_dirichlet(NumericVector eta) {
  int q_minus_1 = eta.size();
  int q = q_minus_1 + 1;
  NumericVector beta_tilde(q);
  
  double denom = 1.0;
  for (int i = 0; i < q_minus_1; ++i) {
    denom += std::exp(eta[i]);
  }
  
  for (int i = 0; i < q_minus_1; ++i) {
    beta_tilde[i] = std::exp(eta[i]) / denom;
  }
  beta_tilde[q - 1] = 1.0 / denom;
  
  return beta_tilde;
}

// [[Rcpp::export]]
double log_jacobian_eta_to_dirichlet(NumericVector eta) {
  int q_minus_1 = eta.size();
  double sum_exp = 0.0;
  for (int i = 0; i < q_minus_1; ++i) {
    sum_exp += std::exp(eta[i]);
  }
  
  // log determinant of the Jacobian of inverse ALR transform
  double log_det = 0.0;
  log_det -= (q_minus_1 + 1) * std::log(1.0 + sum_exp);
  for (int i = 0; i < q_minus_1; ++i) {
    log_det += eta[i];  // since log(exp(eta_i)) = eta_i
  }
  
  return log_det;
}