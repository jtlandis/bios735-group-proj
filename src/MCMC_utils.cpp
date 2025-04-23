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
