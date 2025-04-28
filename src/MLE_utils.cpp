
#include <RcppEigen.h>
#include <Rmath.h>
using namespace Rcpp;

NumericVector get_cov_parts(
  const NumericMatrix& X,
  const NumericVector& gamma,
  // const List& list,
  int q
) {

  Rcpp::List list;
  int n = X.nrow();
  int p = gamma.size();

  // Compute covariate part: exp(X * gamma)
  NumericVector cov_parts(n);
  RObject attr = X.attr(".non_empty");
  if (!attr.isNULL()) {
    List list(attr);
    int t = 0;
    for (int j = 0; j < p; ++j) {
      NumericVector t_index = list[j + q];
      n = t_index.size();
      for (int i = 0; i < n; ++i) {
        t = t_index[i];
        cov_parts[t] += X(t, j + q) * gamma[j];
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double dot = 0.0;
      for (int j = 0; j < p; ++j) {
        dot += X(i, j + q) * gamma[j];
      }
      cov_parts[i] = dot;
    }
  }

  cov_parts = exp(cov_parts);

  return cov_parts;

}



// [[Rcpp::export]]
NumericVector get_mt_cpp(const NumericVector& Y, const NumericMatrix& X,
                         const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  // int p = gamma.size();
  // int t = 0;
  NumericVector mt(n, 0.0);

  // Compute covariate part: exp(X * gamma)
  NumericVector cov_parts = get_cov_parts(X, gamma, q);

  // Compute mt
  double sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);
  for (int t = 0; t < n; ++t) {
    double ar_term = 0.0;
    for (int l = 0; l < q; ++l) {
      ar_term += beta[l] * X(t, l);
    }
    double mt_t = ar_term + (1.0 - sum_beta) * cov_parts[t];
    if (mt_t < 0) {
      mt[t] = 0;
    } else {
      mt[t] = mt_t;
    }
  }
  return mt;
}

// [[Rcpp::export]]
NumericMatrix get_mt_grad_cpp(const NumericVector& Y, const NumericMatrix& X,
                               const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  int p = gamma.size();

  NumericMatrix grad_mt(n, q + p);

  // Compute covariate part: exp(X * gamma)
  NumericVector cov_parts = get_cov_parts(X, gamma, q);

  // Fill in gradients w.r.t. beta
  for (int l = 0; l < q; ++l) {
    for (int t = 0; t < n; ++t) {
      grad_mt(t, l) = X(t, l) - cov_parts[t];
    }
  }

  // Fill in gradients w.r.t. gamma
  double sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);
  sum_beta = 1 - sum_beta;
  RObject attr = X.attr(".non_empty");
  if (attr.isNULL()) {
    for (int t = 0; t < n; ++t) {
      for (int j = 0; j < p; ++j) {
        grad_mt(t, q + j) = sum_beta * cov_parts[t] * X(t, q + j);
      }
    }
  } else {
    List list(attr);
    int t = 0;
    for (int j = 0; j < p; ++j) {
      NumericVector t_index = list[j + q];
      n = t_index.size();
      for (int i = 0; i < n; ++i) {
        t = t_index[i];
        grad_mt(t, q + j) = sum_beta * cov_parts[t] * X(t, j + q);
      }
    }
  }

  return grad_mt;
}


// [[Rcpp::export]]
double loglik_cpp(const NumericVector& Y, const NumericMatrix& X,
                  const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  //int q = beta.size();
  NumericVector mt = get_mt_cpp(Y, X, beta, gamma);

  double ll = 0.0;
  for (int t = 0; t < n; ++t) {
    ll += Y[t] * std::log(mt[t]) - mt[t] - R::lgammafn(Y[t] + 1) ;
  }
  return ll;
}


// [[Rcpp::export]]
NumericVector loglik_grad_cpp(const NumericVector& Y, const NumericMatrix& X,
                               const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  int p = gamma.size();
  int total_params = q + p;

  NumericVector mt = get_mt_cpp(Y, X, beta, gamma);
  NumericMatrix mt_grad = get_mt_grad_cpp(Y, X, beta, gamma);

  NumericVector grad(total_params, 0.0);

  for (int t = 0; t < n; ++t) {
    double weight = Y[t] / mt[t] - 1.0;
    for (int j = 0; j < total_params; ++j) {
      grad[j] += weight * mt_grad(t, j);
    }
  }

  return grad;
}


// [[Rcpp::export]]
Rcpp::NumericVector proj_beta_cpp(Rcpp::NumericVector beta, double epsilon = 1e-4) {
  int q = beta.size();
  Rcpp::NumericVector beta_proj(q);

  // Project to non-negative orthant
  for (int i = 0; i < q; ++i) {
    beta_proj[i] = std::max(beta[i], 0.0);
  }

  double sum_proj = Rcpp::sum(beta_proj);

  if (sum_proj > 1.0) {
    // Sort in decreasing order
    Rcpp::NumericVector u = Rcpp::clone(beta_proj).sort(true);
    Rcpp::NumericVector cssv(q);
    cssv[0] = u[0];
    for (int i = 1; i < q; ++i) {
      cssv[i] = cssv[i - 1] + u[i];
    }

    int rho = 0;
    for (int i = 0; i < q; ++i) {
      if (u[i] > (cssv[i] - 1.0) / (i + 1)) {
        rho = i;
      }
    }
    double theta = (cssv[rho] - 1.0) / (rho + 1);

    for (int i = 0; i < q; ++i) {
      beta_proj[i] = (1 - epsilon) * std::max(beta_proj[i] - theta, 0.0);
    }
  }

  return beta_proj;
}
