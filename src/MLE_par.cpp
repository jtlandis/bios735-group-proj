#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_mt_cpp(const NumericVector& Y, const NumericMatrix& X, 
                         const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  int p = gamma.size();
  NumericVector mt(n, 0.0);
  
  // Compute covariate part: exp(X * gamma)
  NumericVector cov_parts(n);
  for (int i = 0; i < n; ++i) {
    double dot = 0.0;
    for (int j = 0; j < p; ++j) {
      dot += X(i, j) * gamma[j];
    }
    cov_parts[i] = std::exp(dot);
  }
  
  // Compute mt
  double sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);
  for (int t = q; t < n; ++t) {
    double ar_term = 0.0;
    for (int l = 0; l < q; ++l) {
      ar_term += beta[l] * Y[t - l - 1];
    }
    mt[t] = ar_term + (1.0 - sum_beta) * cov_parts[t];
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
  NumericVector cov_parts(n);
  for (int i = 0; i < n; ++i) {
    double dot = 0.0;
    for (int j = 0; j < p; ++j) {
      dot += X(i, j) * gamma[j];
    }
    cov_parts[i] = std::exp(dot);
  }
  
  // Fill in gradients w.r.t. beta
  for (int l = 0; l < q; ++l) {
    for (int t = q; t < n; ++t) {
      grad_mt(t, l) = Y[t - l - 1] - cov_parts[t];
    }
  }
  
  // Fill in gradients w.r.t. gamma
  double sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);
  for (int t = q; t < n; ++t) {
    for (int j = 0; j < p; ++j) {
      grad_mt(t, q + j) = (1.0 - sum_beta) * cov_parts[t] * X(t, j);
    }
  }
  
  return grad_mt;
}


// [[Rcpp::export]]
double loglik_cpp(const NumericVector& Y, const NumericMatrix& X, 
                  const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  NumericVector mt = get_mt_cpp(Y, X, beta, gamma);
  
  double ll = 0.0;
  for (int t = q; t < n; ++t) {
    ll += Y[t] * std::log(mt[t]) - mt[t];
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
  
  for (int t = q; t < n; ++t) {
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


// [[Rcpp::export]]
Rcpp::List proj_grad_descent_cpp(const Rcpp::NumericVector& Y,
                                 const Rcpp::NumericMatrix& X,
                                 Rcpp::NumericVector beta0,
                                 Rcpp::NumericVector gamma0,
                                 double lr = 1e-4,
                                 int maxIter = 100,
                                 double tol = 1e-5,
                                 bool return_allIters = false,
                                 bool verbose = false) {
  
  NumericVector beta = clone(beta0);
  NumericVector gamma = clone(gamma0);
  int q = beta0.size();
  int p = gamma0.size();
  double f = -loglik_cpp(Y, X, beta, gamma);
  double ep = INFINITY;
  int iter = 1;
  
  Rcpp::List beta_iters, gamma_iters;
  if (return_allIters) {
    beta_iters = Rcpp::List(maxIter);
    gamma_iters = Rcpp::List(maxIter);
  }
  
  while (ep > tol && iter <= maxIter) {
    if (verbose) Rcpp::Rcout << "Iteration: " << iter << "\n";
    
    Rcpp::NumericVector grad = -loglik_grad_cpp(Y, X, beta, gamma);
    Rcpp::NumericVector beta_grad = grad[Rcpp::Range(0, q - 1)];
    Rcpp::NumericVector gamma_grad = grad[Rcpp::Range(q, q + p - 1)];
    
    for (int i = 0; i < p; ++i) {
      gamma[i] -= lr * gamma_grad[i];
    }
    
    for (int i = 0; i < q; ++i) {
      beta[i] -= lr * beta_grad[i];
    }
    
    beta = proj_beta_cpp(beta);
    
    double fnew = -loglik_cpp(Y, X, beta, gamma);
    ep = std::abs(fnew - f);
    f = fnew;
    
    if (return_allIters) {
      beta_iters[iter - 1] = Rcpp::clone(beta);
      gamma_iters[iter - 1] = Rcpp::clone(gamma);
    }
    
    if (verbose) Rcpp::Rcout << "ep = " << ep << "\n";
    iter++;
  }
  
  if (verbose) Rcpp::Rcout << "Stopped at iteration " << iter << "\n";
  
  if (return_allIters) {
    return Rcpp::List::create(
      Rcpp::Named("beta") = beta,
      Rcpp::Named("gamma") = gamma,
      Rcpp::Named("beta_iters") = beta_iters,
      Rcpp::Named("gamma_iters") = gamma_iters
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("beta") = beta,
      Rcpp::Named("gamma") = gamma
    );
  }
}