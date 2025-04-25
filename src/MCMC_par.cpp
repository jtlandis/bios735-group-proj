#include <Rcpp.h>
#include "MCMC_utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_fit_par_cpp(const NumericVector& Y, const NumericMatrix& X, 
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
double loglik_par_cpp(const NumericVector& Y, const NumericMatrix& X, 
                      const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  NumericVector mt = get_fit_par_cpp(Y, X, beta, gamma);
  
  double ll = 0.0;
  for (int t = q; t < n; ++t) {
    ll += Y[t] * std::log(mt[t]) - mt[t];
  }
  return ll;
}


// [[Rcpp::export]]
List run_mcmc_par_cpp(const NumericVector& Y,
                      const NumericMatrix& X,
                      int q = 5, 
                      int n_iter = 1000, 
                      Nullable<List> hyperparams = R_NilValue, 
                      double proposal_sd = 0.05,
                      bool verbose = false) {
  // Set hyperparameters
  // int n = Y.size();
  int p = X.ncol();
  
  // Initialize hyperparameters if null
  NumericVector mu(p, 0.0);
  double a_gamma = 1;
  double b_gamma = 1;
  NumericVector alpha(q, 1.0 / q);
  double a_tau = 1.0;
  double b_tau = 1.0;
  
  if (hyperparams.isNotNull()) {
    List h = hyperparams.get();
    if (h.containsElementNamed("mu")) mu = h["mu"];
    if (h.containsElementNamed("a_gamma")) a_gamma = h["a_gamma"];
    if (h.containsElementNamed("b_gamma")) b_gamma = h["b_gamma"];
    if (h.containsElementNamed("alpha")) alpha = h["alpha"];
    if (h.containsElementNamed("a_tau")) a_tau = h["a_tau"];
    if (h.containsElementNamed("b_tau")) b_tau = h["b_tau"];
  }
  
  // Initialize parameters
  NumericVector gamma(p);
  NumericVector sigsq_gamma(p);
  for (int j = 0; j < p; ++j) {
    sigsq_gamma[j] = 1 / R::rgamma(a_gamma, 1/b_gamma);
    gamma[j] = R::rnorm(mu[j], sqrt(sigsq_gamma[j]));
  }
  NumericVector beta_tilde = rdirichlet_cpp(alpha);
  double tau = R::rbeta(a_tau, b_tau);
  NumericVector beta(q);
  for (int i = 0; i < q; ++i) beta[i] = tau * beta_tilde[i];
  NumericMatrix gamma_samples(n_iter, p);
  NumericMatrix sigsq_gamma_samples(n_iter, p);
  NumericMatrix beta_samples(n_iter, q);
  NumericVector tau_samples(n_iter);

  
  // --- MCMC loop ---
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // Propose gamma
    NumericVector gamma_prop(p);
    for (int i = 0; i < p; ++i)
      gamma_prop[i] = R::rnorm(gamma[i], proposal_sd);
    
    double lp_new = loglik_par_cpp(Y, X, beta, gamma_prop) + dnorm_vec_cpp(gamma_prop, mu, sigsq_gamma, true);
    double lp_old = loglik_par_cpp(Y, X, beta, gamma) + dnorm_vec_cpp(gamma, mu, sigsq_gamma, true);
    if (std::log(R::runif(0, 1)) < (lp_new - lp_old)) {
      gamma = gamma_prop;
    }
    
    // Sample variance of gamma (Gibbs)
    for (int j = 0; j < p; ++j) {
      double resid_sq = std::pow(gamma[j] - mu[j], 2.0);
      double shape_post = a_gamma + 0.5;
      double rate_post = b_gamma + 0.5 * resid_sq;
      sigsq_gamma[j] = 1.0 / R::rgamma(shape_post, 1.0 / rate_post);  // Inverse-gamma
    }
    
    // Propose tau (global scale for AR terms) with random walk in inv-logit space
    double eta_curr_tau = std::log(tau / (1.0 - tau)); // current eta = logit(tau)
    double eta_prop_tau = eta_curr_tau + R::rnorm(0, proposal_sd); // symmetric RW proposal
    double tau_prop = 1.0 / (1.0 + std::exp(-eta_prop_tau)); // inverse logit
    
    // double tau_prop = R::rbeta(a_tau, b_tau);
    NumericVector beta_prop(q);
    for (int i = 0; i < q; ++i) beta_prop[i] = tau_prop * beta_tilde[i];
    
    lp_new = loglik_par_cpp(Y, X, beta_prop, gamma) + 
      R::dbeta(tau_prop, a_tau, b_tau, 1) +
      std::log(tau_prop) + std::log(1.0 - tau_prop); // Jacobian
    lp_old = loglik_par_cpp(Y, X, beta, gamma) + 
      R::dbeta(tau, a_tau, b_tau, 1) +
      std::log(tau) + std::log(1.0 - tau); // Jacobian
  
    if (std::log(R::runif(0, 1)) < (lp_new - lp_old)) {
      tau = tau_prop;
      beta = clone(beta_prop);
    }
    
    // Propose beta_tilde via random walk in unconstrained space
    
    // NumericVector beta_tilde_prop = rdirichlet_cpp(alpha);
    NumericVector beta_tilde_prop(q);
    if (q == 1) {
      beta_tilde_prop[0] = 1.0;  // Only one component on the simplex
    } else {
      NumericVector eta_curr_beta = dirichlet_to_eta(beta_tilde);
      NumericVector eta_prop_beta = clone(eta_curr_beta);
      for (int i = 0; i < q - 1; ++i) {
        eta_prop_beta[i] += R::rnorm(0, proposal_sd);
      }
      beta_tilde_prop = eta_to_dirichlet(eta_prop_beta);
    }
    
    for (int i = 0; i < q; ++i) {
      beta_prop[i] = tau * beta_tilde_prop[i];
    }
    
    // compute log-lik + prior + Jacobian (skip Jacobian if q == 1)
    lp_new = loglik_par_cpp(Y, X, beta_prop, gamma) +
      ddirichlet_cpp(beta_tilde_prop, alpha, true) +
      (q > 1 ? log_jacobian_eta_to_dirichlet(dirichlet_to_eta(beta_tilde_prop)) : 0.0);
    
    lp_old = loglik_par_cpp(Y, X, beta, gamma) +
      ddirichlet_cpp(beta_tilde, alpha, true) +
      (q > 1 ? log_jacobian_eta_to_dirichlet(dirichlet_to_eta(beta_tilde)) : 0.0);
    
    // lp_new = loglik_par_cpp(Y, X, beta_prop, gamma) +
    //   ddirichlet_cpp(beta_tilde_prop, alpha, true);
    // lp_old = loglik_par_cpp(Y, X, beta, gamma) +
    //   ddirichlet_cpp(beta_tilde, alpha, true);
    
    if (std::log(R::runif(0, 1)) < (lp_new - lp_old)) {
      beta_tilde = clone(beta_tilde_prop);
      beta = clone(beta_prop);
    }
    
    // Save samples
    for (int j = 0; j < q; ++j) beta_samples(iter, j) = beta[j];
    for (int j = 0; j < p; ++j) gamma_samples(iter, j) = gamma[j];
    for (int j = 0; j < p; ++j) sigsq_gamma_samples(iter, j) = sigsq_gamma[j];
    tau_samples[iter] = tau;

    if (verbose) {
      if ((iter + 1) % 100 == 0) Rcpp::Rcout << "Iteration " << iter + 1 << std::endl;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_samples,
    Rcpp::Named("gamma") = gamma_samples,
    Rcpp::Named("tau") = tau_samples,
    Rcpp::Named("sigsq_gamma") = sigsq_gamma_samples
  );
}