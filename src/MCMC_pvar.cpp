#include <Rcpp.h>
#include "MCMC_utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_fit_item_cpp(const NumericVector& Y, const NumericVector& X, 
                               const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  //int p = gamma.size();
  NumericVector mt(n, 0.0);
  
  // Compute covariate part: exp(X * gamma)
  NumericVector cov_parts(n);
  for (int i = 0; i < n; ++i) {
    double dot = gamma[0] + X(i) * gamma[1];
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
NumericMatrix get_fit_pvar_cpp(const NumericMatrix& Y, 
                               const NumericMatrix& X,
                               const NumericMatrix& Beta,
                               const NumericMatrix& Gamma) {
  int T = Y.nrow();
  int N = Y.ncol();
  int q = Beta.ncol();
  //int p = Gamma.ncol();
  
  NumericMatrix mt(T,N);
  
  for (int i = 0; i < N; ++i) {
    double beta_sum = 0.0;
    for (int l = 0; l < q; ++l) {
      beta_sum += Beta(i, l);
    }
    
    for (int t = q; t < T; ++t) {
      // AR component
      double ar_part = 0.0;
      for (int l = 0; l < q; ++l) {
        ar_part += Beta(i, l) * Y(t - l - 1, i);
      }
      
      // Covariate component
      double cov_dot = Gamma(i, 0) + Gamma(i, 1) * X(t, i);
      double cov_part = std::exp(cov_dot);
      
      mt(t, i) = ar_part + (1.0 - beta_sum) * cov_part;
    }
  }
  return mt;
}

// [[Rcpp::export]]
double loglik_item_cpp(const NumericVector& Y, const NumericVector& X, 
                       const NumericVector& beta, const NumericVector& gamma) {
  int n = Y.size();
  int q = beta.size();
  NumericVector mt = get_fit_item_cpp(Y, X, beta, gamma);
  
  double ll = 0.0;
  for (int t = q; t < n; ++t) {
    ll += Y[t] * std::log(mt[t]) - mt[t];
  }
  return ll;
}


// [[Rcpp::export]]
double loglik_pvar_cpp(const NumericMatrix& Y, const NumericMatrix& X, 
                      const NumericMatrix& Beta, const NumericMatrix& Gamma) {
  int T = Y.nrow();
  int N = Y.ncol();
  int q = Beta.ncol();
  
  NumericMatrix mt = get_fit_pvar_cpp(Y, X, Beta, Gamma);
  
  double ll = 0.0;
  for (int i = 0; i < N; ++i) {
    for (int t = q; t < T; ++t) {
      if (mt(t, i) > 0) {
        ll += Y(t, i) * std::log(mt(t, i)) - mt(t, i);
      }
    }
  }
  return ll;
}


// [[Rcpp::export]]
List run_mcmc_pvar_cpp(const NumericMatrix& Y,
                       const NumericMatrix& X,
                       const IntegerVector& G,
                       int B, // number of brands
                       int q = 5, 
                       int n_iter = 1000, 
                       Nullable<List> hyperparams = R_NilValue, 
                       double proposal_sd = 0.05,
                       bool verbose = false) {
  // Set hyperparameters
  int T = Y.nrow();
  int n_items = Y.ncol();
  int p = 2;
  
  // Initialize hyperparameters if null
  double mu0 = 0;
  double sigsq0 = 1;
  double a_gamma = 1.0;
  double b_gamma = 1.0;
  NumericVector alpha0(q, 1.0 / q);
  double a_tau = 1.0;
  double b_tau = 1.0;
  
  if (hyperparams.isNotNull()) {
    List h = hyperparams.get();
    if (h.containsElementNamed("mu0")) mu0 = h["mu0"];
    if (h.containsElementNamed("sigsq0")) sigsq0 = h["sigsq0"];
    if (h.containsElementNamed("a_gamma")) a_gamma = h["a_gamma"];
    if (h.containsElementNamed("b_gamma")) b_gamma = h["b_gamma"];
    if (h.containsElementNamed("alpha0")) alpha0 = h["alpha0"];
    if (h.containsElementNamed("a_tau")) a_tau = h["a_tau"];
    if (h.containsElementNamed("b_tau")) b_tau = h["b_tau"];
  }
  
  // Initialize parameters: item-level and brand-level
  NumericMatrix Gamma_item(n_items, p);
  NumericMatrix mu_gamma_brand(B, p);
  NumericMatrix sigsq_gamma_brand(B, p);
  
  NumericMatrix Beta_item(n_items, q);
  NumericMatrix alpha_beta_brand(B, q);
  NumericVector tau_brand(B);
  
  // Brand-Level initialization
  for (int b = 0; b < B; ++b) {
    for (int j = 0; j < p; ++j) {
      mu_gamma_brand(b, j) = R::rnorm(mu0, pow(sigsq0, 0.5));
      sigsq_gamma_brand(b, j) = 1 / R::rgamma(a_gamma, 1/b_gamma);
    }
    NumericVector alpha_beta_b = rdirichlet_cpp(alpha0);
    for (int l = 0; l < q; ++l) {
      alpha_beta_brand(b, l) = alpha_beta_b[l];
    }
    tau_brand[b] = R::rbeta(a_tau, b_tau);
  }
  
  // Item-Level initialization
  for (int i = 0; i < n_items; ++i) {
    int b = G[i] - 1;
    
    // gamma_i ~ N(mu_b, diag(sigsq_b))
    for (int j = 0; j < p; ++j) {
      Gamma_item(i, j) = R::rnorm(mu_gamma_brand(b, j), pow(sigsq_gamma_brand(b,j), 0.5));
    }
    
    // beta_tilde_i ~ Dirichlet(alpha_b)
    NumericVector alpha_beta_b = alpha_beta_brand(b, _);
    NumericVector beta_tilde = rdirichlet_cpp(alpha_beta_b);
    for (int l = 0; l < q; ++l) {
      Beta_item(i, l) = tau_brand[b] * beta_tilde[l];
    }
  }
  
  List Gamma_samples(n_iter); // item-level covariate effect
  List mu_samples(n_iter); // brand-level means of gamma
  List sigsq_samples(n_iter); // brand-level variances of gamma
  List Beta_samples(n_iter); // item-level AR effects
  List alpha_samples(n_iter); // brand-level Dirichlet params for AR effects
  List tau_samples(n_iter); // brand-level scaling factors for AR effects
  
  // --- MCMC loop ---
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // --- ITEM-LEVEL SAMPLING ---
    for (int i = 0; i < n_items; ++i) {
      
      // Data for i-th item
      NumericVector Y_i = Y(_, i);
      NumericVector X_i = X(_, i);
      int b = G[i] - 1;
      
      // Current gamma and beta for item i
      NumericVector gamma_i = Gamma_item(i, _);
      NumericVector beta_i = Beta_item(i, _);
      NumericVector beta_tilde_i(q);
      for (int l = 0; l < q; ++l) {
        beta_tilde_i[l] = beta_i[l] / tau_brand[b];
      }
      
      // Propose new gamma_i
      NumericVector mu_b = mu_gamma_brand(b, _);
      NumericVector sigsq_b = sigsq_gamma_brand(b, _);
      NumericVector gamma_prop_i(p);
      for (int j = 0; j < p; ++j) {
        gamma_prop_i[j] = R::rnorm(gamma_i[j], proposal_sd);
      }
      
      // Compute posterior likelihood ratio
      double lp_new_gam = loglik_item_cpp(Y_i, X_i, beta_i, gamma_prop_i) + 
        dnorm_vec_cpp(gamma_prop_i, mu_b, sigsq_b, true);
      double lp_old_gam = loglik_item_cpp(Y_i, X_i, beta_i, gamma_i) +
        dnorm_vec_cpp(gamma_i, mu_b, sigsq_b, true);
      
      // Accept/reject proposed gamma_i
      if (std::log(R::runif(0, 1)) < lp_new_gam - lp_old_gam) {
        for (int j = 0; j < p; ++j) {
          Gamma_item(i, j) = gamma_prop_i[j];
        }
        gamma_i = clone(gamma_prop_i); // updated gamma_i for beta update
      }
      
      // Propose beta_tilde
      NumericVector alpha_b = alpha_beta_brand(b, _);
      double tau_b = tau_brand[b];
      NumericVector beta_tilde_prop_i = rdirichlet_cpp(alpha_b);
      NumericVector beta_prop_i(q);
      for (int l = 0; l < q; ++l) {
        beta_prop_i[l] = tau_b * beta_tilde_prop_i[l];
      }
      
      // Compute posterior likelihood ratio
      double lp_new_beta = loglik_item_cpp(Y_i, X_i, beta_prop_i, gamma_i) +
        ddirichlet_cpp(beta_tilde_prop_i, alpha_b, true);
      double lp_old_beta = loglik_item_cpp(Y_i, X_i, beta_i, gamma_i) +
        ddirichlet_cpp(beta_tilde_i, alpha_b, true);
      
      // Accept/rejected proposed beta_i
      if (std::log(R::runif(0, 1)) < (lp_new_beta - lp_old_beta)) {
        for (int l = 0; l < q; ++l) {
          Beta_item(i, l) = beta_prop_i[l];
        }
        // beta_tilde_i = clone(beta_tilde_prop_i);
        // beta_i = clone(beta_prop_i);
      }
    }
    
    // --- BRAND-LEVEL SAMPLING ---
    for (int b = 0; b < B; ++b) {
      
      // Identify indices for items in brand b
      std::vector<int> brand_indices;
      for (int i = 0; i < n_items; ++i) {
        if (G[i] - 1 == b) {
          brand_indices.push_back(i);
        }
      }
      int n_b = brand_indices.size();
      
      // Data and item-level params for b-th brand
      NumericMatrix Y_b(T, n_b);
      NumericMatrix X_b(T, n_b);
      NumericMatrix Beta_b(n_b, q);
      NumericMatrix Gamma_b(n_b, p);
      for (int j = 0; j < n_b; ++j) {
        int item_idx = brand_indices[j];
        Y_b(_, j) = Y(_, item_idx);
        X_b(_, j) = X(_, item_idx);
        Beta_b(j, _) = Beta_item(item_idx, _);
        Gamma_b(j, _) = Gamma_item(item_idx, _);
      }
      
      // Sample mu_gamma_brand (Gibbs)
      for (int j = 0; j < p; ++j) {
        double sum_gamma = 0.0;
        for (int k = 0; k < n_b; ++k) {
          int i = brand_indices[k];
          sum_gamma += Gamma_item(i, j);
        }
        
        double sigsq_bj = sigsq_gamma_brand(b, j);
        double post_var = 1.0 / (n_b / sigsq_bj + 1.0 / sigsq0);
        double post_mean = post_var * (sum_gamma / sigsq_bj);
        mu_gamma_brand(b, j) = R::rnorm(post_mean, pow(post_var, 0.5));
      }
      
      // Sample sigsq_gamma_brand (Gibbs)
      for (int j = 0; j < p; ++j) {
        double sum_sq = 0.0;
        double mu_bj = mu_gamma_brand(b, j);
        for (int k = 0; k < n_b; ++k) {
          int i = brand_indices[k];
          double diff = Gamma_item(i, j) - mu_bj;
          sum_sq += diff * diff;
        }
        
        double shape_post = a_gamma + n_b / 2.0;
        double rate_post = b_gamma + 0.5 * sum_sq;
        sigsq_gamma_brand(b, j) = 1.0 / R::rgamma(shape_post, 1.0 / rate_post);
      }
      
      // Sample tau_brand (MH)
      double tau_b = tau_brand[b];
      double tau_prop_b = R::rbeta(a_tau, b_tau);
      NumericMatrix Beta_prop_b = clone(Beta_b);
      for (int j = 0; j < n_b; ++j) {
        for (int l = 0; l < q; ++l) {
          Beta_prop_b(j, l) = (tau_prop_b / tau_b) * Beta_b(j,l);
        }
      }
      
      double lp_new_tau = loglik_pvar_cpp(Y_b, X_b, Beta_prop_b, Gamma_b) + R::dbeta(tau_prop_b, a_tau, b_tau, 1);
      double lp_old_tau = loglik_pvar_cpp(Y_b, X_b, Beta_b, Gamma_b) + R::dbeta(tau_b, a_tau, b_tau, 1);
      
      if (std::log(R::runif(0, 1)) < (lp_new_tau - lp_old_tau)) {
        tau_brand[b] = tau_prop_b;
        tau_b = tau_prop_b;
        Beta_b = clone(Beta_prop_b);
      }
      
      // Sample alpha_beta_brand (MH)
      NumericVector alpha_prop_b = rdirichlet_cpp(alpha0);
      NumericVector alpha_b = alpha_beta_brand(b, _);
      double lp_new_alpha = ddirichlet_cpp(alpha_prop_b, alpha0, true); // start with log-prior
      double lp_old_alpha = ddirichlet_cpp(alpha_b, alpha0, true); // start with log-prior
      
      // add log-lik of unscaled beta (beta tilde)
      for (int j = 0; j < n_b; ++j) {
        NumericVector beta_tilde_i(q);
        for (int l = 0; l < q; ++l) {
          beta_tilde_i[l] = Beta_b(j, l) / tau_b;
        }
        lp_new_alpha += ddirichlet_cpp(beta_tilde_i, alpha_prop_b, true);
        lp_old_alpha += ddirichlet_cpp(beta_tilde_i, alpha_b, true);
      }
      
      if (std::log(R::runif(0, 1)) < (lp_new_alpha - lp_old_alpha)) {
        for (int l = 0; l < q; ++l) {
          alpha_beta_brand(b, l) = alpha_prop_b[l];
        }
      }
    }
    
    // Store samples
    Gamma_samples[iter] = clone(Gamma_item);
    Beta_samples[iter] = clone(Beta_item);
    mu_samples[iter] = clone(mu_gamma_brand);
    sigsq_samples[iter] = clone(sigsq_gamma_brand);
    tau_samples[iter] = clone(tau_brand);
    alpha_samples[iter] = clone(alpha_beta_brand);
    
    // Print out iteration after completion
    if (verbose) {
      if ((iter + 1) % 100 == 0) Rcpp::Rcout << "Iteration " << iter + 1 << std::endl;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("Beta_samples") = Beta_samples,           // n_iter list of n_items x q
    Rcpp::Named("Gamma_samples") = Gamma_samples,         // n_iter list of n_items x p
    Rcpp::Named("mu_gamma_samples") = mu_samples,         // n_iter list of B x p
    Rcpp::Named("sigsq_gamma_samples") = sigsq_samples,   // n_iter list of B x p
    Rcpp::Named("tau_samples") = tau_samples,             // n_iter list of B-length vectors
    Rcpp::Named("alpha_samples") = alpha_samples          // n_iter list of B x q
  );
}