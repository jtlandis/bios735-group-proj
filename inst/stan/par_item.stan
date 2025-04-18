// Poisson Autoregression (PAR) Model for a Single Item

// Data block
// -------------
data {
  int<lower=1> T_obs;               // number of timepoints
  int<lower=1> q;                   // AR lag order
  int<lower=1> p;                   // number of covariates

  array[T_obs] int<lower=0> y;     // outcome time series
  matrix[T_obs, p] X;              // covariate matrix

  // prior hyperparameters
  vector[p] mu_gamma;              // prior mean for gamma
  cov_matrix[p] Sigma_gamma;       // prior covariance for gamma
  vector[q] alpha;                 // Dirichlet prior for beta_tilde
  real<lower=0> a_tau;             // Beta prior shape 1 for tau
  real<lower=0> b_tau;             // Beta prior shape 2 for tau
}

// Parameters block
// ------------------
parameters {
  vector[p] gamma;                 // covariate effect
  simplex[q] beta_tilde;          // normalized AR weights
  real<lower=0, upper=1> tau;     // total AR weight
}

// Transformed parameters block
// -----------------------------
transformed parameters {
  vector[q] beta = tau * beta_tilde; // scaled AR weights
}

// Model block
// -------------
model {
  // Priors
  gamma ~ multi_normal(mu_gamma, Sigma_gamma);
  tau ~ beta(a_tau, b_tau);
  beta_tilde ~ dirichlet(alpha);

  // Likelihood
  for (t in (q + 1):T_obs) {
    real ar_part = 0;
    for (l in 1:q) {
      ar_part += beta[l] * y[t - l];
    }
    real cov_part = exp(dot_product(X[t], gamma));
    real m_t = ar_part + (1 - sum(beta)) * cov_part;
    y[t] ~ poisson(m_t);
  }
}

// Generated quantities (optional)
// --------------------------------
generated quantities {
  vector[T_obs - q] log_lik;
  for (t in (q + 1):T_obs) {
    real ar_part = 0;
    for (l in 1:q) {
      ar_part += beta[l] * y[t - l];
    }
    real cov_part = exp(dot_product(X[t], gamma));
    real m_t = ar_part + (1 - sum(beta)) * cov_part;
    log_lik[t - q] = poisson_lpmf(y[t] | m_t);
  }
}
