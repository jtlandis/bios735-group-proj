//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> T_obs; // total number of timepoints
  int<lower=1> q; // AR lag order
  int<lower=1> p; // number of covariates
  int<lower=0> y[T_obs]; // observed counts
  matrix[T_obs, p] X; // covariate matrix
  
  // hyperparams
  vector[p] mu_gamma; // prior mean for gamma
  cov_matrix[p] Sigma_gamma; // prior covariance for gamma
  vector[q] alpha; // Dirichlet prior params for beta_tilde
  real<lower=0> a_tau; // prior shape1 for tau
  real<lower=0> b_tau; // prior shape2 for tau
}

// The parameters accepted by the model
parameters {
  real<lower=0, upper=1> tau;
  simplex[q] beta_tilde; // normalized AR weights
  vector[p] gamma;
}

transformed parameters {
  vector[q] beta = tau * beta_tilde;
}

model {
  // priors
  gamma ~ multi_normal(mu_gamma, Sigma_gamma);
  tau ~ beta(a_tau, b_tau);
  beta_tilde ~ dirichlet(alpha);
  
  // likelihood
  for (t in (q + 1):T_obs) {
    real ar_part = 0;
    for (l in 1:q) {
      ar_part += beta[l] * y[t-l];
    }
    real covariate_part = exp(dot_product(row(X,t), gamma));
    real mt = ar_part + (1 - sum(beta)) * covariate_part;
    y[t] ~ poisson(mt);
  }
}
