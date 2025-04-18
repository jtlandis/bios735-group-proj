// Vector Poisson Autoregressive (PAR) Model with Brand, Item, and Time Effects

data {
  int<lower=1> N;             // total number of observations
  int<lower=1> I;             // number of items
  int<lower=1> B;             // number of brands
  int<lower=1> Q;             // number of lags
  int<lower=1> P;             // number of covariates (e.g., promo)
  int<lower=1> T;             // number of unique time points

  array[N] int<lower=0> y;          // outcome
  array[N] int<lower=1, upper=I> item;
  array[I] int<lower=1, upper=B> brand;
  matrix[N, P] x;             // covariate matrix
  matrix[N, Q] lag_y;         // lagged outcome matrix
  array[N] int<lower=1, upper=T> time; // time index per row
}

parameters {
  // Item-specific effects
  array[I] vector[P] gamma;         // covariate effects
  array[I] simplex[Q] beta_raw;     // normalized AR weights (sum to 1)
  array[I] real<lower=0, upper=1> tau; // total AR weight
  array[I] real eta;                // item-level intercepts

  // Brand-level priors
  array[B] vector[P] mu_g;
  array[B] real mu_brand;           // brand-level intercepts
  real<lower=0> sigma_g;
  real<lower=0> sigma_eta;

  // Shared effects
  real mu;                          // global intercept
  vector[T] f;                      // shared time effect
}

transformed parameters {
  vector[N] m;
  for (n in 1:N) {
    vector[Q] beta_i = tau[item[n]] * beta_raw[item[n]];
    real ar_part = dot_product(beta_i, lag_y[n]);
    real lin_pred = mu + mu_brand[brand[item[n]]] + eta[item[n]] + f[time[n]] + dot_product(gamma[item[n]], x[n]);
    real cov_part = (1 - tau[item[n]]) * exp(lin_pred);
    m[n] = ar_part + cov_part;
  }
}

model {
  // Hyperpriors
  sigma_g ~ normal(0, 1);
  sigma_eta ~ normal(0, 1);

  // Brand-level priors
  for (b in 1:B) {
    mu_g[b] ~ normal(0, 1);
    mu_brand[b] ~ normal(0, 1);
  }

  // Item-level priors
  for (i in 1:I) {
    gamma[i] ~ normal(mu_g[brand[i]], sigma_g);
    eta[i] ~ normal(0, sigma_eta);
    tau[i] ~ beta(2, 2);
    beta_raw[i] ~ dirichlet(rep_vector(1.0, Q));
  }

  // Time effect prior
  f ~ normal(0, 1);

  // Observation model
  y ~ poisson(m);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_lpmf(y[n] | m[n]);
  }
}
