library(MCMCpack)
library(ggplot2)

# Apply Stan to Real Data -------------------------------------------------
set.seed(123)

data <- read.csv('hierarchical_sales_data.csv')

# Data
y <- data$QTY_B1_1
X <- as.matrix(cbind(1,data$PROMO_B1_1))
T_obs <- length(y)
q <- 5
p <- ncol(X)

# Priors
mu_gamma <- rep(0, p)
Sigma_gamma <- diag(p)
alpha <- rep(1, q)
a_tau <- 2
b_tau <- 2

# Stan data
stan_data <- list(
  T_obs = T_obs,
  q = q,
  p = p,
  y = y,
  X = X,
  mu_gamma = mu_gamma,
  Sigma_gamma = Sigma_gamma,
  alpha = alpha,
  a_tau = a_tau,
  b_tau = b_tau
)

# Sample with stan
fit <- stan(
  file = "poisson_ar.stan",
  data = stan_data,
  iter = 1000,
  chains = 4,
  seed = 123
)

print(fit, pars = c("tau", "beta_tilde", "gamma"))


# Assess goodness-of-fit
posterior_samples <- rstan::extract(fit)
n_draws <- length(posterior_samples$tau)
mt_draws <- matrix(NA, nrow = n_draws, ncol = T_obs)

for (s in 1:n_draws) {
  beta <- posterior_samples$tau[s] * posterior_samples$beta_tilde[s, ]
  gamma <- posterior_samples$gamma[s, ]
  
  for (t in (q + 1):T_obs) {
    ar_part <- sum(beta * y[(t - 1):(t - q)])
    cov_part <- exp(X[t, ] %*% gamma)
    mt_draws[s, t] <- ar_part + (1 - sum(beta)) * cov_part
  }
}

mt_mean <- colMeans(mt_draws, na.rm=T)
mt_lower <- apply(mt_draws, 2, quantile, probs = 0.025, na.rm = TRUE)
mt_upper <- apply(mt_draws, 2, quantile, probs = 0.975, na.rm = TRUE)

df_plot <- data.frame(
  time = 1:T_obs,
  y_obs = y,
  y_fit = mt_mean,
  lower = mt_lower,
  upper = mt_upper
)

ggplot(df_plot[(q+1):T_obs,], aes(x = time)) +
  geom_line(aes(y = y_obs), color = "black", size = 1.5, alpha = 0.5) +
  geom_line(aes(y = y_fit), color = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
  labs(y = "Counts", title = "Posterior Predictive Mean and 95% CI") +
  theme_minimal()
