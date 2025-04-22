library(MCMCpack)
library(LaplacesDemon)
library(MASS)
#setwd('~/Dropbox/UNC/Spring2025/BIOS735/project/')


# Functions ---------------------------------------------------------------
####### Helper Functions for forecasting
split_data <- function(Y, X, q = 5, prop.train = 0.75) {
  n <- length(Y)
  T_train <- round(prop.train * n)
  
  # Define indices
  train_inds <- 1:T_train
  test_inds <- (T_train - q + 1):n  # includes tail of training
  
  list(
    Ytrain = Y[train_inds],
    Xtrain = X[train_inds,],
    Ytest = Y[test_inds],
    Xtest = X[test_inds,]
  )
}

# extract estimated fit/one-step forecasts based on posterior samples
get_fit_MCMC <- function(Y, X, samples, burn.in = 0.5) {
  n <- length(Y)
  p <- ncol(X)
  q <- ncol(samples$beta)
  n_iter <- length(samples$tau)
  ss_idx <- round(burn.in * n_iter):n_iter
  
  mt_samples <- array(NA, dim = c(n, length(ss_idx)))
  
  for (s in 1:length(ss_idx)) {
    iter <- ss_idx[s]
    beta <- samples$beta[iter, ]
    gamma <- samples$gamma[iter, ]
    
    m_t <- rep(NA, n)
    for (t in (q + 1):n) {
      ar_term <- sum(beta * Y[(t - 1):(t - q)])
      cov_term <- (1 - sum(beta)) * exp(sum(X[t, ] * gamma))
      m_t[t] <- ar_term + cov_term
    }
    mt_samples[, s] <- m_t
  }
  mt_fit <- apply(mt_samples, 1, mean)
  return(mt_fit)
}

# get posterior means of parameters after MCMC
get_estimates_MCMC <- function(samples, burn.in = 0.5) {
  n_iter <- length(samples$tau)
  ss_idx <- round(burn.in * n_iter):n_iter
  beta_est <- apply(samples$beta[ss_idx,,drop=F], 2, mean)
  gamma_est <- apply(samples$gamma[ss_idx,,drop=F], 2, mean)
  return(list(beta = beta_est, gamma = gamma_est))
}

# get full T-step forecast using posterior mean of parameters
get_forecast_MCMC <- function(Y, X, samples, burn.in = 0.5, reps = 10) {
  estimates <- get_estimates_MCMC(samples, burn.in)
  beta <- estimates$beta
  gamma <- estimates$gamma
  
  n <- length(Y)
  q <- length(beta)
  sum_beta <- sum(beta)
  
  Y_forecast_reps <- array(dim=c(n, reps))
  Y_forecast_reps[1:q,] <- Y[1:q]
  for (s in 1:reps) {
    for (t in (q+1):n) {
      ar_term <- sum(beta * Y_forecast_reps[(t-1):(t-q),s])
      x_term <- exp(sum(gamma * X[t,]))
      mt <- ar_term + (1 - sum_beta) * x_term
      Y_forecast_reps[t,s] <- rpois(1, mt)
    }
  }
  Y_forecast <- round(apply(Y_forecast_reps,1,mean))
  return(Y_forecast)
}

# compare fit/forecast with observed
evaluate <- function(Y, Ypred, q = 5) {
  Y <- Y[-c(1:q)]
  Ypred <- Ypred[-c(1:q)]
  RMSE <- sqrt(mean((Y - Ypred)^2))
  RMSE
}

###### MCMC / Metropolis-within-Gibbs Algorithm
get_loglik <- function(Y, X, beta, gamma) {
  n <- length(Y)
  q <- length(beta)
  loglik <- 0
  sum_beta <- sum(beta)
  for (t in (q+1):n) {
    ar_term <- sum(beta * Y[(t-1):(t-q)])
    x_term <- exp(sum(gamma * X[t,]))
    mt <- ar_term + (1 - sum_beta) * x_term
    loglik <- loglik + Y[t]*log(mt) - mt
  }
  return(loglik)
}

run_MCMC <- function(Y, X, q = 5, n_iter = 1000,
                     hyperparams=NULL, proposal_sd = 0.05) {
  # Set hyperparameters
  n <- length(Y)
  p <- ncol(X)
  if (is.null(hyperparams)) {
    hyperparams <- vector(mode="list", 6)
    names(hyperparams) <- c("mu", "Psi", "nu", "alpha","a","b")
    hyperparams$mu <- rep(0, p)
    hyperparams$Psi <- diag(p)
    hyperparams$nu <- p + 2
    hyperparams$alpha <- rep(1/q, q)
    hyperparams$a <- 1
    hyperparams$b <- 1
  }
  
  # Initialize parameters
  beta_tilde <- rdirichlet(1, hyperparams$alpha)
  tau <- rbeta(1, hyperparams$a, hyperparams$b)
  beta <- tau * beta_tilde
  Sigma <- rinvwishart(hyperparams$nu, hyperparams$Psi)
  gamma <- mvrnorm(1, hyperparams$mu, Sigma)
  samples <- list(beta = matrix(NA, n_iter, q),
                  gamma = matrix(NA, n_iter, p),
                  tau = rep(NA, n_iter),
                  Sigma = array(NA, dim=c(p,p,n_iter)))
  
  # Run MCMC
  for (iter in 1:n_iter) {
    # --- Sample gamma (MH) ---
    gamma_prop <- mvrnorm(1, gamma, proposal_sd * diag(p))
    lp_new <- get_loglik(Y, X, beta, gamma_prop) + dmvn(gamma_prop, hyperparams$mu, Sigma, log = TRUE)
    lp_old <- get_loglik(Y, X, beta, gamma) + dmvn(gamma, hyperparams$mu, Sigma, log = TRUE)
    if (log(runif(1)) < lp_new - lp_old) {
      gamma <- gamma_prop
    }
    
    # --- Sample Sigma (Gibbs) ---
    resid <- gamma - mu
    S <- tcrossprod(resid)
    Sigma <- rinvwishart(hyperparams$nu + 1, hyperparams$Psi + S)
    
    # --- Sample tau (MH) ---
    tau_prop <- rbeta(1, hyperparams$a, hyperparams$b)
    beta_prop <- as.vector(tau_prop * beta_tilde)
    lp_new <- get_loglik(Y, X, beta_prop, gamma) + dbeta(tau_prop, hyperparams$a, hyperparams$b, log=TRUE)
    lp_old <- get_loglik(Y, X, beta, gamma) + dbeta(tau, hyperparams$a, hyperparams$b, log=TRUE)
    if (log(runif(1)) < lp_new - lp_old) {
      tau <- tau_prop
      beta <- beta_prop
    }
    
    # --- Sample beta_tilde (MH) ---
    beta_tilde_prop <- rdirichlet(1, hyperparams$alpha)
    beta_prop <- as.vector(tau * beta_tilde_prop)
    lp_new <- get_loglik(Y, X, beta_prop, gamma) + LaplacesDemon::ddirichlet(beta_tilde_prop, hyperparams$alpha, log=TRUE)
    lp_old <- get_loglik(Y, X, beta, gamma) + LaplacesDemon::ddirichlet(beta_tilde, hyperparams$alpha, log=TRUE)
    if (log(runif(1)) < lp_new - lp_old) {
      beta_tilde <- beta_tilde_prop
      beta <- beta_prop
    }
    
    # --- Store samples ---
    samples$beta[iter,] <- beta
    samples$gamma[iter,] <- gamma
    samples$tau[iter] <- tau
    samples$Sigma[,,iter] <- Sigma
    
    if (iter %% 100 == 0) cat("Iteration", iter, "\n")
  }
  return(samples)
}




# Simulate Data -----------------------------------------------------------
set.seed(12)
n <- 1000
p <- 2
q <- 3
beta_tilde <- rdirichlet(1, alpha = rep(1/q, q))
tau <- rbeta(n = 1, 2, 1)
beta <- beta_tilde * tau
mu <- rep(1,p)
Sigma <- rinvwishart(nu = p+2, S = diag(rep(1,p)))
gamma <- mvrnorm(n = 1, mu, Sigma)
Y0 <- 1
X <- cbind(1, rbinom(n, 1, prob = 0.25))

Y <- rep(NA,n)
mt <- rep(NA,n)
Y[1:q] <- rpois(q, Y0)
for (t in (q+1):n) {
  mt[t] <- sum(beta * Y[(t-1):(t-q)]) + (1-sum(beta)) * exp(sum(gamma*X[t,]))
  Y[t] <- rpois(n = 1, lambda = mt[t])
}
plot(Y,type="l")


# Evaluate on Simulated Data ----------------------------------------------

set.seed(123)
prop.train = 0.75
n_iter <- 2000
burn.in <- 0.5
q <- 3

# Split data
split <- split_data(Y, X, q, prop.train)

# Fit to training
samples <- run_MCMC(split$Ytrain, split$Xtrain, q, n_iter)

# Get estimates
estimates <- get_estimates_MCMC(samples, burn.in = 0.5)

# Get fit on training
train_fit <- get_fit_MCMC(split$Ytrain, split$Xtrain, samples, burn.in=0.5)
plot(split$Ytrain,type="l")
lines(train_fit, col="blue")
evaluate(split$Ytrain, train_fit, q)

# Get one-step ahead forecasts (uses previous test observation to make predictions)
forecast_onestep <- get_fit_MCMC(split$Ytest, split$Xtest, samples)
plot(split$Ytest,type="l")
lines(forecast_onestep, col="blue")
evaluate(split$Ytest, forecast_onestep, q)

# Get full forecast on test set (not using any observations in test set to make prediction)
forecast_full <- get_forecast_MCMC(split$Ytest, split$Xtest, samples, burn.in=0.5)
plot(split$Ytest,type="l")
lines(forecast_full, col="blue")
evaluate(split$Ytest, forecast_full, q)


# Fit to Real Data ------------------------------------------------
data <- read.csv('hierarchical_sales_data.csv')

Y <- data$QTY_B1_1
X <- cbind(1, data$PROMO_B1_1)

set.seed(123)
samples <- run_MCMC(Y, X, q = 3, n_iter = 2000)
estimates <- get_estimates_MCMC(samples, burn.in = 0.5)
fit <- get_fit_MCMC(Y, X, samples, burn.in=0.5)

plot(Y, type = "l")
lines(fit, col = "blue")
evaluate(Y, fit, q)