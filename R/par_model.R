#' Poisson Autoregressive Log-Likelihood
#'
#' Computes the log-likelihood for a Poisson autoregressive model for a single item.
#'
#' @param y Vector of counts (length T)
#' @param x Matrix of covariates of dimension T × p
#' @param beta Vector of length q (autoregressive weights)
#' @param gamma Vector of length p (covariate effects)
#' @param mu Intercept (optional, scalar)
#' @param q Number of lags
#'
#' @return Log-likelihood value (numeric scalar)
#' @export
par_loglik <- function(y, x, beta, gamma, mu = 0, q = length(beta)) {
  T <- length(y)
  stopifnot(nrow(x) == T)
  
  ll <- 0
  for (t in (q + 1):T) {
    y_lags <- y[(t - 1):(t - q)]
    x_t <- x[t, ]
    
    auto_part <- sum(beta * y_lags)
    linear_part <- mu + sum(x_t * gamma)
    mix_weight <- 1 - sum(beta)
    
    m_t <- auto_part + mix_weight * exp(linear_part)
    if (m_t <= 0) m_t <- 1e-10  # safeguard
    
    ll <- ll + y[t] * log(m_t) - m_t - lgamma(y[t] + 1)
  }
  
  return(ll)
}

#' Fit PAR Model Using Optim
#'
#' Fit a non-hierarchical PAR model using maximum likelihood and optim().
#'
#' @param y Vector of counts
#' @param x Matrix of covariates (T × p)
#' @param q Number of lags
#'
#' @return A list with estimated parameters and log-likelihood
#' @export
fit_par_optim <- function(y, x, q = 1) {
  T <- length(y)
  p <- ncol(x)
  
  init <- c(rep(0.05, q), rep(0.01, p))  # beta + gamma
  
  nll <- function(par) {
    beta <- par[1:q]
    gamma <- par[(q + 1):(q + p)]
    -par_loglik(y, x, beta, gamma, mu = 0, q = q)
  }
  
  grad <- function(par) {
    beta <- par[1:q]
    gamma <- par[(q + 1):(q + p)]
    -par_gradient(y, x, beta, gamma, mu = 0, q = q)
  }
  
  fit <- optim(
    par = init,
    fn = nll,
    gr = grad,
    method = "BFGS",
    control = list(maxit = 200)
  )
  
  list(
    beta = fit$par[1:q],
    gamma = fit$par[(q + 1):(q + p)],
    loglik = -fit$value,
    converged = fit$convergence == 0,
    fit = fit
  )
}

#' Gradient of PAR Log-Likelihood
#'
#' Computes the gradient of the log-likelihood for a Poisson autoregressive model.
#'
#' @param y Vector of observed counts (length T)
#' @param x Matrix of covariates (T × p)
#' @param beta Vector of length q
#' @param gamma Vector of length p
#' @param mu Intercept (scalar)
#' @param q Number of lags
#'
#' @return Named numeric vector of length q + p (gradients of beta and gamma)
#' @export
par_gradient <- function(y, x, beta, gamma, mu = 0, q = length(beta)) {
  T <- length(y)
  p <- length(gamma)
  
  grad_beta <- rep(0, q)
  grad_gamma <- rep(0, p)
  
  for (t in (q + 1):T) {
    y_lags <- y[(t - 1):(t - q)]
    x_t <- x[t, ]
    linear_part <- mu + sum(x_t * gamma)
    
    a_t <- sum(beta * y_lags)
    b_t <- (1 - sum(beta)) * exp(linear_part)
    m_t <- a_t + b_t
    
    prefactor <- (y[t] / m_t) - 1
    
    for (l in 1:q) {
      grad_beta[l] <- grad_beta[l] + prefactor * (y[t - l] - exp(linear_part))
    }
    
    for (j in 1:p) {
      grad_gamma[j] <- grad_gamma[j] + prefactor * (1 - sum(beta)) * x_t[j] * exp(linear_part)
    }
  }
  
  grad <- c(beta = grad_beta, gamma = grad_gamma)
  return(grad)
}