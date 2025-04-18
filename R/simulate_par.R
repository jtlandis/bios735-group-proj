#' Simulate a Poisson Autoregressive (PAR) Time Series
#'
#' @param T Number of timepoints
#' @param beta Vector of lag coefficients (length q)
#' @param gamma Vector of covariate coefficients
#' @param x Matrix of covariates of dimension T × p (optional; generated if NULL)
#' @param mu Scalar intercept
#' @param seed Random seed
#'
#' @return A list with y (counts), x (covariates), and the true parameters
#' @export
simulate_par <- function(T = 200, beta = c(0.2), gamma = c(1.5), x = NULL, mu = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  q <- length(beta)
  p <- length(gamma)
  
  if (is.null(x)) {
    x <- matrix(rbinom(T * p, size = 1, prob = 0.2), nrow = T, ncol = p)
  }
  
  y <- rep(0, T)
  y[1:q] <- rpois(q, lambda = 5)  # initialize with arbitrary counts
  
  for (t in (q + 1):T) {
    y_lags <- y[(t - 1):(t - q)]
    linear_part <- mu + sum(x[t, ] * gamma)
    mix_weight <- 1 - sum(beta)
    m_t <- sum(beta * y_lags) + mix_weight * exp(linear_part)
    y[t] <- rpois(1, lambda = m_t)
  }
  
  list(y = y, x = x, beta = beta, gamma = gamma, mu = mu)
}