#' @title Beta Matrix
#' @description Creates a matrix of beta values for a given number of time points (t) and lag values (q).
#' @param betas Numeric vector of beta values.
#' @param t Integer number of time points.
#' @param q Integer number of lag groups.
#' @return Matrix of beta values.
#' @examples
#' beta_mat(betas = c(.5, 1.5), t = 5, q = 2)
beta_mat <- function(betas, t, q) {
  stopifnot(
    rlang::is_integerish(t),
    rlang::is_integerish(q),
    is.numeric(betas),
    length(betas) == q,
    t > 0,
    q < t
  )

  out <- matrix(0, nrow = t, ncol = t)

  shift <- 0:(t - 1)
  along <- seq_len(q)
  for (i in along) {
    index <- seq(i + 1, t) + (t * (0:(t - i - 1)))
    out[index] <- betas[i]
  }
  out
}

#' @title Log-Likelihood of y and x given parameters
#' @param y response vector,
#' @param x indicator vector
#' @param alpha scalar value for intercept effect
#' @param betas beta vector for each lag time point
#' @param gamma scalar value for effect of `x`
#' @param sigma2 sigma squared
#' @return Log-likelihood value.
#' @examples
#' betas <- c(.1, .2, .3)
#' gamma <-  10
#' alpha <- 3
#' t <- 10
#' x <- c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0)
#' A <- diag(t) - beta_mat(betas, t, q = length(betas))
#' A_inv <- solve(A)
#' A_inv
#' sigma2 <- 1
#' V <- sigma2 * (A_inv %*% t(A_inv))
#' mu <- A_inv %*% (alpha + (gamma * x))
#'
#' y <- mvtnorm::rmvnorm(n = 1, mean = mu, sigma = V)
#'
#' ## log-likelihood on values sampled from our given parameters
#' ll(as.vector(y), x, alpha, betas, gamma, sigma2)
#' ## log-likelihood on dummy values (more negative)
#' ll(1:10, x, alpha, betas, gamma, sigma2)
#' ll(rep(1, 10), x, alpha, betas, gamma, sigma2)
ll <- function(
  y,
  x,
  alpha,
  betas,
  gamma,
  sigma2,
  n = 1
) {
  t <- length(y)
  B <- beta_mat(betas, t, q = length(betas))
  A <- diag(t) - B
  #A_inv <- solve(A)
  #mu <- A_inv %*% (alpha + (gamma * x))
  #V <- sigma2 * (A_inv %*% t(A_inv))

  Ay <- A %*% y
  resid <- Ay - alpha - (gamma * x)

  loglik <- -n * t / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * sum(resid^2)

  loglik
}
