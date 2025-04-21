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
  if (length(gamma) == 0) gamma <- 0
  ll <- 0
  for (t in (q + 1):T) {
    y_lags <- y[(t - 1):(t - q)]
    x_t <- x[t, ]

    auto_part <- sum(beta * y_lags)
    linear_part <- mu + sum(x_t * gamma)
    mix_weight <- 1 - sum(beta)

    m_t <- auto_part + mix_weight * exp(linear_part)
    if (m_t <= 0) m_t <- 1e-10 # safeguard

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

  init <- c(rep(0.05, q), rep(0.01, p)) # beta + gamma

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
      grad_gamma[j] <- grad_gamma[j] +
        prefactor * (1 - sum(beta)) * x_t[j] * exp(linear_part)
    }
  }

  grad <- c(beta = grad_beta, gamma = grad_gamma)
  grad[!is.na(grad)]
}



par_bfgs <- function(
  y,
  x,
  beta = numeric(0),
  gamma,
  max_iter = 1000,
  trace_mod = 0
) {
  q <- length(beta)
  p <- length(gamma)

  beta_slice <- seq_len(q)
  gamma_slice <- seq_len(p) + q

  ll <- function(par, y, x) {
    beta <- par[beta_slice]
    gamma <- par[gamma_slice]
    par_loglik(y, x, beta, gamma, mu = 0, q = q)
  }

  grad <- function(par, y, x) {
    beta <- par[beta_slice]
    gamma <- par[gamma_slice]
    par_gradient(y, x, beta, gamma, mu = 0, q = q)
  }

  general_bfgs(
    par = c(beta, gamma),
    f = ll,
    g = grad,
    y = y,
    x = x,
    H = par_hessian(y, x, beta, gamma) |> solve(),
    max_iter = max_iter,
    step_size = 0.1,
    trace_mod = trace_mod
  )
}

# utility function to shift some data
# vector of length `size` by some ammount
# n <= q.
# shift of 0 returns model data
shifter <- function(q, size) {
  seq_ <- seq_len(size - q)
  function(data, n) {
    if (n > q) {
      stop(sprintf("n must be less than or equal to %i", q))
    }
    data[seq_ + q - n]
  }
}

par_hessian <- function(y, x, beta, gamma) {
  q <- length(beta)
  p <- length(gamma)
  T <- length(y)

  shift <- shifter(q, T)
  y_lags_slice <- rev(seq_len(q) - 1L)
  x_t <- x[(q + 1):T, , drop = FALSE]
  m_t <- Auto_t <- numeric(T - q)
  Linear_t <- if (p > 0) as.vector(exp(x_t %*% gamma)) else rep(1, T - q)

  weight_ <- 1 - sum(beta)
  for (t in seq_along(Auto_t)) {
    y_lags_slice <- y_lags_slice + 1L
    y_lags <- y[y_lags_slice]
    Auto_t[t] <- sum(beta * y_lags)
    m_t[t] <- Auto_t[t] + (weight_ * Linear_t[t])
  }
  # create Hessian
  P <- p + q
  H <- matrix(0, P, P)
  param_types <- c(rep("beta", q), rep("gamma", p))
  yt_mt <- shift(y, 0) / m_t
  w <- yt_mt / m_t
  for (i in seq_along(param_types)) {
    i_block <- param_types[i]

    for (j in i:P) {
      j_block <- param_types[j]
      block <- sprintf("%s_%s", i_block, j_block)
      switch(
        block,
        beta_beta = {
          H[i, j] <- H[j, i] <- -sum(
            w * (shift(y, i) - Linear_t) * (shift(y, j) - Linear_t)
          )
        },
        beta_gamma = {
          gamma_j <- j - q
          x_j <- x_t[, gamma_j, drop = FALSE]
          H[i, j] <- H[j, i] <- -sum(
            w *
              (shift(y, i) - Linear_t) *
              weight_ *
              Linear_t *
              x_j +
              ((yt_mt - 1) * Linear_t * x_j)
          )
        },
        gamma_beta = {
          gamma_i <- i - q
          x_i <- x_t[, gamma_i, drop = FALSE]
          H[i, j] <- H[j, i] <- -sum(
            w *
              (shift(y, j) - Linear_t) *
              weight_ *
              Linear_t *
              x_i +
              ((yt_mt - 1) * Linear_t * x_i)
          )
        },
        gamma_gamma = {
          gamma_i <- i - q
          x_i <- x_t[, gamma_i, drop = FALSE]
          gamma_j <- j - q
          x_j <- x_t[, gamma_j, drop = FALSE]
          H[i, j] <- H[j, i] <- -sum(
            x_i *
              x_j *
              (
                w *
                  (weight_^2) *
                  (Linear_t^2) +
                  (Linear_t * (yt_mt - 1) * weight_)
              )
          )
        }
      )
    }
  }

  H
}