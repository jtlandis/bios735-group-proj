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
  p <- length(gamma)
  if (length(gamma) == 0) gamma <- 0
  ll <- 0
  beta_seq <- seq_len(q)
  gamma_seq <- seq_len(p) + q
  for (t in seq_len(T)) {
    y_lags <- x[t, beta_seq]
    x_t <- x[t, gamma_seq]

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
  p <- ncol(x) - q
  beta_seq <- seq_len(q)
  gamma_seq <- seq_len(p) + q
  init <- c(rep(0.05, q), rep(0.01, p)) # beta + gamma

  nll <- function(par) {
    beta <- par[beta_seq]
    gamma <- par[gamma_seq]
    -par_loglik(y, x, beta, gamma, mu = 0, q = q)
  }

  grad <- function(par) {
    beta <- par[beta_seq]
    gamma <- par[gamma_seq]
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
    beta = fit$par[beta_seq],
    gamma = fit$par[gamma_seq],
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

  beta_seq <- seq_len(q)
  gamma_seq <- seq_len(p) + q

  grad_beta <- rep(0, q)
  grad_gamma <- rep(0, p)

  grad <- c(
    beta = grad_beta,
    gamma = grad_gamma
  )

  for (t in seq_len(T)) {
    y_lags <- x[t, beta_seq]
    x_t <- x[t, gamma_seq]
    linear_part <- mu + sum(x_t * gamma)

    a_t <- sum(beta * y_lags)
    b_t <- (1 - sum(beta)) * exp(linear_part)
    m_t <- a_t + b_t

    prefactor <- (y[t] / m_t) - 1

    for (l in beta_seq) {
      grad[l] <- grad[l] + prefactor * (y_lags[l] - exp(linear_part))
    }

    for (j in gamma_seq) {
      grad[j] <- grad[j] +
        prefactor * (1 - sum(beta)) * x_t[j - q] * exp(linear_part)
    }
  }

  grad
}

par_bfgs <- function(
  y,
  x,
  beta = numeric(0),
  gamma,
  max_iter = 1000,
  trace_mod = 0,
  ...
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
    H = diag(q + p), #par_hessian(y, x, beta, gamma) |> solve(),
    max_iter = max_iter,
    trace_mod = trace_mod,
    dynamic_steps = FALSE,
    ...
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

assert_no_dup_time <- function(data, time, .call = rlang::caller_env()) {
  time <- rlang::enquo(time)
  any_dup <- data |>
    summarise(.dup = anyDuplicated(!!time), .groups = "drop")
  if (any(any_dup$.dup)) {
    if (ncol(any_dup) == 1) {
      rlang::abort(
        c(
          sprintf(
            "Duplicate times found in column `%s`",
            rlang::as_label(time)
          ),
          "!" = "only one time point is permitted per group",
          "i" = "consider subsetting your data, or use `dplyr::group_by(data, ...)`"
        ),
        call = .call
      )
    } else {
      w <- filter(any_dup, .dup > 0) |> select(-.dup)
      rlang::abort(
        c(
          "Duplicate times found in groups:",
          {
            w_names <- names(w)
            seq_row <- seq_len(nrow(w))
            l <- vapply(
              seq_row,
              function(i) {
                paste0(w_names, ": ", unlist(w[i, ]), collapse = ", ")
              },
              character(1)
            )
            n <- length(l)
            names(l) <- rep("*", n)
            if (n < 4) {
              l
            } else {
              c(l[1:3], "i" = paste0("... and ", n - 3, " more"))
            }
          },
          "!" = "only one time point is permitted per group",
          "i" = "consider subsetting your data, or use `dplyr::group_by(data, ...)`"
        ),
        call = .call
      )
    }
  }
  invisible(NULL)
}

# assume data is grouped...
# mm_ <- par_model_mat(group_by(filter(data_set_tidy,brand %in% c("B1", "B2"), item %in% c("1","2")), brand, item), QTY ~ PROMO*(brand:item), time = DATE, nlag=4)
# bfgs_cpp(mm_$Y, mm_$X, )

#' Design matrix for analysis
#' @export
par_model_mat <- function(
  data,
  formula,
  time = NULL,
  nlag = 1
) {
  time <- rlang::enexpr(time)
  q <- nlag
  stopifnot(
    rlang::is_formula(formula)
  )
  assert_no_dup_time(data, !!time)
  y_sym <- rlang::f_lhs(formula)
  covar <- rlang::f_rhs(formula)
  seq_ <- seq_len(q)
  lags <- setNames(seq_, sprintf("lag%i", seq_)) |>
    lapply(function(l) expr(lag(!!y_sym, !!l)))
  lag_names <- rlang::syms(names(lags))
  # check if any duplicate times

  data <- data |>
    arrange(!!time) |>
    mutate(!!!lags) |>
    slice(-seq_len(.env$q))

  lag_formula <- switch(
    match(length(lag_names), c(0, 1, 2), nomatch = 4L),
    NULL,
    lag_names[[1]],
    call("+", lag_names[[1]], lag_names[[2]]),
    {
      out <- call("+", lag_names[[1]], lag_names[[2]])
      for (i in 3:length(lag_names)) {
        out <- call("+", out, lag_names[[i]])
      }
      out
    }
  )

  if (!is.null(lag_formula)) {
    covar <- expr(!!lag_formula + !!covar)
  }

  formula2 <- expr(~!!covar)

  X <- eval(
    expr(
      model.matrix(
        object = !!formula2,
        data = data
      )
    )
  )
  #insure
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, c(names(lags), colnames(X)[-which(colnames(X) %in% names(lags))])]
  }
  Y <- pull(data, !!y_sym)
  xnames <- colnames(X)
  list(
    Y = Y,
    X = X,
    beta = setNames(rep(0, q), names(lags)),
    gamma = setNames(rep(0, ncol(X) - q), xnames[-seq_len(q)]),
    .data = data
  )
}
