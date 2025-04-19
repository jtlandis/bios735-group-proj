#' @title Maximum Likelihood Estimation for NAR Model Parameters
#' @name mle_nar_model
#' @description This function estimates the parameters of a Normal Autoregressive (NAR) model using Maximum Likelihood Estimation (MLE).
#' @param y a matrix of dimension t by n
#' @param x a matrix of dimension t by n
#' @param alpha a scalar coefficient for the intercept term
#' @param betas a vector of coefficients for the autoregressive terms
#' @param gamma a scalar coefficient for the promotion
#' @return an estiamte for the respective Parameter
NULL

#' @rdname mle_nar_model
mle_alpha <- function(y, x, betas, gamma) {
  t <- nrow(y)
  n <- ncol(y)
  B <- beta_mat(betas, t)
  A <- diag(t) - B
  Ay <- A %*% y
  gammax <- gamma * x
  sum(Ay - gammax) / (t * n)
}

#' @rdname mle_nar_model
mle_gamma <- function(y, x, betas, alpha) {
  t <- nrow(y)
  B <- beta_mat(betas, t, q = length(betas))
  A <- diag(t) - B
  Ay <- A %*% y
  sum(Ay * x - (alpha * x)) / sum(x * x)
}

#' @rdname mle_nar_model
mle_sigma2 <- function(y, x, betas, alpha, gamma) {
  t <- nrow(y)
  n <- ncol(y)
  B <- beta_mat(betas, t)
  A <- diag(t) - B
  Ay <- A %*% y
  gammax <- gamma * x
  r <- Ay - gammax - alpha
  sum(r * r) / (t * n)
}

shift_by <- function(x, n) {
  N <- length(x)
  c(rep(0, n), x[seq_len(N - n)])
}

shift_mat_row_by <- function(mat, n) {
  index_zero <- seq_len(n)
  index_shift <- seq_len(nrow(mat) - n)
  mat[-index_zero, ] <- mat[index_shift, ]
  mat[index_zero, ] <- 0
  mat
}

align_head_tail <- function(vec_head, vec_tail, n) {
  list(
    vec_head[seq_len(length(vec_head) - n)],
    vec_tail[-seq_len(n)]
  )
}

align_head_tail_mat <- function(vec_head, vec_tail, n) {
  list(
    vec_head[seq_len(nrow(vec_head) - n), ],
    vec_tail[-seq_len(n), ]
  )
}

#' @rdname mle_nar_model
mle_beta <- function(y, x, alpha, betas, gamma, index) {
  n <- ncol(y)
  # the derivative wrt beta_index
  z <- shift_mat_row_by(y, index) * -1
  beta_sums <- vapply(
    seq_along(betas),
    function(i, betas, y, z) {
      beta <- betas[i]
      aligned <- align_head_tail_mat(y, z, i)
      beta * sum(aligned[[1]] * aligned[[2]])
    },
    numeric(1),
    betas = betas,
    y = y,
    z = z
  )
  time_sums <- sum(y * z)
  scale <- beta_sums[index] / betas[index]
  beta_sums <- beta_sums[-index]

  mle <- time_sums - sum(beta_sums) - sum(alpha * z) - sum(gamma * x * z)

  mle / scale
}

mle_params <- function(y, x, alpha, betas, gamma, sigma2) {
  alpha_prior <- alpha
  betas_prior <- betas
  gamma_prior <- gamma
  sigma2_prior <- sigma2

  alpha <- mle_alpha(y, x, betas = betas_prior, gamma = gamma_prior)
  betas <- vapply(
    X = seq_along(betas),
    FUN = mle_beta,
    FUN.VALUE = numeric(1),
    y = y,
    x = x,
    alpha = alpha_prior,
    betas = betas_prior,
    gamma = gamma_prior
  )
  gamma <- mle_gamma(y, x, betas = betas_prior, alpha = alpha_prior)
  sigma2 <- mle_sigma2(
    y,
    x,
    betas = betas_prior,
    alpha = alpha_prior,
    gamma = gamma_prior
  )

  grid <- do.call(
    expand.grid,
    args = c(
      list(
        alpha = seq(alpha_prior, alpha, len = 2)
      ),
      stats::setNames(
        Map(
          function(b, bp) {
            seq(b, bp, len = 2)
          },
          b = betas,
          bp = betas_prior
        ),
        nm = sprintf("beta_%i", seq_along(betas))
      ),
      list(
        gamma = seq(gamma_prior, gamma, len = 2),
        sigma2 = c(sigma2_prior, sigma2, len = 2)
      )
    )
  )
  beta_seq <- grep("beta_", colnames(grid))
  ll_out <- numeric(nrow(grid))
  alphas <- grid$alpha
  gammas <- grid$gamma
  sigmas <- grid$sigma2
  for (i in seq_len(nrow(grid))) {
    ll_out[i] <- ll(
      y = y,
      x = x,
      alpha = alphas[i],
      betas = unlist(grid[i, beta_seq]),
      gamma = gammas[i],
      sigma2 = sigmas[i]
    )
  }
  #out <- dplyr::rowwise(grid) |>
  #  dplyr::mutate(
  #    ll = ll(
  #      y = .env$y,
  #      x = .env$x,
  #      alpha = alpha,
  #      betas = dplyr::c_across(dplyr::starts_with("beta_")),
  #      gamma = gamma,
  #      sigma2 = sigma2
  #    )
  #  ) |>
  #  dplyr::ungroup() |>
  #  dplyr::slice(which.max(ll))

  index <- which.max(ll_out)

  list(
    alpha = alphas[index],
    betas = c(unlist(grid[index, beta_seq])),
    gamma = gammas[index],
    sigma2 = sigmas[index],
    ll = ll_out[index]
  )
}

## --- Work in progres! ---

#' @rdname mle_nar_model
#' @examples
#' betas <- c(.1, .2, .3)
#' gamma <-  10
#' alpha <- 3
#' t <- 100
#' x <- rep(c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0), 10)
#' A <- diag(t) - beta_mat(betas, t, q = length(betas))
#' A_inv <- solve(A)
#' A_inv
#' sigma2 <- 1
#' V <- sigma2 * (A_inv %*% t(A_inv))
#' mu <- A_inv %*% (alpha + (gamma * x))
#'
#' y <- mvtnorm::rmvnorm(n = 1, mean = mu, sigma = V)
#' bad_alphas <- 0
#' bad_betas <- c(1, 1, 1)
#' bad_gamma <- 0
#' bad_sigma2 <- 2
#' #should adjust to a while loop that
#' # compares prior iter.
#' optim_mle_params(
#'  y = y,
#'  x = x,
#'  alpha = bad_alphas,
#'  betas = bad_betas,
#'  gamma = bad_gamma,
#'  sigma2 = bad_sigma2,
#'  times = 150,
#'  trace_mod = 10)
optim_mle_params <- function(
  y,
  x,
  alpha,
  betas,
  gamma,
  sigma2,
  times = 100,
  trace_mod = 1
) {
  prior_params <- list(
    alpha = alpha,
    betas = betas,
    gamma = gamma,
    sigma2 = sigma2
  )

  for (i in seq_len(times)) {
    params <- mle_params(
      y = y,
      x = x,
      alpha = prior_params$alpha,
      betas = prior_params$betas,
      gamma = prior_params$gamma,
      sigma2 = prior_params$sigma2
    )

    if (i %% trace_mod == 0) {
      cat(
        sprintf(
          "step: %i, ll: %0.05f, alpha: %0.05f,  %s, gamma: %0.3f, sigma2: %0.3f\n",
          i,
          params$ll,
          params$alpha,
          paste(
            names(params$betas),
            round(params$betas, 3),
            sep = ": ",
            collapse = ", "
          ),
          params$gamma,
          params$sigma2
        )
      )
    }
    if (identical(params, prior_params)) {
      warning("Step was identical to last")
      return(params)
    }

    prior_params <- params
  }

  params
}
