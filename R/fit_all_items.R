#' Fit PAR Model to All Items in a Dataset
#'
#' Loops through a list of item-specific data and fits the Poisson autoregressive model.
#'
#' @param item_list A list of data objects (each with $y and $x)
#' @param q Number of lags to use in the model
#'
#' @return A list of model fits (each a list of $beta, $gamma, $loglik, etc.)
#' @export
fit_all_items <- function(item_list, q = 1) {
  purrr::map(seq_along(item_list), function(i) {
    item <- item_list[[i]]
    tryCatch(
      {
        fit <- fit_par_optim(
          y = item$y,
          x = item$x,
          q = q
        )
        fit$item_id <- i
        return(fit)
      },
      error = function(e) {
        message("Model failed for item ", i, ": ", e$message)
        return(NULL)
      }
    )
  }) %>%
    purrr::compact()
}

#' Summarize PAR Fit Results Across Items
#'
#' Extracts beta, gamma, loglik, and convergence status from a list of fitted models.
#'
#' @param fits List of model fits (output from `fit_all_items()`)
#'
#' @return A tibble with one row per item
#' @export
summarize_par_fits <- function(fits) {
  purrr::map_dfr(seq_along(fits), function(i) {
    fit <- fits[[i]]
    tibble(
      item_id = i,
      loglik = fit$loglik,
      converged = fit$converged,
      beta = list(fit$beta),
      gamma = list(fit$gamma)
    )
  }) %>%
    tidyr::unnest_wider(beta, names_sep = "_b") %>%
    tidyr::unnest_wider(gamma, names_sep = "_g")
}

#' Fit the Vector PAR model in Stan
#'
#' @param stan_data A list from prepare_data_stan_vectorpar()
#' @param stan_file Path to the Stan file
#' @param iter Number of iterations
#' @param chains Number of chains
#' @param seed Random seed
#'
#' @return A cmdstanr fit object
#' @export
fit_vector_par_stan <- function(
  stan_data,
  stan_file = "inst/stan/vector_par.stan",
  iter = 1000,
  chains = 4,
  seed = 123
) {
  #if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  #  stop("cmdstanr is required. Install it with: install.packages('cmdstanr')")
  #}

  mod <- cmdstanr::cmdstan_model(
    system.file("stan/vector_par.stan", package = "pastasales")
  )
  #mod <- stanmodels$vector_par
  fit <- mod$sample(
    data = stan_data,
    iter_sampling = iter,
    iter_warmup = iter,
    chains = chains,
    seed = seed,
    refresh = 500
  )

  return(fit)
}

#' Prepare and Fit PAR Model for a Single Item using Stan
#'
#' @param df A data.frame with columns `y`, `time`, and covariates (e.g., `promo`)
#' @param q Number of AR lags
#' @param stan_file Path to Stan file (default = "inst/stan/par_item.stan")
#' @param iter Number of iterations (default = 1000)
#' @param chains Number of MCMC chains (default = 4)
#' @param seed Random seed
#'
#' @return A cmdstanr fit object
#' @export
fit_par_item <- function(
  df,
  q = 1,
  stan_file = "inst/stan/par_item.stan",
  iter = 1000,
  chains = 4,
  seed = 123
) {
  stopifnot("y" %in% names(df), "time" %in% names(df))

  # Add AR lags
  for (l in 1:q) {
    df[[paste0("lag", l)]] <- dplyr::lag(df$y, l)
  }
  df <- df %>%
    dplyr::filter(dplyr::if_all(dplyr::starts_with("lag"), ~!is.na(.)))

  y <- df$y
  X <- as.matrix(df |> dplyr::select(dplyr::starts_with("promo")))

  T_obs <- length(y)
  p <- ncol(X)

  stan_data <- list(
    T_obs = T_obs,
    q = q,
    p = p,
    y = y,
    X = X,
    mu_gamma = rep(0, p),
    Sigma_gamma = diag(p),
    alpha = rep(1 / q, q),
    a_tau = 1,
    b_tau = 1
  )

  #mod <- stanmodels$par_item
  mod <- cmdstanr::cmdstan_model(
    system.file("stan/par_item.stan", package = "pastasales")
  )
  fit <- mod$sample(
    data = stan_data,
    iter_sampling = iter,
    iter_warmup = iter,
    chains = chains,
    seed = seed,
    refresh = 500
  )
  return(fit)
}

#' Prepare and Fit PAR Model for a Single Item using Stan (with Item Intercept)
#'
#' @param df A data.frame with columns `y`, `time`, and covariates (e.g., `promo`)
#' @param q Number of AR lags
#' @param stan_file Path to Stan file (default = "inst/stan/fit_par_intercept.stan")
#' @param iter Number of iterations (default = 1000)
#' @param chains Number of MCMC chains (default = 4)
#' @param seed Random seed
#'
#' @return A cmdstanr fit object
#' @export
fit_par_item_int <- function(
  df,
  q = 1,
  stan_file = "inst/stan/par_item_intercept.stan",
  iter = 1000,
  chains = 4,
  seed = 123
) {
  stopifnot("y" %in% names(df), "time" %in% names(df))

  # Add AR lags
  for (l in 1:q) {
    df[[paste0("lag", l)]] <- dplyr::lag(df$y, l)
  }
  df <- df %>%
    dplyr::filter(dplyr::if_all(dplyr::starts_with("lag"), ~!is.na(.)))

  y <- df$y
  X <- as.matrix(df |> dplyr::select(dplyr::starts_with("promo")))

  T_obs <- length(y)
  p <- ncol(X)

  stan_data <- list(
    T_obs = T_obs,
    q = q,
    p = p,
    y = y,
    X = X,
    mu_gamma = rep(0, p),
    Sigma_gamma = diag(p),
    alpha = rep(1 / q, q),
    a_tau = 1,
    b_tau = 1
  )

  #mod <- stanmodels$par_item_intercept
  mod <- cmdstanr::cmdstan_model(
    system.file("stan/par_item_intercept.stan", package = "pastasales")
  )
  fit <- mod$sample(
    data = stan_data,
    iter_sampling = iter,
    iter_warmup = iter,
    chains = chains,
    seed = seed,
    refresh = 500
  )
  return(fit)
}

#' @export
fit_par_bfgs <- function(
  mod,
  global_tol = 0.1,
  verbose = FALSE,
  maxIter = 1000
) {
  res <- bfgs_cpp2(
    mod$Y,
    mod$X,
    mod$beta,
    mod$gamma,
    verbose = verbose,
    maxIter = maxIter
  )
  ll <- res$objective
  eps <- Inf
  # may be not important now that I've fixed the bugs
  # in bfgs_cpp2...
  while (eps > global_tol) {
    if (verbose > 2) {
      cat("Inverse Hessian Reset. Current eps:", eps, "\n")
    }
    res <- bfgs_cpp2(
      mod$Y,
      mod$X,
      res$beta,
      res$gamma,
      verbose = verbose,
      iter = res$iter,
      maxIter = maxIter
    )
    ll_new <- res$objective
    eps <- abs(ll_new - ll)
    ll <- ll_new
  }
  c(
    res[-4],
    list(
      ll = ll,
      eps = eps
    )
  )
}
