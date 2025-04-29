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

#' Fit PAR Model Using BFGS Optimization
#'
#' Fits a PAR model using the BFGS optimization algorithm.
#'
#' @param spec Model specification (output from `make_par_model_spec()`)
#' @param global_tol Tolerance for convergence of the global optimization
#' @param verbose Verbosity level (0 = silent, 1 = progress, 2 = detailed)
#' @param maxIter Maximum number of iterations for the optimization
#'
#' @return A list containing the fitted model parameters
#' @export
fit_par_bfgs <- function(
  spec,
  global_tol = 0.1,
  verbose = FALSE,
  maxIter = 1000
) {
  assert_valid_par_model_spec(spec)
  res <- bfgs_cpp(
    spec$Y,
    spec$X,
    spec$beta,
    spec$gamma,
    verbose = verbose,
    maxIter = maxIter
  )
  ll <- res$objective
  eps <- Inf
  # may be not important now that I've fixed the bugs
  # in bfgs_cpp...
  while (eps > global_tol) {
    if (verbose > 2) {
      cat("Inverse Hessian Reset. Current eps:", eps, "\n")
    }
    res <- bfgs_cpp(
      spec$Y,
      spec$X,
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
