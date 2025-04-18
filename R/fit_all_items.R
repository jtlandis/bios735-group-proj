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
    tryCatch({
      fit <- fit_par_optim(
        y = item$y,
        x = item$x,
        q = q
      )
      fit$item_id <- i
      return(fit)
    }, error = function(e) {
      message("Model failed for item ", i, ": ", e$message)
      return(NULL)
    })
  }) %>% purrr::compact()
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
    unnest_wider(beta, names_sep = "_b") %>%
    unnest_wider(gamma, names_sep = "_g")
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
fit_vector_par_stan <- function(stan_data,
                                stan_file = "inst/stan/vector_par.stan",
                                iter = 1000,
                                chains = 4,
                                seed = 123) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is required. Install it with: install.packages('cmdstanr')")
  }
  
  mod <- cmdstanr::cmdstan_model(stan_file)
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
