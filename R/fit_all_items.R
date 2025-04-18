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