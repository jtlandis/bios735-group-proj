#' Prepare Raw Dataset for Time Series Modeling
#'
#' Adds item/brand indices, arranges by time, and computes lags needed for PAR model.
#'
#' @param data A raw dataframe with columns DATE and multiple columns starting with QTY and PROMO. 
#' @return A list containing time series across all items for sales outcome (Y) and promotions covariate (X), group/brand indices (G), items, and number of items and groups
#' @export
prepare_data_allitems <- function(data) {
  
  # get sales, promotions, and brand for all items
  varnames <- names(data)
  sales <- data[varnames[startsWith(varnames, "QTY")]]
  promos <- data[varnames[startsWith(varnames, "PROMO")]]
  Y <- as.matrix(sales)
  X <- as.matrix(promos)
  item_codes <- substring(names(sales), 5)
  G <- as.numeric(substr(item_codes, 2, 2))
  B <- length(unique(G))
  return(list(Y = Y, X = X, 
              G = G, num_brands = B, 
              items = item_codes, num_items = length(G)))
}


#' Prepare and Lag Tidy Dataset for PAR Modeling
#'
#' Adds item/brand indices, arranges by time, and computes lags needed for PAR model.
#'
#' @param data A tidy dataframe with columns DATE, QTY, PROMO, item, brand
#' @param q Number of lags to compute (default = 1)
#' @param covariates A character vector of covariate column names to include in X
#'
#' @return A long-format data.frame with lagged values, numeric item/brand IDs, and covariate matrix
#' @export
prepare_item_data <- function(data, q = 1, covariates = "PROMO") {
  stopifnot("DATE" %in% names(data), "QTY" %in% names(data), "item" %in% names(data))
  
  data_proc <- data %>%
    mutate(
      item_id = as.integer(as.factor(item)),
      brand_id = as.integer(as.factor(brand)),
      DATE = as.Date(DATE)
    ) %>%
    arrange(item_id, DATE)
  
  # Add lag columns (lag1, lag2, ..., lagq)
  for (l in 1:q) {
    data_proc <- data_proc %>%
      group_by(item_id) %>%
      mutate(!!paste0("lag", l) := lag(QTY, l)) %>%
      ungroup()
  }
  
  # Drop rows missing any lags
  data_proc <- data_proc %>%
    filter(if_all(starts_with("lag"), ~ !is.na(.)))
  
  # Add covariate matrix X
  data_proc <- data_proc %>%
    mutate(across(all_of(covariates), as.numeric))  # ensure numeric
  
  return(data_proc)
}

#' Split Preprocessed Data into Per-Item Lists
#'
#' After using `prepare_item_data()`, split into list of per-item data with y and X matrix.
#'
#' @param data A preprocessed dataframe (output of prepare_item_data)
#' @param covariates Character vector of covariate column names
#'
#' @return A list of item-specific lists each containing y and x
#' @export
split_by_item <- function(data, covariates = "PROMO") {
  data %>%
    group_by(item_id) %>%
    group_split() %>%
    purrr::map(~ list(
      y = .x$QTY,
      x = as.matrix(.x[, covariates, drop = FALSE])
    ))
}

