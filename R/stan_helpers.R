#' Prepare Stan Data for Vector PAR Model
#'
#' @param data Long-format tibble from simulate_vector_par()
#' @param q Number of lags
#'
#' @return A list ready to pass to Stan
#' @export
prepare_data_stan_vectorpar <- function(data, q = 1) {
  stopifnot(all(c("item_id", "brand_id", "time", "y", "promo") %in% names(data)))
  
  item_ids <- sort(unique(data$item_id))
  brand_ids <- sort(unique(data$brand_id))
  p <- 1 # promo only for now
  q_seq <- seq_len(q)
  data <- data %>%
    group_by(item_id) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(
      !!!lapply(
        stats::setNames(q_seq, paste0("lag", q_seq)),
        function(x) rlang::expr(lag(y, !!x))
      )
    ) %>%
    ungroup()

  # Drop rows with NA lags
  data <- data %>%
    dplyr::filter(if_all(starts_with("lag"), ~!is.na(.)))

  N <- nrow(data)
  I <- length(item_ids)
  B <- length(brand_ids)
  brand <- data %>%
    distinct(item_id, brand_id) %>%
    arrange(item_id) %>%
    pull(brand_id)
  T <- max(data$time)
  Q <- q
  P <- p
  
  list(
    N = N,
    I = I,
    B = B,
    Q = Q,
    P = P,
    T = max(data$time),
    y = data$y,
    item = data$item_id,
    brand = brand,
    x = matrix(data$promo, ncol = P),
    lag_y = as.matrix(data %>% select(starts_with("lag"))),
    time = data$time
  )
}