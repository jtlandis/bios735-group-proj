#' Simulate a Poisson Autoregressive (PAR) Time Series
#'
#' @param T Number of timepoints
#' @param beta Vector of lag coefficients (length q)
#' @param gamma Vector of covariate coefficients
#' @param x Matrix of covariates of dimension T × p (optional; generated if NULL)
#' @param mu Scalar intercept
#' @param seed Random seed
#'
#' @return A list with y (counts), x (covariates), and the true parameters
#' @export
simulate_par <- function(T = 200, beta = c(0.2), gamma = c(1.5), x = NULL, mu = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  q <- length(beta)
  p <- length(gamma)
  
  if (is.null(x)) {
    x <- matrix(rbinom(T * p, size = 1, prob = 0.2), nrow = T, ncol = p)
  }
  
  y <- rep(0, T)
  y[1:q] <- rpois(q, lambda = 5)  # initialize with arbitrary counts
  
  for (t in (q + 1):T) {
    y_lags <- y[(t - 1):(t - q)]
    linear_part <- mu + sum(x[t, ] * gamma)
    mix_weight <- 1 - sum(beta)
    m_t <- sum(beta * y_lags) + mix_weight * exp(linear_part)
    y[t] <- rpois(1, lambda = m_t)
  }
  
  list(y = y, x = x, beta = beta, gamma = gamma, mu = mu)
}

#' Simulate Hierarchical Vector Poisson Autoregressive (PAR) Time Series
#'
#' @param n_items Number of items
#' @param n_brands Number of brands
#' @param T Number of timepoints
#' @param q Number of lags
#' @param p Number of covariates (default = 1)
#' @param promo_prob Probability of promotion being on
#' @param sigma_mu SD for brand-level mean of gamma
#' @param rho_g SD for item-level gamma noise
#' @param sigma_eta SD for item-level intercepts
#' @param seed Random seed
#'
#' @return A tibble with simulated long-format data
#' @export
simulate_vector_par <- function(n_items = 30, n_brands = 3, T = 100, q = 1, p = 1,
                                promo_prob = 0.2,
                                sigma_mu = 1, rho_g = 0.5, sigma_eta = 0.5,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  item_ids <- 1:n_items
  brand_ids <- sample(1:n_brands, n_items, replace = TRUE)
  
  mu_g_list <- lapply(1:n_brands, function(b) rnorm(p, mean = 0, sd = sigma_mu))
  gamma_list <- lapply(1:n_items, function(i) {
    brand <- brand_ids[i]
    rnorm(p, mean = mu_g_list[[brand]], sd = rho_g)
  })
  
  beta_list <- lapply(1:n_items, function(i) {
    beta_raw <- rgamma(q, shape = 2)
    beta <- beta_raw / (1.5 * sum(beta_raw))
    return(beta)
  })
  
  mu_brand <- rnorm(n_brands, 0, 1)
  eta_item <- rnorm(n_items, 0, sigma_eta)
  f_time <- rnorm(T, 0, 1)
  
  covariate_array <- array(NA, dim = c(n_items, T, p))
  for (i in 1:n_items) {
    for (j in 1:p) {
      covariate_array[i, , j] <- rbinom(T, size = 1, prob = promo_prob)
    }
  }
  
  y_array <- matrix(0, n_items, T)
  y_array[, 1:q] <- matrix(rpois(n_items * q, lambda = 5), n_items, q)
  
  for (t in (q + 1):T) {
    for (i in 1:n_items) {
      beta <- beta_list[[i]]
      gamma <- gamma_list[[i]]
      x_t <- covariate_array[i, t, ]
      
      ar_term <- sum(beta * y_array[i, (t - 1):(t - q)])
      lin_pred <- mu_brand[brand_ids[i]] + eta_item[i] + f_time[t] + sum(gamma * x_t)
      mix_weight <- 1 - sum(beta)
      
      mu_t <- ar_term + mix_weight * exp(lin_pred)
      y_array[i, t] <- rpois(1, lambda = mu_t)
    }
  }
  
  long_df <- purrr::map_dfr(1:n_items, function(i) {
    tibble(
      item_id = i,
      brand_id = brand_ids[i],
      time = 1:T,
      y = y_array[i, ],
      promo = covariate_array[i, , 1]
    )
  })
  
  attr(long_df, "true_params") <- list(
    gamma = gamma_list,
    beta = beta_list,
    brand_ids = brand_ids,
    mu_brand = mu_brand,
    eta = eta_item,
    f = f_time
  )
  
  return(long_df)
}