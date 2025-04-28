#' Poisson Vector Autoregressive (PVAR) Log-Likelihood
#'
#' Computes the log-likelihood for a Poisson vector autoregressive model for N items
#'
#' @param Y Matrix of counts (T x N)
#' @param X Matrix of single covariate (T x N)
#' @param Beta Matrix of autoregressive weights (N x q)
#' @param Gamma Matrix of intercepts and covariate effect (N x 2)
#' @param q Number of lags
#'
#' @return Log-likelihood value (numeric scalar)
#' @export
pvar_loglik <- function(Y, X, Beta, Gamma, q = ncol(Beta)) {
  T <- nrow(Y)
  N <- ncol(Y)
  stopifnot(nrow(X) == T)
  stopifnot(ncol(X) == N)
  Mt <- matrix(nrow = T, ncol = N)
  ll <- 0
  
  # loop over all items
  for (i in 1:N) {
    
    # loop over timepoints after lag
    for (t in (q+1):T) {
      y_lags <- Y[(t - 1):(t - q), i]
      x_t <- X[t, i]
      auto_part <- sum(Beta[i, ] * y_lags)
      linear_part <- Gamma[i, 1] + Gamma[i, 2] * x_t
      mix_weight <- 1 - sum(Beta[i, ])
      m_t <- auto_part + mix_weight * exp(linear_part)
      if (m_t <= 0) m_t <- 1e-10 # safeguard
      ll <- ll + Y[t, i] * log(m_t) - m_t - lgamma(Y[t, i] + 1)
    }
  }
  
  return(ll)
}


#' Fit PVAR Model Using MCMC
#'
#' Fit a hierarchical multi-item Poisson Vector AR model using MCMC (Metropolis-within-Gibbs sampler)
#'
#' @param Y Matrix of counts (T x N)
#' @param X Matrix of covariates (T × N)
#' @param G Length N vector indexing group with values from 1 to total number of groups
#' @param q Number of lags
#' @param mcmc_iter Number of iterations for MCMC. Default 2000
#' @param burn_in Proportion of burn-in samples for computing parameter estimates (default 0.5)
#' @param hyperparams Default NULL. Can be specified as named list of hyperparameters
#' @param proposal_sd Default 0.05. Standard deviation for Metropolis-Hasting proposal density
#' @param return_mcmc Default TRUE. Set FALSE to only return point estimates without MCMC samples
#' @param verbose Default FALSE. Set TRUE to print iterations while running
#'
#' @return A list with estimated parameters, log-likelihood, and (optionally) MCMC samples
#' @export
fit_pvar_mcmc <- function(Y, X, G, q = 1, mcmc_iter = 2000, burn_in = 0.5, hyperparams = NULL, proposal_sd = 0.05, return_mcmc = TRUE, verbose = FALSE) {
  B <- length(unique(G))
  
  # run MCMC
  mcmc <- run_mcmc_pvar_cpp(Y, X, G, B = B, q = q, n_iter = mcmc_iter, hyperparams = hyperparams, proposal_sd = proposal_sd, verbose = verbose)

  # get parameter estimates
  ss_idxs <- round(burn_in * mcmc_iter):mcmc_iter
  Gamma_array <- simplify2array(mcmc$Gamma_samples)
  Beta_array <- simplify2array(mcmc$Beta_samples)
  Gamma_est <- apply(Gamma_array[, , ss_idxs,drop=F], c(1,2), mean)
  Beta_est <- apply(Beta_array[, , ss_idxs,drop=F], c(1,2), mean)
  
  # ll <- pvar_loglik(Y, X, Beta, Gamma) # slower
  ll <- loglik_pvar_cpp(Y, X, Beta_est, Gamma_est) - sum(lgamma(Y[(q+1):nrow(Y),] + 1)) # faster log-lik computation
  
  if (return_mcmc == T) res <- list(Beta = Beta_est, Gamma = Gamma_est, loglik = ll, mcmc_samps = mcmc)
  if (return_mcmc == F) res <- list(Beta = Beta_est, Gamma = Gamma_est, loglik = ll)
  
  res
}


#' Fit PAR Model for all items Using MCMC
#'
#' Fit a hierarchical multi-item Poisson Vector AR model using MCMC (Metropolis-within-Gibbs sampler)
#'
#' @param Y Matrix of counts (T x N)
#' @param X Matrix of covariates (T × N)
#' @param q Number of lags
#' @param mcmc_iter Number of iterations for MCMC. Default 2000
#' @param burn_in Proportion of burn-in samples for computing parameter estimates (default 0.5)
#' @param hyperparams Default NULL. Can be specified as named list of hyperparameters
#' @param proposal_sd Default 0.05. Standard deviation for Metropolis-Hasting proposal density
#' @param return_mcmc Default TRUE. Set FALSE to only return point estimates without MCMC samples
#' @param verbose Default FALSE. Set TRUE to print iterations while running
#'
#' @return A list with estimated parameters, log-likelihood, and (optionally) list of MCMC samples for each item
#' @export
fit_par_mcmc_allitems <- function(Y, X, q = 1, mcmc_iter = 2000, burn_in = 0.5, hyperparams = NULL, proposal_sd = 0.05, return_mcmc = TRUE, verbose = FALSE) {
  N <- ncol(Y)
  mcmcs <- vector(mode = "list", N)
  Gamma_est <- matrix(nrow = N, ncol = 2)
  Beta_est <- matrix(nrow = N, ncol = q)
  ss_idxs <- round(burn_in * mcmc_iter):mcmc_iter
  
  # loop over all items and do MCMC
  for (i in 1:N) {
    if (verbose) print(paste0('Item = ', i))
    mcmc_i <- run_mcmc_par_cpp(Y[,i], cbind(1, X[,i]), q, n_iter = mcmc_iter, hyperparams = hyperparams, proposal_sd = proposal_sd, verbose = F)
    Beta_est[i,] <- colMeans(mcmc_i$beta[ss_idxs,,drop=F])
    Gamma_est[i,] <- colMeans(mcmc_i$gamma[ss_idxs,,drop=F])
    mcmcs[[i]] <- mcmc_i
  }
  
  ll <- loglik_pvar_cpp(Y, X, Beta_est, Gamma_est) - sum(lgamma(Y[(q+1):nrow(Y),] + 1))

  if (return_mcmc == T) res <- list(Beta = Beta_est, Gamma = Gamma_est, loglik = ll, mcmc_samps = mcmcs)
  if (return_mcmc == F) res <- list(Beta = Beta_est, Gamma = Gamma_est, loglik = ll)
  
  res
}



