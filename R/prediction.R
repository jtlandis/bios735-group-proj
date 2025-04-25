#' Bayesian PAR/PVAR model training
#'
#' Fits the Bayesian PAR (single item at a time) and/or PVAR (multi-item/hierarchical) model(s) on training set, which is defined as first timepoints up to a specified cutoff
#'
#' @param data Preprocessed list of multi-item time series matrices, obtained from 'prepare_data_allitems' function on raw dataset
#' @param ptrain Proportion between 0-1 specifying the % of timepoints to be used for model training. Default of 0.8
#' @param model_type String: Specify 'par' to fit PAR model, 'pvar' to fit PVAR model, or 'both' to fit both. If else, fits PAR model for each item. Default 'par'
#' @param q Integer lag for autoregressive terms. Default of 1
#' @param mcmc_iter Integer number of MCMC iterations for PAR/PVAR model fit. Default of 2000
#' @param verbose Set TRUE to return progress of model fitting/MCMC. Default TRUE
#'
#' @return List containing model fits and MCMC samples for specified model(s)
#' @export
train_bayes_pvar <- function(data, ptrain = 0.8, model_type = 'par', q = 1, mcmc_iter = 2000, verbose = TRUE) {
  T <- nrow(data$Y)
  train_idx <- 1:round(T * ptrain)
  Ytrain <- data$Y[train_idx, ]
  Xtrain <- data$X[train_idx, ]
  G <- data$G
  
  if (!(model_type %in% c('par','pvar','both'))) {
    print("Warning: Must specify model_type as 'par', 'pvar', or 'both'")
    print("Defaulting to par model")
    model_type <- 'par'
  }
  
  # run MCMC for PAR model if specified
  if (model_type %in% c('par', 'both')) {
    print('Running MCMC for PAR model on training set')
    par_mcmc <- fit_par_mcmc_allitems(Ytrain, Xtrain, q, mcmc_iter = mcmc_iter, verbose = verbose)
  }
  
  # run MCMC for PVAR model if specified
  if (model_type %in% c('pvar', 'both')) {
    print('Running MCMC for PVAR model on training set')
    pvar_mcmc <- fit_pvar_mcmc(Ytrain, Xtrain, G, q, mcmc_iter = mcmc_iter, verbose = verbose)
  }
  
  # save results
  if (model_type == 'par') {
    res <- list(PAR = par_mcmc, train_idx = train_idx)
  } else if (model_type == 'pvar') {
    res <- list(PVAR = pvar_mcmc, train_idx = train_idx)
  } else if (model_type == 'both') {
    res <- list(PAR = par_mcmc, PVAR = pvar_mcmc, train_idx = train_idx)
  } else {
    res <- NULL
  }
  
  return(res)
}

#' Compute Fitted Conditional Means from PAR/PVAR Model
#'
#' Given outcome, covariates, and parameter matrices (Beta, Gamma), compute the conditional means \eqn{m_t} for each time point and item.
#'
#' @param Y Matrix of observed outcomes (T × N), where T = # of timepoints and N = # of items
#' @param X Matrix of covariates (T × N), matching the layout of Y (e.g. PROMO variable for each item)
#' @param Beta Matrix of autoregressive weights (N × q), where q is the number of lags
#' @param Gamma Matrix of covariate effects (N × 2), for intercept + promo
#' @param q Number of lags
#'
#' @return Matrix of estimated conditional means \eqn{\hat{Y}_t = m_t} (T × N)
#' @export
get_fit_pvar <- function(Y, X, Beta, Gamma, q = ncol(Beta)) {
  T <- nrow(Y)
  N <- ncol(Y)
  Mt <- matrix(NA, T, N)
  
  for (i in 1:N) {
    for (t in (q + 1):T) {
      y_lags <- Y[(t - 1):(t - q), i]
      x_t <- X[t, i]
      auto_part <- sum(Beta[i, ] * y_lags)
      linear_part <- Gamma[i, 1] + Gamma[i, 2] * x_t
      mix_weight <- 1 - sum(Beta[i, ])
      Mt[t, i] <- auto_part + mix_weight * exp(linear_part)
    }
  }
  return(Mt)
}

#' 1-Step Forecast Using Posterior-based Prediction
#'
#' Given outcome, covariates, and parameter matrices (Beta, Gamma), compute the conditional means \eqn{m_t} for each time point and item.
#'
#' @param Y Matrix of outcomes (T × N), where T = # of timepoints and N = # of items
#' @param X Matrix of covariates (T × N), matching the layout of Y (e.g. PROMO variable for each item)
#' @param mcmc_samples List of MCMC samples from PAR/PVAR model for all items. Can be returned from fit_par_mcmc_all_items / fit_pvar_mcmc / train_bayes_model
#' @param model_type Specify 'par' or 'pvar'
#' @param burn_in MCMC burn-in proportion. Default 0.5
#' @param verbose Set TRUE to output progress. Default FALSE
#'
#' @return Matrix of predicted conditional means (T × N), based on posterior predictive distribution from MCMC
#' @export
onestep_forecast_mcmc <- function(Y, X, mcmc_samples, model_type, burn_in = 0.5, verbose = TRUE) {
  n_items <- ncol(Y)
  Mt_1step_est <- 0*Y
  
  # par
  if (model_type == "par") {
    mcmc_iter <- length(mcmc_samples[[1]]$tau)
    ss_idx <- round(burn_in * mcmc_iter):mcmc_iter
    n_samps <- length(ss_idx)
    Beta_samps <- simplify2array(lapply(mcmc_samples, function(x) x$beta[ss_idx,,drop=F]))
    Gamma_samps <- simplify2array(lapply(mcmc_samples, function(x) x$gamma[ss_idx,,drop=F]))

    for (s in 1:n_samps) {
      if (verbose == T & s %% 100 == 0) print(paste0('Iteration = ', s, ' / ', n_samps))
      Beta_s <- Beta_samps[s,,]
      if (is.null(dim(Beta_s))) {
        Beta_s <- as.matrix(Beta_s)
      } else {
        Beta_s <- t(Beta_s)
      }
      Gamma_s <- t(Gamma_samps[s,,])
      Mt_s <- get_fit_pvar_cpp(Y, X, Beta_s, Gamma_s)
      Mt_1step_est <- Mt_1step_est + (1/n_samps)*Mt_s
    }
  }
  
  if (model_type == "pvar") {
    mcmc_iter <- length(mcmc_samples$tau_samples)
    ss_idx <- round(burn_in * mcmc_iter):mcmc_iter
    n_samps <- length(ss_idx)
    Beta_samps <- simplify2array(mcmc_samples$Beta_samples)[,,ss_idx,drop=F]
    Gamma_samps <- simplify2array(mcmc_samples$Gamma_samples)[,,ss_idx,drop=F]
    
    for (s in 1:n_samps) {
      if (verbose == T & s %% 100 == 0) print(paste0('Iteration = ', s, ' / ', n_samps))
      Beta_s <- Beta_samps[,,s]
      if (is.null(dim(Beta_s))) {
        Beta_s <- as.matrix(Beta_s)
      } else {
        Beta_s <- Beta_s
      }
      Gamma_s <- Gamma_samps[,,s]
      Mt_s <- get_fit_pvar_cpp(Y, X, Beta_s, Gamma_s)
      Mt_1step_est <- Mt_1step_est + (1/n_samps)*Mt_s
    }
  }
  
  return(Mt_1step_est)
}






