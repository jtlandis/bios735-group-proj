% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{onestep_forecast_mcmc}
\alias{onestep_forecast_mcmc}
\title{1-Step Forecast Using Posterior-based Prediction}
\usage{
onestep_forecast_mcmc(
  Y,
  X,
  mcmc_samples,
  model_type,
  burn_in = 0.5,
  verbose = TRUE
)
}
\arguments{
\item{Y}{Matrix of outcomes (T × N), where T = # of timepoints and N = # of items}

\item{X}{Matrix of covariates (T × N), matching the layout of Y (e.g. PROMO variable for each item)}

\item{mcmc_samples}{List of MCMC samples from PAR/PVAR model for all items. Can be returned from fit_par_mcmc_all_items / fit_pvar_mcmc / train_bayes_model}

\item{model_type}{Specify 'par' or 'pvar'}

\item{burn_in}{MCMC burn-in proportion. Default 0.5}

\item{verbose}{Set TRUE to output progress. Default FALSE}
}
\value{
Matrix of predicted conditional means (T × N), based on posterior predictive distribution from MCMC
}
\description{
Given outcome, covariates, and parameter matrices (Beta, Gamma), compute the conditional means \eqn{m_t} for each time point and item.
}
