% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_par.R
\name{simulate_par_intercept}
\alias{simulate_par_intercept}
\title{Simulate Time Series from PAR Model with Item Intercept}
\usage{
simulate_par_intercept(
  T = 100,
  q = 1,
  gamma = 1,
  beta_lags = 0.5,
  eta = 0,
  promo_prob = 0.2,
  seed = NULL
)
}
\arguments{
\item{T}{Number of time points}

\item{q}{Number of AR lags}

\item{gamma}{Covariate coefficient vector}

\item{beta_lags}{AR coefficients (length q, should sum <= 1)}

\item{eta}{Item-specific intercept}

\item{promo_prob}{Probability of promotion at each time}

\item{seed}{Optional random seed}
}
\value{
A data.frame with y, time, promo
}
\description{
Simulate Time Series from PAR Model with Item Intercept
}
