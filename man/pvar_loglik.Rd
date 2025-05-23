% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvar_model.R
\name{pvar_loglik}
\alias{pvar_loglik}
\title{Poisson Vector Autoregressive (PVAR) Log-Likelihood}
\usage{
pvar_loglik(Y, X, Beta, Gamma, q = ncol(Beta))
}
\arguments{
\item{Y}{Matrix of counts (T x N)}

\item{X}{Matrix of single covariate (T x N)}

\item{Beta}{Matrix of autoregressive weights (N x q)}

\item{Gamma}{Matrix of intercepts and covariate effect (N x 2)}

\item{q}{Number of lags}
}
\value{
Log-likelihood value (numeric scalar)
}
\description{
Computes the log-likelihood for a Poisson vector autoregressive model for N items
}
