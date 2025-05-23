% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_all_items.R
\name{summarize_par_fits}
\alias{summarize_par_fits}
\title{Summarize PAR Fit Results Across Items}
\usage{
summarize_par_fits(fits)
}
\arguments{
\item{fits}{List of model fits (output from \code{fit_all_items()})}
}
\value{
A tibble with one row per item
}
\description{
Extracts beta, gamma, loglik, and convergence status from a list of fitted models.
}
