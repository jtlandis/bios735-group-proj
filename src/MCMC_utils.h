// MCMC_utils.h
#ifndef MCMC_UTILS_H
#define MCMC_UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

double dnorm_vec_cpp(NumericVector x, NumericVector mu, NumericVector sigsq, bool log = true);
NumericVector rdirichlet_cpp(NumericVector alpha_m);
double ddirichlet_cpp(NumericVector x, NumericVector alpha, bool log = true);
NumericVector dirichlet_to_eta(NumericVector beta_tilde);
NumericVector eta_to_dirichlet(NumericVector eta);
double log_jacobian_eta_to_dirichlet(NumericVector eta);

#endif