// MCMC_utils.h
#ifndef MCMC_UTILS_H
#define MCMC_UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

double dnorm_vec_cpp(NumericVector x, NumericVector mu, NumericVector sigsq, bool log = true);
NumericVector rdirichlet_cpp(NumericVector alpha_m);
double ddirichlet_cpp(NumericVector x, NumericVector alpha, bool log = true);

#endif