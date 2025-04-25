// MLE_utils.h
#ifndef MLE_UTILS_H
#define MLE_UTILS_H

#include <RcppEigen.h>
#include <Rmath.h>

Rcpp::NumericVector get_cov_parts(
  const Rcpp::NumericMatrix& X,
  const Rcpp::NumericVector& gamma,
  // const List& list,
  int q
);
Rcpp::NumericVector get_mt_cpp(const Rcpp::NumericVector& Y, const Rcpp::NumericMatrix& X,
                               const Rcpp::NumericVector& beta, const Rcpp::NumericVector& gamma);
Rcpp::NumericMatrix get_mt_grad_cpp(const Rcpp::NumericVector& Y, const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericVector& beta, const Rcpp::NumericVector& gamma);
double loglik_cpp(const Rcpp::NumericVector& Y, const Rcpp::NumericMatrix& X,
                  const Rcpp::NumericVector& beta, const Rcpp::NumericVector& gamma);
Rcpp::NumericVector loglik_grad_cpp(const Rcpp::NumericVector& Y, const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericVector& beta, const Rcpp::NumericVector& gamma);
Rcpp::NumericVector proj_beta_cpp(Rcpp::NumericVector beta, double epsilon = 1e-4);

#endif
