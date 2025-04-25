
#include <RcppEigen.h>
#include <Rmath.h>
#include "MLE_utils.h"
using namespace Rcpp;

Rcpp::NumericVector grad_direc(
  const Eigen::MatrixXd& H,
  Rcpp::NumericVector& grad
) {
  int n = grad.size();
  Rcpp::NumericVector direc(n);
  for (int i = 0; i < n; ++i) {
    double dot = 0.0;
    for (int j = 0; j < n; ++j) {
      dot += H(i, j) * grad[j];
    }
    direc[i] = -dot;
  }
  return direc;
}

// x * t(y)
Eigen::MatrixXd simple_dot(const NumericVector& x, const NumericVector& y, double div = 1.0) {
  int n = x.size();
  Eigen::MatrixXd result(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      result(i, j) = x[i] * y[j] / div;
    }
  }

  return result;
}

Eigen::MatrixXd update_H(
  Eigen::MatrixXd& H,
  const Rcpp::NumericVector& step,
  const Rcpp::NumericVector& diff) {
  int n = step.size();
  double yts = 0;
  for (int i = 0; i < n; i++) {
    yts += diff[i] * step[i];
  }

  Eigen::MatrixXd left = -simple_dot(step, diff, yts);
  Eigen::MatrixXd right = -simple_dot(diff, step, yts);
  Eigen::MatrixXd sst = simple_dot(step, step, yts);

  // simulate I - M
  for (int i = 0; i < n; i++) {
    left(i, i) += 1;
    right(i, i) += 1;
  }

  Eigen::MatrixXd result = left * H * right;

  return result + sst;

}

double armijo_line_search(const Rcpp::NumericVector& Y,
                          const Rcpp::NumericMatrix& X,
                          const Rcpp::NumericVector& beta,
                          const Rcpp::NumericVector& gamma,
                          const Rcpp::NumericVector& direc,
                          IntegerVector& beta_slice,
                          IntegerVector& gamma_slice,
                          double& f,
                          double alpha = 1,
                          double c = 1e-4,
                          double tau = 0.5) {
  // double f_old = -loglik_cpp(Y, X, beta, gamma);
  Rcpp::NumericVector grad = -loglik_grad_cpp(Y, X, beta, gamma);
  double grad_dot_direc = sum(grad * direc);
  Rcpp::NumericVector beta_new = clone(beta);
  Rcpp::NumericVector gamma_new = clone(gamma);
  Rcpp::NumericVector grand_beta = direc[beta_slice];
  Rcpp::NumericVector grand_gamma = direc[gamma_slice];
  while (true) {
    beta_new = beta + alpha * grand_beta;
    gamma_new = gamma + alpha * grand_gamma;
    double f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
    // Rcpp::NumericVector grad_new = -loglik_grad_cpp(Y, X, beta_new, gamma_new);
    // double grad_new_dot_direc = Rcpp::sum(grad_new * direc);

    if (f_new <= f + c * alpha * grad_dot_direc) {
      // if(std::abs(grad_new_dot_direc) < std::abs(c2 * grad_dot_direc)) {
        // both conditions are met
        f = f_new;
        break;
    //   }
    //   alpha /=  0.8;
    // } else {

    }
    alpha *= tau;

  }

  return alpha;
}

// [[Rcpp::export]]
Rcpp::List bfgs_cpp(
  const Rcpp::NumericVector& Y,
  const Rcpp::NumericMatrix& X,
  Rcpp::NumericVector beta0,
  Rcpp::NumericVector gamma0,
  int maxIter = 100,
  double tol = 1e-5,
  int verbose = 0,
  int iter = 0
) {

  NumericVector beta = clone(beta0);
  NumericVector gamma = clone(gamma0);
  int q = beta0.size();
  int p = gamma0.size();
  double objective = -loglik_cpp(Y, X, beta, gamma);
  double old_objective = objective;
  double ep = INFINITY;
  // int alpha_index = 0;
  IntegerVector beta_slice(q);
  NumericVector beta_proxy(q);
  for (int i = 0; i < q; i++) {
    beta_slice[i] = i;
  }
  IntegerVector gamma_slice(p);
  NumericVector gamma_proxy(q);
  for (int i = 0; i < p; i++) {
    gamma_slice[i] = i + q;
  }
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(p + q, p + q);
  // NumericMatrix H = NumericMatrix::diag(p + q, 1.0);
  Rcpp::NumericVector grad = -loglik_grad_cpp(Y, X, beta, gamma);
  double step = 0;
  while (ep > tol && iter <= maxIter) {

    Rcpp::NumericVector direc = grad_direc(H, grad);
    // step = line_search_cpp3(Y, X, beta, gamma, direc, alpha_index, objective, beta_slice, gamma_slice);
    step = armijo_line_search(Y, X, beta, gamma, direc, beta_slice, gamma_slice, objective);
    // if (verbose) Rcout << " direction: " << direc << "\n";
    direc = direc * step;
    beta_proxy = direc[beta_slice];
    beta += beta_proxy;
    gamma_proxy = direc[gamma_slice];
    gamma += gamma_proxy;
    Rcpp::NumericVector grad_new = -loglik_grad_cpp(Y, X, beta, gamma);
    Rcpp::NumericVector grad_diff = grad_new - grad;
    H = update_H(H, direc, grad_diff);
    //double f_new = -loglik_cpp(Y, X, beta, gamma);
    ep = std::abs(objective - old_objective);
    old_objective = objective;
    grad = grad_new;
    if (verbose) {
      Rcout << "Iteration: " << iter << ", Step size: " << step << ", eps: " << ep << ", ll: " << objective;
      if (verbose > 1) {
        Rcout << ", beta: " << beta << ", gamma: " << gamma;
      }
      Rcout << "\n";
    }

    iter++;
  }


  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("iter") = iter,
    Rcpp::Named("objective") = -objective
  );
}









// [[Rcpp::export]]
Rcpp::List proj_grad_descent_cpp(const Rcpp::NumericVector& Y,
                                 const Rcpp::NumericMatrix& X,
                                 Rcpp::NumericVector beta0,
                                 Rcpp::NumericVector gamma0,
                                 double lr = 1e-4,
                                 int maxIter = 100,
                                 double tol = 1e-5,
                                 bool return_allIters = false,
                                 bool verbose = false) {

  NumericVector beta = clone(beta0);
  NumericVector gamma = clone(gamma0);
  int q = beta0.size();
  int p = gamma0.size();
  double f = -loglik_cpp(Y, X, beta, gamma);
  double ep = INFINITY;
  int iter = 1;

  Rcpp::List beta_iters, gamma_iters;
  if (return_allIters) {
    beta_iters = Rcpp::List(maxIter);
    gamma_iters = Rcpp::List(maxIter);
  }

  while (ep > tol && iter <= maxIter) {
    if (verbose) Rcpp::Rcout << "Iteration: " << iter << "\n";

    Rcpp::NumericVector grad = -loglik_grad_cpp(Y, X, beta, gamma);
    Rcpp::NumericVector beta_grad = grad[Rcpp::Range(0, q - 1)];
    Rcpp::NumericVector gamma_grad = grad[Rcpp::Range(q, q + p - 1)];

    for (int i = 0; i < p; ++i) {
      gamma[i] -= lr * gamma_grad[i];
    }

    for (int i = 0; i < q; ++i) {
      beta[i] -= lr * beta_grad[i];
    }

    beta = proj_beta_cpp(beta);

    double fnew = -loglik_cpp(Y, X, beta, gamma);
    ep = std::abs(fnew - f);
    f = fnew;

    if (return_allIters) {
      beta_iters[iter - 1] = Rcpp::clone(beta);
      gamma_iters[iter - 1] = Rcpp::clone(gamma);
    }

    if (verbose) Rcpp::Rcout << "ep = " << ep << "\n";
    iter++;
  }

  if (verbose) Rcpp::Rcout << "Stopped at iteration " << iter << "\n";

  if (return_allIters) {
    return Rcpp::List::create(
      Rcpp::Named("beta") = beta,
      Rcpp::Named("gamma") = gamma,
      Rcpp::Named("beta_iters") = beta_iters,
      Rcpp::Named("gamma_iters") = gamma_iters,
      Rcpp::Named("epsilon") = ep
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("beta") = beta,
      Rcpp::Named("gamma") = gamma,
      Rcpp::Named("epsilon") = ep
    );
  }
}
