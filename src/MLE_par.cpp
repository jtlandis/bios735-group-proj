
#include <RcppEigen.h>
#include <Rmath.h>
#include "MLE_utils.h"
using namespace Rcpp;

Rcpp::NumericVector grad_direc(
  const Rcpp::NumericMatrix& H,
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

Rcpp::NumericVector grad_direc2(
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

  // assume x and y are the same size
  // skips some steps
NumericMatrix simple_dot(const NumericVector& x, const NumericVector& y, double div) {
  int n = x.size();
  NumericMatrix result(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      result(i, j) = - x[i] * y[j] / div;
    }
  }

  return result;
}

// x * t(y)
Eigen::MatrixXd simple_dot2(const NumericVector& x, const NumericVector& y, double div = 1.0) {
  int n = x.size();
  Eigen::MatrixXd result(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      result(i, j) = x[i] * y[j] / div;
    }
  }

  return result;
}

void update_H(
  NumericMatrix& H,
  const Rcpp::NumericVector& step,
  const Rcpp::NumericVector& diff) {
  int n = step.size();
  double yts = 0;
  for (int i = 0; i < n; i++) {
    yts += diff[i] * step[i];
  }
  NumericMatrix left = simple_dot(step, diff, yts);
  NumericMatrix right = simple_dot(diff, step, yts);
  NumericMatrix sst = simple_dot(step, step, yts);

  // simulate I - M
  for (int i = 0; i < n; i++) {
    left(i, i) += 1;
    right(i, i) += 1;
  }

  NumericMatrix H_new(n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      NumericVector row = left(i, _);
      NumericVector col = H(_, j);
      H_new(i, j) = sum(row * col);
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      NumericVector row = H_new(i, _);
      NumericVector col = right(_, j);
      H(i, j) = sum(row * col) + sst(i, j);
    }
  }

}

Eigen::MatrixXd update_H2(
  Eigen::MatrixXd& H,
  const Rcpp::NumericVector& step,
  const Rcpp::NumericVector& diff) {
  int n = step.size();
  double yts = 0;
  for (int i = 0; i < n; i++) {
    yts += diff[i] * step[i];
  }

  Eigen::MatrixXd left = -simple_dot2(step, diff, yts);
  Eigen::MatrixXd right = -simple_dot2(diff, step, yts);
  Eigen::MatrixXd sst = simple_dot2(step, step, yts);

  // simulate I - M
  for (int i = 0; i < n; i++) {
    left(i, i) += 1;
    right(i, i) += 1;
  }

  Eigen::MatrixXd result = left * H * right;

  return result + sst;

}



double line_search_cpp2(const Rcpp::NumericVector& Y,
                       const Rcpp::NumericMatrix& X,
                       const Rcpp::NumericVector& beta,
                       const Rcpp::NumericVector& gamma,
                       const Rcpp::NumericVector& direc,
                       double alpha = 1) {
  const double alphas[18] = {
    -0.000000001,
    -0.00000001,
    -0.0000001,
    -0.000001,
    -0.00001,
    -0.0001,
    -0.001,
    -0.01,
    -0.1,
     0.1,
     0.01,
     0.001,
     0.0001,
     0.00001,
     0.000001,
     0.0000001,
     0.00000001,
     0.000000001
  };

  int index = 0;
  double objective = INFINITY;
  int n = 18;
  int q = beta.size();
  Rcpp::NumericVector beta_new = clone(beta);
  Rcpp::NumericVector gamma_new = clone(gamma);
  Rcpp::Range beta_slice = Rcpp::Range(0, q - 1);
  Rcpp::Range gamma_slice = Rcpp::Range(q, direc.size() - 1);
  Rcpp::NumericVector grand_beta = direc[beta_slice];
  Rcpp::NumericVector grand_gamma = direc[gamma_slice];

  for (int i = 0; i<n; i++) {
    double step = alphas[i];
    beta_new = beta + step * grand_beta;
    gamma_new = gamma + step * grand_gamma;
    double f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
    if (f_new < objective) {
      index = i;
      objective = f_new;
    }
  }
  return alphas[index];
}

// attempt at optimizing line search such that we
// calculate `loglik_cpp()` as little as possible
//
// keep a reference "index" into a constant set of
// scalars.
// &objective is the loglik_cpp() of a prior iteration
//
// immediately check if steping in direction (with current
// index scalar) is lower than prior loglik_cpp().
//
// if so --> current scalar works but maybe we could step
//           further... thus try larger magnitude scalar of direction
//           until it is no longer lower (or nan is reached)
// if not --> current scalar is too large...
//            thus try smaller scalars (till end of list)
//
// regardless of these branches, the oposite direction is
// always checked (starting from the smallest magnitude).
// if at any point a smaller objective could not be found
// then this "reverse" search is stopped immediately
double line_search_cpp3(const Rcpp::NumericVector& Y,
                       const Rcpp::NumericMatrix& X,
                       const Rcpp::NumericVector& beta,
                       const Rcpp::NumericVector& gamma,
                       const Rcpp::NumericVector& direc,
                       int& index,
                       double& objective) {
  const double alphas[20] = {
     0.5,
     0.1,
     0.01,
     0.001,
     0.0001,
     0.00001,
     0.000001,
     0.0000001,
     0.00000001,
     0.000000001,
    -0.000000001,
    -0.00000001,
    -0.0000001,
    -0.000001,
    -0.00001,
    -0.0001,
    -0.001,
    -0.01,
    -0.1,
    -0.5
  };

  //double objective = INFINITY;
  // int n = 20;
  int q = beta.size();
  Rcpp::NumericVector beta_new = clone(beta);
  Rcpp::NumericVector gamma_new = clone(gamma);
  Rcpp::Range beta_slice = Rcpp::Range(0, q - 1);
  Rcpp::Range gamma_slice = Rcpp::Range(q, direc.size() - 1);
  Rcpp::NumericVector grand_beta = direc[beta_slice];
  Rcpp::NumericVector grand_gamma = direc[gamma_slice];

  double step = alphas[index];
  beta_new = beta + step * grand_beta;
  gamma_new = gamma + step * grand_gamma;
  double f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
  // Rcout << "f: " << f_new << " obj: " << objective << " step: " << step << "\n";
  if (f_new < objective && !std::isnan(f_new)) {

    // f_new is smaller than the objective and is a number
    objective = f_new;
    // now try larger values
    if (index < 10) {
      // in the positive region
      // loop through starting from index
      // and checking larger scalings until
      // f_new is no longer smaller.
      for(int i = index - 1; i >= 1; i--) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
        } else {
          // we have increased enough
          break;
        }
      }
      // for completeness, check i == 0
      if (index == 1) {
        step = alphas[0];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = 0;
        }
      }
      // for verbosity, check other direction...
      // starting with the smallest value
      for(int i = 10; i < 20; i++) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
        } else {
          // failed to find better value
          break;
        }
      }
    } else {
      // in the negatives
      // try larger negatives - i.e. increase index
      for(int i = index + 1; i < 20; i++) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
        } else {
          // we have decreased enough
          break;
        }
      }

      // for verbosity, check other direction...
      // starting with the smallest value
      for(int i = 9; i > 0; i--) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
        } else {
          // we have decreased enough
          break;
        }
      }
      // for completeness, check i == 0
      if (index == 1) {
        step = alphas[0];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = 0;
        }
      }

    }
  } else {
    // Rcout << "f_new was either nan or greater than objective\n";
    // f_new is larger than the objective or f_new is not a number
    // *** try a smaller alpha ***
    //
    // if we are looking in the positive direction
    if (index < 10) {
      // loop through up until index 10 - i.e. when it becomes
      // negative
      for(int i = index + 1; i < 10; i++) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        // Rcout << "f: " << f_new << " obj: " << objective << " step: " << step << "\n";
        // once we find an a value is just small enough - stop
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
          break;
        }
      }
      // for verbosity, check other direction...
      // starting with the smallest value
      // and keep inceased so long as objective
      // decreases
      for(int i = 10; i < 20; i++) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
        } else {
          // failed to find better value
          break;
        }
      }

    } else {
      // we are in the negatives
      for(int i = index - 1; i > 9; i--) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
          // we have decreased enough
          break;
        }
      }

      // for verbosity, check other direction...
      // starting with the smallest value
      for(int i = 9; i > 0; i--) {
        step = alphas[i];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = i;
        } else {
          // we have decreased enough
          break;
        }
      }
      // for completeness, check i == 0
      if (index == 1) {
        step = alphas[0];
        beta_new = beta + step * grand_beta;
        gamma_new = gamma + step * grand_gamma;
        f_new = -loglik_cpp(Y, X, beta_new, gamma_new);
        if (f_new < objective && !std::isnan(f_new)) {
          objective = f_new;
          index = 0;
        }
      }
    }

  }
  return alphas[index];
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
  bool verbose = false
) {

  NumericVector beta = clone(beta0);
  NumericVector gamma = clone(gamma0);
  int q = beta0.size();
  int p = gamma0.size();
  double f = -loglik_cpp(Y, X, beta, gamma);
  double ep = INFINITY;
  int iter = 1;
  Rcpp::Range beta_slice = Rcpp::Range(0, q - 1);
  Rcpp::Range gamma_slice = Rcpp::Range(q, p + q - 1);
  NumericMatrix H = NumericMatrix::diag(p + q, 1.0);
  Rcpp::NumericVector grad = -loglik_grad_cpp(Y, X, beta, gamma);
  double step = 0;
  while (ep > tol && iter <= maxIter) {

    Rcpp::NumericVector direc = grad_direc(H, grad);

    step = line_search_cpp2(Y, X, beta, gamma, direc);
    if (verbose) Rcout << "Step size: " << step << " direction: " << direc << "\n";
    direc = direc * step;
    beta += direc[beta_slice];
    gamma += direc[gamma_slice];
    Rcpp::NumericVector grad_new = -loglik_grad_cpp(Y, X, beta, gamma);
    Rcpp::NumericVector grad_diff = grad_new - grad;
    update_H(H, direc, grad_diff);
    double f_new = -loglik_cpp(Y, X, beta, gamma);
    ep = std::abs(f_new - f);
    f = f_new;
    grad = grad_new;
    if (verbose) Rcout << "Iteration: " << iter << ", eps: " << ep << ", ll: " << f << ", beta: " << beta << ", gamma: " << gamma << "\n";
    iter++;
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma
  );
}

// [[Rcpp::export]]
Rcpp::List bfgs_cpp2(
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

    Rcpp::NumericVector direc = grad_direc2(H, grad);
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
    H = update_H2(H, direc, grad_diff);
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
