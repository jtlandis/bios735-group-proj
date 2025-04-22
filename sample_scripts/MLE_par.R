library(Rcpp)
# setwd('~/Dropbox/UNC/Spring2025/BIOS735/project/')
sourceCpp('MLE_par.cpp')

### Functions for MLE for single item Poisson AR
get_mt <- function(Y, X, beta, gamma) {
  n <- length(Y)
  q <- length(beta)
  mt <- rep(NA, n)
  sum_beta <- sum(beta)
  cov_parts <- exp(X %*% as.matrix(gamma))
  for (t in (q+1):n) {
    mt[t] <- t(rev(beta)) %*% Y[(t-q):(t-1)] + (1-sum_beta)*cov_parts[t]
  }
  return(mt)
}
get_mt_grad <- function(Y, X, beta, gamma) {
  n <- length(Y)
  q <- length(beta)
  p <- length(gamma)
  
  grad_mt <- matrix(nrow = n, ncol = q + p)
  colnames(grad_mt) <- c(paste0("beta",1:q), paste0("gamma",1:p))
  
  # get derivatives wrt beta
  cov_parts <- exp(X %*% as.matrix(gamma))
  for (l in 1:q) {
    grad_mt[(q+1):n, l] <- Y[(q+1):n - l] - cov_parts[(q+1):n]
  }
  
  # get derivatives wrt gamma
  grad_mt[(q+1):n, -c(1:q)] <- (1 - sum(beta)) * cov_parts[(q+1):n] * X[(q+1):n,]
  
  return(grad_mt)
}
loglik <- function(Y, X, beta, gamma) {
  n <- length(Y)
  q <- length(beta)
  mt <- get_mt(Y, X, beta, gamma)
  ll <- sum(Y[(q+1):n] * log(mt[(q+1):n])) - sum(mt[(q+1):n])
  return(ll)
}
loglik_grad <- function(Y, X, beta, gamma) {
  n <- length(Y)
  q <- length(beta)
  mt <- get_mt(Y, X, beta, gamma)
  mt_grad <- get_mt_grad(Y, X, beta, gamma)
  ll_grad <- t((Y[(q+1):n] / mt[(q+1):n] - 1)) %*% mt_grad[(q+1):n,]
  return(ll_grad)
}
proj_beta <- function(beta, epsilon = 1e-4) {
  beta_proj <- pmax(beta, 0)
  sum_proj <- sum(beta_proj)
  if (sum_proj > 1) {
    u <- sort(beta_proj, decreasing = T)
    cssv <- cumsum(u)
    rho <- which(u > (cssv - 1) / seq_along(u)) |> tail(1)
    theta <- (cssv[rho] - 1) / rho
    beta_proj <- (1-epsilon) * pmax(beta_proj - theta, 0)
  }
  return(beta_proj)
}
proj_grad_descent <- function(Y, X, beta0, gamma0, lr=1e-5, maxIter=100, tol=1e-5, return_allIters = F, verbose = F) {
  q <- length(beta0)
  ep <- Inf
  iter <- 1
  beta <- beta0
  gamma <- gamma0
  f <- -loglik(Y, X, beta, gamma)
  
  if (return_allIters == T) {
    beta_iters <- vector(mode = "list", maxIter)
    gamma_iters <- vector(mode = "list", maxIter)
  }
  
  while (ep > tol & iter <= maxIter) {
    if (verbose == T) print(paste0('Iteration: ', iter))
    fgrad <- -loglik_grad(Y, X, beta, gamma)
    beta <- proj_beta(beta - lr * fgrad[1:q])
    gamma <- gamma - lr * fgrad[-c(1:q)]
    fnew <- -loglik(Y, X, beta, gamma)
    ep <- abs(fnew - f)
    f <- fnew
    if (verbose == T) print(paste0('ep = ', ep))
    if (return_allIters == T) {
      beta_iters[[iter]] <- beta
      gamma_iters[[iter]] <- gamma
    }
    iter <- iter + 1
  }
  
  if (verbose == T) print(paste0('Gradient descent stopped at iter = ', iter))
  
  if (return_allIters) {
    return(list(beta = beta, gamma = gamma, beta_iters = beta_iters, gamma_iters = gamma_iters))
  } else {
    return(list(beta = beta, gamma = gamma))
  }
}



### Test on data
data <- read.csv('hierarchical_sales_data.csv')
Y <- data$QTY_B1_1
X <- cbind(1, data$PROMO_B1_1) # include extra col of 1 for intercept

# get estimates using gradient descent
q <- 3
p <- 2
beta0 <- rep(0, q) # initial value for AR effect
gamma0 <- rep(0, p) # initial value for cov effect

estimates_R <- proj_grad_descent(Y, X, beta0, gamma0, lr = 1e-5, maxIter = 1000, verbose = T, return_allIters = T)
estimates_cpp <- proj_grad_descent_cpp(Y, X, beta0, gamma0, lr = 1e-5, maxIter = 1000, verbose = T, return_allIters = T)

# convergence plots
beta1_R <- unlist(lapply(estimates_R$beta_iters, function(x) x[1]))
beta1_cpp <-unlist(lapply(estimates_cpp$beta_iters, function(x) x[1]))
plot(beta1_R)
plot(beta1_cpp)

# get model fit
fit_R <- get_mt(Y, X, estimates_R$beta, estimates_R$gamma)
fit_cpp <- get_mt(Y, X, estimates_cpp$beta, estimates_cpp$gamma)

plot(Y, type = "l", col="black")
lines(fit_R, col="red")
lines(fit_cpp, col="blue")


