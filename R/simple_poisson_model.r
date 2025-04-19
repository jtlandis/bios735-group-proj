# log(\lambda_{it}) = \alpha + \sum_{l=1}^q \beta_l y_{i,t-1} + \gamma \x_{it}

sample_pois <- function(x, alpha, betas, gamma) {
  n <- ncol(x)
  t <- nrow(x)
  y <- matrix(NA, nrow = t, ncol = n)
  q <- length(betas)
  int <- alpha + (gamma * x)

  # get inital points
  y[1, ] <- rpois(n = n, lambda = exp(int[1, ]))
  # build up to q + 1
  for (i in 2:(q + 1)) {
    b <- betas[seq_len(i - 1)]
    slice_i <- (i - 1):1
    y[i, ] <- rpois(
      n = n,
      lambda = exp(int[i, ] + colSums(b * y[slice_i, , drop = FALSE]))
    )
  }

  if (t > q + 1) {
    for (i in (q + 2):t) {
      slice_i <- slice_i + 1
      y[i, ] <- rpois(
        n = n,
        lambda = exp(int[i, ] + colSums(betas * y[slice_i, , drop = FALSE]))
      )
    }
  }
  y
}

ll_pois <- function(y, x, alpha, betas, gamma) {
  Ay <- beta_mat(betas, t = nrow(y)) %*% y
  lambdas <- exp(alpha + Ay + (gamma * x))

  sum(y * log(lambdas) - lambdas - log(factorial(y)))
}

ll_pois_prop <- function(y, x, alpha, betas, gamma) {
  Ay <- beta_mat(betas, t = nrow(y)) %*% y
  lambdas <- exp(alpha + Ay + (gamma * x))

  sum(y * log(lambdas) - lambdas)
}

lambdas <- function(y, x, alpha, betas, gamma) {
  exp(alpha + (beta_mat(betas, t = nrow(y)) %*% y) + (gamma * x))
}

dll_pois <- function(y, x, alpha, betas, gamma) {
  lambdas <- exp(alpha + (beta_mat(betas, t = nrow(y)) %*% y) + (gamma * x))
  diff <- y - lambdas
  list(
    alpha = sum(diff),
    beta = vapply(
      seq_along(betas),
      function(i, d) {
        sum(d * shift_mat_row_by(y, i))
      },
      numeric(1),
      d = diff
    ),
    gamma = sum(diff * x)
  )
}

ddll_pois <- function(y, x, alpha, betas, gamma) {
  lambdas <- exp(alpha + (beta_mat(betas, t = nrow(y)) %*% y) + (gamma * x))

  d_beta <- lapply(
    seq_along(betas),
    function(i, y) {
      shift_mat_row_by(y, i)
    },
    y = y
  )
  data <- c(
    list(
      alpha = 1
    ),
    setNames(d_beta, sprintf("beta_%i", seq_along(betas))),
    list(gamma = x)
  )
  dim <- length(betas) + 2
  out <- matrix(0, nrow = dim, ncol = dim)
  for (k in seq_len(length(out)) - 1) {
    i <- k %% dim + 1
    j <- k %/% dim + 1
    out[i, j] <- sum(lambdas * data[[i]] * data[[j]])
  }
  -out
}

bfgs <- function(y, x, alpha, betas, gamma, max_iter = 1000, trace_mod = 1) {
  q <- length(betas)
  if (q == 0) {
    beta_slice <- 0
  } else {
    beta_slice <- 2:(q + 1)
  }
  gamma_slice <- q + 2
  ident <- diag(q + 2)
  f <- function(params) {
    ll_pois(y, x, params[1], params[beta_slice], params[gamma_slice])
  }
  g <- function(params) {
    unlist(dll_pois(y, x, params[1], params[beta_slice], params[gamma_slice]))
  }

  H <- ddll_pois(y, x, alpha, betas, gamma) |> solve()

  p <- function(params) {
    -H %*% g(params)
  }

  step_by <- function(params, by) {
    p(params) * by
  }

  update <- function(s, y) {
    syt <- s %*% t(y)
    yts <- as.vector(y %*% s)
    yst <- y %*% t(s)
    sst <- s %*% t(s)

    (ident - (syt / yts)) %*% H %*% (ident - (yst / yts)) + (sst / yts)
  }

  params <- c(alpha, betas, gamma)
  grade <- g(params)
  LL <- f(params)
  count_cache <- count <- 0L
  step_size <- 0.1
  diffs <- numeric(max_iter)
  diff_slice <- 0:2
  while (count < max_iter) {
    count <- count + 1L
    s <- step_by(params, step_size)
    params_ <- params + s
    grade_ <- g(params_)
    H_ <- update(s, grade_ - grade)
    LL_ <- f(params_)
    LL_diff <- LL_ - LL
    diffs[count] <- LL_diff
    if (count %% trace_mod == 0) {
      cat(
        sprintf(
          "Iteration %i: LL = %.3f, LL_old: %.3f, params: %s, step size: %.6f\n",
          count,
          LL_,
          LL,
          paste(round(params_, 3), collapse = ", "),
          step_size
        )
      )
    }

    if (abs(LL_ - LL) < 1e-6) {
      cat("INFO: Convergence reached\n")
      return(params_)
    }

    # if we every increase in LL try again
    # with a smaller step step_size
    if (LL_ < LL) {
      if (count - count_cache > 1) {
        cat("INFO: Decreasing step size\n")
        step_size <- step_size / 2
        #count <- count - 1L
        count_cache <- count
        next
      }
    } else {
      if (count - count_cache > 2) {
        .diff <- diff(diffs[(count - 2):count])
        slp <- .diff / mean(.diff) - 1
        if (all(slp > -0.2) && all(slp < 0.2)) {
          cat("INFO: Increasing step size\n")
          step_size <- step_size * 1.8
          count_cache <- count
        }
      }
    }

    LL <- LL_
    grade <- grade_
    params <- params_
    H <- H_
  }

  params
}
