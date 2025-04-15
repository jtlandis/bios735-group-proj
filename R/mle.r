mle_alpha <- function(y, x, betas, gamma, n = 1) {
  t <- length(y)
  B <- beta_mat(betas, t, q = length(betas))
  A <- diag(t) - B
  Ay <- A %*% y
  gammax <- gamma * x
  sum(Ay - gammax) / (t * n)
}

mle_gamma <- function(y, x, betas, alpha, n = 1) {
  t <- length(y)
  B <- beta_mat(betas, t, q = length(betas))
  A <- diag(t) - B
  Ay <- A %*% y
  ((x %*% Ay) - sum(alpha * x)) / sum(x * x)
}

mle_sigma2 <- function(y, x, betas, alpha, gamma, n = 1) {
  t <- length(y)
  B <- beta_mat(betas, t, q = length(betas))
  A <- diag(t) - B
  Ay <- A %*% y
  gammax <- gamma * x
  r <- Ay - gammax - alpha
  sum(r * r) / (t * n)
}

shift_by <- function(x, n) {
  N <- length(x)
  c(rep(0, n), x[seq_len(N - n)])
}

align_head_tail <- function(vec_head, vec_tail, n) {
  list(
    vec_head[seq_len(length(vec_head) - n)],
    vec_tail[-seq_len(n)]
  )
}

mle_beta <- function(y, x, alpha, betas, gamma, n = 1, index) {
  # the derivative wrt beta_index
  z <- shift_by(y, index) * -1
  beta_sums <- vapply(
    seq_along(betas),
    function(i, betas, y, z) {
      beta <- betas[i]
      aligned <- align_head_tail(y, z, i)
      sum(aligned[[1]] * aligned[[2]])
    },
    numeric(1),
    betas = betas,
    y = y,
    z = z
  )
  time_sums <- sum(y * z)
  scale <- -1 * beta_sums[index] / betas[index]
  beta_sums <- beta_sums[-index]

  mle <- time_sums - sum(beta_sums) - sum(alpha * z) - sum(gamma * x * z)

  mle / scale
}
