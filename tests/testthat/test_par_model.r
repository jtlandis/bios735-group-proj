simple_data <- subset(data_set_tidy, brand = "B1")

test_that("par_model: Intercept only", {
  mod <- expect_no_error(
    par_model_mat(simple_data, QTY ~ 1, nlag = 0, time = DATE)
  )
  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )
  res <- expect_no_error(bfgs_cpp2(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  # note this is barely better
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})

test_that("par_model: lags only", {
  mod <- expect_no_error(
    par_model_mat(simple_data, QTY ~ -1, nlag = 3, time = DATE)
  )
  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )
  res <- expect_no_error(bfgs_cpp2(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})

test_that("par_model: lags and covar", {
  mod <- expect_no_error(
    par_model_mat(simple_data, QTY ~ PROMO, nlag = 3, time = DATE)
  )
  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )
  res <- expect_no_error(bfgs_cpp2(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})
