simple_data <- subset(data_set_tidy, brand == "B1")

test_that("par_model: Intercept only", {
  mod <- expect_no_error(
    par_model_mat(simple_data, QTY ~ 1, nlag = 0, time = DATE, groups = item)
  )
  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )
  res <- expect_no_error(bfgs_cpp(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  # note this is barely better
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})

test_that("par_model: lags only", {
  mod <- expect_no_error(
    par_model_mat(simple_data, QTY ~ -1, nlag = 3, time = DATE, groups = item)
  )
  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )
  res <- expect_no_error(bfgs_cpp(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})

test_that("par_model: lags and covar", {
  mod <- expect_no_error(
    par_model_mat(
      simple_data,
      QTY ~ PROMO,
      nlag = 3,
      time = DATE,
      groups = item
    )
  )
  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )
  res <- expect_no_error(bfgs_cpp(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})

test_that("par_model: estimate covariates by ave items", {
  mod <- data_set_tidy |>
    filter(item %in% c("1", "5")) |>
    par_model_mat(
      formula = QTY ~ PROMO + brand,
      time = DATE,
      nlag = 2,
      groups = c(brand, item)
    )

  ll_old <- expect_no_error(
    loglik_cpp(mod$Y, mod$X, mod$beta, mod$gamma)
  )

  res <- expect_no_error(bfgs_cpp(mod$Y, mod$X, mod$beta, mod$gamma))
  expect_snapshot(print(res))
  expect_true(loglik_cpp(mod$Y, mod$X, res$beta, res$gamma) > ll_old)
})
