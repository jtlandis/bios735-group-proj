test_that("model_mat errors", {
  expect_error(
    par_model_mat(
      data = data_set_tidy,
      formula = QTY ~ PROMO + brand,
      time = DATE,
      nlag = 2
    )
  )

  expect_error(
    par_model_mat(
      data = data_set_tidy,
      formula = QTY ~ PROMO + brand,
      time = DATE,
      nlag = 2,
      groups = item
    )
  )

  expect_no_error(
    par_model_mat(
      data = data_set_tidy,
      formula = QTY ~ PROMO + brand,
      time = DATE,
      nlag = 2,
      groups = c(brand, item)
    )
  )
})

test_that("model_mat returns correct size w/ groups", {
  .dummy_data <- data.frame(
    x = rep(1:5, 3),
    y = rpois(15, 4),
    g = rep(1:3, each = 5)
  )

  mod <- par_model_mat(
    data = .dummy_data,
    y ~ x,
    time = x,
    groups = g,
    nlag = 0
  )

  expect_equal(nrow(mod$X), 15)
  expect_equal(length(mod$Y), 15)

  mod <- par_model_mat(
    data = .dummy_data,
    y ~ x,
    time = x,
    groups = g,
    nlag = 1
  )
  # minus 1 data point per group
  expect_equal(nrow(mod$X), 12)
  expect_equal(length(mod$Y), 12)

  mod <- par_model_mat(
    data = .dummy_data,
    y ~ x,
    time = x,
    groups = g,
    nlag = 2
  )
  # minus 1 data point per group
  expect_equal(nrow(mod$X), 9)
  expect_equal(length(mod$Y), 9)
})
