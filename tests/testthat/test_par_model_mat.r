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
