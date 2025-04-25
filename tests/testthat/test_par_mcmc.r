
# mcmc is random samples so its hard to know
# if it is working _exactly_ the same.
# we compare with original bfgs values
test_that("par_mcmc works", {
  mcmc <- fit_par_mcmc(
    y = data_set_raw$QTY_B1_1,
    x = cbind(1, data_set_raw$PROMO_B1_1)
  )
  expect_equal(mcmc$beta, 0.409, tolerance = 1e-2)
  expect_equal(mcmc$gamma, c(1.407, 1.527), tolerance = 1e-2)
})
