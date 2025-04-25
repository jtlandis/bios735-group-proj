
# mcmc is random samples so its hard to know
# if it is working _exactly_ the same.
# we compare with original bfgs values
test_that("par_mcmc works - Brand1 item 1", {
  mcmc <- fit_par_mcmc(
    y = data_set_raw$QTY_B1_1,
    x = cbind(1, data_set_raw$PROMO_B1_1)
  )
  expect_equal(mcmc$beta, 0.409, tolerance = 1e-2)
  expect_equal(mcmc$gamma, c(1.407, 1.527), tolerance = 1e-2)
})

#test_that("par_mcmc works - Brand2 item 1", {
#  mcmc <- fit_par_mcmc(
#    y = data_set_raw$QTY_B2_1,
#    x = cbind(1, data_set_raw$PROMO_B2_1)
#  )
#  expect_equal(mcmc$beta, -0.0064, tolerance = 1e-2)
#  expect_equal(mcmc$gamma, c(-0.780, 1.352), tolerance = 1e-2)
#})

#test_that("par_mcmc works - Brand3 item 1", {
#  mcmc <- fit_par_mcmc(
#    y = data_set_raw$QTY_B3_1,
#    x = cbind(1, data_set_raw$PROMO_B3_1)
#  )
#  expect_equal(mcmc$beta, 0.0982, tolerance = 1e-2)
#  expect_equal(mcmc$gamma, c(-0.669, 1.385), tolerance = 1e-2)
#})

#test_that("par_mcmc works - Brand4 item 1", {
#  mcmc <- fit_par_mcmc(
#    y = data_set_raw$QTY_B4_1,
#    x = cbind(1, data_set_raw$PROMO_B4_1)
#  )
#  expect_equal(mcmc$beta, 0.330, tolerance = 1e-2)
#  expect_equal(mcmc$gamma, c(-1.350, 1.928), tolerance = 1e-2)
#})
