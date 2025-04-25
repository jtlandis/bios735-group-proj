devtools::load_all()
mod <- par_model_mat(
  data_set_tidy,
  QTY ~ PROMO * brand * item,
  time = DATE,
  nlag = 4,
  groups = c(brand, item)
)
.call <- rlang::expr(
  par_model_mat(
    data_set_tidy,
    QTY ~ PROMO * brand * item,
    time = DATE,
    nlag = 4,
    groups = c(brand, item)
  )
)
temp_file <- tempfile()
time_expr <- function(x) {
  .then <- Sys.time()
  x <- force(x)
  print(.then - Sys.time())
  x
}
sink(temp_file, split = TRUE)
res <- time_expr({
  fit_par_bfgs(mod, global_tol = 1e-3, verbose = 3)
})
sink()
walk_promo_brand_item <- readLines(temp_file)
mod_promo_brand_item <- list(mod_call = .call, mod_fit = res)
usethis::use_data(walk_promo_brand_item, overwrite = TRUE)
usethis::use_data(mod_promo_brand_item, overwrite = TRUE)
