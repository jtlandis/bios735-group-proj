devtools::load_all()
mod <- data_set_tidy |>
  par_model_mat(
    QTY ~ PROMO * brand * item,
    time = DATE,
    nlag = 4,
    groups = c(brand, item)
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
mod$beta <- res$beta
mod$gamma <- res$gamma
sink()
walk_promo_brand_item <- readLines(temp_file)
mod_promo_brand_item <- mod
usethis::use_data(walk_promo_brand_item, overwrite = TRUE)
usethis::use_data(mod_promo_brand_item, overwrite = TRUE)
