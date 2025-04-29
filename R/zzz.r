# To remove warnings of "unused global variable"

utils::globalVariables(
  c("run_mcmc_par_cpp", "run_mcmc_pvar_cpp", "loglik_pvar_cpp")
)

`:=` <- rlang::`:=`
QTY <- DATE <- brand_id <- item_id <- brand <- item <- y <- time <- NULL
...i <- ...pseudo_col <- ..i <- .dup <- .env <- i <- ll <- NULL
