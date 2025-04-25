# bios735-group-proj

```{r}
library(pastasales)
library(dplyr)

model_spec <- data_set_tidy |>
  # smaller subset for faster result
  filter(
    brand %in% c("B1", "B2"),
    item  %in% c("1", "2")
  ) |>
  # prepares data for fit function
  par_model_mat(
    # model design
    QTY ~ PROMO + brand*item,
    # inform which column is the time variable
    # so that we can insure data is in correct
    # order
    time = DATE,
    # this function will error if it detects
    # duplicate time points per time series
    groups = c(brand, item),
    nlag = 2
  )

model_spec

res <- fit_par_mcmc(model_spec$Y, model_spec$X, q = length(model_spec$beta))
res$beta
res$gamma
res2 <- fit_par_bfgs(model_spec)
res2$beta
res2$gamma

```
