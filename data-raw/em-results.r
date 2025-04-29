devtools::load_all()
library(dplyr)
library(tidyr)
## long-running process
### discovered a few bugs when any data
### has missing values... - did not fix this yet...
model_em_all <- data_set_tidy |>
  mutate(
    data_set_tidy_index = seq_len(n()),
    brand = as.factor(brand),
    item = as.integer(item)
  ) |>
  group_by(brand, item) |>
  summarise(
    model = par_model_mat(
      pick(DATE, QTY, PROMO, data_set_tidy_index),
      QTY ~ PROMO,
      nlag = 5
    ) |>
      list(),
    em = {
      print(cur_group())
      par_em_effective(
        model[[1]],
        col = PROMO,
        subset = PROMO == 1,
        tol = 1e-6,
        show_plots = TRUE
      ) |>
        list()
    }
  ) |>
  rowwise() |>
  mutate(
    pi = em$pi,
    ll = em$ll,
    ll_old = em$ll_original
  )

model_em_all <- dplyr::transmute(
  model_em_all,
  item = item,
  pi = pi,
  em_spec = list(em$spec),
  original_spec = list(em$original_spec)
)

usethis::use_data(model_em_all, overwrite = TRUE)
