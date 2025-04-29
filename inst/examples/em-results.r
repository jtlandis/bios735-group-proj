library(pastasales)
library(dplyr)
library(tidyr)
library(ggplot2)
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
  original_spec = list(em$original_spec),
  em_loglik = ll,
  original_loglik = ll_old
)

em_data <- model_em_all |>
  mutate(
    deviance_em = get_deviance(em_spec),
    deviance_orig = get_deviance(original_spec),
    loglik_em = get_loglik(em_spec),
    loglik_orig = get_loglik(original_spec)
  )

panel_a <- filter(
  em_data,
  (brand == "B1" & item == "1") |
    (brand == "B2" & item == "34") |
    (brand == "B3" & item == "3") |
    (brand == "B4" & item == "5")
) |>
  reframe(
    brand = brand,
    item = item,
    pi = pi,
    Date = em_spec$.data$DATE,
    Quantity = em_spec$.data$QTY,
    PROMO = em_spec$.data$PROMO,
    z = em_spec$.data$z,
    em_mt = get_par_mt(em_spec),
    orig_mt = get_par_mt(original_spec)
  )
gg_panel_a <- panel_a |>
  ggplot(aes(Date)) +
  geom_line(aes(y = Quantity, color = factor("data"))) +
  geom_line(aes(y = orig_mt, color = factor("model")), linetype = "dotted") +
  geom_line(aes(y = em_mt, color = factor("em")), linetype = "dotted") +
  geom_rug(data = ~subset(.x, PROMO == 1)) +
  geom_rug(data = ~subset(.x, z == 1), color = "yellow") +
  geom_text(
    aes(label = label, x = x, y = y),
    data = ~summarise(
      group_by(.x, brand),
      label = unique(paste("Item:", item, "\npi estimate:", round(pi, 2))),
      x = min(Date),
      y = max(Quantity)
    ),
    hjust = 0,
    vjust = 1
  ) +
  #facet_grid(rows = vars(brand), scales = "free_y") +
  facet_wrap(~brand, nrow = 1, scales = "free_y") +
  scale_color_manual(
    values = c("data" = "grey", "model" = "blue", "em" = "red")
  ) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(color = NULL)

panel_b <- em_data |>
  pivot_longer(
    cols = starts_with("deviance"),
    names_to = "fit",
    names_pattern = "deviance_(.*)",
    names_transform = ~case_when(.x == "em" ~ "EM", TRUE ~ "Original"),
    values_to = "deviance"
  )
gg_panel_b <- panel_b |>
  mutate(fit = factor(fit, levels = c("Original", "EM"))) |>
  ggplot(aes(x = fit, y = deviance)) +
  geom_point(aes(color = pi)) +
  geom_line(aes(group = interaction(brand, item), color = pi)) +
  geom_text(
    data = ~summarise(group_by(.x, brand), mpi = mean(pi), y = 20000),
    mapping = aes(label = paste("Mean pi:", round(mpi, 2)), y = y),
    x = 2
  ) +
  facet_wrap(~brand, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Model Fit", y = "Deviance") +
  scale_color_viridis_c() +
  theme(legend.title = element_text(vjust = .8))

cowplot::plot_grid(gg_panel_a, gg_panel_b, nrow = 2, labels = c("A", "B"))
