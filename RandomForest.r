# Load necessary libraries
library(randomForest)
library(ranger)
library(dplyr)
library(lubridate)
library(ggplot2)
library(caret)
library(tidyverse)
library(slider)

#########################################################
#     One-step prediction                               #
#########################################################
# ---------- 0. read, parse date ------------------------------------------------
sales <- read_csv("hierarchical_sales_data.csv") %>%
  mutate(DATE = ymd(DATE)) %>% 
  arrange(DATE)                     # ensure chronological order

# ---------- 1. lag features for the focal item ---------------------------------
sales <- sales %>%
  mutate(
    lag1_QTY_B1_1 = lag(QTY_B1_1, 1),           # yesterday’s sales
    lag7_QTY_B1_1 = lag(QTY_B1_1, 7)            # same weekday last week
  )

# ---------- 2. total Brand-1 sales each day ------------------------------------
sales <- sales %>%
  mutate(
    total_B1     = select(., starts_with("QTY_B1_")) %>% rowSums(),
    roll7_totalB1 = slide_dbl(total_B1, mean, .before = 6, .complete = TRUE)  # 7-day MA
  )

# ---------- 3. calendar + promotion variables ---------------------------------
sales <- sales %>%
  mutate(
    dow     = wday(DATE, label = TRUE),          # factor
    month   = month(DATE),
    PROMO_B1_1 = factor(PROMO_B1_1, levels = c(0,1))
  )

# ---------- 4. drop rows with NA lags (first week) -----------------------------
sales_model <- sales %>% drop_na()

# ---------- 5. train / test split + Random-Forest ------------------------------
set.seed(123)
n<-nrow(sales)
train_idx <- 1:round(n*0.8)
train <- sales_model[train_idx, ]
test  <- sales_model[-train_idx, ]

# ---------- 5½. hyper-parameter tuning  ------------------------------

library(ranger)
library(furrr)        # parallel mapping
library(tictoc)       # simple timing
#plan(multisession, workers = parallel::detectCores() - 1)

set.seed(123)

param_grid <- expand.grid(
  mtry           = c(2, 3, 4, 5),          # √p ± 1
  min.node.size  = c(3, 5, 10),            # leaf sizes
  sample.fraction= c(.6, .8, 1.0)          # sub-bagging vs bootstrap
)

tic("grid search")
tuning_out <- future_map_dfr(
  1:nrow(param_grid),
  function(i) {
    p <- param_grid[i, ]
    fit <- ranger(
      QTY_B1_1 ~ PROMO_B1_1 + dow + month + lag1_QTY_B1_1 + lag7_QTY_B1_1 +
                   total_B1 + roll7_totalB1,
      data            = train,    # <-- all rows (OOB handles CV)
      num.trees       = 500,
      mtry            = p$mtry,
      min.node.size   = p$min.node.size,
      sample.fraction = p$sample.fraction,
      importance      = "impurity",
      oob.error       = TRUE
    )
    tibble(
      mtry            = p$mtry,
      min.node.size   = p$min.node.size,
      sample.fraction = p$sample.fraction,
      RMSE_oob        = sqrt(fit$prediction.error)
    )
  }
)
toc()


## 1) heatmap of OOB-RMSE by (mtry, min.node.size)
tuning_out2 <- tuning_out %>%                       # keep lowest RMSE per tile
  group_by(mtry, min.node.size) %>% 
  summarise(RMSE_oob = min(RMSE_oob), .groups = "drop")

tuning <- ggplot(tuning_out2, aes(factor(mtry), factor(min.node.size), fill = RMSE_oob)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(RMSE_oob, 2)), size = 3) +
  scale_fill_viridis_c(option = "E") +
  labs(title = "Best OOB-RMSE per (mtry, min.node.size)",
       x = "mtry", y = "min.node.size", fill = "RMSE") +
  theme_minimal(base_size = 12)
# ggsave("bios735/hyperpara_tuning.jpg",tuning, width = 10, height = 6, dpi = 300)

## 2) profile of best RMSE versus sample.fraction
best_by_sf <- tuning_out %>%
  group_by(sample.fraction) %>%
  slice_min(RMSE_oob, n = 1)

tuning2<-ggplot(best_by_sf, aes(sample.fraction, RMSE_oob)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Best RMSE per bootstrap fraction",
       x = "sample.fraction",
       y = "OOB-RMSE") +
  theme_minimal(base_size = 12)
# ggsave("bios735/hyperpara_tuning2.jpg",tuning2, width = 10, height = 6, dpi = 300)
# ──────────  pick the best setting and continue  ────────────────────────────
best <- tuning_out %>% slice_min(RMSE_oob, n = 1)

mtry_best           <- best$mtry
min_node_best       <- best$min.node.size
sample_fraction_best<- best$sample.fraction

cat("Best grid point:",
    "mtry =", mtry_best,
    "min.node.size =", min_node_best,
    "sample.fraction =", sample_fraction_best,
    "OOB-RMSE =", round(best$RMSE_oob,3), "\n")

# now proceed with the TIME-RESPECTING train/test split


rf_fit <- ranger(
  QTY_B1_1 ~ PROMO_B1_1 + dow + month + lag1_QTY_B1_1 + lag7_QTY_B1_1 +
              total_B1 + roll7_totalB1,
  data       = train,
  num.trees  = 500,
  min.node.size = 10,
  mtry = 4,
  sample.fraction = 1,
  importance = "impurity"   # Gini importance
)

# ---------- 6. evaluate --------------------------------------------------------
pred  <- predict(rf_fit, data = test)$predictions
mse <- mean((pred - test$QTY_B1_1)^2)
rmse  <- sqrt(mean((pred - test$QTY_B1_1)^2))
cat("Test MSE:", round(mse, 3), "\n")
cat("Test RMSE:", round(rmse, 3), "\n")
cat("Variable importance:\n")
print(sort(rf_fit$variable.importance, decreasing = TRUE))

imp_df <- enframe(rf_fit$variable.importance,
                  name  = "Variable",
                  value = "Gini") %>% 
  arrange(desc(Gini))              # high to low

# Optional: keep only top N variables to avoid clutter
# imp_df <- imp_df %>% slice_max(Gini, n = 20)

plot_gini<- ggplot(imp_df, aes(x = reorder(Variable, Gini), y = Gini)) +
  geom_col(fill = "dodgerblue4") +
  coord_flip() +
  labs(title = "Gini variable importance (random forest)",
       x = NULL, y = "Importance") +
  theme_minimal(base_size = 12)
# ggsave("bios735/rf_gini_importance.jpg", plot_gini, width = 10, height = 6, dpi = 300)

# ---------------7. Time–series plot: Actual vs Predicted on the test set -----------------------
plot_df <- test %>%                                   # rows already in chrono order
  mutate(
    Actual    = QTY_B1_1,
    Predicted = pred
  ) %>% 
  select(DATE, Actual, Predicted) %>% 
  pivot_longer(-DATE, names_to = "Series", values_to = "Value")

plot_pred <- ggplot(plot_df, aes(DATE, Value, colour = Series)) +
  geom_line(size = 0.9) +
  labs(title = "Random-Forest predictions vs. actual sales",
       x = NULL, y = "Units sold (Brand 1 Item 1)", colour = "") +
  scale_colour_manual(values = c(Actual = "grey",
                                 Predicted = "red")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# ggsave("bios735/rf_prediction.jpg", plot_pred, width = 10, height = 6, dpi = 300)


################################################
#     H-step forecast  (H = length(test))      #
################################################

library(dplyr)
library(tidyr)
library(ranger)

H <- nrow(test)

# ----------------------- 1.  Build a "rolling" data frame that holds both history + future -----------------------
rolling <- bind_rows(train, test %>% mutate(QTY_B1_1 = NA_real_))  # NA so we can fill

# -----------------------2. Predict step-by-step, updating lag features each time -----------------------
start_ix <- nrow(train)

for (h in 1:H) {
  idx      <- start_ix + h
  lag1_ix  <- idx - 1
  lag7_ix  <- idx - 7

  # refresh lag features for *this* row
  rolling$lag1_QTY_B1_1[idx] <- rolling$QTY_B1_1[lag1_ix]
  rolling$lag7_QTY_B1_1[idx] <- if (lag7_ix > 0) rolling$QTY_B1_1[lag7_ix] else NA

  # total_B1 and roll7_totalB1 — keep the last observed item share
  last_share <- tail(with(rolling[1:(idx-1), ], QTY_B1_1 / total_B1), 1)
  rolling$total_B1[idx] <- rolling$lag1_QTY_B1_1[idx] / last_share
  if (idx >= 7)
    rolling$roll7_totalB1[idx] <- mean(rolling$total_B1[(idx-6):idx])

  # run the forest on the single-row data frame
  pred <- predict(rf_fit, data = rolling[idx, ])$predictions
  rolling$QTY_B1_1[idx] <- pred      # store prediction so next step can use it
}

# ----------------------- 3.  Collect predictions & actuals for plotting / metrics -----------------------
forecast_df <- rolling[(start_ix+1):(start_ix+H), ] %>%
  select(DATE, Predicted = QTY_B1_1) %>%
  mutate(Actual = test$QTY_B1_1)

# RMSE over the *recursive* horizon
mse_H <- mean((forecast_df$Predicted - forecast_df$Actual)^2, na.rm = TRUE)
cat("Recursive H-step MSE =", round(mse_H, 3), "\n")
rmse_H <- sqrt(mean((forecast_df$Predicted - forecast_df$Actual)^2, na.rm = TRUE))
cat("Recursive H-step RMSE =", round(rmse_H, 3), "\n")

# ----------------------- 4.Prediction vs. Actual – full test horizon -----------------------                             
plot_df <- forecast_df %>%
             pivot_longer(cols = c(Predicted, Actual),
                          names_to = "Series", values_to = "Units")

plot_pred_h <- ggplot(plot_df, aes(DATE, Units, colour = Series)) +
  geom_line(size = 0.9) +
  labs(title = paste0(" H -step forecast vs. truth"),
       x = NULL, y = "Units sold (Brand-1 Item-1)", colour = "") +
  scale_colour_manual(values = c(Actual = "grey", Predicted = "red")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
# ggsave("bios735/rf_performance_hstep.jpg", plot_pred_h, width = 10, height = 6, dpi = 300)
