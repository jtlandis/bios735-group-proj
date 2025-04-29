# Load necessary libraries
library(randomForest)
library(ranger)
library(dplyr)
library(lubridate)
library(ggplot2)
library(caret)
library(tidyverse)
library(pastasales)

# ---------- 0. read, parse date ------------------------------------------------
sales <- data_set_raw %>%
  mutate(DATE = ymd(DATE)) %>%
  arrange(DATE) # ensure chronological order

# ---------- 1. lag features for the focal item ---------------------------------
sales <- sales %>%
  mutate(
    lag1_QTY_B1_1 = lag(QTY_B1_1, 1), # yesterday’s sales
    lag7_QTY_B1_1 = lag(QTY_B1_1, 7) # same weekday last week
  )

# ---------- 2. total Brand-1 sales each day ------------------------------------
sales <- sales %>%
  mutate(
    total_B1 = select(., starts_with("QTY_B1_")) %>% rowSums(),
    roll7_totalB1 = slide_dbl(total_B1, mean, .before = 6, .complete = TRUE) # 7-day MA
  )

# ---------- 3. calendar + promotion variables ---------------------------------
sales <- sales %>%
  mutate(
    dow = wday(DATE, label = TRUE), # factor
    month = month(DATE),
    PROMO_B1_1 = factor(PROMO_B1_1, levels = c(0, 1))
  )

# ---------- 4. drop rows with NA lags (first week) -----------------------------
sales_model <- sales %>% drop_na()

# ---------- 5. train / test split + Random-Forest ------------------------------
set.seed(123)
n <- nrow(data_rf)
train_idx <- 1:round(n * 0.8)
train <- sales_model[train_idx, ]
test <- sales_model[-train_idx, ]

rf_fit <- ranger(
  QTY_B1_1 ~
    PROMO_B1_1 +
      dow +
      month +
      lag1_QTY_B1_1 +
      lag7_QTY_B1_1 +
      total_B1 +
      roll7_totalB1,
  data = train,
  num.trees = 500,
  importance = "impurity" # Gini importance
)

# ---------- 6. evaluate --------------------------------------------------------
pred <- predict(rf_fit, data = test)$predictions
rmse <- sqrt(mean((pred - test$QTY_B1_1)^2))
cat("Test RMSE:", round(rmse, 3), "\n")
cat("Variable importance:\n")
print(sort(rf_fit$variable.importance, decreasing = TRUE))

imp_df <- enframe(
  rf_fit$variable.importance,
  name = "Variable",
  value = "Gini"
) %>%
  arrange(desc(Gini)) # high to low

# Optional: keep only top N variables to avoid clutter
# imp_df <- imp_df %>% slice_max(Gini, n = 20)

ggplot(imp_df, aes(x = reorder(Variable, Gini), y = Gini)) +
  geom_col(fill = "dodgerblue4") +
  coord_flip() +
  labs(
    title = "Gini variable importance (random forest)",
    x = NULL,
    y = "Importance"
  ) +
  theme_minimal(base_size = 12)
# ggsave("rf_gini_importance.jpg", width = 10, height = 6, dpi = 300)
# -----------------------------------------------------------
# 7.  Time–series plot: Actual vs Predicted on the test set
# -----------------------------------------------------------

plot_df <- test %>% # rows already in chrono order
  mutate(
    Actual = QTY_B1_1,
    Predicted = pred
  ) %>%
  select(DATE, Actual, Predicted) %>%
  pivot_longer(-DATE, names_to = "Series", values_to = "Value")

ggplot(plot_df, aes(DATE, Value, colour = Series)) +
  geom_line(size = 0.9) +
  labs(
    title = "Random-Forest predictions vs. actual sales",
    x = NULL,
    y = "Units sold (Brand 1 Item 1)",
    colour = ""
  ) +
  scale_colour_manual(values = c(Actual = "grey", Predicted = "red")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# ggsave("ver2nr_rf_tuned_performance.jpg", width = 10, height = 6, dpi = 300)
