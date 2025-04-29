library(pastasales)
library(rlist)

data <- prepare_data_allitems(pastasales::data_set_raw)
brand <- data$G
B <- data$num_brands
n_items <- data$num_items

#### Run MCMC for Bayesian PAR and PVAR models on training data (80% of timepoints)
# May take around 20min to run (even using Rcpp)
# Set lag (q) equal to 5
set.seed(1234)
ptrain = 0.8
nMC <- 10000
q <- 5
trained_models <- train_bayes_pvar(data, ptrain = ptrain, model_type = 'both', q = q, mcmc_iter = nMC, verbose = TRUE)
# list.save(trained_models, 'par_pvar_q5_p80_10000.rdata')
# trained_models <- list.load('par_pvar_q5_p80_10000.rdata')


#### Make Traceplots corresponding to random item in each brand
make_traceplot <- function(MCMC_df, param_name = "beta1", param_label = beta[1], model = "PAR Model", alpha = 0.75) {
  param_label_expr <- substitute(param_label)
  MCMC_df %>%
    ggplot(aes(x = MCMC_iter, y = !!sym(param_name), group = item))+
    geom_line(aes(color = brand), alpha = alpha)+
    scale_color_brewer(palette = 'Set1')+
    theme_bw(base_size=16)+
    ggtitle(
      substitute("Traceplot for " * x * " for random items in each brand (" * m * ")", 
                 list(x = param_label_expr, m = model))
    ) +
    ylab(param_label_expr) +
    xlab('MCMC Iteration')+
    labs(color = "Brand")
}

set.seed(12)
nitems_plot <- 2
rand_item_idxs <- sapply(1:B, function(b) which(brand == b) %>% sample(nitems_plot, replace = FALSE))
brand_ind <- brand[rand_item_idxs]
n_items_sel <- length(rand_item_idxs)

burn_in <- 0.5
ss_idxs <- (round(burn_in * nMC)):nMC


# Get post burn-in samples from PAR model for all items 
# PAR results are saved as list across items, each element is a list including samples for each model parameter
PAR_Beta_samps <- simplify2array(lapply(trained_models$PAR$mcmc_samps, function(x) x$beta[ss_idxs,]))
PAR_Gamma_samps <- simplify2array(lapply(trained_models$PAR$mcmc_samps, function(x) x$gamma[ss_idxs,]))

# Put PAR MCMC samples for all items in dataframe
MCMC_df_PAR <- data.frame()
for (idx in 1:n_items_sel) {
  i <- rand_item_idxs[idx]
  b <- brand_ind[idx]
  
  beta_df <- as.data.frame(PAR_Beta_samps[ , , i])
  colnames(beta_df) <- paste0("beta", 1:ncol(beta_df))
  
  gamma_df <- as.data.frame(PAR_Gamma_samps[ , , i])
  colnames(gamma_df) <- paste0("gamma", c(1:ncol(gamma_df))-1)
  
  MCMC_df_b <- data.frame(
    MCMC_iter = ss_idxs,
    item = paste0("Item ", i),
    brand = paste0("Brand ", b),
    item_idx = rand_item_idxs[b],
    beta_df,
    gamma_df
  )
  MCMC_df_PAR <- rbind(MCMC_df_PAR, MCMC_df_b)
}

# Traceplots for PAR model
make_traceplot(MCMC_df_PAR, param_name = "beta1", param_label = beta[1], model = "PAR Model")
make_traceplot(MCMC_df_PAR, param_name = "beta2", param_label = beta[2], model = "PAR Model")
make_traceplot(MCMC_df_PAR, param_name = "gamma0", param_label = gamma[0], model = "PAR Model")
make_traceplot(MCMC_df_PAR, param_name = "gamma1", param_label = gamma[1], model = "PAR Model")


# Get post burn-in samples from PVAR model
PVAR_Beta_samps <- simplify2array(trained_models$PVAR$mcmc_samps$Beta_samples)[, , ss_idxs]
PVAR_Gamma_samps <- simplify2array(trained_models$PVAR$mcmc_samps$Gamma_samples)[, , ss_idxs]

# Put PVAR MCMC samples in dataframe
MCMC_df_PVAR <- data.frame()
for (idx in 1:n_items_sel) {
  i <- rand_item_idxs[idx]
  b <- brand_ind[idx]
  
  beta_df <- as.data.frame(t(PVAR_Beta_samps[i , , ])) # different indexing from PAR model
  colnames(beta_df) <- paste0("beta", 1:ncol(beta_df))
  
  gamma_df <- as.data.frame(t(PVAR_Gamma_samps[i , , ]))
  colnames(gamma_df) <- paste0("gamma", c(1:ncol(gamma_df))-1)
  
  MCMC_df_b <- data.frame(
    MCMC_iter = ss_idxs,
    item = paste0("Item ", i),
    brand = paste0("Brand ", b),
    item_idx = rand_item_idxs[b],
    beta_df,
    gamma_df
  )
  MCMC_df_PVAR <- rbind(MCMC_df_PVAR, MCMC_df_b)
}

# Traceplots for PVAR model
make_traceplot(MCMC_df_PVAR, param_name = "beta1", param_label = beta[1], model = "PVAR Model")
make_traceplot(MCMC_df_PVAR, param_name = "beta2", param_label = beta[2], model = "PVAR Model")
make_traceplot(MCMC_df_PVAR, param_name = "gamma0", param_label = gamma[0], model = "PVAR Model")
make_traceplot(MCMC_df_PVAR, param_name = "gamma1", param_label = gamma[1], model = "PVAR Model")



#### Plot estimates for all items
make_estplot <- function(Est_df, param_name = "beta1", param_label = beta[1]) {
  param_label_expr <- substitute(param_label)
  Est_df %>%
    ggplot(aes(x = brand, y = !!sym(param_name), fill = brand))+
    geom_boxplot()+
    scale_fill_brewer(palette = 'Set1')+
    facet_wrap(vars(model))+
    theme_bw(base_size = 16)+
    ggtitle(
      substitute("Estimates for " * x * " for all items",
                 list(x = param_label_expr))
    )+
    ylab(param_label_expr)+
    labs(fill = "Brand")+
    theme(axis.title.x = element_blank())
}

PAR_Beta <- trained_models$PAR$Beta
PAR_tau <- rowSums(PAR_Beta)
PAR_Gamma <- trained_models$PAR$Gamma
PVAR_Beta <- trained_models$PVAR$Beta
PVAR_tau <- rowSums(PVAR_Beta)
PVAR_Gamma <- trained_models$PVAR$Gamma

# Put estimates in dataframe
PAR_df <- data.frame(model = "PAR", item = 1:n_items, brand = paste0("Brand ", brand), beta = PAR_Beta, gamma = PAR_Gamma, tau = PAR_tau)
PVAR_df <- data.frame(model = "PVAR", item = 1:n_items, brand = paste0("Brand ",brand), beta = PVAR_Beta, gamma = PVAR_Gamma, tau = PVAR_tau)
Est_df <- rbind(PAR_df, PVAR_df) %>%
  mutate(item = factor(item),
         brand = factor(brand))
names(Est_df) <- c("model", "item", "brand", paste0("beta", 1:q), paste0("gamma",0:1), "tau")

# Make plot for AR Terms
Est_AR_df <- Est_df %>% pivot_longer(cols = starts_with("beta"), names_to = "term", values_to = "est")
Est_AR_df %>%
  ggplot(aes(x = brand, y = est, fill = brand, group = interaction(brand, term)))+
  geom_boxplot(position = position_dodge(width = 0.8))+
  facet_wrap(vars(model))+
  theme_bw(base_size = 16)+
  scale_fill_brewer(palette = 'Set1')+
  ggtitle("Estimates for AR(5) Coefficients for all items")+
  ylab('AR Coefficients')+
  labs(fill = "Brand")+
  theme(axis.title.x = element_blank())
  
# Make plot for intercept and promotion effect estimates
make_estplot(Est_df, param_name = "gamma0", param_label = gamma[0])
make_estplot(Est_df, param_name = "gamma1", param_label = gamma[1])+ylim(c(0,10))

# Plot the post-burn-in samples for brand-level mean promo effect (PVAR)
Mu_samps <- simplify2array(trained_models$PVAR$mcmc_samps$mu_gamma_samples)[, , ss_idxs]
Mu_df <- data.frame()
for (b in 1:B) {
  Mu_b <- Mu_samps[b, , ]
  Mu_df_b <- data.frame(MCMC_iter = ss_idxs, brand = paste0("Brand ", b), 
                        mu0 = Mu_samps[b, 1, ], mu1 = Mu_samps[b, 2, ])
  Mu_df <- rbind(Mu_df, Mu_df_b)
}
Mu_df$brand <- factor(Mu_df$brand)
Mu_df %>%
  ggplot(aes(x = brand, y = mu1, fill = brand))+
  geom_boxplot()+
  theme_bw(base_size = 16)+
  scale_fill_brewer(palette = 'Set1')+
  ggtitle("PVAR MCMC Samples for brand-level mean promo effect")+
  ylab(expression(mu[1])) +
  labs(fill = "Brand")+
  theme(axis.title.x = element_blank())+
  ylim(c(0,10))




##### Prediction/Forecasting
n_timepoints <- nrow(data$Y)
train_cutoff <- round(n_timepoints * ptrain)

# Include q timepoints before cutoff in test set
Ytest <- data$Y[(train_cutoff - q):n_timepoints, ]
Xtest <- data$X[(train_cutoff - q):n_timepoints, ]
DATE_Test <- data_set_raw$DATE[(train_cutoff - q):n_timepoints]

# Get one step forecast (touching test set)
Mtest_1step_PAR <- get_fit_pvar(Ytest, Xtest, PAR_Beta, PAR_Gamma)
Mtest_1step_PVAR <- get_fit_pvar(Ytest, Xtest, PVAR_Beta, PVAR_Gamma)

# Get H step forecast (not touching test set)
Mtest_PAR_temp <- Mtest_1step_PAR
Mtest_PAR_temp[1:q,] <- Ytest[1:q,]
Mtest_PVAR_temp <- Mtest_1step_PVAR
Mtest_PVAR_temp[1:q,] <- Ytest[1:q,]

Mtest_Hstep_PAR <- get_fit_pvar(Mtest_PAR_temp, Xtest, PAR_Beta, PAR_Gamma)
Mtest_Hstep_PVAR <- get_fit_pvar(Mtest_PVAR_temp, Xtest, PVAR_Beta, PVAR_Gamma)

# Combine into df
Pred_df <- data.frame()
for (i in 1:n_items) {
  Pred_df_i <- data.frame(item = i, brand = brand[i], 
                          DATE = DATE_Test, QTY = Ytest[, i], PROMO = Xtest[, i],
                          QTY_1step_PAR = Mtest_1step_PAR[, i], QTY_1step_PVAR = Mtest_1step_PVAR[, i],
                          QTY_Hstep_PAR = Mtest_Hstep_PAR[, i], QTY_Hstep_PVAR = Mtest_Hstep_PVAR[, i])
  Pred_df <- rbind(Pred_df, Pred_df_i)
}
Pred_df$brand <- factor(paste0("Brand ", Pred_df$brand))
Pred_df$item <- factor(Pred_df$item)

# Convert to long format for 1-step and H-step
Pred_df_long_1step <- Pred_df %>%
  pivot_longer(cols = c("QTY", "QTY_1step_PAR", "QTY_1step_PVAR"),
               names_to = "Model", values_to = "Sales") %>%
  mutate(Forecast = factor(case_when(Model == "QTY" ~ "True Value",
                                     Model == "QTY_1step_PAR" ~ "PAR",
                                     Model == "QTY_1step_PVAR" ~ "PVAR")))

Pred_df_long_Hstep <- Pred_df %>%
  pivot_longer(cols = c("QTY", "QTY_Hstep_PAR", "QTY_Hstep_PVAR"),
               names_to = "Model", values_to = "Sales") %>%
  mutate(Forecast = factor(case_when(Model == "QTY" ~ "True Value",
                                     Model == "QTY_Hstep_PAR" ~ "PAR",
                                     Model == "QTY_Hstep_PVAR" ~ "PVAR")))


## Plot 1-step and H-step forecasts
plot_forecast <- function(df, item_idx = 1, step = 'One', brand = 1) {
  df %>%
    dplyr::filter(item == item_idx) %>%
    ggplot(aes(x = DATE, y = Sales, color = Forecast, alpha = Forecast))+
    geom_line(size=0.7)+
    scale_color_manual(values = c("red","blue","black"))+
    scale_alpha_manual(values = c(1,1,0.5))+
    theme_bw(base_size=16)+
    ggtitle(paste0(step,'-Step Forecasts on Test Set for Brand ',brand,' Item'))+
    labs(color = "Type", alpha = "Type")
}

# brand 1 item 1 (1-step)
plot_forecast(Pred_df_long_1step, item_idx = 1, step = 'One', brand = 1)

# brand 1 item 1 (H-step)
plot_forecast(Pred_df_long_Hstep, item_idx = 1, step = 'H', brand = 1)

# brand 4 item 1/idx = 109 (1-step)
plot_forecast(Pred_df_long_1step, item_idx = 109, step = 'One', brand = 4)

#brand 4 item 1/idx = 109 (H-step)
plot_forecast(Pred_df_long_Hstep, item_idx = 109, step = 'H', brand = 4)


## Plot brand-level mean forecasts
Pred_df_avg <- Pred_df_long %>%
  group_by(DATE, brand, Forecast) %>%
  summarise(Sales = mean(Sales), .groups = 'drop')

Pred_df_avg %>%
  ggplot(aes(x = DATE, y = Sales, color = Forecast, alpha = Forecast)) +
  geom_line(size = 0.7) +
  facet_wrap(vars(brand)) +
  scale_color_manual(values = c("red", "blue", "black")) +
  scale_alpha_manual(values = c(1,1,0.5))+
  theme_bw(base_size = 16) +
  labs(color = "Type", alpha = "Type")+
  scale_y_log10()+
  ggtitle('Brand-Level Mean One-Step Forecasts on Test Set (log-scale)')
  
# Get deviance and MSE for all items
get_deviance <- function(y, mu) {
  log_term <- y * log(y / mu)
  log_term[y == 0] <- 0
  D <- 2 * sum(log_term - (y - mu), na.rm=T)
  D
}
get_mse <- function(y, mu) {
  n_nona <- sum(!is.na(mu))
  sum((y-mu)^2, na.rm=T)/n_nona
}

Dev_1step_PAR <- sapply(1:n_items, function(i) get_deviance(Ytest[,i], Mtest_1step_PAR[,i]))
Dev_1step_PVAR <- sapply(1:n_items, function(i) get_deviance(Ytest[,i], Mtest_1step_PVAR[,i]))
Dev_Hstep_PAR <- sapply(1:n_items, function(i) get_deviance(Ytest[,i], Mtest_Hstep_PAR[,i]))
Dev_Hstep_PVAR <- sapply(1:n_items, function(i) get_deviance(Ytest[,i], Mtest_Hstep_PVAR[,i]))

MSE_1step_PAR <- sapply(1:n_items, function(i) get_mse(Ytest[,i], Mtest_1step_PAR[,i]))
MSE_1step_PVAR <- sapply(1:n_items, function(i) get_mse(Ytest[,i], Mtest_1step_PVAR[,i]))
MSE_Hstep_PAR <- sapply(1:n_items, function(i) get_mse(Ytest[,i], Mtest_Hstep_PAR[,i]))
MSE_Hstep_PVAR <- sapply(1:n_items, function(i) get_mse(Ytest[,i], Mtest_Hstep_PVAR[,i]))


# put in df
MSE_df <- data.frame(item = 1:n_items, brand = factor(paste0("Brand ", brand)),
                     MSE_1step_PAR = MSE_1step_PAR, MSE_1step_PVAR = MSE_1step_PVAR,
                     MSE_Hstep_PAR = MSE_Hstep_PAR, MSE_Hstep_PVAR = MSE_Hstep_PVAR)

MSE_df_long_1step <- MSE_df %>% 
  pivot_longer(cols = c("MSE_1step_PAR", "MSE_1step_PVAR"), names_to = "Model", values_to = "MSE") %>%
  mutate(Model = factor(case_when(Model == "MSE_1step_PAR" ~ "PAR",
                                  TRUE ~ "PVAR")))
MSE_df_long_Hstep <- MSE_df %>% pivot_longer(cols = c("MSE_Hstep_PAR", "MSE_Hstep_PVAR"), names_to = "Model", values_to = "MSE") %>%
  mutate(Model = factor(case_when(Model == "MSE_Hstep_PAR" ~ "PAR",
                                  TRUE ~ "PVAR")))

# plot PAR/PVAR MSE for all items for 1-step case and H-step case
MSE_df_long_1step %>%
  ggplot(aes(x = brand, y = MSE, fill = brand))+
  facet_wrap(vars(Model))+
  geom_boxplot()+
  scale_fill_brewer(palette = 'Set1')+
  theme_bw(base_size = 16) +
  scale_y_log10()+
  labs(fill = 'Brand')+
  ggtitle('MSE for One-step forecasts across items (log-scale)')+
  xlab(NULL)


MSE_df_long_Hstep %>%
  ggplot(aes(x = brand, y = MSE, fill = brand))+
  facet_wrap(vars(Model))+
  geom_boxplot()+
  scale_fill_brewer(palette = 'Set1')+
  theme_bw(base_size = 16) +
  scale_y_log10()+
  labs(fill = 'Brand')+
  ggtitle('MSE for H-step forecasts across items (log-scale)')+
  xlab(NULL)

