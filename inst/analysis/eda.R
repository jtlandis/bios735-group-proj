library(pastasales)
library(ggplot2)
library(tidyverse)
library(pheatmap)



# Exploratory Plots -------------------------------------------------------

# Plot Mean time series for each brand
data_eda <- pastasales::data_set_tidy %>%
  mutate(item_id = factor(paste0(brand, "_", item))) %>%
  mutate(
    brand = factor(paste0("Brand ", as.numeric(factor(brand))), levels = paste0("Brand ", 1:4))
  )

data_eda_brand <- data_eda %>%
  group_by(brand, DATE, .add = TRUE) %>%
  summarise(QTY_brand_mean = mean(QTY),
            PROMO_brand_prop = mean(PROMO),
            .groups = "drop")

data_eda_brand %>%
  ggplot(aes(x = DATE, y = QTY_brand_mean, color = brand))+
  geom_line(aes(alpha = PROMO_brand_prop), linewidth = 0.5) +  # color by promo status
  facet_wrap(vars(brand))+
  theme_bw(base_size = 16)+
  scale_x_date(date_breaks = "6 month", 
               date_labels = "%m/%y")+
  scale_y_log10(limits = c(0.25, 100)) +
  scale_color_brewer(palette = 'Set1')+
  scale_alpha(range = c(0.15, 1))+
  #scale_color_gradient(low = 'blue', high = 'orange')+
  ylab('Mean Sales (log-scale)')+
  xlab('Date')+
  theme(legend.position = "right",
        panel.grid.minor = element_blank())+
  labs(alpha = "Promotion")+
  guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  

### Make correlation heatmap across items, showing within and across brand variations
data_eda_items <- prepare_data_allitems(pastasales::data_set_raw)
brand <- data_eda_items$G
B <- data_eda_items$num_brands

# get correlation matrix
corr_data_items <- cor(data_eda_items$Y)

# make gaps to group by brand
brand_sizes <- table(brand)
gaps <- cumsum(brand_sizes)
gaps <- gaps[-length(gaps)]

# make col/rownames to only show 1 label per brand
custom_names <- rep("", length(brand))
for (b in 1:B) {
  inds_b <- which(brand == b)
  custom_names[round(mean(inds_b))] <- paste0("Brand ", b) # assign middle index for name
}
rownames(corr_data_items) <- custom_names
colnames(corr_data_items) <- custom_names

pheatmap(
  corr_data_items,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps,
  gaps_col = gaps,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_legend = FALSE,
  main = "Item-Item Correlations Grouped by Brand"
)