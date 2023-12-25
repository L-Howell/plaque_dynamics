#libraries ####
library(svDialogs)
library(ggplot2)
library(tidyr)
library(plyr)
library(gganimate)
library(transformr)
library(dplyr)
library(tidyverse)
library(broom)
library(magick)
library(patchwork)
library(lubridate)
library(zoo)
library(extrafont)
library(showtext)
library(rstatix)
library(changepoint)
library(ggsignif)
library(purrr)
library(RColorBrewer)
library(FSA)
library(rstatix)


#Statistics ####

#Figure 2 t0 for each reporter
# Conduct the pairwise t-test
test_result <- results_df %>% 
  filter(cell_type=="BSC") %>% 
  pairwise_t_test(formula = t0 ~ reporter,
                  p.adjust.method = "none")  # Adjust this if you want p-value adjustment for multiple comparisons
# Print the result
print(test_result)

test_result <- results_df %>% 
  filter(reporter=="late_front") %>% 
  pairwise_t_test(formula = t0 ~ cell_type,
                  p.adjust.method = "none")  # Adjust this if you want p-value adjustment for multiple comparisons
print(test_result)

test_result <- results_df %>% 
  filter(reporter=="early_front") %>% 
  pairwise_t_test(formula = growth_rate ~ cell_type,
                  p.adjust.method = "none")  # Adjust this if you want p-value adjustment for multiple comparisons
print(test_result)

test_result <- results_df %>% 
  filter(reporter=="late_front") %>% 
  pairwise_t_test(formula = growth_rate ~ cell_type,
                  p.adjust.method = "none")  # Adjust this if you want p-value adjustment for multiple comparisons

# Print the result
print(test_result)

#Comparing growth rate between cell types

# Subset the data for 'early_front' and 'late_front'
early_front_df <- subset(results_df, reporter == "early_front")
late_front_df <- subset(results_df, reporter == "late_front")
kruskal_late <- kruskal.test(t0 ~ cell_type, data = late_front_df)
kruskal_early <- kruskal.test(t0 ~ cell_type, data = early_front_df)
# Dunn test for 'early_front' across cell types
dunn_early_front <- dunn.test(x=early_front_df$growth_rate, 
                              g=early_front_df$cell_type,
                              method="bonferroni")

# Dunn test for 'late_front' across cell types
dunn_late_front <- dunn.test(x=late_front_df$growth_rate, 
                             g=late_front_df$cell_type,
                             method="bonferroni")

# Now, for the comparisons within each cell type for different reporters
bsc_df <- subset(results_df, cell_type == "BSC")
ha_cat_df <- subset(results_df, cell_type == "HA_CAT")

# Dunn test within 'BSC' for different reporters
dunn_bsc <- dunn.test(x=bsc_df$growth_rate, 
                      g=bsc_df$reporter,
                      method="bonferroni")

# Dunn test within 'HA_CAT' for different reporters
dunn_ha_cat <- dunn.test(x=ha_cat_df$growth_rate, 
                         g=ha_cat_df$reporter,
                         method="bonferroni")

# Print out the results
print(dunn_early_front)
print(dunn_late_front)
print(dunn_bsc)
print(dunn_ha_cat)


#Comparing t0 between cell types

# Subset the data for 'early_front' and 'late_front'
early_front_df <- subset(results_df, reporter == "early_front")
late_front_df <- subset(results_df, reporter == "late_front")
kruskal_late <- kruskal.test(t0 ~ cell_type, data = late_front_df)
kruskal_early <- kruskal.test(t0 ~ cell_type, data = early_front_df)
# Dunn test for 'early_front' across cell types
dunn_early_front <- dunn.test(x=early_front_df$t0, 
                              g=early_front_df$cell_type,
                              method="bonferroni")

# Dunn test for 'late_front' across cell types
dunn_late_front <- dunn.test(x=late_front_df$t0, 
                             g=late_front_df$cell_type,
                             method="bonferroni")

# Now, for the comparisons within each cell type for different reporters
bsc_df <- subset(results_df, cell_type == "BSC")
ha_cat_df <- subset(results_df, cell_type == "HA_CAT")

# Dunn test within 'BSC' for different reporters
dunn_bsc <- dunn.test(x=bsc_df$t0, 
                      g=bsc_df$reporter,
                      method="bonferroni")

# Dunn test within 'HA_CAT' for different reporters
dunn_ha_cat <- dunn.test(x=ha_cat_df$t0, 
                         g=ha_cat_df$reporter,
                         method="bonferroni")

# Print out the results
print(dunn_early_front)
print(dunn_late_front)
print(dunn_bsc)
print(dunn_ha_cat)


#Multivariate regression model ####
# Filter the data for the BSC cell type
bsc_data <- variable_importance %>% 
  dplyr::filter(cell_type == "BSC")

# Fit a multiple linear regression model
bsc_model <- lm(rolled_wavefront_micron ~ t0 + growth_rate, data = bsc_data)

# Summarize the model to get the R^2 value and other statistics
bsc_model_summary <- summary(bsc_model)

# Print the summary
print(bsc_model_summary)

# To extract the R-squared value specifically
r_squared <- bsc_model_summary$r.squared
print(paste("R-squared for the BSC model is:", r_squared))

