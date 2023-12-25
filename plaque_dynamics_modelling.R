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
library(FSA)


final_df <- read.csv("tidy_results.csv")

#Optimising parameter combinations ####
#WARNING: Moderately computationally intensive process. Check, and consider carefully, the number of combinations generated in param_grid

# Requires a df callled manual_t0s with the columns 'Plaque_ID', 'reporter', 'Start_hour'.
# manual_t0 values should be set based on manual observation of the timepoint at which a plaque spreads from a single infected cell
# to multiple infected cells. Take manual observations for a subset of your data (at least 5 plaques per condition) and assess the accuracy
# of the automated process vs your manual observations

# Define your specific cell type
cell_type_to_optimize <- "BSC"  # replace with the cell type you want to optimize for

# Filter the data for the specific cell type
filtered_final_df <- final_df %>% filter(cell_type == cell_type_to_optimize)
filtered_final_df$Plaque_ID <- as.character(filtered_final_df$Plaque_ID)
manual_t0s$Plaque_ID <- as.character(manual_t0s$Plaque_ID)
# Define the ranges for each variable (these should be adjusted to your needs)
window_size_range <- seq(5, 50, by = 5)
window_size_2_range <- seq(5, 40, by = 5)
factor_threshold_range <- seq(1.1, 2, by = 0.1)
min_increase_range <- seq(0, 60, by = 5)

# Create a grid of all possible combinations of these variables
param_grid <- expand.grid(
  window_size = window_size_range,
  window_size_2 = window_size_2_range,
  factor_threshold = factor_threshold_range,
  min_increase = min_increase_range
)


# Initialize a data frame to store results
results <- data.frame(
  window_size = numeric(),
  window_size_2 = numeric(),
  factor_threshold = numeric(),
  min_increase = numeric(),
  reporter = character(),
  Plaque_ID = numeric(),
  sum_accuracy = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each parameter set in param_grid
for(i in 1:nrow(param_grid)) {
  # Extract parameters for this iteration
  window_size <- param_grid[i, "window_size"]
  window_size_2 <- param_grid[i, "window_size_2"]
  factor_threshold <- param_grid[i, "factor_threshold"]
  min_increase <- param_grid[i, "min_increase"]
  
  # The list to collect accuracies for each reporter and plaque combination
  t0_results <- data.frame(Plaque_ID = integer(), reporter = character(), t0 = numeric(), stringsAsFactors = FALSE)
  
  
  # Loop over each Plaque_ID 
  for(plaque_id in unique(filtered_final_df$Plaque_ID)) {
    for(reporter_type in c('early_front', 'late_front')) {
      # Filter for current Plaque_ID and reporter
      plaque_data <- filtered_final_df %>%
        filter(Plaque_ID == plaque_id, reporter == reporter_type) %>%
        arrange(time_hours)
      
      # Check if the data frame is empty after filtering
      if(nrow(plaque_data) == 0) {
        warning(sprintf("No data for Plaque_ID '%s' with reporter '%s'.", plaque_id, reporter_type))
        next # Skip this iteration
      }
      
      # Check for NA or NaN values in time_hours
      if(any(is.na(plaque_data$time_hours) | is.nan(plaque_data$time_hours))) {
        warning(sprintf("NA or NaN values in 'time_hours' for Plaque_ID '%s' with reporter '%s'.", plaque_id, reporter_type))
        next # Skip this iteration
      }
      # Ensure time_hours is in increasing order
      if(!is.strictly.increasing(plaque_data$time_hours)) {
        warning(sprintf("'time_hours' is not strictly increasing for Plaque_ID '%s' with reporter '%s'.", plaque_id, reporter_type))
        next # Skip this iteration
      }
      # Add NA at the start to align the length of First_Derivative with the original data
      plaque_data$First_Derivative <- c(NA, diff(plaque_data$rolled_wavefront_micron) / diff(plaque_data$time_hours))
      
      # Apply LOESS smoothing to the first derivative (excluding the first NA value)
      loess_fit <- loess(First_Derivative ~ time_hours, data = plaque_data[-1, ], span = 0.1)
      
      # Pad the predicted values with NA at the start to align with the original data
      plaque_data$Smoothed_First_Derivative_LOESS <- c(NA, predict(loess_fit))
      
      # Remove NA values from the Smoothed_First_Derivative_LOESS column only
      non_na_smoothed_values <- plaque_data$Smoothed_First_Derivative_LOESS[-1]
      non_na_times <- plaque_data$time_hours[-1]
      
      # Proceed with changepoint detection on non-NA data
      cpt <- cpt.meanvar(non_na_smoothed_values, method = "PELT", penalty="Hannan-Quinn", test.stat = "Normal", minseglen = 5)
     
      # Settings for the analysis
      valid_changepoints <- sapply(cpt@cpts, function(cp) {
        # Ensure we don't go before the start or beyond the end of our data
        start_index <- max(cp - window_size, 1)
        end_index <- min(cp + window_size, nrow(plaque_data))
        end_index_2 <- min(cp + window_size_2, nrow(plaque_data))
        # Calculate means before and after the changepoint
        pre_cp_mean <- mean(plaque_data$rolled_wavefront_micron[start_index:cp], na.rm = TRUE)
        post_cp_mean <- mean(plaque_data$rolled_wavefront_micron[(cp + 1):end_index], na.rm = TRUE)
        
        # Check the conditions - Uncomment print() calls if you need to check function outputs mid-loop.
        is_increase <- post_cp_mean > pre_cp_mean * factor_threshold
        #print(paste("increases:", toString(is_increase)))
        is_significant_increase <- (post_cp_mean - pre_cp_mean) > min_increase
        #print(paste("significant increases:", toString(is_significant_increase)))
        no_negative_derivative <- !any(plaque_data$Smoothed_First_Derivative_LOESS[cp:end_index_2] < 0)
        #print(paste("positive leading first derivative:", toString(no_negative_derivative)))
        
        # Only return TRUE if all conditions are met
        is_valid_cp <- is_increase & is_significant_increase & no_negative_derivative
        #print(paste("valid cpt:", toString(is_valid_cp)))
        
        return(is_valid_cp)
      })
      
      # Retain changepoints where the condition is TRUE (no negative derivatives in the window and a significant mean increase)
      valid_changepoint_indexes <- cpt@cpts[valid_changepoints] +1
      # Extract the time (in hours) at which each changepoint occurs
      valid_changepoint_times <- plaque_data$time_hours[valid_changepoint_indexes]
      
      # If there are valid changepoints, proceed
      if(length(valid_changepoint_times) > 0) {
        t0 <- min(valid_changepoint_times)
        t0_results <- rbind(t0_results, data.frame(Plaque_ID = plaque_id, reporter = reporter_type, t0 = t0, valid_cpt_found = TRUE))
      } else {
        # No valid changepoints found, so append NA
        t0_results <- rbind(t0_results, data.frame(Plaque_ID = plaque_id, reporter = reporter_type, t0 = NA, valid_cpt_found = FALSE))
      }
    }
  }
  
  # Now check if all combinations had valid changepoints
  if(all(t0_results$valid_cpt_found)) {
    # Calculate the accuracy for each Plaque_ID and reporter
    t0_accuracy <- t0_results %>%
      left_join(manual_t0s, by = c("Plaque_ID", "reporter")) %>%
      mutate(accuracy = if_else(is.na(t0), NA_real_, abs(t0 - Start_hour)))
    
    # Calculate the summary of t0_accuracy for this parameter set
    sum_accuracy <- sum(t0_accuracy$accuracy, na.rm = TRUE)
  } else {
    # At least one combination had no valid changepoints, set sum_accuracy to NA
    sum_accuracy <- NA
  }
  new_row <- data.frame(
    window_size = window_size, 
    window_size_2 = window_size_2, 
    factor_threshold = factor_threshold, 
    min_increase = min_increase, 
    sum_accuracy = sum_accuracy
  )
  results <- rbind(results, new_row)
}
# Filter out the NA sum_accuracy rows and order the results by sum_accuracy to find the best set
results <- results[!is.na(results$sum_accuracy),]
BSC_results <- results[order(results$sum_accuracy),]

# The first row now contains the parameter set with the lowest (i.e. least deviation from manual_t0 values) sum_accuracy
optimal_parameters <- results[1,]


#Modelling ####

# Settings for the changepoint analysis of each cell_type
#Window size determines the size of the region around each changepoint in which to search for an increase in mean
#window size 2 determiens the size of the region leading each changepoint in which to search for a negative first derivative
#factor threshold determines how much larger the leading mean needs to be than the trailing mean for a valid changepoint
#min increase determines how much absolute difference there must be between the leading mean and the trailing mean for a valid changepoint
settings <- list(
  "HA_CAT" = list(window_size = 10, window_size_2 = 5, factor_threshold = 1.1, min_increase = 5),
  "BSC" = list(window_size = 25, window_size_2 = 5, factor_threshold = 1.1, min_increase = 50)
)

# Initialize a data frame to store the results
results_df <- data.frame(Plaque_ID = character(),
                         reporter = character(),
                         t0 = numeric(),
                         growth_rate = numeric(),
                         r_squared = numeric(),
                         starting_size = numeric(),
                         cell_type = character(),
                         stringsAsFactors = FALSE)

# Loop over unique combinations of Plaque_ID and reporter
for(plaque_id in unique(final_df$Plaque_ID)) {
  for(reporter_type in c('early_front', 'late_front')) {
    # Filter for current Plaque_ID and reporter
    plaque_data <- final_df %>%
      dplyr::filter(Plaque_ID == plaque_id, reporter == reporter_type) %>%
      arrange(time_hours)
    
    # Determine which settings to use based on the cell_type in plaque_data
    current_cell_type <- unique(plaque_data$cell_type)
    if (current_cell_type %in% names(settings)) {
      # Use settings specific to the current cell_type
      current_settings <- settings[[current_cell_type]]
    } else {
      next # Skip this iteration if the cell_type is not recognized
    }
    
    # Apply the settings
    window_size <- current_settings$window_size
    window_size_2 <- current_settings$window_size_2
    factor_threshold <- current_settings$factor_threshold
    min_increase <- current_settings$min_increase
    # Add NA at the start to align the length of First_Derivative with the original data
    plaque_data$First_Derivative <- c(NA, diff(plaque_data$rolled_wavefront_micron) / diff(plaque_data$time_hours))
    
    # Apply LOESS smoothing to the first derivative (excluding the first NA value)
    loess_fit <- loess(First_Derivative ~ time_hours, data = plaque_data[-1, ], span = 0.1)
    
    # Pad the predicted values with NA at the start to align with the original data
    plaque_data$Smoothed_First_Derivative_LOESS <- c(NA, predict(loess_fit))
    
    # Remove NA values from the Smoothed_First_Derivative_LOESS column only
    non_na_smoothed_values <- plaque_data$Smoothed_First_Derivative_LOESS[-1]
    non_na_times <- plaque_data$time_hours[-1]
    # Proceed with changepoint detection on non-NA data
    cpt <- cpt.meanvar(non_na_smoothed_values, method = "PELT", penalty="Hannan-Quinn", test.stat = "Normal", minseglen = 5)
    # Settings for the analysis
    
    valid_changepoints <- sapply(cpt@cpts, function(cp) {
      # Ensure we don't go before the start or beyond the end of our data
      start_index <- max(cp - window_size, 1)
      end_index <- min(cp + window_size, nrow(plaque_data))
      end_index_2 <- min(cp + window_size_2, nrow(plaque_data))
      # Calculate means before and after the changepoint
      pre_cp_mean <- mean(plaque_data$rolled_wavefront_micron[start_index:cp], na.rm = TRUE)
      post_cp_mean <- mean(plaque_data$rolled_wavefront_micron[(cp + 1):end_index], na.rm = TRUE)
      
      # Check the conditions
      is_increase <- post_cp_mean > pre_cp_mean * factor_threshold
      #print(paste("increases:", toString(is_increase)))
      is_significant_increase <- (post_cp_mean - pre_cp_mean) > min_increase
      #print(paste("significant increases:", toString(is_significant_increase)))
      no_negative_derivative <- !any(plaque_data$Smoothed_First_Derivative_LOESS[cp:end_index_2] < 0)
      #print(paste("positive leading first derivative:", toString(no_negative_derivative)))
      
      # Only return TRUE if all conditions are met
      is_valid_cp <- is_increase & is_significant_increase & no_negative_derivative
      #print(paste("valid cpt:", toString(is_valid_cp)))
      
      return(is_valid_cp)
    })
    
    # Now, we only keep the changepoints where the condition is TRUE (no negative derivatives in the window and a significant mean increase)
    valid_changepoint_indexes <- cpt@cpts[valid_changepoints] +1
    # Now you can extract the time for the valid changepoints if needed
    valid_changepoint_times <- plaque_data$time_hours[valid_changepoint_indexes]
    
    # If there are valid changepoints, proceed
    if(length(valid_changepoint_times) > 0) {
      # Take the earliest changepoint time
      t0 <- min(valid_changepoint_times)
      
      # Filter data for time points after t0 and fit linear model
      lm_data <- plaque_data %>%
        dplyr::filter(time_hours >= t0&time_hours<=45)
      lm_fit <- lm(rolled_wavefront_micron ~ time_hours, data = lm_data)
      
      # Extract metrics
      growth_rate <- coef(lm_fit)["time_hours"]
      r_squared <- summary(lm_fit)$r.squared
      starting_size <- median(plaque_data$rolled_wavefront_micron[plaque_data$time_hours < t0])
      cell_type <- unique(plaque_data$cell_type)
      
      # Append to results dataframe
      results_df <- rbind(results_df, data.frame(Plaque_ID = plaque_id,
                                                 reporter = reporter_type,
                                                 t0 = t0,
                                                 growth_rate = growth_rate,
                                                 r_squared = r_squared,
                                                 starting_size = starting_size,
                                                 cell_type = cell_type))
    }else { # Warning message if no valid changepoints are found
      warning(sprintf("No valid changepoints found for Plaque_ID '%s' with reporter '%s'.", plaque_id, reporter_type))
    }
  }
}

results_df_summary <- results_df %>% 
  group_by(cell_type, reporter) %>% 
  dplyr::summarise(mean_t0=mean(t0, na.rm=TRUE),
                   mean_t0_initial=mean(t0, na.rm=TRUE),
                   mean_growthRate=mean(growth_rate, na.rm=TRUE),
                   sd_t0=sd(t0, na.rm=TRUE),
                   sd_t0_initial=sd(t0, na.rm=TRUE),
                   sd_growthRate=sd(growth_rate, na.rm=TRUE))

manual_t0s<- read.csv("manual_t0s.csv") %>% 
  dplyr::select(Plaque_ID, reporter, Start_hour)

t0_accuracy <- results_df %>% 
  left_join(by=c("Plaque_ID", "reporter"), manual_t0s) %>% 
  dplyr::mutate(accuracy=t0-Start_hour)

t0_accuracy_summary <- t0_accuracy %>%
  group_by(reporter, cell_type) %>% 
  dplyr::summarise(mean_accuracy=mean(accuracy),
                   sum_accuracy=sum(abs(accuracy)),
                   mean_startHours=mean(Start_hour),
                   mean_t0=mean(t0),
                   sd_accuracy=sd(accuracy),
                   sd_startHours=sd(Start_hour),
                   sd_t0=sd(t0))

write.csv(results_df, "modeling_results.csv")
write.csv(t0_accuracy, "modeling_accuracy.csv")
