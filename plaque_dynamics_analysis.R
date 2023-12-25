#RStudio v1.4.1103 ####
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
# Uses data drawn from the Radial_Profiles_and_Plots FIJI macro ####
# this script compares a single channel's data from many plaques 

# setting directories ####

#prompts user for input directory.
#creates an output directory based on the location of the input directory
input_directory <- choose.dir(default = getwd(), caption = "Select the folder that contains all your plaque data (e.g. plaque_1")
setwd(input_directory)
dir.create(file.path("R_output_channel_comparisons/"))
output_directory <- paste0(file.path(input_directory, "R_output_channel_comparisons/"))

#user input ####
#asks the user to provide a name for their dataset, and to the name of the channel to analyse.
#assumes "channel_1' is found in the names of the .csv's that contain channel 1 data.
name_1 <- "TXRED"
name_2 <- "FITC"
channel_1 <- "TXRED"
channel_2 <- "FITC"
#specify the time between frames, and the start time post infection
#required for converting[Time] units from slice# to HPI
#if the named objects already exist, the associated line is not executed.

interval <- 15

start_time <- 8

#creates a list of all .csv's in input_directory
setwd(input_directory)
temp <- list.files(path = input_directory, pattern = ".tiff_binary", recursive = T)

#filters "temp" to those that contain "channel_1" in their name.
filenames_1 <- Filter(function(x) grepl(channel_1, x), temp)
filenames_2 <- Filter(function(x) grepl(channel_2, x), temp)

mydata_1 <- data.frame(lapply(filenames_1[1], read.csv, colClasses=c("NULL","numeric", "numeric", "numeric")))
mydata_2 <- data.frame(lapply(filenames_2[1], read.csv, colClasses=c("NULL","numeric", "numeric", "numeric")))

#Reads in the contents of "mydata" and merges into a single dataframe.

mydata_assembled_1 <- do.call("cbind", lapply(filenames_1[-1], read.csv, colClasses=c("NULL","NULL","NULL", "numeric")))  
mydata_assembled_2 <- do.call("cbind", lapply(filenames_2[-1], read.csv, colClasses=c("NULL","NULL","NULL", "numeric")))  

mydata_assembled_for_real_1 <- cbind(mydata_1, mydata_assembled_1)
mydata_assembled_for_real_2 <- cbind(mydata_2, mydata_assembled_2)

columns <- ncol(mydata_assembled_for_real_1)

#renames columns. 

names(mydata_assembled_for_real_1)[1:(columns)] <- c("Time", "Location", paste0("Plaque_", 1:(columns-2)))
names(mydata_assembled_for_real_2)[1:(columns)] <- c("Time", "Location", paste0("Plaque_", 1:(columns-2)))
#converts Time column to HPI based on user inputted interval and start_time.

mydata_assembled_for_real_1 <- mydata_assembled_for_real_1 %>%
  mutate(time_in_hours= (Time*(interval/60))+(start_time-(interval/60)))
mydata_assembled_for_real_2 <- mydata_assembled_for_real_2 %>%
  mutate(time_in_hours= (Time*(interval/60))+(start_time-(interval/60)))

#pivots data to make plaque ID and intensity single columns.
data_to_plot_1 <- mydata_assembled_for_real_1 %>% 
  pivot_longer(-c("Time", "Location", "time_in_hours"), names_to = "Plaque_ID", values_to = "Intensity") %>% 
  group_by(Plaque_ID) %>%
  mutate(time_integer= as.integer(time_in_hours)) %>% #integer value time column for gganimate transition state mapping
  ungroup()%>% 
  group_by(Time, Plaque_ID) 
data_to_plot_2 <- mydata_assembled_for_real_2 %>% 
  pivot_longer(-c("Time", "Location", "time_in_hours"), names_to = "Plaque_ID", values_to = "Intensity") %>% 
  group_by(Plaque_ID) %>%
  mutate(time_integer= as.integer(time_in_hours)) %>% #integer value time column for gganimate transition state mapping
  ungroup()%>% 
  group_by(Time, Plaque_ID)

#Early reporter front calculation ####

wavefront_early <- data_to_plot_1 %>% 
  unite("filter", c("Plaque_ID","Time"), remove = FALSE) %>%
  #filter(Time<=80) %>% 
  #filter(filter != "Plaque_2_65") %>% #optional removal of specific slices of data
  #filter(filter != "Plaque_2_66") %>%
  #filter(filter != "Plaque_2_67") %>% 
  # filter(filter != "Plaque_2_68") %>% 
  # filter(filter != "Plaque_2_69") %>% 
  # filter(filter != "Plaque_2_70") %>% 
  group_by(Plaque_ID) %>%
  dplyr::select(-filter) %>% 
  filter(Intensity < 249) %>%
  ungroup() %>% 
  group_by(Time, Plaque_ID) %>% 
  slice_max(Location, n = 10, with_ties=FALSE)%>% 
  group_by(mean = (row_number()-1) %/% 10) %>%  #calculates the mean of every group of 4 rows
  mutate(plaque_front_location = median(Location)) %>% 
  slice(which(row_number() %% 10 == 1)) %>% 
  group_by(Plaque_ID) %>%
  dplyr::select(Plaque_ID, Time, plaque_front_location) 

# wavefront_late <- data_to_plot_1 %>% 
#     unite("filter", c("Plaque_ID","Time"), remove = FALSE) %>%
#     filter(Time>80) %>% 
#   #filter(filter != "Plaque_2_65") %>% #optional removal of specific slices of data
#   #filter(filter != "Plaque_2_66") %>%
#   #filter(filter != "Plaque_2_67") %>% 
#   # filter(filter != "Plaque_2_68") %>% 
#   # filter(filter != "Plaque_2_69") %>% 
#   # filter(filter != "Plaque_2_70") %>% 
#     group_by(Plaque_ID) %>%
#     select(-filter) %>% 
#     filter(Intensity < 250) %>%
#     ungroup() %>% 
#     group_by(Time, Plaque_ID) %>% 
#     slice_max(Location, n = 60, with_ties=FALSE)%>% 
#     group_by(mean = (row_number()-1) %/% 60) %>%  #calculates the mean of every group of 4 rows
#     mutate(plaque_front_location = mean(Location)) %>% 
#     slice(which(row_number() %% 60 == 1)) %>% 
#     group_by(Plaque_ID) %>%
#     select(Plaque_ID, Time, plaque_front_location) 

wavefront_Early <- wavefront_early %>% 
  #rbind(wavefront_early, wavefront_late)%>%
  group_by(Plaque_ID) %>% 
  filter(plaque_front_location-lag(plaque_front_location) < 30) %>%
  filter(plaque_front_location-lag(plaque_front_location, n=2) < 30) %>% #filters out errors introduced by imaging artifacts that jump the front far outwards
  dplyr::select(Plaque_ID, Time, plaque_front_location) %>% 
  group_by(Plaque_ID, Time) %>% 
  dplyr::rename(early_front = plaque_front_location)

#Late reporter front calculation ####
wavefront_early <- data_to_plot_2 %>% 
  unite("filter", c("Plaque_ID","Time"), remove = FALSE) %>%
  #filter(Time<=80) %>% 
  #filter(filter != "Plaque_2_65") %>% #optional removal of specific slices of data
  #filter(filter != "Plaque_2_66") %>%
  #filter(filter != "Plaque_2_67") %>% 
  # filter(filter != "Plaque_2_68") %>% 
  # filter(filter != "Plaque_2_69") %>% 
  # filter(filter != "Plaque_2_70") %>% 
  group_by(Plaque_ID) %>%
  dplyr::select(-filter) %>% 
  filter(Intensity < 252) %>%
  ungroup() %>% 
  group_by(Time, Plaque_ID) %>% 
  slice_max(Location, n = 10, with_ties=FALSE)%>% 
  group_by(mean = (row_number()-1) %/% 10) %>%  #calculates the mean of every group of 4 rows
  mutate(plaque_front_location = median(Location)) %>% 
  slice(which(row_number() %% 10 == 1)) %>% 
  group_by(Plaque_ID) %>%
  dplyr::select(Plaque_ID, Time, plaque_front_location) 

# wavefront_late <- data_to_plot_2 %>% 
#   unite("filter", c("Plaque_ID","Time"), remove = FALSE) %>%
#   filter(Time>80) %>% 
#   #filter(filter != "Plaque_2_65") %>% #optional removal of specific slices of data
#   #filter(filter != "Plaque_2_66") %>%
#   #filter(filter != "Plaque_2_67") %>% 
#   # filter(filter != "Plaque_2_68") %>% 
#   # filter(filter != "Plaque_2_69") %>% 
#   # filter(filter != "Plaque_2_70") %>% 
#   group_by(Plaque_ID) %>%
#   select(-filter) %>% 
#   filter(Intensity < 250) %>%
#   ungroup() %>% 
#   group_by(Time, Plaque_ID) %>% 
#   slice_max(Location, n = 60, with_ties=FALSE)%>% 
#   group_by(mean = (row_number()-1) %/% 60) %>%  #calculates the mean of every group of 4 rows
#   mutate(plaque_front_location = mean(Location)) %>% 
#   slice(which(row_number() %% 60 == 1)) %>% 
#   group_by(Plaque_ID) %>%
#   select(Plaque_ID, Time, plaque_front_location) 

wavefront_Late <- wavefront_early %>% 
  #rbind(wavefront_early, wavefront_late)%>%
  group_by(Plaque_ID) %>% 
  filter(plaque_front_location-lag(plaque_front_location) < 30) %>%
  filter(plaque_front_location-lag(plaque_front_location, n=2) < 30) %>% #filters out errors introduced by imaging artifacts that jump the front far outwards
  dplyr::select(Plaque_ID, Time, plaque_front_location) %>% 
  group_by(Plaque_ID, Time) %>% 
  dplyr::rename(late_front = plaque_front_location)

#merging DFs ####
final_df <- wavefront_Early %>%
  dplyr::select(Plaque_ID, Time, early_front) %>% 
  left_join(wavefront_Late %>% 
              dplyr::select(Plaque_ID, Time, late_front),
            by = c("Plaque_ID", "Time")) %>% 
  mutate(cell_type = case_when(
    Plaque_ID %in% paste0("Plaque_", 1:9) ~ "BSC",
    Plaque_ID %in% paste0("Plaque_", 10:19) ~ "HA_CAT",
    TRUE ~ NA_character_ # For any other cases 
  )) %>% 
  ungroup()
# Create a dataframe with each Plaque_ID and all Time values (1:289) and assign cell_type
expanded_df <- distinct(final_df, Plaque_ID) %>%
  rowwise() %>%
  mutate(cell_type = case_when(
    Plaque_ID %in% paste0("Plaque_", 1:9) ~ "BSC",
    Plaque_ID %in% paste0("Plaque_", 10:19) ~ "HA_CAT",
    TRUE ~ NA_character_  # For any other cases
  )) %>%
  do(data.frame(Plaque_ID = .$Plaque_ID, Time = 1:289, cell_type = .$cell_type)) %>%
  ungroup()

# Join with the original data
final_df <- expanded_df %>%
  left_join(final_df, by = c("Plaque_ID", "Time", "cell_type"))

# Pivot to create wavefront column
final_df <- final_df %>%
  pivot_longer(cols=4:5, names_to = "reporter", values_to="wavefront") %>% 
  group_by(cell_type, Plaque_ID, reporter) %>% 
  dplyr::mutate(rolled_wavefront = rollmean(wavefront, k=5, align="center",fill=NA))

# Fill in the rolled_wavefront column
final_df <- final_df %>%
  arrange(Plaque_ID, Time) %>%
  tidyr::fill(rolled_wavefront, .direction = "down") %>%
  tidyr::fill(rolled_wavefront, .direction = "up") %>%
  ungroup()

# Calculate other necessary columns
final_df <- final_df %>%
  dplyr::mutate(rolled_wavefront_micron = rolled_wavefront * 1.62,
                time_hours = ((Time * interval) / 60) + start_time) %>% 
  group_by(Plaque_ID,reporter) %>%
  #dplyr::mutate(rolled_wavefront_micron2 = rolled_wavefront_micron - min(rolled_wavefront_micron))  %>%
  dplyr::mutate(min_val = ifelse(reporter == "late_front", min(rolled_wavefront_micron, na.rm=TRUE), NA_real_)) %>%
  ungroup() %>% 
  tidyr::fill(min_val, .direction = "updown") %>%
  dplyr::mutate(rolled_wavefront_micron = rolled_wavefront_micron - min_val) %>%
  dplyr::select(-min_val)

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



#Statistics ####
library(rstatix)
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

#Plots ####

# multiple plaques plots ####
p <- final_df %>%
  filter(reporter=="early_front") %>% 
  filter(cell_type=="BSC") %>% 
  ggplot() +
  geom_path(aes(x=time_hours, y=rolled_wavefront_micron, colour = Plaque_ID), size=0.25) +
  labs(x = "Time (HPI)", y = "Plaque radius (µm)", colour="Plaque \nNumber") +
  scale_color_brewer(palette = "Set1",labels = labels)+
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))

ggsave("BSC_multiple_plaques_radius.pdf",width=55, height=45, units="mm", p)

p <- final_df %>%
  filter(reporter=="early_front") %>% 
  filter(cell_type=="BSC") %>% 
  ggplot() +
  geom_path(aes(x=time_hours, y=pi*rolled_wavefront_micron^2, colour = Plaque_ID), size=0.25) +
  labs(x = "Time (HPI)", y = expression("Plaque area (µm"^2*")"), colour="Plaque \nNumber") +
  scale_color_brewer(palette = "Set1",labels = labels)+
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"),
        legend.box.margin = margin(0, 0, 0, -13))

ggsave("BSC_multiple_plaques_area.pdf",width=65, height=45, units="mm", p)



#Means plots ####
#Font loading
showtext_auto()
font_add("Arial", "C:/Windows/Fonts/arial.ttf")
summary_df <- final_df %>% 
  group_by(time_hours, cell_type, reporter) %>% 
  dplyr::summarise(mean_rolled_wavefront_micron=mean(rolled_wavefront_micron, na.rm=TRUE))

p <- summary_df %>% 
  dplyr::filter(cell_type=="BSC") %>% 
  ggplot() +
  geom_line(aes(x=time_hours, y=pi*mean_rolled_wavefront_micron^2, group = reporter, colour=reporter)) +
  scale_color_manual(values=c("red3","green3"))+
  labs(x = "Time (HPI)", y = expression("Plaque Area (µm"^2*")")) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))

ggsave("BSC_mean_area.pdf",width=60, height=50, units="mm", p)


p <- summary_df %>% 
  dplyr::filter(cell_type=="BSC") %>% 
  ggplot() +
  geom_line(aes(x=time_hours, y=mean_rolled_wavefront_micron, group = reporter, colour=reporter)) +
  scale_color_manual(values=c("red3","green3"))+
  labs(x = "Time (HPI)", y = "Plaque edge distance \nfrom center (µm)") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))

ggsave("BSC_mean_location.pdf",width=60, height=50, units="mm", p)


p <- summary_df %>%
  ggplot(aes(x = time_hours, y = mean_rolled_wavefront_micron, color = reporter, linetype = cell_type)) + 
  geom_path(size = 0.5) +
  labs(title = "", y = "Mean plaque radius (µm)", x = "Time (HPI)") +
  scale_color_manual(values = c("early_front" = "red2", "late_front" = "green3"),
                     name = "Reporter",  # Setting the title for color legend
                     labels = c("early_front" = "Early", "late_front" = "Late")) +
  scale_linetype_discrete(name = "Cell Type",  # Setting the title for linetype legend
                          labels = c("BSC" = "BS-C-1", "HA_CAT" = "HaCaT")) +  # Setting custom linetype names
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    text = element_text(family = "Arial")
  )

ggsave("mean_wavefront_BSC_vs_HACAT.pdf",width=75, height=65, units="mm", p)
#Model fitting example plot for Fig. 2
# Specify the desired plaque and reporter
desired_plaque <- "Plaque_7"  #Plaque_9 looks good

# Get data for the desired plaque (both reporters)
subset_data_raw <- final_df[final_df$Plaque_ID == desired_plaque, ]

# Fetch the model fitting results for this plaque from results_df
subset_data_results <- results_df[results_df$Plaque_ID == desired_plaque, ]

# Plotting
p <- ggplot(subset_data_raw, aes(x = time_hours, y = rolled_wavefront_micron, color = reporter)) +
  geom_path(aes(color = reporter), size = 0.5) +  # Plot the raw data points
  scale_color_manual(values = c("early_front" = "red2", "late_front" = "green3")) +
  geom_vline(data = subset_data_results, aes(xintercept = t0, color = reporter), linetype = "dashed") +
  geom_smooth(data = subset(subset_data_raw, time_hours >= min(subset_data_results$t0) & time_hours <= 60), 
              method = "lm", se = FALSE, aes(group = reporter), color = "black", linetype = "dashed", fullrange = TRUE, size=0.5) +
  #geom_hline(data = subset_data_results, aes(yintercept = starting_size),color = "purple3", linetype = "dashed")+
  labs(title = "",
       x = "Time (HPI)", 
       y = "Plaque radius (µm)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_text(size = 9),
        text = element_text(family = "Arial"))+
  coord_cartesian(xlim=c(-5,80))

ggsave("model_fitting.pdf",width=90, height=65, units="mm", p)

#Velocity plots####
t0_accuracy_join <- t0_accuracy %>% 
  filter(reporter=="early_front") %>% 
  dplyr::select(t0, Plaque_ID)
calculate_smoothed_derivative <- function(df) {
  df <- df %>% arrange(time_hours)
  df$First_Derivative <- c(NA, diff(df$rolled_wavefront_micron) / diff(df$time_hours))
  loess_fit <- loess(First_Derivative ~ time_hours, data = df[-1, ], span = 0.2)
  df$Smoothed_First_Derivative_LOESS <- c(NA, predict(loess_fit))
  return(df)
}
plaque_derivatives <- final_df %>%
  filter(cell_type=="BSC") %>% 
  filter(reporter == "early_front") %>%
  group_by(Plaque_ID) %>%
  do(calculate_smoothed_derivative(.)) %>% 
  left_join(by="Plaque_ID", t0_accuracy_join) %>% 
  group_by(Plaque_ID) %>% 
  filter(time_hours >= t0) #i.e. once the infection is spreading how does early_front velocity behave

number_of_plaques <- length(unique(plaque_derivatives$Plaque_ID))
labels <- as.character(1:number_of_plaques)  # Creates character vector c("1", "2", ..., "20")

p <- ggplot(plaque_derivatives, aes(x=time_hours, y=Smoothed_First_Derivative_LOESS, colour = Plaque_ID)) +
  geom_path(size=0.25) +
  labs(x = "Time (HPI)", y = "Velocity (µm/hr)", colour="Plaque \nNumber") +
  scale_color_brewer(palette = "Set1",labels = labels)+
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"),
        legend.box.margin = margin(0, 0, 0, -13))

ggsave("BSC_velocity_all_plaques.pdf",width=65, height=45, units="mm", p)

#EL delay plot 1 ####
# Step 1: Separate the data into early and late fronts
early_front <- final_df %>% 
  #dplyr::filter(cell_type=="HA_CAT") %>% 
  dplyr::filter(reporter == "early_front")%>% 
  dplyr::select(time_hours, rolled_wavefront_micron, Plaque_ID, cell_type) 

late_front <- final_df %>%
  #dplyr::filter(cell_type=="HA_CAT") %>% 
  dplyr::filter(reporter == "late_front") %>% 
  dplyr::select(time_hours, rolled_wavefront_micron, Plaque_ID, cell_type) 

# Step 2: Join the early and late data frames to calculate the distance differences
distance_difference <- early_front %>%
  inner_join(late_front, by = c("Plaque_ID", "time_hours"), suffix = c("_early", "_late")) %>%
  mutate(distance_diff = rolled_wavefront_micron_early - rolled_wavefront_micron_late)

late_front <- late_front %>%
  group_by(Plaque_ID) %>%
  arrange(time_hours) %>%
  mutate(First_Derivative = c(NA, diff(rolled_wavefront_micron) / diff(time_hours)))

# Step 3: Calculate the smoothed first derivative for the late front
late_front <- late_front %>%
  group_by(Plaque_ID, cell_type) %>%
  nest() %>%
  mutate(
    model = map(data, ~ loess(First_Derivative ~ time_hours,
                              data = .x %>% dplyr::filter(!is.na(First_Derivative)), 
                              span = 0.35, 
                              family = "symmetric")),  # Inserted family argument here
    predictions = map2(model, data, ~ predict(.x, newdata = .y))
  ) %>%
  dplyr::select(Plaque_ID, data, predictions, cell_type) %>%
  unnest(c(data, predictions)) %>%
  rename(Smoothed_First_Derivative = predictions) %>%
  ungroup()


# Step 4: Combine the smoothed derivative with the distance difference
combined_data <- distance_difference %>%
  left_join(late_front %>% dplyr::select(Plaque_ID, time_hours, Smoothed_First_Derivative), by = c("Plaque_ID", "time_hours")) %>%
  dplyr::mutate(Time_Difference = distance_diff / Smoothed_First_Derivative) %>% 
  dplyr::mutate(Time_Difference = ifelse(Time_Difference < 0, 0, Time_Difference))

threshold <- 0.8  # Set a threshold value for the derivative that you consider significant. i.e. a minimum velocity that is still considered moving
max_time_delay <- 24  # Set a reasonable maximum delay, e.g. the maximum reasonable replication time for your virus

combined_data <- combined_data %>%
  mutate(
    Time_Delay_Hours = ifelse(
      abs(Smoothed_First_Derivative) < threshold | is.na(Smoothed_First_Derivative), 
      NA, 
      distance_diff / Smoothed_First_Derivative
    ),
    Time_Delay_Hours = ifelse(Time_Delay_Hours < 0, 0, Time_Delay_Hours),  # Correct negative values
    Time_Delay_Hours = pmin(Time_Delay_Hours, max_time_delay)  # Cap the time delay at max_time_delay hours
  )

#Spatial plot
p <- combined_data %>% 
  dplyr::filter(cell_type_early=="BSC") %>% 
  left_join(by=c("Plaque_ID"), t0_accuracy_join) %>%
  group_by(Plaque_ID) %>% 
  dplyr::mutate(time_hours_corrected=time_hours-t0) %>% #synchronise by plaque stage
  dplyr::filter(time_hours_corrected>=0) %>% #remove pre-green start data
  group_by(time_hours_corrected) %>%
  ggplot(aes(x = time_hours_corrected, y = distance_diff, colour = Plaque_ID)) +
  geom_path(size=0.25) +
  labs(x = "Time (hrs)", y = "Spatial \ndelay (µm)", color="Plaque \nNumber") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",labels = labels)+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))

ggsave("BSC_spatial_delay.pdf",width=75, height=75/2, units="mm", p)
#Temporal plot #
p <- combined_data %>% 
  dplyr::filter(cell_type_early=="BSC") %>% 
  left_join(by=c("Plaque_ID"), t0_accuracy_join) %>%
  group_by(Plaque_ID) %>% 
  dplyr::mutate(time_hours_corrected=time_hours-t0) %>% #synchronise by plaque stage
  dplyr::filter(time_hours_corrected>=0) %>% #remove pre-green start data
  group_by(time_hours_corrected) %>%
  ggplot(aes(x = time_hours_corrected, y = Time_Delay_Hours, colour = Plaque_ID)) +
  geom_path(size=0.25) +
  labs(x = "Time (hrs)", y = "Temporal \ndelay (hrs)", color="Plaque \nNumber") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",labels = labels)+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))

ggsave("BSC_temporal_delay.pdf",width=75, height=75/2, units="mm", p)

#Cell type comparison temporal delay
p <- combined_data %>% 
  left_join(by=c("Plaque_ID"), t0_accuracy_join) %>%
  group_by(Plaque_ID) %>% 
  dplyr::mutate(time_hours_corrected=time_hours-t0) %>% #synchronise by plaque stage
  dplyr::filter(time_hours_corrected>=0) %>% #remove pre-green start data
  group_by(time_hours_corrected, cell_type_late) %>%
  dplyr::summarise(mean_delay=mean(Time_Delay_Hours, na.rm=TRUE)) %>% 
  dplyr::filter(time_hours_corrected<69) %>% 
  ggplot(aes(x = time_hours_corrected, y = mean_delay, colour=cell_type_late)) +
  geom_path(size=0.125) +
  labs(x = "Time (hrs post-t0)", y = "Temporal \ndelay (hrs)", color="Cell type") +
  theme_minimal() +
  scale_color_manual(values = c("BSC" = "blue2", "HA_CAT" = "orange2"),
                     labels = c("BS-C-1", "HaCaT")) +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))
# geom_smooth(data = subset(subset_data_raw, time_hours >= min(subset_data_results$t0) & time_hours <= 60), 
#             method = "lm", se = FALSE, aes(group = reporter), color = "black", linetype = "dashed", fullrange = TRUE, size=0.5) +

ggsave("BSC_vs_HACAT_temporal_delay.pdf",width=70, height=65/2, units="mm", p)

#Cell type comparison spatial delay
p <- combined_data %>% 
  left_join(by=c("Plaque_ID"), t0_accuracy_join) %>%
  group_by(Plaque_ID) %>% 
  dplyr::mutate(time_hours_corrected=time_hours-t0) %>% #synchronise by plaque stage
  dplyr::filter(time_hours_corrected>=0) %>% #remove pre-green start data
  group_by(time_hours_corrected, cell_type_late) %>%
  dplyr::summarise(mean_delay=mean(distance_diff, na.rm=TRUE)) %>% 
  dplyr::filter(time_hours_corrected<69) %>% 
  ggplot(aes(x = time_hours_corrected, y = mean_delay, colour=cell_type_late)) +
  geom_path(size=0.125) +
  labs(x = "Time (hrs)", y = "Spatial \ndelay (hrs)", color="Cell type") +
  theme_minimal() +
  scale_color_manual(values = c("BSC" = "blue2", "HA_CAT" = "orange2"),
                     labels = c("BS-C-1", "HaCaT")) +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))

ggsave("BSC_vs_HACAT_spatial_delay.pdf",width=70, height=65/2, units="mm", p)

#Supplementary Fig. 1 ####

# Unique Plaque_IDs
unique_plaque_ids <- unique(final_df$Plaque_ID)
unique_plaque_ids_sorted <- unique_plaque_ids[order(as.numeric(gsub("Plaque_", "", unique_plaque_ids)))]

# Function to create a plot for a single Plaque_ID
create_plot <- function(plaque_id) {
  # Subset the data
  subset_data_raw <- final_df[final_df$Plaque_ID == plaque_id, ]
  subset_data_results <- results_df[results_df$Plaque_ID == plaque_id, ]

  # Create the plot
  p <- ggplot(subset_data_raw, aes(x = time_hours, y = rolled_wavefront_micron, color = reporter)) +
    geom_path(aes(color = reporter), size = 0.5) +
    scale_color_manual(values = c("early_front" = "red2", "late_front" = "green3")) +
    geom_vline(data = subset_data_results, aes(xintercept = t0, color = reporter), linetype = "dashed") +
    geom_smooth(data = subset(subset_data_raw, time_hours >= min(subset_data_results$t0) & time_hours <= 60), 
                method = "lm", se = FALSE, aes(group = reporter), color = "black", linetype = "dashed", fullrange = TRUE, size = 0.5) +
    labs(title = paste(plaque_id),
         x = "Time (HPI)", 
         y = "Plaque radius (µm)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7),
          plot.title =element_text(size = 8)) +
    coord_cartesian(xlim = c(-5, 80))
  
  # Fit linear models and calculate R^2 for each reporter
  vertical_offset <- -9.5
  for(reporter in c("early_front", "late_front")) {
    lm_model <- lm(rolled_wavefront_micron ~ time_hours, data = subset(subset_data_raw, reporter == reporter))
    r_squared <- summary(lm_model)$r.squared
    r_squared_label <- paste("R^2 == ", format(r_squared, digits = 2))
    print(r_squared_label)
    p <- p + annotate("text", x = Inf, y = Inf, label = r_squared_label,
                      hjust = 1.1, vjust = 1.1 - vertical_offset, 
                      color = ifelse(reporter == "early_front", "red2", "green3"), 
                      size = 2, parse = TRUE)
    vertical_offset <- vertical_offset - 1  # Adjust this for spacing
  }
  
  return(p)
}

# Generate plots for each Plaque_ID
plots <- lapply(unique_plaque_ids_sorted, create_plot)

# Assemble the plots into a grid
grid_plot <- do.call(grid.arrange, c(plots, ncol = 4))

ggsave("plaque_plots_grid.pdf", grid_plot, width = 200, height = 200, unit="mm")
#Inferred MOI from single cell data ####
#Single-cell level delay between Early and Late taken from 'Single-cell analysis of VACV infection reveals pathogen-driven 
#timing of early and late phases and host-limited dynamics of virus production'
# Create the dataframe
data <- data.frame(
  MOI = c(1, 10, 50, 100),
  start_early = c(5.452111, 3.499631, 2.430895, 2.08513)
)

# Perform the logarithmic regression
model <- lm(start_early ~ log(MOI), data = data)

# Display the summary of the model to get the regression equation details
summary(model)

#Estimate local MOI from a known EL_delay
# Coefficients from the regression equation
a <- 5.3632
b <- -0.7380

# Function to calculate MOI based on known EL_delay
calculate_MOI <- function(EL_delay) {
  exp((EL_delay - a) / b)
}

#BS-C-1 early phase max delay: 9.034860
#BS-C-1 min delay: 2.187430
known_EL_delay <- 2.187430 # Replace with the actual EL_delay value
estimated_MOI <- calculate_MOI(known_EL_delay)

# Print the estimated MOI
print(estimated_MOI)

#Boxplots ####
p <- results_df %>%
  group_by(Plaque_ID, cell_type, reporter) %>%
  ggplot(aes(x = interaction(cell_type, reporter), y = t0, fill = reporter)) +
  geom_boxplot(width = 0.8, size = 0.1, position = position_dodge(width = 1), outlier.size = 0.5, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("early_front" = "red2", "late_front" = "green3")) +
  labs(title = "", y = expression("Mean" ~italic(t)[0] ~ " (HPI)"), x = "") +
  scale_x_discrete(labels = c("BSC.early_front" = "BS-C-1 Early", 
                              "BSC.late_front" = "BS-C-1 Late", 
                              "HA_CAT.early_front" = "HaCaT Early", 
                              "HA_CAT.late_front" = "HaCaT Late")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),  # Adjusted for readability
    axis.title = element_text(size = 7),
    text = element_text(family = "Arial")
  ) +
  coord_cartesian(ylim = c(3, 37))+ 
  geom_signif(
    comparisons = list(
      c("BSC.early_front", "HA_CAT.early_front"),
      c("BSC.late_front", "HA_CAT.late_front"),
      c("HA_CAT.early_front", "HA_CAT.late_front"),
      c("BSC.early_front", "BSC.late_front")
    ),
    annotations = c("****", "****", "*", "****"),
    map_signif_level=TRUE,
    y_position = c(23, 26, 29, 32),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

ggsave("mean_t0_cell_type.pdf",width=40, height=65, units="mm", p)

p <- results_df %>%
  group_by(Plaque_ID, cell_type, reporter) %>%
  ggplot(aes(x = interaction(cell_type, reporter), y = growth_rate, fill = reporter)) +
  geom_boxplot(width = 0.8, size = 0.1, position = position_dodge(width = 1), outlier.size = 0.5, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("early_front" = "red2", "late_front" = "green3")) +
  labs(title = "", y = "Mean growth rate (µm/hr)", x = "") +
  scale_x_discrete(labels = c("BSC.early_front" = "BS-C-1 Early", 
                              "BSC.late_front" = "BS-C-1 Late", 
                              "HA_CAT.early_front" = "HaCaT Early", 
                              "HA_CAT.late_front" = "HaCaT Late")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),  # Adjusted for readability
    axis.title = element_text(size = 7),
    text = element_text(family = "Arial")
  ) +
  coord_cartesian(ylim = c(3, 23))+ 
  geom_signif(
    comparisons = list(
      c("BSC.early_front", "HA_CAT.early_front"),
      c("BSC.late_front", "HA_CAT.late_front")
    ),
    annotations = c("**", "**"),
    map_signif_level=TRUE,
    y_position = c(17.5, 17.5),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

ggsave("mean_growth_rate_by_cell_type_reporter.pdf",width=40, height=65, units="mm", p)



#Mean final size by cell type
results_df_join<- results_df %>% 
  dplyr::filter(reporter=="late_front")

variable_importance <- final_df %>% 
  group_by(cell_type) %>% 
  dplyr::filter(reporter=="late_front") %>% 
  dplyr::filter(time_hours==max(time_hours)) %>% 
  ungroup() %>% 
  dplyr::select(Plaque_ID, rolled_wavefront_micron) %>% 
  left_join(by=c("Plaque_ID"),results_df_join)%>%
  # Group by cell_type to apply normalization within each cell type
  group_by(cell_type) %>%
  # Apply min-max normalization to each variable
  mutate(
    t0_normalized = (t0 - min(t0, na.rm = TRUE)) / (max(t0, na.rm = TRUE) - min(t0, na.rm = TRUE)),
    growth_rate_normalized = (growth_rate - min(growth_rate, na.rm = TRUE)) / (max(growth_rate, na.rm = TRUE) - min(growth_rate, na.rm = TRUE)),
    rolled_wavefront_micron_normalized = (rolled_wavefront_micron - min(rolled_wavefront_micron, na.rm = TRUE)) / (max(rolled_wavefront_micron, na.rm = TRUE) - min(rolled_wavefront_micron, na.rm = TRUE))
  ) %>%
  ungroup() # Ungroup if further operations require non-grouped data



# Calculate Pearson's correlation coefficient squared (R^2) for each cell type
# Calculate R squared values for each cell type and each relationship
r_squared_t0 <- variable_importance %>%
  group_by(cell_type) %>%
  summarize(r_squared = round(cor(rolled_wavefront_micron, t0, use = "complete.obs")^2, digits = 2)) %>%
  ungroup()

r_squared_growth_rate <- variable_importance %>%
  group_by(cell_type) %>%
  summarize(r_squared = round(cor(rolled_wavefront_micron, growth_rate, use = "complete.obs")^2, digits = 2)) %>%
  ungroup()

# Define color mapping for cell types
colors <- c("BSC" = "blue2", "HA_CAT" = "orange2")

r_squared_t0 <- r_squared_t0 %>%
  mutate(y_position = 1 - (row_number() - 1) / 12)

# Now create the plot
plot_t0 <- ggplot(variable_importance, aes(x = t0_normalized, y = rolled_wavefront_micron_normalized, color = cell_type)) +
  geom_point(alpha = 0.5, size=0.5) +
  geom_smooth(method = "lm", aes(fill = cell_type), fullrange = TRUE, size=0.5, alpha = 0.2) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title = element_text(size = 7),
    text = element_text(family = "Arial")
  ) +
  labs(title = "",
       y = "Normalized final plaque radius (µm)", 
       x = expression("Normalized "~italic(t)[0] ~ " (HPI)")) +
  geom_text(
    data = r_squared_t0, 
    aes(x = 0.95, y = y_position, label = sprintf("R² = %.2f", r_squared), color = cell_type), 
    hjust = 1.2, vjust = 0, size = 2
  )

r_squared_growth_rate <- r_squared_growth_rate %>%
  mutate(y_position = 1.2 - (row_number() - 1) / 12)

# Now create the plot
plot_growth_rate <- ggplot(variable_importance, aes(x = growth_rate_normalized, y = rolled_wavefront_micron_normalized, color = cell_type)) +
  geom_point(alpha = 0.5, size=0.5) +
  geom_smooth(method = "lm", aes(fill = cell_type), fullrange = TRUE, size=0.5, alpha = 0.2) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title = element_text(size = 7),
    text = element_text(family = "Arial")
  ) +
  labs(title = "",
       y = "Normalized final plaque radius (µm)", 
       x = "Normalized growth rate (µm/hr)") +
  geom_text(
    data = r_squared_growth_rate, 
    aes(x = 0.95, y = y_position, label = sprintf("R² = %.2f", r_squared), color = cell_type), 
    hjust = 1.2, vjust = 0, size = 2
  )

ggsave("plot_t0_combined.pdf",width=50, height=55, units="mm", plot_t0)
ggsave("plot_growth_rate_combined.pdf", width=50, height=55, units="mm", plot_growth_rate)


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

#Establishment frequency ####
#Statistical testing
# Create the data frame
establishment_frequency <- data.frame(
  cell_type = c("BSC", "HA_CAT"),
  established = c(12 + 15, 21),
  failed = c(18 + 29, 111)
)
# Convert the data frame to a contingency table
contingency_table <- as.matrix(establishment_frequency[, -1])  # Exclude the first column (cell_type) and convert to matrix
rownames(contingency_table) <- establishment_frequency$cell_type

# Print the contingency table
print(contingency_table)

# Perform the chi-squared test
chi_squared_test <- chisq.test(contingency_table)

# Print the results of the chi-squared test
print(chi_squared_test)

established <- contingency_table[,1] #subset established count
failed <- contingency_table[,2] #subset failed count
ratio <- established / (failed+established) # Calculate establishment frequency

# Create a dataframe
results <- data.frame(
  cell_type = rownames(contingency_table),  # Extract cell types from the row names of the matrix
  ratio = ratio  # Use the calculated ratio
)


p_establishment <- results %>%
  ggplot(aes(x = cell_type, y = ratio*100, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 1)) +  # Use geom_bar with stat = "identity"
  scale_fill_manual(values = c("BSC" = "blue2", "HA_CAT" = "orange2")) +
  labs(title = "", y = "Establishment frequency (%)", x = "") +
  scale_x_discrete(labels = c("BSC" = "BS-C-1", 
                              "HA_CAT" = "HaCaT")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.title = element_text(size = 7),
    text = element_text(family = "Arial")
  ) +
  coord_cartesian(ylim = c(0, NA)) +  # Adjust as needed
  geom_signif(
    comparisons = list(c("BSC", "HA_CAT")),
    annotations = c("**"),  # Use the actual significance level from your chi-squared test
    map_signif_level = TRUE,
    y_position = c(45),  # Adjust based on your plot
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )
ggsave("establishment_frequency.pdf",width=45, height=60, units="mm", p_establishment)

#Optimising parameter combinations ####
#WARNING: Moderately computationally intensive process. Check, and consider carefully, the number of combinations generated in param_grid
# Assuming 'final_df' and 'manual_t0s' are already defined and available
# manual_t0 values should be set based on manual observation of the timepoint at which a plaque spreads from a single infected cell
#to multiple infected cells. Take manual observations for a subset of your data (at least 5 plaques per condition) and assess the accuracy
#of the automated process vs your manual observations
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
    
    
    # Your existing loop starts here
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
  
  # The first row now contains the parameter set with the lowest sum_accuracy
  optimal_parameters <- results[1,]