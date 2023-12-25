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
#Name the fluorescence channels to analyse
channel_1 <- "TXRED"
channel_2 <- "FITC"

#specify the time between frames in minutes, and the start time post infection in hours
#required for converting[Time] units from slice to HPI
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

wavefront_Early <- wavefront_early %>% 
  group_by(Plaque_ID) %>% 
  filter(plaque_front_location-lag(plaque_front_location) < 30) %>%
  filter(plaque_front_location-lag(plaque_front_location, n=2) < 30) %>% #filters out errors introduced by imaging artifacts that jump the front far outwards
  dplyr::select(Plaque_ID, Time, plaque_front_location) %>% 
  group_by(Plaque_ID, Time) %>% 
  dplyr::rename(early_front = plaque_front_location)

#Late reporter front calculation ####
wavefront_late <- data_to_plot_2 %>% 
  unite("filter", c("Plaque_ID","Time"), remove = FALSE) %>%
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

wavefront_Late <- wavefront_late %>% 
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
  dplyr::mutate(min_val = ifelse(reporter == "late_front", min(rolled_wavefront_micron, na.rm=TRUE), NA_real_)) %>%
  ungroup() %>% 
  tidyr::fill(min_val, .direction = "updown") %>%
  dplyr::mutate(rolled_wavefront_micron = rolled_wavefront_micron - min_val) %>%
  dplyr::select(-min_val)

write.csv(final_df, "tidy_results.csv")
