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
library(FSA)

#Load font
showtext_auto()
font_add("Arial", "C:/Windows/Fonts/arial.ttf")

#Figure 1 ####
#Panel B
raw_data_temp <- data_to_plot_1 %>% 
  filter(Plaque_ID=="Plaque_7") %>% 
  filter(Time==280)

p <- raw_data_temp %>% 
  mutate(Intensity = max(Intensity) - Intensity + min(Intensity)) %>% 
  ggplot() +
  geom_line(aes(x=Location, y=Intensity), colour = "black", size=0.25) +
  labs(x = "Radial position", y = expression("Probability")) +
  theme_minimal()+
  theme(legend.position = "right",
        axis.text = element_text(size = 4),
        axis.title = element_text(size = 5),
        text=element_text(family ="Arial")) 

ggsave("radial_quantitation_example.pdf",width=35, height=35, units="mm", p)

final_df %>% 
  group_by(Plaque_ID, cell_type) %>% 
  dplyr::summarise(n=n())
#Panel C
p <- final_df %>% 
  filter(cell_type=="BSC") %>% 
  filter(Plaque_ID=="Plaque_7") %>% 
  #filter(reporter=="early_front") %>% 
  ggplot() +
  geom_line(aes(x=time_hours, y=pi*rolled_wavefront_micron^2, colour = reporter), size=0.25) +
  labs(x = "Time (HPI)", y = expression("Plaque Area (µm"^2*")")) +
  scale_color_manual(values=c("red3","green3"))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial")) 

ggsave("Single_plaque_area.pdf",width=50, height=50, units="mm", p)

#Figure 2 ####
#Panel A
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

#Panel B
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

#Panel C
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

#Panel D
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

#Figure 3####
#Panel A 
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

#Panel B
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

#Panel C
results_df_summary <- results_df %>% 
  group_by(cell_type,reporter) %>% 
  dplyr::summarise(mean_growth_rate=mean(growth_rate))
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

#Panel D
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

gsave("BSC_vs_HACAT_temporal_delay.pdf",width=70, height=65/2, units="mm", p)

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
  labs(x = "Time (hrs)", y = "Spatial \ndelay (µm)", color="Cell type") +
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

#Panel E
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


#Panel F
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

#Supplementary Figure 1 ####

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

#spare ####
#Plotting wavefronts with plaques synchronised by t0.
p <- combined_data %>% 
  left_join(by=c("Plaque_ID"), t0_accuracy_join) %>%
  group_by(Plaque_ID) %>% 
  dplyr::mutate(time_hours_corrected=time_hours-t0) %>% #synchronise by plaque stage
  dplyr::filter(time_hours_corrected>=0) %>% #remove pre-green start data
  group_by(time_hours_corrected, cell_type_late) %>%
  ggplot(aes(x = time_hours_corrected, y = rolled_wavefront_micron_late, group=Plaque_ID, colour=cell_type_late)) +
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
