#!/usr/bin/env Rscript

#### DOC ####
' Plot windspeed plots from a given .RData file containing HYSPLIT trajectory output files.

Usage:
  wind_speed_plots.R <RDATA>...

  Arguments:
   RDATA                 One or more .RData files output by Hysplit_wind_analysis.R [https://github.com/etd530/Hysplit_R_interface]

  Options:
    -v, --verbose         Print the progression of the program execution to the terminal (Standard Error).
    -h, --help            Show this message and exit.

' -> doc

#### LIBS ####
library(docopt)
library(ggplot2)
library(viridis)

#### ARGS ####
args <- docopt(doc, version = 'wind_speed_plots.R 0.1')
files_list <- unlist(args$RDATA)

#### MAIN ####
pdf("windspeed_plots.pdf", height = 7, width = 14)
for (file in files_list) {
  load(file)
  trajs <- trajs[[1]]
  
  hours_since_start_vec <- c()
  wind_speed_mean_vec <- c()
  wind_speed_se_vec <- c()
  start_height_vec <- c()
  for (start_height in unique(trajs$start_height)) {
    for (hour_since_start in unique(trajs$hour.inc)) {
      df_subset <- trajs[trajs$start_height == start_height & trajs$hour.inc == hour_since_start,]
      hours_since_start_vec[length(hours_since_start_vec) + 1] = hour_since_start
      wind_speed_mean_vec[length(wind_speed_mean_vec) + 1] = mean(df_subset$wind_speed)
      wind_speed_se_vec[length(wind_speed_se_vec) + 1] = sd(df_subset$wind_speed)/sqrt(nrow(df_subset))
      start_height_vec[length(start_height_vec) + 1] = start_height
    }
  }
  
  wind_speed_mean_plus_SE = wind_speed_mean_vec + wind_speed_se_vec
  wind_speed_mean_minus_SE = wind_speed_mean_vec - wind_speed_se_vec
  
  wind_speeds_df <- data.frame(hours_since_start_vec, wind_speed_mean_vec, 
                               wind_speed_se_vec, start_height_vec, 
                               wind_speed_mean_plus_SE, wind_speed_mean_minus_SE)
  
  minDate = min(trajs$date2[trajs$hour.inc == 0])
  maxDate = max(trajs$date2[trajs$hour.inc == 0])
  direction = "backwards"
  duration = max(abs(trajs$hour.inc))
  
  # WARNING: LEVELS FOR HEIGHTS FACTOR ARE CURRENTLY HARDCODED!!!!!
  windspeed_plot <- ggplot(data = wind_speeds_df, aes(x = abs(hours_since_start_vec), y = wind_speed_mean_vec)) +
    ggtitle(paste0("Wind speed profile from ", minDate, "\nto ", maxDate, ", ", direction, " ", abs(duration), " hours")) +
    scale_color_viridis(begin = 1, end = 0, discrete = T, alpha = 1) +
    scale_fill_viridis(begin = 1, end = 0, discrete = T, alpha = 0.25) + 
    ylab("Wind speed (m/s)") +
    xlab("Hours before observation") +
    theme(plot.title = element_text(hjust = 0.5)) + labs(color = "Start height") +
    geom_line(data = wind_speeds_df, aes(x = abs(hours_since_start_vec), y = wind_speed_mean_vec, 
                                         group = start_height_vec, 
                                         color = factor(start_height_vec, levels = c("500", "1000", "2000"), ordered = T))) +
    geom_line(data = wind_speeds_df, aes(x = abs(hours_since_start_vec), y = wind_speed_mean_plus_SE, 
                                         group = start_height_vec, 
                                         color = factor(start_height_vec, levels = c("500", "1000", "2000"), ordered = T))) +
    geom_line(data = wind_speeds_df, aes(x = abs(hours_since_start_vec), y = wind_speed_mean_minus_SE, 
                                         group = start_height_vec, 
                                         color = factor(start_height_vec, levels = c("500", "1000", "2000"), ordered = T))) +
    geom_ribbon(aes(ymin = wind_speed_mean_minus_SE, ymax = wind_speed_mean_plus_SE, 
                    group = start_height_vec, fill = factor(start_height_vec, levels = c("500", "1000", "2000"), ordered = T)),
                    alpha = 0.5) +
    guides(fill = "none") +
    xlim(0, max(abs(duration)))
  
  plot(windspeed_plot)
}
dev.off()