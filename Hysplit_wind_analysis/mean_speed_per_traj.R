#!/opt/R/4.2.3/bin/Rscript

#### DOC ####
' Given an .RData file containing HYSPLIT trajectory output files, compute the mean speed of each trajectory by taking the lenght of all the segements in the trajectory and the duration of the trajectory.

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
library(raster)

#### ARGS ####
args <- docopt(doc, version = 'wind_speed_plots.R 0.1')
files_list <- unlist(args$RDATA)

#### MAIN ####
sink("mean_traj_speeds.csv")
print("file,starting_altitude,mean_speed(m/s),se_speed(m/s)")
for (file in files_list) {
  load(file)
  trajs <- trajs[[1]]
  trajs_lengths <- c()
  duration <- max(abs(trajs$hour.inc))
  for (start_height in unique(trajs$start_height)) {
    for (start_date in unique(trajs$date)) {
      df_subset <- trajs[trajs$start_height == start_height & trajs$date == start_date,]
      traj_len = 0 # variable to add the length of all segment to get trajectory length
      for (i in 1:nrow(df_subset)){
        if (i != 1){
          lat1 = df_subset$lat[i]
          lon1 = df_subset$lon[i]
          
          lat2 = df_subset$lat[i-1]
          lon2 = df_subset$lon[i-1]
          
          segment_len = pointDistance(cbind(lon1, lat1), c(lon2, lat2), lonlat = T)
          traj_len = traj_len + segment_len
        }
      }
      trajs_lengths[length(trajs_lengths) + 1] <- traj_len
    }
    trajs_mean_speeds <- trajs_lengths/(3600*duration)
    trajs_mean_speeds_mean <- mean(trajs_mean_speeds)
    trajs_mean_speeds_SE <- sd(trajs_mean_speeds)/sqrt(length(trajs_mean_speeds))
    print(paste(file, start_height, trajs_mean_speeds_mean, trajs_mean_speeds_SE, sep = ","))
  }
}
sink()
