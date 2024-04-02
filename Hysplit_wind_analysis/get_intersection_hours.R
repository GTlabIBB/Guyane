#!/opt/R/4.2.3/bin/Rscript

#### DOC ####
' Given an .shp file containing traejctory intersections, compute the mean time to intersection of the trajectories in each file by starting altitudes.

Usage:
  wind_speed_plots.R <SHAPEFILE>...

  Arguments:
   SHAPEFILE                 One or more .shp files output by Rdata_to_shapefile.R

  Options:
    -v, --verbose         Print the progression of the program execution to the terminal (Standard Error).
    -h, --help            Show this message and exit.

' -> doc

# Options:
#--shapefile = Shapefile with the intersection points to work on

#### LIBS ####
library(terra)
library(lubridate)
# library(R.utils)
library(docopt)
library(stringr)

#### ARGS ####
args <- docopt(doc, version = 'get_intersection_hours.R 0.1')
files_list <- unlist(args$SHAPEFILE)

# # SET COMMAND LINE ENVIRONMENTS
# args <- R.utils::commandArgs(asValues=TRUE)
# str(args)
# # Turn arguments into R variables
# args <- R.utils::commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]
# keys <- attachLocally(args)
# cat("\n Command-line arguments attached to global environment:\n\n");
# 
# str(mget(keys, envir=globalenv()))

#### VARS ####
# shapefile = "intersections/intersections_2013-10-27_2013-10-28_200_053.934W_05.746N_africa_atlantic_coast.shp"

#### MAIN ####
sink("intersection_times.csv")
print("time_interval,starting_altitude,mean_intersection_time(hours),se_intersection_time(hours),proportion_of_arrivals")
for (shapefile in files_list){
  intersections = vect(shapefile)
  time_interval = str_match(shapefile, "[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{4}-[0-9]{2}-[0-9]{2}")[,1]
  start_day = as.integer(strsplit(strsplit(time_interval, "_")[[1]][1], "-")[[1]][3])
  end_day = as.integer(strsplit(strsplit(time_interval, "_")[[1]][2], "-")[[1]][3])
  duration_days = end_day - start_day
  total_trajs = ifelse(duration_days == 1, 25, 49) ### WARNING: This is hardcoded for the Guyane case for now!!!! needs fixing!!!
  for (height in unique(intersections$start_heig)){
    mask <- intersections$start_heig == height
    int_dates = lubridate::as_datetime(intersections$int_date[mask], tz = "Brazil/East")
    start_dates = lubridate::as_datetime(intersections$start_date[mask], tz = "Brazil/East")
    times_to_arrival = as.integer(difftime(int_dates, start_dates, units = "hours"))
    mean_arrival_time = mean(times_to_arrival)
    se_arrival_time = sd(times_to_arrival)/sqrt(length(times_to_arrival))
    number_of_arrivals = length(times_to_arrival)
    proportion_of_arrivals = round(100*(number_of_arrivals/total_trajs), digits = 2)
    print(paste(time_interval, height, mean_arrival_time, se_arrival_time, proportion_of_arrivals, sep = ","))
  }
  # print("Trajectories file:")
  # print(shapefile)
  # print("Intersection times (mean +- SE) at 500m:")
  # 
  # print(mean_arrival_time)
  # print(se_arrival_time)
  # hist(times_to_arrival)
}
sink()