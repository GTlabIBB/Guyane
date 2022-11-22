#!/usr/bin/Rscript
#### HYSPLIT Program ####

#### Load packages ####
library(splitr)       # to work with Hysplit (to download files mostly)
library(opentraj)     # to work with Hysplit (does the calculations and plotting)
library(lubridate)    # for parsing dates
library(ggplot2)      # plotting library
library(raster)       # needed for the lines that change the raster values, from maxValue until setValues, and for pointDistnace()
library(geosphere)    # needed for bearing()
library(viridis)      # colorblind-friendly color palettes
library(plyr)

#### Variables for the runs ####
dayList <- c(21:31)                          # put the days of the month here, without caring about short months
monthList <- c(10)                      # months go here
yearList <- c(2013)                    # years
dayblocks <- list(c(21:22), c(22:23), c(23:24), c(24:25), c(25:26),
                  c(26:27), c(27:28), c(28:29), c(29:30), c(30:31))       # Set the blocks of days you want to run together to plot in the same map
coord <- list(c(5.745974, -53.934047))       # coordinates
height <- c(500, 1000, 2000)                               # height of the winds at starting point
duration <- -200                             # how long forwards or backwards should the run go
times <- list(c("06:00", "06:00"))          # first and last hour on which the trajectories should start (put the same to run just at one hour)
hourInt <- 1                                # at which intervals should you start new trajectories (every 2 hours, etc.)
TZ <- "Brazil/East"                         # time zone to use
bb <- as.matrix(cbind(c(-60.0, -12.0), c(10.0, 45.0))) #limits of the map for the plots
colnames(bb) <- c("min", "max")
rownames(bb) <- c("x", "y")

##### FUNCTIONS #####
#modified version of PlotTrajFreq, so that we can change the scale of the plot and the color scale
plotRaster=function (spGridDf, background = T, overlay = NA, overlay.color = "white", 
                     pdf = F, file.name = "output", bb,...) 
{
  if (pdf == T) {
    pdf(file.name, paper = "USr", height = 0, width = 0)
  }
  oldpar <- par(no.readonly = TRUE)
  par(mar = c(0, 0, 0, 0) + 2)
  plot.add <- F
  extra.args <- list(...)
  if (!"main" %in% names(extra.args)) {
    extra.args$main <- NULL
  }
  if (background == T) {
    bb
    PlotBgMap(spGridDf, xlim = bb[1, ], ylim = bb[2, ], 
              axes = TRUE)
    grid(col = "white")
    plot.add <- T
  }
  
  grays <- colorRampPalette(c("yellow", "orange", "orangered", "red"))(12) #names are color range, number is how many colors to generate
  
  grays[length(grays)+1] <- "#FFFFFF00"
  grays[length(grays)+1] <- "#000000"
  
  #if you change the number of colors in the previous line you must change breaks and legend accordingly
  image(spGridDf, col = grays, breaks = (c(0, 0.02, 0.04, 0.06, 
                                           0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.99, 1)), add = plot.add)
  legend("topleft", legend = c("0.00 - 0.02", "0.02 - 0.04", "0.04 - 0.06", 
                               "0.06 - 0.08", "0.08 - 0.10", "0.10 - 0.12", "0.12 - 0.14", 
                               "0.14 - 0.16", "0.16 - 0.18", "0.18 - 0.2" ,"0.2 - 0.22", "0.22 - 0.24"), fill = grays)
  do.call(title, extra.args)
  if (!missing(overlay)) {
    plot(overlay, add = T, col = "black", border = "black")
  }
  par(oldpar)
  if (pdf == T) {
    dev.off()
  }
}

# modified version of ProcTraj from opentraj
ProcTrajMod = function (lat = 51.5, lon = -45.1, hour.interval = 1, name = "london", 
                        start.hour = "00:00", end.hour = "23:00", met, out, hours = 12, 
                        height = 100, hy.path, ID = 1, dates, script.name = "test", 
                        add.new.column = F, new.column.name, new.column.value, tz = "GMT", 
                        clean.files = TRUE) 
{
  wd <- getwd()
  script.extension <- ".sh"
  OS <- "unix"
  if (.Platform$OS.type == "windows") {
    script.extension <- ".bat"
    OS <- "windows"
  }
  hy.split.wd <- file.path(hy.path, "working")
  hy.split.wd <- normalizePath(hy.split.wd)
  setwd(hy.split.wd)
  folder.name = paste("process_", ID, sep = "")
  process.working.dir <- file.path(hy.split.wd, folder.name)
  dir.create(process.working.dir, showWarnings = FALSE)
  process.working.dir <- normalizePath(process.working.dir)
  setwd(process.working.dir)
  hy.split.exec.dir <- file.path(hy.path, "exec", "hyts_std")
  bdyfiles.path <- file.path(hy.path, "bdyfiles")
  symb.link.files <- list.files(path = bdyfiles.path)
  for (i in 1:length(symb.link.files)) {
    from <- normalizePath(file.path(bdyfiles.path, symb.link.files[[i]]))
    to <- file.path(process.working.dir, symb.link.files[[i]])
    file.copy(from, to)
  }
  control.file.number <- 1
  script.name <- paste(script.name, "_", ID, script.extension, 
                       sep = "")
  dates.and.times <- laply(.data = dates, .fun = function(d) {
    start.day <- paste(d, start.hour, sep = " ")
    end.day <- paste(d, end.hour, sep = " ")
    posix.date <- seq(as.POSIXct(start.day, tz), as.POSIXct(end.day, 
                                                            tz), by = paste(hour.interval, "hour", sep = " "))
    as.character(posix.date)
  })
  dates.and.times <- unique(dates.and.times)
  hour.interval <- paste(hour.interval, "hour", sep = " ")
  for (i in 1:length(dates.and.times)) {
    control.file <- "CONTROL"
    date <- as.POSIXct(dates.and.times[i], tz = tz)
    control.file.extension <- paste(as.character(ID), "_", 
                                    control.file.number, sep = "")
    control.file <- paste(control.file, control.file.extension, 
                          sep = ".")
    year <- format(date, "%y")
    Year <- format(date, "%Y")
    month <- format(date, "%m")
    day <- format(date, "%d")
    hour <- format(date, "%H")
    script.file <- file(script.name, "w")
    if (OS == "unix") {
      cat("#!/bin/bash", file = script.file, sep = "\n")
    }
    line <- paste("echo", year, month, day, hour, ">", control.file, 
                  sep = " ")
    cat(line, file = script.file, sep = "\n")
    line <- paste("echo 1 >>", control.file, sep = " ")
    cat(line, file = script.file, sep = "\n")
    line <- paste("echo", lat, lon, height, ">>", control.file, 
                  sep = " ")
    cat(line, file = script.file, sep = "\n")
    line <- paste("echo", hours, ">>", control.file, sep = " ")
    cat(line, file = script.file, sep = "\n")
    line <- paste("echo 0 >> ", control.file, "\n", "echo 10000.0 >> ", 
                  control.file, "\n", "echo 3 >> ", control.file, 
                  "\n", sep = "")
    cat(line, file = script.file, sep = "")
    months <- as.numeric(unique(format(date, "%m")))
    months <- c(months, months + 1:2)
    months <- months - 1
    months <- months[months <= 12]
    if (length(months) == 2) {
      months <- c(min(months) - 1, months)
    }
    for (i in 1:3) {
      AddMetFiles(months[i], Year, met, script.file, control.file)
    }
    line <- paste("echo ./ >>", control.file, sep = " ")
    cat(line, file = script.file, sep = "\n")
    line <- paste("echo tdump", "_", ID, "_", year, month, 
                  day, hour, " >> ", control.file, sep = "")
    cat(line, file = script.file, sep = "\n")
    line <- paste(hy.split.exec.dir, control.file.extension, 
                  sep = " ")
    cat(line, file = script.file, sep = "\n")
    close(script.file)
    if (OS == "unix") {
      system(paste0("sh ", script.name))
    }
    else {
      system(paste0(script.name))
    }
    control.file.number <- control.file.number + 1
  }
  traj <- ReadFiles(process.working.dir, ID, dates.and.times, 
                    tz)
  if (add.new.column == T) {
    if (!missing(new.column.name) & !missing(new.column.value)) {
      traj[new.column.name] <- new.column.value
    }
    else {
      stop("Parameters 'new.column.name' and 'new.column.value' are not defined.")
    }
  }
  if (!missing(out)) {
    file.name <- paste(out, name, Year, ".RData", sep = "")
    save(traj, file = file.name)
  }
  setwd(hy.split.wd)
  if (clean.files == T) {
    unlink(folder.name, recursive = TRUE)
  }
  setwd(wd)
  traj
}

##### END FUNCTIONS #####

# Creating the list with the days of interest ####
dateList <- as.vector(5) # just creating a vector
for(i in 1:length(yearList)) {
  for (j in 1:length(monthList)) {
    for (k in 1:length(dayList)) {
      dateList[length(dayList)*length(monthList)*(i-1)+length(dayList)*(j-1)+k] <- as.character(paste(yearList[i], monthList[j], dayList[k], sep = "-"))
    }
  }
} # this loop generates the dates

dateList <- as.Date(dateList) # change from strings to Date objects
dateList <- na.omit(dateList) # remove NAs (i.e. remove impossible dates such as February 31)

# Getting the meteorological files ####
# ProcTraj() requires one meteorological file before and after the limits of what you want to plot,
# so we download the files for the months previous and next to the ones of interest

# Generate list with the dates of the previous month
prevDates <- seq.Date(if (monthList[1]==01) as.Date(paste(yearList[1]-1, "12", "01", sep="-"), "%Y-%m-%d", tz = TZ)
                      else as.Date(paste(yearList[1], monthList[1]-1, "01", sep="-"), "%Y-%m-%d", tz = TZ), # initial date
                      if (monthList[1]==01) as.Date(paste(yearList[length(yearList)]-1, "12", "01", sep="-"), "%Y-%m-%d", tz = TZ)
                      else as.Date(paste(yearList[length(yearList)], monthList[1]-1, "01", sep="-"), "%Y-%m-%d", tz = TZ), # final date
                      by = "year",                                   # interval
                      length.out = NULL)                             # period length

# Get meteorological files:
# I think with the loop is better because if you put the whole list directly in the "days" argument it seems to go over all the months in the middle
# and download unnecessary files
for(i in 1:length(prevDates)) {
  get_met_reanalysis(days = prevDates[i], duration = 12, direction = "forward",
                     path_met_files = "C:/hysplit/working")
}

# Generate list with dates of next month
postDates <- seq.Date(if (monthList[length(monthList)]==12) as.Date(paste(yearList[1]+1, "01", "01", sep = "-"), "%Y-%m-%d", tz = TZ)
                      else as.Date(paste(yearList[1], monthList[length(monthList)]+1, "01", sep = "-"), "%Y-%m-%d", tz = TZ), # initial date
                      if (monthList[length(monthList)]==12) as.Date(paste(yearList[length(yearList)]+1, "01","01", sep = "-"), "%Y-%m-%d", tz = TZ)
                      else as.Date(paste(yearList[length(yearList)], monthList[length(monthList)]+1,"01", sep = "-"), "%Y-%m-%d", tz = TZ), # final date
                      by = "year",                                   # interval
                      length.out = NULL)                             # period length

# Get the files; as before the loop should save time
for(i in 1:length(postDates)) {
  get_met_reanalysis(days = postDates[i], duration = 12, direction = "forward",
                     path_met_files = "C:/hysplit/working")
}

# Get the files for the dates of interest; since the files are by month it does not matter here if the list has more days than needed, it will skip files already downloaded
for(i in 1:length(dateList)) {
  get_met_reanalysis(days = dateList[i], duration = 48, direction = "backward",
                     path_met_files = "C:/hysplit/working")
}


# Calculate trajectories ####
pdf("./Guiana_winds_21-31_daily.pdf")
for (n in coord){
  for (i in monthList){
    for (j in dayblocks){
      if (exists("merged_trajs")){rm(merged_trajs)}
      if (exists("merged_trajlines_df")) {rm(merged_trajlines_df)}
      for(h in height){
        if (exists("traj")){rm(traj)}
        for (dayNum in 1:length(j)) {
          
          # Set starting and ending hours depending on if its the first, last or a middle day
          if (dayNum == 1) {
            startHour = times[[1]][1]
            endHour = "23:00"
          } else if (dayNum == length(j)){
            startHour = "00:00"
            endHour = times[[1]][2]
          } else {
            startHour = "00:00"
            endHour = "23:00"
          }
          ###Calculate the trajectories
          CurrentTraj <- ProcTrajMod(lat = n[1], lon = n[2],
                                     hour.interval = hourInt, name = "traj", start.hour = startHour, end.hour = endHour,
                                     met = "C:/hysplit/working/", out = "C:/hysplit/working/Out_files/", hours = duration, height = h, 
                                     hy.path = "C:/hysplit/", dates = dateList[day(dateList) %in% j][dayNum], tz = TZ)
          
          #Bind trajectories for previous days to current one
          if (exists("traj")) {
            traj <- rbind(traj, CurrentTraj)
          } else {
            traj <- CurrentTraj
          }
        }
        
        # get the SpatialLinesDf
        traj_lines<-Df2SpLines(traj, crs = "+proj=longlat +datum=NAD27")
        traj_lines_df<-Df2SpLinesDf(traj_lines, traj, add.distance = T, add.azimuth = T)
        
        # add starting height to the SpatialLinesDf
        traj_lines_df$start_height <- h
        
        # combine all SpatialLinesDf of the same run and different heights in a single file
        if (exists("merged_trajlines_df")) {
          merged_trajlines_df <- rbind(merged_trajlines_df, traj_lines_df)
        } else {
          merged_trajlines_df <- traj_lines_df
        }
        
        yearString <- if (length(yearList)==1) {yearList[1]} else {paste0(yearList[1],"-",yearList[length(yearList)])}
        
        #add the starting height of the trajectories as a variable
        traj$start_height <- h
        
        #join all trajectories of the same days and different starting heights in a single df
        if (exists("merged_trajs")) {
          merged_trajs <- rbind(merged_trajs, traj)
        } else {
          merged_trajs <- traj
        }
      }
      
      rm(traj)
      
      #Plot the trajectories in a map
      setAlpha = ifelse(merged_trajlines_df$day == 28 & merged_trajlines_df$hour == 6, 1, 0.25)
      PlotBgMap(merged_trajlines_df, xlim = bb[1, ], ylim = bb[2, ], axes = TRUE)
      plot(merged_trajlines_df,
           col = ifelse(merged_trajlines_df$start_height == 500, 
                        viridis(n=1, alpha = setAlpha, begin = 1), 
                        ifelse(merged_trajlines_df$start_height == 1000,
                               viridis(n=1, alpha = setAlpha, begin = 0.5),
                               viridis(n=1, alpha = setAlpha, begin=0))),
           add = T,
      )
      
      title(main = paste0(month.name[i]," ", j[1], " ", times[[1]][1], " to ", 
                          month.name[i], " ", j[length(j)], " ", times[[1]][2], " ", 
                          yearString),
            outer = T, line = -1.6
      )
      
      legend(bb[1,1], bb[2,2], legend = c("500m", "1000m", "2000m"), bg = "transparent",
             fill = c(viridis(n=1, begin = 1), viridis(n=1, begin = 0.5), 
                      viridis(n=1, begin = 0)), title = "Starting altitude (m AGL)")
      
      # Calculate distance from origin and angle relative to origin for each trajectory position
      merged_trajs$dist <- pointDistance(cbind(merged_trajs$lon, merged_trajs$lat), c(coord[[1]][2], coord[[1]][1]), lonlat = T)
      merged_trajs$angle <- bearing(c(coord[[1]][2], coord[[1]][1]), cbind(merged_trajs$lon, merged_trajs$lat))
      
      #we take out values for the starting points, as they are still at the origin so distance is zero and it makes no sense to calculate an angle
      merged_trajs$dist[merged_trajs$hour.inc==0] <- NA
      merged_trajs$angle[merged_trajs$hour.inc==0] <- NA
      
      #add 360 to negative azimuths
      for (ang in 1:length(merged_trajs$angle)) {
        if (!is.na(merged_trajs$angle[ang]) & (merged_trajs$angle[ang] < 0)) {
          merged_trajs$angle[ang] <- merged_trajs$angle[ang] + 360
        }
      }
      
      # turn starting height into factor
      merged_trajs$start_height <- factor(merged_trajs$start_height, levels = c("2000", "1000", "500"))
      
      # Build windrose plot for -200h
      
      windrose = ggplot(data=merged_trajs[merged_trajs$hour.inc == -200,], aes(x=angle, y=stat(count/sum(count)), group=start_height,
                                                                               fill=start_height)) + 
        geom_histogram(aes(y = stat(count/sum(count))), bins = 360) +
        coord_polar(start = 0, clip = "off") +
        ggtitle(paste0("Wind directions ", month.name[i]," ", j[1], " ", times[[1]][1], " to ", 
                       month.name[i], " ", j[length(j)], " ", times[[1]][2], " ", 
                       yearString, " (backwards 200h)")) +
        scale_fill_viridis(discrete = T, alpha = 1) +
        scale_x_continuous(breaks =c(0, 90, 180, 270) , limits = c(0, 360), labels = c("N", "E", "S", "W")) +
        theme(plot.title = element_text(hjust = 0.5))
      print(windrose)
    }
  }
}

dev.off()

#####__END__######