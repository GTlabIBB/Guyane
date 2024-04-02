#!/usr/bin/env Rscript
library(terra)
#library(rgeos)
#library(sf)
#library(raster)
library(R.utils)
#library(tibble)
#library(renv)
# 

# Options:
#--WD = Working Directory
#--intersect_files = Comma Separated values of the Shapefiles to intersect with trajectories
#--traj_Rdata = Output RData to load
#--n_lin = Number of lines to vectorize
#--write_shp = Write trajectories as ShapeFile
#--hours_backwards = Number of hours to reconstruct trajectories backwards
#"C:\Program Files\R\R-4.3.0\bin\Rscript.exe" Rdata_to_shp.R --WD "working_directory" --med_line1 "med_line1" --med_line2 "med_line2" --traj_Rdata "traj_Rdata" --n_lin "100"
# MEDITERRANEAN WINDS Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/Mediterranean_winds/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/Mediterranean_winds/med_south_10_deg.shp,/home/rogerlm/Hysplit/Mediterranean_winds/sahel.shp,/home/rogerlm/Hysplit/Mediterranean_winds/med_south_v2.shp --traj_Rdata "/home/eric/mediterranean_migratory_routes/mediterranean_wind_dispersal_1980-2022_Feb-Jun_500-1000-2000m_00-23h_by3hours_10KmRes_72h_backwards_lat_38.0444515708082_lon_15.4349519308499.pdf.RData" --n_lin 5000 --write_shp TRUE

#WD = "/home/rogerlm/Hysplit/Mediterranean_winds/"
#intersect_files = "/home/rogerlm/Hysplit/vcardui_hawaii/asia_east.shp"
#traj_Rdata = "/home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_21.1596368004464_lon_-156.710092919908.pdf.RData"
#hours_backwards = "-24,-48,-72"
#write_shp = T

# HAWAII1 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_18.9125336412606_lon_-155.679264300863.pdf.RData --write_shp TRUE
# HAWAII2 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_19.5234360095577_lon_-154.815816294495.pdf.RData --write_shp TRUE
# HAWAII3 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_19.7284219581233_lon_-156.064903333596.pdf.RData --write_shp TRUE
# HAWAII4 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_20.1277021421099_lon_-155.556125259818.pdf.RData --write_shp TRUE
# HAWAII5 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_20.7710483296044_lon_-155.978632268912.pdf.RData --write_shp TRUE
# HAWAII6 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_20.9251564707922_lon_-156.697557807946.pdf.RData --write_shp TRUE
# HAWAII7 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_21.0991733192148_lon_-157.310303575061.pdf.RData --write_shp TRUE
# HAWAII8 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_21.1596368004464_lon_-156.710092919908.pdf.RData --write_shp TRUE
# HAWAII9 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_21.3078460654165_lon_-157.65059499687.pdf.RData --write_shp TRUE
# HAWAII10 Rscript /home/rogerlm/Hysplit/Rdata_to_shp.R --WD /home/rogerlm/Hysplit/vcardui_hawaii/ --hours_backwards -24,-48,-72 --intersect_files /home/rogerlm/Hysplit/vcardui_hawaii/america_west.shp --traj_Rdata /home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_21.5747702142681_lon_-158.283618150984.pdf.RData --write_shp TRUE

#WD = "/home/rogerlm/Hysplit/vcardui_hawaii/"
#intersect_files = "/home/rogerlm/Hysplit/vcardui_hawaii/asia_east.shp"
#traj_Rdata = "/home/eric/vcardui_hawaii/hawaii_vcardui_1980-2022_by3hours_Apr-Jun_500-1000-2000m_00-23h_10KmRes_72h_backwards_lat_21.1596368004464_lon_-156.710092919908.pdf.RData"
#hours_backwards = "-24,-48,-72"
#write_shp = T

## ERIC GUYANE

#WD = "/home/eric/guyane/"
#intersect_files = "/home/eric/guyane/africa_atlantic_coast.shp"
#traj_Rdata = "/home/eric/guyane/guyane_vcardui_2013_byhour_Oct_500-1000-2000m_00-23h_10KmRes_200h_backwards_daily_21-22.pdf.RData"
#hours_backwards = 200
#write_shp = T

## SET COMMAND LINE ENVIRONMENTS
args <- commandArgs(asValues=TRUE)
str(args)
# Turn arguments into R variables
args <- commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]
keys <- attachLocally(args)
cat("\n Command-line arguments attached to global environment:\n\n");

str(mget(keys, envir=globalenv()))

setwd = WD
dir.create(file.path(WD, "test"), showWarnings = FALSE)
dir.create(file.path(WD, "intersections"), showWarnings = FALSE)
dir.create(file.path(WD, "trajectories"), showWarnings = FALSE)

extractorRData <- function(file, object) {
    #' Function for extracting an object from a .RData file created by R's save() command
    #' Inputs: RData file, object name
    E <- new.env()
    load(file=file, envir=E)
    return(get(object, envir=E, inherits=F))
  }
traj <- extractorRData(traj_Rdata, "trajs")[[1]]
#load(traj_Rdata)
cat("\n")
cat("Data loaded: \n")
write_shp = as.logical(write_shp)
hours = max(abs(traj$hour.inc))+1
hours_backwards = if(exists("hours_backwards")){
  -abs(as.numeric(strsplit(as.character(hours_backwards),",")[[1]]))
}else{
  -as.numeric(hours-1)
}
hours_backwards = if(min(hours_backwards) < 0){
  sort(hours_backwards,decreasing = F)
}else{
  sort(hours_backwards,decreasing = T)
}
hours_backwards = hours_backwards[which(abs(hours_backwards) <= hours)]

cat(" Hours backwards:",hours_backwards,"\n")
if (exists("n_lin")){
  nlin = as.numeric(n_lin)
}else{
  nlin = nrow(traj)/hours
}
if (nlin > nrow(traj)/hours){
  nlin = nrow(traj)/hours
}
if (nlin%%1 != 0){
  nlin = round(nlin)-1
}
cat(" Trajectories:",nlin," \n")
if (exists("intersect_files")){
intersect_files = strsplit(intersect_files,",")[[1]]
cat(" Intersect files: ",intersect_files,"\n\n")
}
traj = traj[1:(nlin*hours),] 
traj$unique = paste0(traj$date,"_",traj$start_height)
cat("Step 0: Variables read \n")
cat(" Input file has",nrow(traj)/hours,"trajectories. From theese,",nlin,"will be processed","\n")
cat(" Trajectories last",hours-1,"hours","\n")

dates = data.frame(date = unique(traj$unique),id = 1:length(unique(traj$unique)))
cat(" Unique identifiers to each line created \n\n")
#fun5 = function(x,dates){
#which(x == dates$date)
#}

#traj$id = vapply(traj$unique,fun5,dates=dates,FUN.VALUE=1)

# ######################  METHOD 1  ######################
# lines = 20
# traj$geom = c("LINESTRING(")
# for (i in (1:nrow(traj))[1:(lines*73)]){
#   id = traj$id[i]
#   n = which(traj$id == id)[1]
#   traj$geom[n] = paste0(traj$geom[n]," ",traj$lat[i]," ",traj$lon[i],",")
#   cat(i,"\n")
# }
# 
# traj$geom = paste0(substr(traj$geom,start=1,stop=nchar(traj$geom)-1),")")
# traj2 = traj[(c(0,73*(1:(nrow(traj)/73)))+1)[1:lines],]
# traj3 = st_as_sf(traj2,wkt="geom")
# 

######################  METHOD 2  ######################

cat("Step 1: Making trajectories points and lines","\n")
p = traj[which(traj$unique == dates$date[which(dates$id==1)]),]  # & traj$hour.inc >= hours_backwards[1]
## from the sp vignette:
l <- as.data.frame(cbind(p$lon, p$lat))
colnames(l) = c("lon","lat")
S0 <- terra::vect(x=p,geom=c("lon","lat"),keepgeom=T)
S0$ID=1
crs(S0) = "epsg:4326"
SL0 = as.lines(S0)

#Extract IDs for file names with the same number of characters
lat = if(substr(l[1,2],1,1) == "-"){
  if(nchar(abs(as.numeric(l[1,2]))) - (nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2]))) < 4){ #returns units +1
      paste0(rep(0,times = 3-(nchar(abs(as.numeric(l[1,2]))) - (nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2]))))),abs(l[1,2]),rep(0,times=3-(nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2])))),"S")
  }else{
      paste0(rep(0,times = 7-nchar(abs(as.numeric(l[1,2])))),abs(l[1,2]),rep(0,times=3-(nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2])))),"S")
  }
}else{
   if(nchar(abs(as.numeric(l[1,2]))) - (nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2]))) < 4){
      paste0(rep(0,times = 3-(nchar(abs(as.numeric(l[1,2]))) - (nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2]))))),abs(l[1,2]),rep(0,times=3-(nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2])))),"N")
  }else{
      paste0(rep(0,times = 7-nchar(abs(as.numeric(l[1,2])))),abs(l[1,2]),rep(0,times=3-(nchar(l[1,2]) - unlist(gregexpr('\\.',l[1,2])))),"N")
  }
}

lon = if(substr(l[1,1],1,1) == "-"){
  if(nchar(abs(as.numeric(l[1,1]))) - (nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1]))) < 5){ #returns units +1
      paste0(rep(0,times = 4-(nchar(abs(as.numeric(l[1,1]))) - (nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1]))))),abs(l[1,1]),rep(0,times=3-(nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1])))),"W")
  }else{
      paste0(rep(0,times = 7-nchar(abs(as.numeric(l[1,1])))),abs(l[1,1]),rep(0,times=3-(nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1])))),"W")
  }
}else{
   if(nchar(abs(as.numeric(l[1,1]))) - (nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1]))) < 5){
      paste0(rep(0,times = 4-(nchar(abs(as.numeric(l[1,1]))) - (nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1]))))),abs(l[1,1]),rep(0,times=3-(nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1])))),"E")
  }else{
      paste0(rep(0,times = 7-nchar(abs(as.numeric(l[1,1])))),abs(l[1,1]),rep(0,times=3-(nchar(l[1,1]) - unlist(gregexpr('\\.',l[1,1])))),"E")
  }
}

for (h in hours_backwards){
  SL0_path = gsub("-","_",paste0(WD,"trajectories/trajectories_lines_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",abs(h),"_",lon,"_",lat,".shp"))
  S0_path = gsub("-","_",paste0(WD,"trajectories/trajectories_points_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",abs(h),"_",lon,"_",lat,".shp"))

  if(file.exists(SL0_path)){

    SL0 = vect(SL0_path)
    S0 = vect(S0_path)

    cat("File exist: Loading \n",SL0_path,"\n",S0_path,"\n\n")
  }else{


    pb = txtProgressBar(min = 0, max = length((2:(length(unique(traj$unique))))), initial = 0) 
    for (i in (2:(length(unique(traj$unique))))){

      i = as.numeric(i)
      p = traj[which(traj$unique == dates$date[which(dates$id==i)]),] # & traj$hour.inc >= hours_backwards[1]
      S0tmp = vect(x=p, geom = c("lon","lat"),keepgeom=T)
      S0tmp$ID = i
      S0 = rbind(S0, S0tmp)
      SL0 = rbind(SL0,as.lines(S0tmp))
      setTxtProgressBar(pb,i)
    }
  close(pb)

    h = abs(h)


    

    if (exists("write_shp")){

      if (write_shp==T){

        cat("Writting trajectories shapefile: \n")
        writeVector(SL0,gsub("-","_",paste0(WD,"trajectories/trajectories_lines_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",h,"_",lat,"_",lon,".shp")),overwrite = T)
        writeVector(S0,gsub("-","_",paste0(WD,"trajectories/trajectories_points_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",h,"_",lat,"_",lon,".shp")),overwrite = T)
        cat(gsub("-","_",paste0("FILE WRITTEN: trajectories_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",h,"_",lat,"_",lon,".shp")),"\n\n")
      }
    }
  }

# trjectories written in file in here. Now intersection points
  cat("Step2: Searching intersection points in ",SL0_path,"\n")
  f12 = function(df){
  length(df[which(df == intrs0$ID)])
  }

  if (exists("intersect_files")){

    for (f in intersect_files){
        name0 = strsplit(f,"/")
        name = (name0[[1]][length(name0[[1]])])
      if (exists(paste0(WD,"intersections/intersections_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",abs(h),"_",lon,"_",lat,name))){
        cat("File_exists:",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",abs(h),"_",lat,"_",lon,name, "\n")
      }else{
        int = aggregate(vect(f))
        crs(SL0) = crs(int)
        intrs0 = intersect(SL0,int)
        intrs0 = intrs0[,which(names(intrs0) %in% names(SL0))]
        # Chech n and i meaning
        for(i in which(relate(int,SL0,relation="intersects"))){
          cat(i,"\n")
          n = which(i == which(relate(int,SL0,relation="intersects")))
          t = nearest(intrs0[n],S0[which(S0$ID == i)])
          intrs0$int_date[n] = as.character(S0$date2[which(S0$ID == i & geom(S0)[,c("x","y")] %in% geom(t)[,c("x","y")] )][1])
          intrs0$start_height[n] = S0$start_height[which(S0$ID == i & geom(S0)[,c("x","y")] %in% geom(t)[,c("x","y")] )][1]
          intrs0$int_wind_speed[n] = S0$wind_speed[which(S0$ID == i & geom(S0)[,c("x","y")] %in% geom(t)[,c("x","y")] )][1]
          intrs0$start_date[n] = as.character(S0$date[which(S0$ID == i & geom(S0)[,c("x","y")] %in% geom(t)[,c("x","y")] )][1])
          intrs0$height[n] = S0$height[which(S0$ID == i & geom(S0)[,c("x","y")] %in% geom(t)[,c("x","y")] )][1]
          intrs0$ID[n] = S0$ID[which(S0$ID == i & geom(S0)[,c("x","y")] %in% geom(t)[,c("x","y")] )][1]
        
        }
      
        intrs0$REP = unname(unlist(as.vector(data.frame(lapply(intrs0$ID,f12))[1,])))
        writeVector(intrs0,paste0(WD,"intersections/intersections_",substr(min(traj$date),1,10),"_",substr(max(traj$date),1,10),"_",abs(h),"_",lon,"_",lat,"_",name),overwrite=T)
        cat("\n WRITTING FILES:",h,"\n")
      }
    }
  }
}
 
cat("\nEND OFF PROCESSING\n")
