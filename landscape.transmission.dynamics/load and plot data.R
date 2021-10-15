# load packages
library(dplyr)
library(mgcv) 
library(tidyverse)
library(ggplot2)
library(rgdal)
library(raster)
library(ggspatial)
library(rasterVis)
library(gridExtra)
library(geosphere)
library(sf) 
library(sp)
library(leaflet)
library(spatstat)
library(plotKML)
library(sampSurf)
library(quickPlot)
library(exactextractr)
library(Rcpp)

# load spatial data

## define area of interest
area_of_interest <- extent(matrix(c(315000,335000, 4300000,433000), nrow=2,byrow=TRUE)) # Bounds the spatial datasets to speed up computation time: (somewhat arbitrarily at the moment)

topography_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release1/UER_dem_filled_1m_v2.tif"
topography_adjusted <- paste("/vsicurl/",topography_raw,sep="")
topography <- raster(topography_adjusted, progress='text')
topography_cropped <- crop(topography, area_of_interest, filename=tempfile(),progress="text")

landcover_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release1/UER_landcover_1m_v4.tif"
landcover_adjusted <- paste("/vsicurl/",landcover_raw,sep="")
landcover <- raster(landcover_adjusted, progress='text')
landcover_cropped <- crop(landcover, area_of_interest, filename=tempfile(), progress="text")

surface_water_cover_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release1/UER_surface_water_1m_v3.tif"
surface_water_cover_adjusted <- paste("/vsicurl/",surface_water_cover_raw,sep="")
surface_water_cover <- raster(surface_water_cover_adjusted, progress='text')
surface_water_cover_cropped <- crop(surface_water_cover, area_of_interest, filename=tempfile(),progress="text")

summer_solar_radiation_adjusted <- "/vsicurl/https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release2/UER_srad_bareearth_day172_3m_v1.tif"
summer_solar_radiation <- raster(summer_solar_radiation_adjusted, progress='text')
summer_solar_radiation_cropped <- crop(summer_solar_radiation, area_of_interest, filename=tempfile(),progress="text")

# visualize spatial data
clearPlot()
Plot(topography)
Plot(landcover)
Plot(surface_water_cover)
Plot(summer_solar_radiation)

# load gpx data and population data

#all_transects <- c('RG' , 'GL' , 'WM', 'TR' , 'NP' , 'CL' , 'LG' , 'MB' , 'HM'  , 'CT' , 'CB' , 'SH' , 'DC' , 'VB' , 'BG' , 'CM' , 'BC' , 'WS' , 'UL' , 'RL' , 'ER' , 'ME' , 'OBJ' , 'TC') ## transects to consider
all_transects<-c("UL")

## function to extract and clean gpx data
get_gpx_tracks <- function(transects = all_transects) {
  gpx_tracks<-data.frame("transect"=character(),"longitude"=numeric(),"latitude"=numeric())
  filenames <- paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/cleaned/",transects, ".gpx") 
  for (transect in transects) {
    file_name<- paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/cleaned/",transect, "_cleaned.gpx") 
    GPX_points <- readGPX(file_name)$tracks[[1]][[1]] %>% dplyr::select(lon, lat) %>% dplyr::rename(Latitude = lat, Longitude = lon)
    new_gpx_tracks<-data.frame("transect"=transect,"longitude"= GPX_points$Longitude,"latitude"=GPX_points$Latitude)
    gpx_tracks<-rbind(gpx_tracks,new_gpx_tracks)
  }
  gpx_tracks
}

gpx_tracks <- get_gpx_tracks() ## retrive data
flax_pops <- read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.csv") ## load flax population data

# extract landscape data and plot flax populations

## initiate map
map<-leaflet() %>% 
  addTiles() %>% 
  addMeasure(primaryLengthUnit = "meters")

## plotting option
plot<-T

## loop to extract data

for(transect in all_transects) ### for each transect
{
  sub_gpx_tracks<-gpx_tracks[which(gpx_tracks$transect==transect),] #### subset gpx data
  flax_sub_pop<-flax_pops[which(flax_pops$transect=="UL"),] #### subset flax population data
  
  line<-st_linestring(cbind(sub_gpx_tracks$longitude,sub_gpx_tracks$latitude)) #### draw line through all gps points
  n.points<-floor(4*as.numeric(st_length(st_sfc(line, crs = "+proj=longlat +datum=WGS84")))) ### set number of points to sample such that sampling will occur every 0.25m (4* distance of gps track in meters)
  points<-st_line_sample(line,n=n.points,type="regular") #### sample points at regular interval along track
  points<-st_sfc(points,crs = "+proj=longlat +datum=WGS84") #### add crs to pio ts
  
  new_line<-st_linestring(st_coordinates(points)[,c("X","Y")]) #### convert points to line
  new_line<-st_sfc(new_line, crs = "+proj=longlat +datum=WGS84") #### add crs
  
  if(plot) {map<-addPolylines(map=map,data = st_transform(new_line,crs="+proj=longlat +datum=WGS84"),col="black",weight=1)} #### add track to map
  
  track_points<-st_coordinates(new_line)[,c("X","Y")] #### extract raw track coordinates
  transect.coded<-data.frame(track_points,"class"="nfz","chunk"=NA) #### set up data frame for flax populations, default class to 'no flax zone' = 'nfz
  
  #### subdivide gps data into 'chunks' of no flax/healthy flax/diseased flax
  
  chunk<-1 #### set starting chunk index to 1
  
  for(i in 1:nrow(flax_sub_pop)) ##### for each flax population
  {
    start.cords<-as.numeric(flax_sub_pop[i,c("start.long","start.lat")]) ##### get starting coordinates
    start.point<-st_sfc(st_point(start.cords),crs = "+proj=longlat +datum=WGS84") ##### set crs
    start.dists<-c(st_distance(st_cast(points,"POINT"),start.point)) ##### get distances between gps coordinates along transect line and population start point
    start.index<-which.min(start.dists) ##### get index of minimun distance
    
    end.cords<-as.numeric(flax_sub_pop[i,c("end.long","end.lat")]) ##### get end coordinates
    end.point<-st_sfc(st_point(end.cords),crs = "+proj=longlat +datum=WGS84") ##### set crs
    end.dists<-c(st_distance(st_cast(points,"POINT"),end.point)) ##### get distances between gpx coordinates along transect line and population end point
    end.index<-which.min(end.dists) ##### get index of minimum distance
    
    if(flax_sub_pop[i,"incidence"]==0) {new.status<-"fz"} ##### set class to 'flax zone' if there is no disease
    if(flax_sub_pop[i,"incidence"]==1) {new.status<-"dfz"} ##### set class to 'diseased flax zone' if there is disease
    
    transect.coded[start.index:(end.index-40),"class"]<-new.status ##### set class of chunk containing flax pop, backing up 10m from recorded end pont (per sampling protocol)
    transect.coded[which(is.na(transect.coded[1:start.index,"chunk"])),"chunk"]<-chunk ##### set chunk index of chunk with class nfz preceeding flax population
    transect.coded[start.index:(end.index-40),"chunk"]<-chunk+1 ##### set chunk index of flax population , backing up 10m from recorded end pont (per sampling protocol)
    chunk<-chunk+ 2 ##### set next hunk index
    if(i==nrow(flax_sub_pop)) {transect.coded[which(is.na(transect.coded$chunk)),"chunk"]<-chunk} ##### after the last flax population, set chunk of remaining nfz if it exists
  }
  
  #### for each chunk, extract landscape data and plot 
  for(chunk in 1:max(transect.coded$chunk))
  {
    sub.transect.coded<-transect.coded[which(transect.coded$chunk==chunk),] ##### subset data
  
    if(sub.transect.coded[1,"class"]=="nfz") {col<-"grey"} ##### set plot color based on class
    if(sub.transect.coded[1,"class"]=="fz") {col<-"green"}
    if(sub.transect.coded[1,"class"]=="dfz") {col<-"yellow"}
    
    
    line<-st_linestring(cbind(sub.transect.coded$X,sub.transect.coded$Y)) ##### draw line through all gps points
    line<-st_sfc(line, crs = "+proj=longlat +datum=WGS84") ##### convert to coordinates
    line_mod <- st_transform(line, crs = 3501) ##### set crs
    ribbon <- st_buffer(line_mod, dist = 2.5, endCapStyle = "FLAT",joinStyle = "ROUND") ##### draw polygon around line
  
    if(plot) {map<-addPolygons(map=map,data = st_transform(ribbon,crs="+proj=longlat +datum=WGS84"),col=col)} ##### add polygon to map
    
    sample_points<-st_sample(ribbon,size=as.numeric(floor(st_area(ribbon))),type="regular") ##### regularly sample 1 point per m2 from polygon 
    
    #if(plot) {map<-addCircleMarkers(map=map,data=st_transform(sample_points,crs="+proj=longlat +datum=WGS84"),col="black",radius = .1,weight=0)}
    
    unlist(raster::extract(x=topography,y=as(sample_points,"Spatial"),method="bilinear"))->x
    length(x)
    dim(st_coordinates(sample_points))
    }
}


raster::extract(x=topography,y=as(ribbon,"Spatial"),method="bilinear")->x


