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

# load gpx data

## transects to consider
#all_transects <- c('RG' , 'GL' , 'WM', 'TR' , 'NP' , 'CL' , 'LG' , 'MB' , 'HM'  , 'CT' , 'CB' , 'SH' , 'DC' , 'VB' , 'BG' , 'CM' , 'BC' , 'WS' , 'UL' , 'RL' , 'ER' , 'ME' , 'OBJ' , 'TC')
all_transects<-c("UL")

## function to extract and clean data
get_gpx_tracks <- function(names = all_transects) {
  gpx_tracks<-data.frame("transect"=character(),"longitude"=numeric(),"latitude"=numeric())
  filenames <- paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/cleaned/",names, ".gpx") 
  for (name in names) {
    file_name<- paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/cleaned/",name, "_cleaned.gpx") 
    GPX_points <- readGPX(file_name)$tracks[[1]][[1]] %>% dplyr::select(lon, lat) %>% dplyr::rename(Latitude = lat, Longitude = lon)
    new_gpx_tracks<-data.frame("transect"=name,"longitude"= GPX_points$Longitude,"latitude"=GPX_points$Latitude)
    gpx_tracks<-rbind(gpx_tracks,new_gpx_tracks)
  }
  gpx_tracks
}

gpx_tracks <- get_gpx_tracks()

# visualize gpx data

line<-st_sfc(st_linestring(cbind(gpx_tracks$longitude,gpx_tracks$latitude))) ## draw line through all gps points
n.points<-floor(as.numeric(st_length(st_sfc(line, crs = "+proj=longlat +datum=WGS84")))) ## set number of points to sample as distance of gps track in meters
points<-st_line_sample(line,n=n.points,type="regular") #3 sample points at a reular interval along track

line<-st_sfc(line, crs = "+proj=longlat +datum=WGS84") # convert to coordinates
points<-st_sfc(points, crs = "+proj=longlat +datum=WGS84")  # convert to coordinates

plot(line) ## plot original track
plot(points,col="red",add=T) ## plot sampled points
plot(st_sfc(st_cast(points,"LINESTRING"),crs = "+proj=longlat +datum=WGS84"),col="blue",add=T) ## plot track made from sampled points


## proof of concept illustrated with fewer points
#line<-st_sfc(st_linestring(cbind(gpx_tracks$longitude,gpx_tracks$latitude)))
#points<-st_line_sample(line,density = 100,type="regular")
#plot(st_sfc(line, crs = "+proj=longlat +datum=WGS84"))
#plot(st_sfc(points, crs = "+proj=longlat +datum=WGS84"),col="red",add=T)
#plot(st_sfc(st_cast(points,"LINESTRING"),crs = "+proj=longlat +datum=WGS84"),col="blue",add=T)



test.point<-st_sfc(st_point(cbind(gpx_tracks$longitude[1],gpx_tracks$latitude[1])),crs = "+proj=longlat +datum=WGS84")
st_distance(st_cast(points,"POINT"),test.point)

x<-st_distance(test.point,st_cast(points,"POINT"))
