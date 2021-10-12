par(bg = 'grey')
options(digits=10)
library(plyr)
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

all_transects <- c('RG' , 'GL' , 'WM', 'TR' , 'NP' , 'CL' , 'LG' , 'MB' , 'HM'  , 'CT' , 'CB' , 'SH' , 'DC' , 'VB' , 'BG' , 'CM' , 'BC' , 'WS' , 'UL' , 'RL' , 'ER' , 'ME' , 'OBJ' , 'TC')

# Reads in transect file: 
flax_GPS_raw <- read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.csv")

# Bounds the spatial datasets to speed up computation time: (somewhat arbitrarily at the moment)
area_of_interest <- extent(matrix(c(315000,335000, 4300000,433000), nrow=2,byrow=TRUE))

topography_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release1/UER_dem_filled_1m_v2.tif"
topography_adjusted <- paste("/vsicurl/",topography_raw,sep="")
topography <- raster(topography_adjusted, progress='text')
#topography_cropped <- crop(topography, area_of_interest, filename=tempfile(),progress="text")


landcover_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release1/UER_landcover_1m_v4.tif"
landcover_adjusted <- paste("/vsicurl/",landcover_raw,sep="")
landcover <- raster(landcover_adjusted, progress='text')
# landcover_cropped <- crop(landcover, area_of_interest, filename=tempfile(), progress="text")


surface_water_cover_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release1/UER_surface_water_1m_v3.tif"
surface_water_cover_adjusted <- paste("/vsicurl/",surface_water_cover_raw,sep="")
surface_water_cover <- raster(surface_water_cover_adjusted, progress='text')
#surface_water_cover_cropped <- crop(surface_water_cover, area_of_interest, filename=tempfile(),progress="text")


summer_solar_radiation_adjusted <- "/vsicurl/https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release2/UER_srad_bareearth_day172_3m_v1.tif"
summer_solar_radiation <- raster(summer_solar_radiation_adjusted, progress='text')

clearPlot()
Plot(topography)
Plot(landcover)
Plot(surface_water_cover)
Plot(summer_solar_radiation)




