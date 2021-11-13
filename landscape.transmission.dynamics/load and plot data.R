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
topography_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release3/UG_dem_1m_v1.tif"
topography_adjusted <- paste("/vsicurl/",topography_raw,sep="")
topography <- raster(topography_adjusted, progress='text')

landcover_raw <- "https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release3/UG_landcover_1m_v4.tif"
landcover_adjusted <- paste("/vsicurl/",landcover_raw,sep="")
landcover <- raster(landcover_adjusted, progress='text')

# visualize spatial data
clearPlot()
Plot(topography)
Plot(landcover)

# load gpx data and population data

#transects <- c('RG' , 'GL' , 'WM', 'TR' , 'NP' , 'CL' , 'LG' , 'MB' , 'HM'  , 'CT' , 'CB' , 'SH' , 'DC' , 'VB' , 'BG' , 'CM' , 'BC' , 'WS' , 'UL' , 'RL' , 'ER' , 'ME' , 'OBJ' , 'TC') ## transects to consider
transects<-c("UL","RL","SH","TC","TR","VB","WM","WS")

## function to extract and clean gpx data
get_gpx_tracks <- function(transect = transects) {
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

## setup data objects
all_transects<-data.frame("transect"=character(),"chunk"=integer(),"flax.presence"=integer(),"incidence"=numeric(),"elevation"=numeric(),"landcover"=numeric())
all_populations<-data.frame("transect"=character(),"chunk"=integer(),"density"=numeric(),"num.H"=numeric(),"num.D"=numeric(),"incidence"=numeric(),"nearest.pop.dist"=numeric(),"nearest.D.pop.dist"=numeric(),"elevation"=numeric(),"mode.landcover"=numeric(),"p.landcover.1"=numeric(),"p.landcover.2"=numeric(),"p.landcover.3"=numeric(),"p.landcover.4"=numeric(),"p.landcover.5"=numeric(),"p.landcover.6"=numeric())

## loop to extract data
for(transect in transects) ### for each transect
{
  print(paste0("starting transect ",transect))
  
  sub_gpx_tracks<-gpx_tracks[which(gpx_tracks$transect==transect),] #### subset gpx data
  flax_sub_pop<-cbind(flax_pops[which(flax_pops$transect==transect),],"chunk"=NA) #### subset flax population data
  
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
    
    if(isTRUE(flax_sub_pop[i,"incidence"]==0)) {new.status<-"fz"} ##### set class to 'flax zone' if there is no disease
    if(isTRUE(flax_sub_pop[i,"incidence"]==1)) {new.status<-"dfz"} ##### set class to 'diseased flax zone' if there is disease
    
    transect.coded[start.index:(end.index-40),"class"]<-new.status ##### set class of chunk containing flax pop, backing up 10m from recorded end pont (per sampling protocol)
    transect.coded[which(is.na(transect.coded[1:start.index,"chunk"])),"chunk"]<-chunk ##### set chunk index of chunk with class nfz preceeding flax population
    transect.coded[start.index:(end.index-40),"chunk"]<-chunk+1 ##### set chunk index of flax population , backing up 10m from recorded end pont (per sampling protocol)
    flax_sub_pop[i,"chunk"]<-chunk+1
    chunk<-chunk+ 2 ##### set next hunk index
    if(i==nrow(flax_sub_pop)) {transect.coded[which(is.na(transect.coded$chunk)),"chunk"]<-chunk} ##### after the last flax population, set chunk of remaining nfz if it exists
    
  }
  
  print(paste0("finished gps processing"))
  
  
  #### for each chunk, extract landscape data and plot 
  for(chunk in unique(transect.coded$chunk))
  {
    sub.transect.coded<-transect.coded[which(transect.coded$chunk==chunk),] ##### subset data
  
    if(sub.transect.coded[1,"class"]=="nfz") {col<-"grey"; flax.presence<-0; incidence<-0} ##### set plot color based on class
    if(sub.transect.coded[1,"class"]=="fz") {col<-"green"; flax.presence<-1; incidence<-0}
    if(sub.transect.coded[1,"class"]=="dfz") {col<-"yellow"; flax.presence<-1; incidence<-1}
    
    
    line<-st_linestring(cbind(sub.transect.coded$X,sub.transect.coded$Y)) ##### draw line through all gps points
    line<-st_sfc(line, crs = "+proj=longlat +datum=WGS84") ##### convert to coordinates
    line_mod <- st_transform(line, crs = 3501) ##### set crs
    ribbon <- st_buffer(line_mod, dist = 2.5, endCapStyle = "FLAT",joinStyle = "ROUND") ##### draw polygon around line
  
    if(plot) {map<-addPolygons(map=map,data = st_transform(ribbon,crs="+proj=longlat +datum=WGS84"),col=col)} ##### add polygon to map
    
    sample_points<-st_sample(ribbon,size=as.numeric(floor(st_area(ribbon))),type="regular") ##### regularly sample 1 point per m2 from polygon 
    sample_points<-st_transform(sample_points,crs=proj4string(topography)) ##### align crs of sample_points and data
    
    #if(plot) {map<-addCircleMarkers(map=map,data=st_transform(sample_points,crs="+proj=longlat +datum=WGS84"),col="black",radius = .1,weight=0)}
    
    if(!file.exists(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/raster.extracted/",transect,".elevation.RDS")) | !file.exists(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/raster.extracted/",transect,".landcover.RDS")) )
    {
      elevation_data<-unlist(raster::extract(x=topography,y=as(sample_points,"Spatial"),method="bilinear",na.rm=T)) ##### extract data
      landcover_data<-unlist(raster::extract(x=landcover,y=as(sample_points,"Spatial"),method="simple",na.rm=T))
      saveRDS(elevation_data,file=paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/raster.extracted/",transect,".chunk.",chunk,".elevation.RDS"))
      saveRDS(landcover_data,file=paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/raster.extracted/",transect,".chunk.",chunk,".landcover.RDS"))
    }
    
    elevation_data<-readRDS(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/raster.extracted/",transect,".chunk.",chunk,".elevation.RDS"))
    landcover_data<-readRDS(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/raster.extracted/",transect,".chunk.",chunk,".elevation.RDS"))
    
    ##### store data
    
    ###### make new data tables
    new_all_transects<-data.frame("transect"=transect,
                                  "chunk"=chunk,
                                  "flax.presence"=flax.presence,
                                  "incidence"=incidence,
                                  "elevation"=elevation_data,
                                  "landcover"=landcover_data)
    all_transects<-rbind(all_transects,new_all_transects)
    
    num.H<-if(chunk %in% flax_sub_pop$chunk) {flax_sub_pop[which(flax_sub_pop$chunk==chunk),"num.H"]} else {0} ###### number of H plants
    num.D<-if(chunk %in% flax_sub_pop$chunk) {flax_sub_pop[which(flax_sub_pop$chunk==chunk),"num.H"]} else {0}  ###### number of D plants
    prevalence<-if(chunk %in% flax_sub_pop$chunk) {num.D/(num.H+num.D)} ###### prevalence
    incidence<-if(chunk %in% flax_sub_pop$chunk) {flax_sub_pop[which(flax_sub_pop$chunk==chunk),"incidence"]} else {0} ###### incidence
    
    ###### get distance to nearest diseased and healthy populations by making use of point indicies. Points are spaced 1m apart, so the difference in indicies represents distance (in m) along transect between populations
    chunk.indicies<-which(transect.coded$chunk==chunk)
    chunk.endpoint.1<-min(chunk.indicies)
    chunk.endpoint.2<-max(chunk.indicies)
    
    flax.indicies.all<-which(transect.coded$class %in% c("dfz","fz"))
    flax.indicies<-flax.indicies.all[which(!(flax.indicies.all %in% chunk.indicies))]
    disease.indicies.all<-which(transect.coded$class=="dfz")
    disease.indicies<-disease.indicies.all[which(!(disease.indicies.all %in% chunk.indicies))]
    
    nearest.pop.dist<-min(min(abs(chunk.endpoint.1-flax.indicies)),min(abs(chunk.endpoint.2-flax.indicies)))
    nearest.D.pop.dist<-min(min(abs(chunk.endpoint.1-disease.indicies)),min(abs(chunk.endpoint.2-disease.indicies)))
    
    ###### extract landcover varriables
    landcover.table<-table(landcover_data)
    mode.landcover<-as.numeric(names(landcover.table)[which.max(landcover.table)])
    p.landcover.1<-if(1 %in% names(landcover.table)) {landcover.table[which(names(landcover.table)==1)]/sum(landcover.table)} else{0}
    p.landcover.2<-if(2 %in% names(landcover.table)) {landcover.table[which(names(landcover.table)==2)]/sum(landcover.table)} else{0}
    p.landcover.3<-if(3 %in% names(landcover.table)) {landcover.table[which(names(landcover.table)==3)]/sum(landcover.table)} else{0}
    p.landcover.4<-if(4 %in% names(landcover.table)) {landcover.table[which(names(landcover.table)==4)]/sum(landcover.table)} else{0}
    p.landcover.5<-if(5 %in% names(landcover.table)) {landcover.table[which(names(landcover.table)==5)]/sum(landcover.table)} else{0}
    p.landcover.6<-if(6 %in% names(landcover.table)) {landcover.table[which(names(landcover.table)==6)]/sum(landcover.table)} else{0}
    
    ###### insert new data tables into main tables
    new_all_populations<-data.frame("transect"=transect,
                                    "chunk"=chunk,
                                    "density"=num.H+num.D,
                                    "num.H"=num.H,
                                    "num.D"=num.D,
                                    "incidence"=incidence,
                                    "nearest.pop.dist"=nearest.pop.dist,
                                    "nearest.D.pop.dist"=nearest.D.pop.dist,
                                    "elevation"=mean(elevation_data),
                                    "mode.landcover"=mode.landcover,
                                    "p.landcover.1"=p.landcover.1,
                                    "p.landcover.2"=p.landcover.2,
                                    "p.landcover.3"=p.landcover.3,
                                    "p.landcover.4"=p.landcover.4,
                                    "p.landcover.5"=p.landcover.5,
                                    "p.landcover.6"=p.landcover.6)
    all_populations<-rbind(all_populations,new_all_populations)
    
    print(paste0("finished data extract chunk ",chunk," of ",max(transect.coded$chunk)))
  }
  print(paste0("finished transect ",transect))
}

map

head(all_transects)
head(all_populations)

saveRDS(all_transects,file="~/Documents/GitHub/flax.rust/landscape.transmission.dynamics/summarized.data/all_transects.RDS")
saveRDS(all_populations,file="~/Documents/GitHub/flax.rust/landscape.transmission.dynamics/summarized.data/all_populations.RDS")


