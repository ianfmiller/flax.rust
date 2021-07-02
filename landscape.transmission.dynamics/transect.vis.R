library(leaflet)
library(plotKML)

pop.data<-read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/prelim.landscape.transects.csv")
path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/RG.gpx")$tracks[[1]][[1]]

RG<-pop.data[which(pop.data$transect=="RG"),]

map<-leaflet() %>% 
  setView(lng = mean(c(RG$start.long,RG$end.long)), lat = mean(c(RG$start.lat,RG$end.lat)), zoom = 13) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap)

map<-addPolylines(map,lng = path.data$lon,lat=path.data$lat,colo="black",weight=5)

for(i in 1:nrow(RG))
{
  map<-addPolylines(map,lng=c(RG$start.long[i],RG$end.long[i]),lat=c(RG$start.lat[i],RG$end.lat[i]),color = "blue",opacity = .75)
}

for(i in 1:nrow(RG))
{
  map<-addCircles(map,lng=RG$start.long[i],lat=RG$start.lat[i],radius=5*(RG$num.H[i]+RG$num.D[i]),fillOpacity=0,color = switch(RG$incidence[i]+1,"blue","yellow"))
}


map
