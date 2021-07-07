library(leaflet)
library(plotKML)

pop.data<-read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.csv")
RG.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/RG.gpx")$tracks[[1]][[1]]
TR.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/TR.gpx")$tracks[[1]][[1]]
WM.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/WM.gpx")$tracks[[1]][[1]]


RG<-pop.data[which(pop.data$transect=="RG"),]
TR<-pop.data[which(pop.data$transect=="TR"),]
WM<-pop.data[which(pop.data$transect=="WM"),]

map<-leaflet() %>% 
  setView(lng = mean(c(pop.data$start.long,pop.data$end.long)), lat = mean(c(pop.data$start.lat,pop.data$end.lat)), zoom = 13) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap)

map<-addPolylines(map,lng = RG.path.data$lon,lat=RG.path.data$lat,colo="black",weight=5)
map<-addPolylines(map,lng = TR.path.data$lon,lat=TR.path.data$lat,colo="black",weight=5)
map<-addPolylines(map,lng = WM.path.data$lon,lat=WM.path.data$lat,colo="black",weight=5)


for(i in 1:nrow(RG))
{
  map<-addPolylines(map,lng=c(RG$start.long[i],RG$end.long[i]),lat=c(RG$start.lat[i],RG$end.lat[i]),color = "blue",opacity = .75)
}

for(i in 1:nrow(TR))
{
  map<-addPolylines(map,lng=c(TR$start.long[i],TR$end.long[i]),lat=c(TR$start.lat[i],TR$end.lat[i]),color = "blue",opacity = .75)
}

for(i in 1:nrow(WM))
{
  map<-addPolylines(map,lng=c(WM$start.long[i],WM$end.long[i]),lat=c(WM$start.lat[i],WM$end.lat[i]),color = "blue",opacity = .75)
}


for(i in 1:nrow(RG))
{
  map<-addCircles(map,lng=RG$start.long[i],lat=RG$start.lat[i],radius=5*(RG$num.H[i]+RG$num.D[i]),fillOpacity=0,color = switch(RG$incidence[i]+1,"blue","yellow"))
}

for(i in 1:nrow(TR))
{
  map<-addCircles(map,lng=TR$start.long[i],lat=TR$start.lat[i],radius=5*(TR$num.H[i]+TR$num.D[i]),fillOpacity=0,color = switch(TR$incidence[i]+1,"blue","yellow"))
}
 
for(i in 1:nrow(WM))
{
  map<-addCircles(map,lng=WM$start.long[i],lat=WM$start.lat[i],radius=5*(WM$num.H[i]+WM$num.D[i]),fillOpacity=0,color = switch(WM$incidence[i]+1,"blue","yellow"))
}


map
