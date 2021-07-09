library(leaflet)
library(plotKML)

data<-read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.csv")
metadata<-read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.metadata.csv")

#for (i in 1:nrow(metadata))
#{
#  assign(paste0(metadata[i,"tag"],".path.data"),readGPX(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/",metadata[i,"tag"],".gpx"))$tracks[[1]][[1]])
#}

RG.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/RG.gpx")$tracks[[1]][[1]]
TR.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/TR.gpx")$tracks[[1]][[1]]
WM.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/WM.gpx")$tracks[[1]][[1]]
NP.path.data<-readGPX("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/NP.gpx")$tracks[[1]][[1]]


RG<-data[which(data$transect=="RG"),]
TR<-data[which(data$transect=="TR"),]
WM<-data[which(data$transect=="WM"),]
NP<-data[which(data$transect=="NP"),]


map<-leaflet() %>% 
  setView(lng = mean(c(data$start.long,data$end.long)), lat = mean(c(data$start.lat,data$end.lat)), zoom = 13) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap)

map<-addPolylines(map,lng = RG.path.data$lon,lat=RG.path.data$lat,colo="black",weight=5)
map<-addPolylines(map,lng = TR.path.data$lon,lat=TR.path.data$lat,colo="black",weight=5)
map<-addPolylines(map,lng = WM.path.data$lon,lat=WM.path.data$lat,colo="black",weight=5)
map<-addPolylines(map,lng = NP.path.data$lon,lat=NP.path.data$lat,colo="black",weight=5)


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


for(i in 1:nrow(NP))
{
  map<-addPolylines(map,lng=c(NP$start.long[i],NP$end.long[i]),lat=c(NP$start.lat[i],NP$end.lat[i]),color = "blue",opacity = .75)
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

for(i in 1:nrow(NP))
{
  map<-addCircles(map,lng=NP$start.long[i],lat=NP$start.lat[i],radius=5*(NP$num.H[i]+NP$num.D[i]),fillOpacity=0,color = switch(NP$incidence[i]+1,"blue","yellow"))
}


map
