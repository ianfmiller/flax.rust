library(leaflet)
library(plotKML)

data<-read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.csv")
metadata<-read.csv("~/Documents/GitHub/flax.rust/data/landscape.transect.data/landscape.transects.metadata.csv")

#for (i in 1:nrow(metadata))
#{
#  assign(paste0(metadata[i,"tag"],".path.data"),readGPX(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/",metadata[i,"tag"],".gpx"))$tracks[[1]][[1]])
#}

transects<-c("RG","TR","WM","NP","OBJ","BG","CM","BC","CL","LG","HM","CT","DC","UL","ER","RL")

for(transect in transects)
{
  assign(paste0(transect,".path.data"),readGPX(paste0("~/Documents/GitHub/flax.rust/data/landscape.transect.data/gpx/",transect,".gpx"))$tracks[[1]][[1]])
}

for(transect in transects)
{
  assign(paste0(transect),data[which(data[,"transect"]==transect),])
}

map<-leaflet() %>% 
  setView(lng = mean(c(data$start.long,data$end.long)), lat = mean(c(data$start.lat,data$end.lat)), zoom = 10) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap)

for(transect in transects)
{
  dat.obj<-eval(parse(text=paste0(transect,".path.data")))
  map<-addPolylines(map,lng = dat.obj$lon,lat=dat.obj$lat,colo="black",weight=5)
}

for(transect in transects)
{
  dat.obj<-eval(parse(text=paste0(transect)))
  for (i in 1:nrow(dat.obj))
  {
    map<-addPolylines(map,lng=c(dat.obj$start.long[i],dat.obj$end.long[i]),lat=c(dat.obj$start.lat[i],dat.obj$end.lat[i]),color = "blue",opacity = .75)
  }
}

for(transect in transects)
{
  dat.obj<-eval(parse(text=paste0(transect)))
  for (i in 1:nrow(dat.obj))
  {
    if(nrow(dat.obj)>=1) 
    {
      map<-addCircles(map,lng=dat.obj$start.long[i],lat=dat.obj$start.lat[i],radius=5*(dat.obj$num.H[i]+dat.obj$num.D[i]),fillOpacity=0,color = switch(dat.obj$incidence[i]+1,"blue","yellow"))
    }
  }
}

map
