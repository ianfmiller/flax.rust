# load + prep data
## load data
epi<-read.csv("~/Documents/GitHub/flax.rust/data/Epidemiology.csv")
within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
plant.locs<-read.csv("~/Documents/GitHub/flax.rust/data/plant locs.csv")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")

## prep data
epi<-epi[which(epi$Year==2020),]
epi$Date.First.Observed.Diseased<-as.Date(epi$Date.First.Observed.Diseased,tryFormats = "%m/%d/%Y")
within.host<-within.host[which(within.host$Year==2020),]
within.host$Date<-as.Date(within.host$Date,tryFormats = "%m/%d/%Y")

# visualize
par(mfrow=c(1,2))
for (site in unique(epi$Site))
{
  sub.data.1<-epi[which(epi$Site==site),]
  sub.data.1<-sub.data.1[order(sub.data.1$Date.First.Observed.Diseased),]
  plant.loc.data<-corrected.plant.locs[which(corrected.plant.locs$Site==site),]
  plot(plant.loc.data$X+plant.loc.data$x,plant.loc.data$Y+plant.loc.data$y,xlim=c(0,10),ylim=c(0,20),main=site,xlab="X",ylab="Y")
  new.infs<-0
  for (date in unique(sub.data.1$Date.First.Observed.Diseased))
  {
    plot(plant.loc.data$X+plant.loc.data$x,plant.loc.data$Y+plant.loc.data$y,xlim=c(0,10),ylim=c(0,20),main=paste0(site,date),xlab="X",ylab="Y")
    print(as.Date(date,origin = "1970-01-01"))
    sub.data.2<-sub.data.1[which(sub.data.1$Date.First.Observed.Diseased==date),]
    points(sub.data.2$X+sub.data.2$x,sub.data.2$Y+sub.data.2$y,col="red")
    sub.data.2
    #Sys.sleep(2)
    new.infs<-c(new.infs,new.infs[length(new.infs)]+dim(sub.data.2)[1])
  }
}

epi<-read.csv("~/Documents/GitHub/flax.rust/data/Epidemiology.csv")
epi<-epi[which(epi$Year==2020),]
plant.locs<-read.csv("~/Documents/GitHub/flax.rust/data/plant locs.csv")
plant.locs<-plant.locs[which(plant.locs$Year==2020),]


for(tag in epi$Tag)
{
  epi.index<-which(epi$Tag==tag)
  epiX<-epi[epi.index,"X"]
  epix<-epi[epi.index,"x"]
  epiY<-epi[epi.index,"Y"]
  epiy<-epi[epi.index,"y"]
  
  loc.index<-which(plant.locs$tag==tag)
  locX<-plant.locs[loc.index,"X"]
  locx<-plant.locs[loc.index,"x"]
  locY<-plant.locs[loc.index,"Y"]
  locy<-plant.locs[loc.index,"y"]
  
  if (any((!(epiX==locX)),(!(epix==locx)),(!(epiY==locY)),(!(epiy==locy)))) {print(tag)}
}
