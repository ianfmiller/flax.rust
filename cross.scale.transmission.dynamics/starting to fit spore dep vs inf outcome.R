# load + prep data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")

# visualize
par(mfrow=c(1,1))
for (site in unique(epi$Site))
{
  sub.data.1<-corrected.epi[which(corrected.epi$Site==site),]
  sub.data.1<-sub.data.1[order(sub.data.1$Date.First.Observed.Diseased),]
  plant.loc.data<-corrected.locs[which(corrected.locs$Site==site),]
  plot(plant.loc.data$X+plant.loc.data$x,plant.loc.data$Y+plant.loc.data$y,xlim=c(0,10),ylim=c(0,20),main=site,xlab="X",ylab="Y",pch=16,cex=.9)
  new.infs<-0
  for (date in unique(sub.data.1$Date.First.Observed.Diseased))
  {
    #plot(plant.loc.data$X+plant.loc.data$x,plant.loc.data$Y+plant.loc.data$y,xlim=c(0,10),ylim=c(0,20),main=paste0(site,as.Date(date,origin = "1970-01-01")),xlab="X",ylab="Y")
    sub.data.2<-sub.data.1[which(sub.data.1$Date.First.Observed.Diseased==date),]
    points(sub.data.2$X+sub.data.2$x,sub.data.2$Y+sub.data.2$y,col="red")
    Sys.sleep(2)
    new.infs<-c(new.infs,new.infs[length(new.infs)]+dim(sub.data.2)[1])
  }
}

epi<-read.csv("~/Documents/GitHub/flax.rust/data/Epidemiology.csv")
epi<-epi[which(epi$Year==2020),]
plant.locs<-read.csv("~/Documents/GitHub/flax.rust/data/plant locs.csv")
plant.locs<-plant.locs[which(plant.locs$Year==2020),]
