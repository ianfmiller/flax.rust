library(mgcv)
set.seed(89757038)

# load and prep data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/infection intensity data prep.R")
delta.infection.intensity<-subset(delta.infection.intensity,time<=8)

# visualize data
layout(matrix(c(1,2,3,3),2,2,byrow = T))
par(mar=c(5,6,2,5))

## histograms
hist(delta.infection.intensity$infection.intensity,main="",breaks=100,xlab="infection intensity",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.infection.intensity$infection.intensity.next-delta.infection.intensity$infection.intensity)/delta.infection.intensity$time,main="",breaks=100,xlab="change in infection intensity per day",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories
par(mar=c(2,6,0,0.5))
plot(c(min(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y")),max(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y"))),c(0,max(infection.intensity$infection.intensity)),type="n",xlab="",ylab="infection intensity",cex.lab=2,cex.axis=2,ylim=c(0,750))
mtext("C",side=3,adj=1,line=-3,cex=2)

i<-0
plot.cols<-sample(rainbow(length(unique(infection.intensity$Tag))))

for (tag in unique(infection.intensity$Tag))
{
  sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
  points(sub.infection.intensity.1$Date,sub.infection.intensity.1$infection.intensity,col=plot.cols[i],type="l",lwd=2)
  i<-i+1
}

par(fig=c(.04,.34,.24,.49),new=T)
plot(c(min(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y")),max(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y"))),c(0,max(infection.intensity$infection.intensity)),type="n",xlab="",ylab="",cex.lab=1,cex.axis=1)
rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4],col="white")
i<-0
plot.cols<-sample(rainbow(length(unique(infection.intensity$Tag))))

for (tag in unique(infection.intensity$Tag))
{
  sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
  points(sub.infection.intensity.1$Date,sub.infection.intensity.1$infection.intensity,col=plot.cols[i],type="l",lwd=2)
  i<-i+1
}
    