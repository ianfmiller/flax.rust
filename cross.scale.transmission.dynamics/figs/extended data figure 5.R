library(mgcv)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.plant.heights<-subset(delta.plant.heights,time<=8)

# visualize data

layout(matrix(c(1,2,3,3),2,2,byrow = T))
par(mar=c(5,6,2,5))

## histograms
hist(delta.plant.heights$height,main="",breaks=100,xlab="plant height (cm)",cex.lab=2,cex.axis=2,cex.main=2,panel.first=grid())
box()
mtext("A",side=3,adj=.95,line=-3,cex=2)
hist((delta.plant.heights$height.next-delta.plant.heights$height)/delta.plant.heights$time,main="",breaks=100,xlab="change in plant height per day (cm)",cex.lab=2,cex.axis=2,cex.main=2,panel.first=grid())
mtext("B",side=3,adj=.95,line=-3,cex=2)
box()

## plot trajectories
par(mar=c(3,6,0,5))
plot(c(min(plant.heights$date),max(plant.heights$date)),c(0,max(plant.heights$max.height)),type="n",xlab="",ylab="plant height (cm)",cex.lab=2,cex.axis=2,panel.first=grid())
mtext("C",side=3,adj=.975,line=-3,cex=2)
i<-0

set.seed(8427602)
plot.cols<-sample(rainbow(180))

for (i in 1:length(unique(plant.heights$tag)))
{
  tag<-unique(plant.heights$tag)[i]
  sub.heights<-plant.heights[which(plant.heights$tag==tag),]
  sub.heights<-sub.heights[order(sub.heights$date),]
  points(sub.heights$date,sub.heights$max.height,col=plot.cols[i],type="l",lwd=.5)
  
}

#export at 1156 x 738