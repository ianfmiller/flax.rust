library(mgcv)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")

delta.pustules<-subset(delta.pustules,time<=8)

# visualize data

layout(matrix(c(1,2,3,3),2,2,byrow = T))
par(mar=c(5,6,2,5))

## histograms
hist(delta.pustules$area,main="",breaks=100,xlab="",cex.lab=2,cex.axis=2,cex.main=2)
mtext(expression('pustule area ('*mm^2*')'),side=1,line = 3.5,cex=2*.83)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.pustules$area.next-delta.pustules$area)/delta.pustules$time,main="",breaks=100,xlab="",cex.lab=2,cex.axis=2,cex.main=2)
mtext(expression('change in pustule area per day ('*mm^2*')'),side=1,line = 3.5,cex=2*.83)

mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories
par(mar=c(2,6,0,0.5))
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="",ylab=expression('pustule area ('*mm^2*')'),cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3,cex=2)
i<-0

plot.cols<-sample(rainbow(2007))

for (tag in unique(pustules$tag))
{
  sub.pustules1<-pustules[which(pustules$tag==tag),]
  
  for (color in unique(sub.pustules1$color))
  {
    sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
    {
      sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
      
      for(pustule.number in unique(sub.pustules3$pustule.number))
      {
        i<-i+1
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        points(sub.pustules4$date,sub.pustules4$area,col=plot.cols[i],type="l",lwd=.5)
      }
    }
  }
}