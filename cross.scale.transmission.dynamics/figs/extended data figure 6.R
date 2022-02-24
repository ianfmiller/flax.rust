library(viridis)

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
diseased.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/diseased.focal.plants.RDS")
site<-"GM"
date0<-"2020-07-01"
date1<-"2020-07-08"

wind.data<-all.weath[which(all.weath$site==site),]
wind.data<-wind.data[which(wind.data$date>as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")),]
wind.data<-wind.data[which(wind.data$date<=as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")),]

xxvals<-seq(-5,5,.05)
yyvals<-seq(-5,5,.05)

I<-10
half.height<-5

preds<-c()

for(xx in xxvals)
{
  for(yy in yyvals)
  {
    preds<-c(preds,predict.kernel.tilted.plume(I=I,H=half.height/100,k=5.739690e-07,Ws=4.451030e-02,A=7.777373e-02,xtarget=xx,ytarget=yy,wind.data=wind.data))
  }
}
out.mat.1<-matrix(preds,length(xxvals),length(yyvals),byrow = T)

I<-1000
half.height<-5

preds<-c()
for(xx in xxvals)
{
  for(yy in yyvals)
  {
    preds<-c(preds,predict.kernel.tilted.plume(I=I,H=half.height/100,k=5.739690e-07,Ws=4.451030e-02,A=7.777373e-02,xtarget=xx,ytarget=yy,wind.data=wind.data))
  }
}
out.mat.2<-matrix(preds,length(xxvals),length(yyvals),byrow = T)

I<-10
half.height<-25

preds<-c()
for(xx in xxvals)
{
  for(yy in yyvals)
  {
    preds<-c(preds,predict.kernel.tilted.plume(I=I,H=half.height/100,k=5.739690e-07,Ws=4.451030e-02,A=7.777373e-02,xtarget=xx,ytarget=yy,wind.data=wind.data))
  }
}
out.mat.3<-matrix(preds,length(xxvals),length(yyvals),byrow = T)

I<-1000
half.height<-25

preds<-c()
for(xx in xxvals)
{
  for(yy in yyvals)
  {
    preds<-c(preds,predict.kernel.tilted.plume(I=I,H=half.height/100,k=5.739690e-07,Ws=4.451030e-02,A=7.777373e-02,xtarget=xx,ytarget=yy,wind.data=wind.data))
  }
}
out.mat.4<-matrix(preds,length(xxvals),length(yyvals),byrow = T)

layout(matrix(c(1,1,1,1,2,2,2,2,9,1,1,1,1,2,2,2,2,5,3,3,3,3,4,4,4,4,5,3,3,3,3,4,4,4,4,10,6,6,6,6,7,7,7,7,8,6,6,6,6,7,7,7,7,8),6,9,byrow = T))
par(mar=c(5,5,4,1))

image(xxvals,yyvals,log10(out.mat.1),zlim=c(-6.5,2.5),col=heat.colors(100,rev=T),xlab="X (m)",ylab="Y (m)",cex.lab=1.75,cex.axis=1.5,main="infection intensity = 10, plant height = 10 cm",cex.main=1.5)
points(0,0,pch=16,cex=1)
mtext("A",adj=1,cex=1.25,line=1)

image(xxvals,yyvals,log10(out.mat.2),zlim=c(-6.5,2.5),col=heat.colors(100,rev=T),xlab="X (m)",ylab="Y (m)",cex.lab=1.75,cex.axis=1.5,main="infection intensity = 1000, plant height = 10 cm",cex.main=1.5)
points(0,0,pch=16,cex=1)
mtext("B",adj=1,cex=1.25,line=1)

image(xxvals,yyvals,log10(out.mat.3),zlim=c(-6.5,2.5),col=heat.colors(100,rev=T),xlab="X (m)",ylab="Y ",cex.lab=1.75,cex.axis=1.5,main="infection intensity = 10, plant height = 50 cm",cex.main=1.5)
points(0,0,pch=16,cex=1)
mtext("C",adj=1,cex=1.25,line=1)

image(xxvals,yyvals,log10(out.mat.4),zlim=c(-6.5,2.5),col=heat.colors(100,rev=T),xlab="X (m)",ylab="Y (m)",cex.lab=1.75,cex.axis=1.5,main="infection intensity = 1000, plant height = 50 cm",cex.main=1.5)
points(0,0,pch=16,cex=1)
mtext("D",adj=1,cex=1.25,line=1)

par(mar=c(0,1,0,7))
adj<-diff(seq(-6.5,2.5,length.out = 100))[1]/2
plot(0,0,type="n",xlim=c(0,1),ylim=c(-6.5-adj,2.5+adj),axes=F,xlab="",ylab="")
for(i in 1:100)
{
  ii<-seq(-6.5,2.5,length.out=100)[i]
  rect(0,ii-adj,1,ii+adj,col=rev(heat.colors(100))[i],border = NA)
}
rect(0,-6.5-adj,1,2.5+adj)
mtext(expression(log[10]*' predicted spore deposition'),line=4,side=4,cex = 1.25)
axis(4,cex.axis=2)

par(mar=c(5,5,4,1))
image(xxvals,yyvals,log10(abs(out.mat.3-out.mat.1)),zlim=c(-8,3),col=rev(magma(100)),xlab="X (m)",ylab="Y (m)",cex.lab=1.75,cex.axis=1.5,main="infection intensity = 10",cex.main=1.5)
contour(xxvals,yyvals,out.mat.3-out.mat.1,add=T,col="blue",levels=c(0),lwd=5,drawlabels = F)
points(0,0,pch=16,cex=1)
mtext("E",adj=1,cex=1.25,line=1)

par(mar=c(5,5,4,1))
image(xxvals,yyvals,log10(abs(out.mat.4-out.mat.2)),zlim=c(-8,3),col=rev(magma(100)),xlab="X (m)",ylab="Y (m)",cex.lab=1.75,cex.axis=1.5,main="infection intensity = 1000",cex.main=1.5)
contour(xxvals,yyvals,out.mat.4-out.mat.2,add=T,col="blue",levels=c(0),lwd=5,drawlabels = F)
points(0,0,pch=16,cex=1)
mtext("F",adj=1,cex=1.25,line=1)

par(mar=c(0,1,0,7))
adj<-diff(seq(-8,3,length.out = 100))[1]/2
plot(0,0,type="n",xlim=c(0,1),ylim=c(-8-adj,3+adj),axes=F,xlab="",ylab="")
for(i in 1:100)
{
  ii<-seq(-8,3,length.out=100)[i]
  rect(0,ii-adj,1,ii+adj,col=rev(magma(100))[i],border = NA)
}
rect(0,-8-adj,1,3+adj)
axis(4,cex.axis=2)
mtext(expression(log[10]*' absolute difference in'),line=3.5,side=4,cex = 1.25)
mtext("predicted spore deposition",line=5,side=4,cex=1.25)


#export ata 1100 x 1400

