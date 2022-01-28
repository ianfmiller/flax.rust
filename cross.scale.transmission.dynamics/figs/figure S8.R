source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
diseased.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/diseased.focal.plants.RDS")
site<-"GM"
date0<-"2020-07-01"
date1<-"2020-07-08"

wind.data<-all.weath[which(all.weath$site==site),]
wind.data<-wind.data[which(wind.data$date>as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")),]
wind.data<-wind.data[which(wind.data$date<=as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")),]

layout(matrix(c(1,1,1,1,2,2,2,2,7,1,1,1,1,2,2,2,2,5,3,3,3,3,4,4,4,4,5,3,3,3,3,4,4,4,4,6),4,9,byrow = T))
par(mar=c(4,5,4,4))

xxvals<-seq(-1,1,.05)
yyvals<-seq(-1,1,.05)

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
out.mat<-matrix(preds,length(xxvals),length(yyvals),byrow = T)
image(xxvals,yyvals,log10(out.mat),zlim=c(-6,3),col=heat.colors(100,rev=T),xlab="X",ylab="Y",cex.lab=1.5,cex.axis=1.5,main="infection intensity = 10, plant height = 10 cm",cex.main=1.5)
contour(xxvals,yyvals,log10(out.mat),add=T)
points(0,0,pch=16,cex=3)
mtext("A",adj=1,cex=1.25,font=2,line=1)

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
out.mat<-matrix(preds,length(xxvals),length(yyvals),byrow = T)
image(xxvals,yyvals,log10(out.mat),zlim=c(-6,3),col=heat.colors(100,rev=T),xlab="X",ylab="Y",cex.lab=1.5,cex.axis=1.5,main="infection intensity = 1000, plant height = 10 cm",cex.main=1.5)
contour(xxvals,yyvals,log10(out.mat),add=T)
points(0,0,pch=16,cex=3)
mtext("B",adj=1,cex=1.25,font=2,line=1)

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
out.mat<-matrix(preds,length(xxvals),length(yyvals),byrow = T)
image(xxvals,yyvals,log10(out.mat),zlim=c(-6,3),col=heat.colors(100,rev=T),xlab="X",ylab="Y",cex.lab=1.5,cex.axis=1.5,main="infection intensity = 10, plant height = 50 cm",cex.main=1.5)
contour(xxvals,yyvals,log10(out.mat),add=T)
points(0,0,pch=16,cex=3)
mtext("C",adj=1,cex=1.25,font=2,line=1)

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
out.mat<-matrix(preds,length(xxvals),length(yyvals),byrow = T)
image(xxvals,yyvals,log10(out.mat),zlim=c(-6,3),col=heat.colors(100,rev=T),xlab="X",ylab="Y",cex.lab=1.5,cex.axis=1.5,main="infection intensity = 1000, plant height = 50 cm",cex.main=1.5)
contour(xxvals,yyvals,log10(out.mat),add=T) 
points(0,0,pch=16,cex=3)
mtext("D",adj=1,cex=1.25,font=2,line=1)

par(mar=c(0,1,0,7))
adj<-diff(seq(-6,3,length.out = 100))[1]/2
plot(0,0,type="n",xlim=c(0,1),ylim=c(-6-adj,3+adj),axes=F,xlab="",ylab="")
for(i in 1:100)
{
  ii<-seq(-6,3,length.out=100)[i]
  rect(0,ii-adj,1,ii+adj,col=rev(heat.colors(100))[i],border = NA)
}
rect(0,-6-adj,1,3+adj)
mtext(expression(log[10]*' predicted spore deposition'),line=4,side=4,cex = 2)
axis(4,cex.axis=2)

par(mar=c(0,0,2,0))
plot(0,0,type="n",axes=F,xlab="",ylab="")
legend("topleft",legend="source\nplant",pch=16,cex=2,pt.cex=3,bty="n")
