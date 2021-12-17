source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
diseased.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/diseased.focal.plants.RDS")
site<-"BT"
tag<-86
date0<-"2020-06-24"
date1<-"2020-07-01"

wind.data<-all.weath[which(all.weath$site==site),]
wind.data<-wind.data[which(wind.data$date>as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")),]
wind.data<-wind.data[which(wind.data$date<=as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")),]

diseased.focal.plants.index<-intersect(which(diseased.focal.plants$Date==date0),which(diseased.focal.plants$Tag==tag))
I<-diseased.focal.plants[diseased.focal.plants.index,"inf.intens"]
half.height<-.5*diseased.focal.plants[diseased.focal.plants.index,"max.height"]/100

preds<-c()
xxvals<-seq(-2,2,.2)
yyvals<-seq(-2,2,.2)
for(xx in xxvals)
{
  for(yy in yyvals)
  {
    preds<-c(preds,predict.kernel.tilted.plume(I=I,H=half.height,k=1.506668e-06,Ws=1.016603e+00,A=1.488799e-01,xtarget=xx,ytarget=yy,wind.data=wind.data))
  }
}
filled.contour(x=seq(-2,2,.01),y=seq(-2,2,.01),z=res.mat,xlim=c(-2,2),ylim=c(-2,2),key.title = mtext(expression(log[10]*'spores per mm'^2),cex=1.5),cex.lab=1.5,xlab="X (meters)",ylab="Y (meters)",)
out.mat<-matrix(preds,21,21,byrow = T)
filled.contour(log10(out.mat),x=xxvals,y=yyvals,labcex=1,zlim=c(-10,3),xlab="X (meters)",ylab="Y (meters)")
