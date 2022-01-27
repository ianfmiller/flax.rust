library(mgcv)
library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/transmission data set building.R")
transmission.data<-transmission.data[which(transmission.data$status==0),]
transmission.data$site<-as.factor(transmission.data$site)
for(i in 1:nrow(transmission.data))
{
  if(is.na(transmission.data[i,"tag"])) {transmission.data[i,"tag"]<-paste0(transmission.data[i,"site"],"X",transmission.data[i,"X"]+transmission.data[i,"x"],"Y",transmission.data[i,"Y"]+transmission.data[i,"y"])}
}
transmission.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")

site.cols<-viridis_pal(alpha=.5)(20)[c(20,15,6,1)]
weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2,1,1)])

layout(matrix(c(rep(16,10),1,1,1,1,2,3,3,3,3,4,4,4,4,rep(17,10),rep(10,10),1,1,1,1,2,3,3,3,3,4,4,4,4,rep(11,10),rep(10,10),12,5,5,5,5,6,6,6,6,7,7,7,7,rep(11,10),rep(10,10),12,5,5,5,5,6,6,6,6,7,7,7,7,rep(11,10),rep(10,10),13,14,14,8,8,8,8,9,9,9,9,15,15,rep(11,10),rep(18,10),13,14,14,8,8,8,8,9,9,9,9,15,15,rep(19,10)),6,33,byrow=T))

par(mar=c(4,4,3,1))
plot(transmission.model,select = 1,scheme = 2,zlim=c(-10,10),too.far=1,hcolors=rev(heat.colors(101)),contour.col="black",cex.lab=1.5,cex.axis=1.5,xlab="",ylab="",main="",rug=F)
mtext(expression(log[10]*' predicted spore deposition'),1,line = 2.25,cex=1)
mtext("plant height (cm)",2,line=2.25,cex=1)
points(log10(transmission.data$tot.spore.deposition),transmission.data$height.cm,pch=16,cex=.1)

par(mar=c(4,1,2.5,1.5))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-10-.1,10+.1),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-10,10,length.out=101)[i]
  rect(0,ii-.1,1,ii+.1,col=rev(heat.colors(101))[i],border = NA)
}
rect(0,-10-.1,1,10+.1)
mtext("B",cex=1.25,font=2)
mtext('te(log10 predicted spore\ndeposition, plant height)',side=2,line=0,cex = .6)
axis(4,cex.axis=1,tck=-.5,padj=-1)

par(mar=c(4,4,3,1))
plot(transmission.model,select = 2,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("mean temperature (°C)",1,line = 2.25,cex=1)
mtext("s(mean temperature)",2,line=2.25,cex=1)
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 3,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("max. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(max. temperature)",2,line=2.25,cex=1)
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 4,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("min. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(min. temperature)",2,line=2.25,cex=1)
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 5,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext(expression('mean abs. humidity ('*g/m^3*')'),1,line = 2.25,cex=1)
mtext("s(mean abs. humidity)",2,line=2.25,cex=1)
grid()
mtext("F",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 6,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("total rainfall (mm)",1,line = 2.25,cex=1)
mtext("s(total rainfall)",2,line=2.25,cex=1)
grid()
mtext("G",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 7,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(tag)",2,line=2.25,cex=1)
grid()
mtext("H",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 8,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",col=site.cols[c(2,1,3,4)],cex=2,pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(site)",2,line=2.25,cex=1)
grid()
mtext("I",adj=1,cex=1.25,font=2)
legend("topleft",legend=c("CC","BT","GM","HM"),pch=16,col=site.cols,cex=1,bty="n",pt.cex = 2)

transmission.data.shuffled= transmission.data[sample(1:nrow(transmission.data)), ]
site.indicies<-c(2,1,3,4)[as.numeric(transmission.data.shuffled$site)]
par(mar=c(5,6,5,2))
plot(jitter(log10(transmission.data.shuffled$tot.spore.deposition)),jitter(transmission.data.shuffled$status.next),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,cex.lab=2,col=site.cols[site.indicies],pch=16,cex=2,panel.first = grid())
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
mtext("A",side=3,adj=1,cex=2)
legend("topleft",legend=c("CC","BT","GM","HM"),col=c(site.cols),pt.cex=c(2,2,2,2,2),pch=16,cex=2,bty="n")

dummy.tot.spore.deposition.levels<-10^seq(-6,2,.1)

date0<-as.POSIXct("2020-07-07 00:00:00",tz="UTC")
date1<-as.POSIXct("2020-07-14 00:00:00",tz="UTC")
site<-"GM"
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")

weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #### pull out temp data for site

#### subst temp rh data to relevant window
temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs

#### subset weather data to relevant window
weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data

#### calculate environmental variable metrics
new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
abs.hum<-13.24732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+temp.rh.sub$temp.c)
new.mean.abs.hum<-mean(abs.hum,na.rm=T)
new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)*(12*24)


new.data<-transmission.data[1:length(dummy.tot.spore.deposition.levels),]
new.data$tot.spore.deposition<-dummy.tot.spore.deposition.levels
new.data$mean.temp<-new.mean.temp
new.data$max.temp<-new.max.temp
new.data$min.temp<-new.min.temp
new.data$mean.abs.hum<-new.mean.abs.hum
new.data$mean.daily.rain<-new.mean.temp
new.data$site<-site

new.data.1<-new.data
new.data.1$height.cm<-5
new.data.2<-new.data
new.data.2$height.cm<-10
new.data.3<-new.data
new.data.3$height.cm<-25
new.data.4<-new.data
new.data.4$height.cm<-50


exclude.vec<-c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(mean.daily.rain)")

plot(0,0,type="n",xlab=expression('predicted '*log[10]*' spore deposition'),ylab="odds of infection",cex.lab=2,cex.axis=2,ylim=c(0,.85),xlim=c(-6,2))
points(log10(new.data.1$tot.spore.deposition),predict(transmission.model,newdata = new.data.1,type = "response",exclude = exclude.vec),type="l",col="khaki1",lwd=4,lty=1)
points(log10(new.data.2$tot.spore.deposition),predict(transmission.model,newdata = new.data.2,type = "response",exclude = exclude.vec),type="l",col="khaki2",lwd=4,lty=1)
points(log10(new.data.3$tot.spore.deposition),predict(transmission.model,newdata = new.data.3,type = "response",exclude = exclude.vec),type="l",col="khaki3",lwd=4,lty=1)
points(log10(new.data.4$tot.spore.deposition),predict(transmission.model,newdata = new.data.4,type = "response",exclude = exclude.vec),type="l",col="khaki4",lwd=4,lty=1)
legend("topleft",legend=c("50cm","25cm","10cm","5cm"),col=c("khaki4","khaki3","khaki2","khaki1"),lty=c(1),bty="n",cex=2,lwd=4,seg.len=5)
mtext("J",side=3,adj=1,cex=2)

