library(mgcv)
library(viridis)
# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/transmission data set building.R")
transmission.data<-transmission.data[which(transmission.data$status==0),]
transmission.data$site<-as.factor(transmission.data$site)
for(i in 1:nrow(transmission.data))
{
  if(is.na(transmission.data[i,"tag"])) {transmission.data[i,"tag"]<-paste0(transmission.data[i,"site"],"X",transmission.data[i,"X"]+transmission.data[i,"x"],"Y",transmission.data[i,"Y"]+transmission.data[i,"y"])}
}
transmission.data$tag<-as.factor(transmission.data$tag)

## visualize data

layout(matrix(c(1,1,2,3,4,5),3,2,byrow = T))
par(mar=c(6,6,6,6))
plot(jitter(log10(transmission.data$tot.spore.deposition)),jitter(transmission.data$status.next),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",main="all",axes=F,cex.lab=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(transmission.data[which(transmission.data$site=="CC"),"spore deposition"]+1e-10)),jitter(transmission.data[which(transmission.data$site=="CC"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[2],xlim=c(min(log10(transmission.data$tot.spore.deposition+1e-10)),max(log10(transmission.data$tot.spore.deposition+1e-10))),main="CC",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(transmission.data[which(transmission.data$site=="BT"),"spore deposition"]+1e-10)),jitter(transmission.data[which(transmission.data$site=="BT"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[3],xlim=c(min(log10(transmission.data$tot.spore.deposition+1e-10)),max(log10(transmission.data$tot.spore.deposition+1e-10))),main="BT",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(transmission.data[which(transmission.data$site=="GM"),"spore deposition"]+1e-10)),jitter(transmission.data[which(transmission.data$site=="GM"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[4],xlim=c(min(log10(transmission.data$tot.spore.deposition+1e-10)),max(log10(transmission.data$tot.spore.deposition+1e-10))),main="GM",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(transmission.data[which(transmission.data$site=="HM"),"spore deposition"]+1e-10)),jitter(transmission.data[which(transmission.data$site=="HM"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[5],xlim=c(min(log10(transmission.data$tot.spore.deposition+1e-10)),max(log10(transmission.data$tot.spore.deposition+1e-10))),main="HM",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()


if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS"))
{
  mod<-bam(status.next~
             s(log10(tot.spore.deposition),height.cm)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(max.abs.hum)+
             s(min.abs.hum)+
             s(mean.solar)+
             s(mean.daily.rain)+
             s(mean.wetness)+
             +offset(log(time))+
             s(tag,bs="re")+
             s(site,bs="re"),
           family=binomial(link="cloglog"),
           select = T,
           method="fREML",
           data=transmission.data,
           control = list(nthreads=4),
           discrete = T)
  
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")
}

## load model

transmission.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")

## model checking
#par(mfrow=c(2,2))
#gam.check(transmission.model) #Indicates that k should be higher all smooths aside from the tensor. Increasing k to significantly higher values is extremelly computationally expensive and leads to overfitting
#concurvity(transmission.model,full=F) #Concurvity is expected between all weather variables. The non-pessimistic estimations of concurvity (under $estimate and $observed) indicate that it is not a serious issue in the model fit.

layout(matrix(c(1,1,1,1,2,3,3,3,3,4,4,4,4,5,5,5,5,18,14,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,15,16,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,17),3,18,byrow = T))

par(mar=c(4,7,3,1.5))
plot(transmission.model,select = 1,scheme = 2,zlim=c(-8,8),too.far=.5,hcolors=rev(heat.colors(101)),contour.col="black",cex.lab=1.5,cex.axis=1.5,xlab=expression(log[10]*' predicted spore deposition'),ylab="plant height (cm)",main="",rug=F)
points(log10(transmission.data$tot.spore.deposition),transmission.data$height.cm,pch=".")
par(mar=c(4,1,3,4))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-8-0.08,8+0.08),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-8,8,length.out=101)[i]
  rect(0,ii-0.08,1,ii+0.08,col=rev(heat.colors(101))[i],border = NA)
}
rect(0,-8-0.08,1,8+0.08)
mtext("A",cex=1.25,font=2)
#mtext(expression('te('*log[10]*' predicted spore deposition, plant height)'),side=4,line=.5)
axis(2,cex.axis=1.5)

par(mar=c(4,4.5,3,4))
plot(transmission.model,select = 2,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean temperature (°C)",ylab="s(mean temperature)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 3,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="max. temperature (°C)",ylab="s(max. temperature)")
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 4,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="min. temperature (°C)",ylab="s(min. temperature)")
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 5,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean abs. humidity ('*g/m^3*')'),ylab="s(mean abs. humidity)")
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 6,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('max. abs. humidity ('*g/m^3*')'),ylab="s(max. abs. humidity)")
grid()
mtext("F",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 7,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('min abs. humidity ('*g/m^3*')'),ylab="s(min. abs. humidity)")
grid()
mtext("G",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean solar radiation ('*W/m^2*')'),ylab="s(mean solalr radiation")
grid()
mtext("H",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 9,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="total rainfall (mm)",ylab="s(total rainfall)")
grid()
mtext("I",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 10,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean leaf wetness (%)",ylab="s(mean leaf wetness)")
grid()
mtext("J",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 11,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(tag)")
grid()
mtext("K",adj=1,cex=1.25,font=2)
plot(transmission.model,select = 12,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(site)")
grid()
mtext("L",adj=1,cex=1.25,font=2)



par(mar=c(6,8,6,0))
layout(matrix(c(1,1,1,1,2),1,5))
plot(log10(transmission.data$tot.spore.deposition),jitter(transmission.data$status.next*1,factor=.5),xlab=expression('predicted '*log[10]*' spore depositin'),ylab="",axes=F,cex.lab=2,ylim=c(-0.1,1.1),xlim=c(-6,3))
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2,line=3,tick = F)
axis(2,at=c(0,.5,1),cex.axis=2,col="red",col.axis="red")
mtext("odds of infection",side=2,line=2.5,cex=2,col="red")
axis(1,cex.axis=1.5)

dummy.tot.spore.deposition.levels<-10^seq(-6,2.3,.1)

new.data.1<-transmission.data[1:length(dummy.tot.spore.deposition.levels),]
new.data.1$tot.spore.deposition<-dummy.tot.spore.deposition.levels
new.data.1$height.cm<-5

new.data.2<-transmission.data[1:length(dummy.tot.spore.deposition.levels),]
new.data.2$tot.spore.deposition<-dummy.tot.spore.deposition.levels
new.data.2$height.cm<-10

new.data.3<-transmission.data[1:length(dummy.tot.spore.deposition.levels),]
new.data.3$tot.spore.deposition<-dummy.tot.spore.deposition.levels
new.data.3$height.cm<-25

new.data.4<-transmission.data[1:length(dummy.tot.spore.deposition.levels),]
new.data.4$tot.spore.deposition<-dummy.tot.spore.deposition.levels
new.data.4$height.cm<-50

points(log10(new.data.1$tot.spore.deposition),predict(transmission.model,newdata = new.data.1,type = "response",exclude = c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(max.abs.hum)","s(min.abs.hum)","s(mean.solar)","s(mean.daily.rain)","s(mean.wetness)")),type="l",col="red",lwd=4,lty=1)
points(log10(new.data.2$tot.spore.deposition),predict(transmission.model,newdata = new.data.2,type = "response",exclude = c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(max.abs.hum)","s(min.abs.hum)","s(mean.solar)","s(mean.daily.rain)","s(mean.wetness)")),type="l",col="red",lwd=4,lty=2)
points(log10(new.data.3$tot.spore.deposition),predict(transmission.model,newdata = new.data.3,type = "response",exclude = c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(max.abs.hum)","s(min.abs.hum)","s(mean.solar)","s(mean.daily.rain)","s(mean.wetness)")),type="l",col="red",lwd=4,lty=3)
points(log10(new.data.4$tot.spore.deposition),predict(transmission.model,newdata = new.data.4,type = "response",exclude = c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(max.abs.hum)","s(min.abs.hum)","s(mean.solar)","s(mean.daily.rain)","s(mean.wetness)")),type="l",col="red",lwd=4,lty=4)
par(mar=c(0,0,0,0))
plot(0,0,type="n",axes=F,xlab = "",ylab="")
legend("left",legend=c("50cm","25cm","10cm","5cm"),col="red",lty=c(4,3,2,1),bty="n",cex=2,lwd=4)
