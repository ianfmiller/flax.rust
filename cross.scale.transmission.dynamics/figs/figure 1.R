library(imager)
library(lubridate)
library(viridis)
site.cols<-viridis_pal(alpha=1)(20)[c(20,15,6,1)]

layout(matrix(c(1,1,1,1,2,2,2,2,3,3,4,4,1,1,1,1,2,2,2,2,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18),4,12,byrow = T))
par(mar=c(.5,.5,.5,.5))

panel.A<-load.image("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/A.jpg")
plot(panel.A,axes=F)
text((par()$usr[2]-par()$usr[1])*.95+par()$usr[1],
     (par()$usr[4]-par()$usr[3])*.975+par()$usr[3],
     labels="A",
     cex=2/par()$cex)

panel.B<-load.image("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/B.jpg")
plot(panel.B,axes=F)
text((par()$usr[2]-par()$usr[1])*.95+par()$usr[1],
     (par()$usr[4]-par()$usr[3])*.975+par()$usr[3],
     labels="B",
     cex=2/par()$cex)



par(mar=c(.5,.5,.5,3.25),xpd=T)

panel.C<-load.image("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/C.jpg")
plot(panel.C,axes=F)
rect(0,dim(panel.C)[2],dim(panel.C)[1],0,border=site.cols[4],density = 0,lwd=5)
text(250,200,"HM",cex=3,col="white")
text((par()$usr[2]-par()$usr[1])*1.1+par()$usr[1],
     (par()$usr[4]-par()$usr[3])*.975+par()$usr[3],
     labels="C",
     cex=2/par()$cex)
panel.D<-load.image("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/D.jpg")
plot(panel.D,axes=F)
rect(0,dim(panel.D)[2],dim(panel.D)[1],0,border=site.cols[3],density = 0,lwd=5)
text(250,200,"GM",cex=3,col="white")
text((par()$usr[2]-par()$usr[1])*1.1+par()$usr[1],
     (par()$usr[4]-par()$usr[3])*.975+par()$usr[3],
     labels="D",
     cex=2/par()$cex)
panel.E<-load.image("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/E.jpg")
plot(panel.E,axes=F)
rect(0,dim(panel.E)[2],dim(panel.E)[1],0,border=site.cols[2],density = 0,lwd=5)
text(250,200,"BT",cex=3,col="white")
text((par()$usr[2]-par()$usr[1])*1.1+par()$usr[1],
     (par()$usr[4]-par()$usr[3])*.975+par()$usr[3],
     labels="E",
     cex=2/par()$cex)
panel.F<-load.image("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/F.jpg")
plot(panel.F,axes=F)
rect(0,dim(panel.F)[2],dim(panel.F)[1],0,border=site.cols[1],density = 0,lwd=5)
text(250,200,"CC",cex=3,col="white")
text((par()$usr[2]-par()$usr[1])*1.1+par()$usr[1],
     (par()$usr[4]-par()$usr[3])*.975+par()$usr[3],
     labels="F",
     cex=2/par()$cex)

par(mar=c(4,4,3,3),xpd=F)
site.cols<-viridis_pal(alpha=.75)(20)[c(20,15,6,1)]
source("~/Documents/Github/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")

cc.temp.rh<-all.temp.rh[which(all.temp.rh$site=="CC"),]
cc.temp.rh.dates<-seq.POSIXt(as.POSIXct("2020-06-16 00:00:00",tz="UTC"),as.POSIXct("2020-07-26 00:00:00",tz="UTC"),by="1 day")

bt.temp.rh<-all.temp.rh[which(all.temp.rh$site=="BT"),]
bt.temp.rh.dates<-seq.POSIXt(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day")

gm.temp.rh<-all.temp.rh[which(all.temp.rh$site=="GM"),]
gm.temp.rh.dates<-seq.POSIXt(as.POSIXct("2020-06-17 00:00:00",tz="UTC"),as.POSIXct("2020-07-27 00:00:00",tz="UTC"),by="1 day")

hm.temp.rh<-all.temp.rh[which(all.temp.rh$site=="HM"),]
hm.temp.rh.dates<-seq.POSIXt(as.POSIXct("2020-06-19 00:00:00",tz="UTC"),as.POSIXct("2020-07-09 00:00:00",tz="UTC"),by="1 day")

pull.daily.means.temp.rh<-function(temp.rh.data,dates)
{
  out.data<-data.frame(date=POSIXct(),mean.temp=numeric(),max.temp=numeric(),min.temp=numeric(),mean.abs.hum=numeric())
  for(i in 1:(length(dates)-1))
  {
    date0<-dates[i]
    date1<-dates[i+1]
    sub.temp.rh<-temp.rh.data[intersect(which(temp.rh.data$date.time>=date0),which(temp.rh.data$date.time<date1)),]
    mean.temp<-mean(sub.temp.rh$temp.c,na.rm=T)
    max.temp<-max(sub.temp.rh$temp.c,na.rm=T)
    min.temp<-min(sub.temp.rh$temp.c,na.rm=T)
    abs.hum<-13.24732*exp((17.67*sub.temp.rh$temp.c)/(sub.temp.rh$temp.c+243.5))*sub.temp.rh$rh/(273.15+sub.temp.rh$temp.c)
    mean.abs.hum<-mean(abs.hum,na.rm=T)
    out.data<-rbind(out.data,data.frame(date=date0,mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,mean.abs.hum=mean.abs.hum))
  }
  out.data
}

cc.temp.rh.daily.means<-pull.daily.means.temp.rh(cc.temp.rh,cc.temp.rh.dates)
bt.temp.rh.daily.means<-pull.daily.means.temp.rh(bt.temp.rh,bt.temp.rh.dates)
gm.temp.rh.daily.means<-pull.daily.means.temp.rh(gm.temp.rh,gm.temp.rh.dates)
hm.temp.rh.daily.means<-pull.daily.means.temp.rh(hm.temp.rh,hm.temp.rh.dates)

cc.weath<-all.weath[which(all.weath$site=="CC"),]
cc.weath.dates<-seq.POSIXt(as.POSIXct("2020-06-23 00:00:00",tz="UTC"),as.POSIXct("2020-07-26 00:00:00",tz="UTC"),by="1 day")

bt.weath<-all.weath[which(all.weath$site=="BT"),]
bt.weath.dates<-seq.POSIXt(as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day")

gm.weath<-all.weath[which(all.weath$site=="GM"),]
gm.weath.dates<-seq.POSIXt(as.POSIXct("2020-06-24 00:00:00",tz="UTC"),as.POSIXct("2020-07-27 00:00:00",tz="UTC"),by="1 day")

hm.weath<-all.weath[which(all.weath$site=="HM"),]
hm.weath.dates<-seq.POSIXt(as.POSIXct("2020-06-26 00:00:00",tz="UTC"),as.POSIXct("2020-07-27 00:00:00",tz="UTC"),by="1 day")

pull.daily.means.weath<-function(weath.data,dates)
{
  out.data<-data.frame(date=POSIXct(),mean.daily.rain=numeric())
  for(i in 1:(length(dates)-1))
  {
    date0<-dates[i]
    date1<-dates[i+1]
    sub.weath<-weath.data[intersect(which(weath.data$date>=date0),which(weath.data$date<date1)),]
    new.mean.daily.rain<-mean(sub.weath$rain,na.rm=T)*(12*24)
    out.data<-rbind(out.data,data.frame(date=date0,mean.daily.rain=new.mean.daily.rain))
  }
  out.data
}

cc.weath.daily.means<-pull.daily.means.weath(cc.weath,cc.weath.dates)
bt.weath.daily.means<-pull.daily.means.weath(bt.weath,bt.weath.dates)
gm.weath.daily.means<-pull.daily.means.weath(gm.weath,gm.weath.dates)
hm.weath.daily.means<-pull.daily.means.weath(hm.weath,hm.weath.dates)


plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(5,19),xlab="",ylab="",axes=F)
mtext("mean temperature (°C)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(cc.temp.rh.daily.means$date,cc.temp.rh.daily.means$mean.temp,type="l",col=site.cols[1],lwd=2)
points(bt.temp.rh.daily.means$date,bt.temp.rh.daily.means$mean.temp,type="l",col=site.cols[2],lwd=2)
points(gm.temp.rh.daily.means$date,gm.temp.rh.daily.means$mean.temp,type="l",col=site.cols[3],lwd=2)
points(hm.temp.rh.daily.means$date,hm.temp.rh.daily.means$mean.temp,type="l",col=site.cols[4],lwd=2)
mtext("G",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(11,36),xlab="",ylab="",axes=F)
mtext("max. temperature (°C)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(cc.temp.rh.daily.means$date,cc.temp.rh.daily.means$max.temp,type="l",col=site.cols[1],lwd=2)
points(bt.temp.rh.daily.means$date,bt.temp.rh.daily.means$max.temp,type="l",col=site.cols[2],lwd=2)
points(gm.temp.rh.daily.means$date,gm.temp.rh.daily.means$max.temp,type="l",col=site.cols[3],lwd=2)
points(hm.temp.rh.daily.means$date,hm.temp.rh.daily.means$max.temp,type="l",col=site.cols[4],lwd=2)
mtext("H",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(-6,12),xlab="",ylab="",axes=F)
mtext("min. temperature (°C)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(cc.temp.rh.daily.means$date,cc.temp.rh.daily.means$min.temp,type="l",col=site.cols[1],lwd=2)
points(bt.temp.rh.daily.means$date,bt.temp.rh.daily.means$min.temp,type="l",col=site.cols[2],lwd=2)
points(gm.temp.rh.daily.means$date,gm.temp.rh.daily.means$min.temp,type="l",col=site.cols[3],lwd=2)
points(hm.temp.rh.daily.means$date,hm.temp.rh.daily.means$min.temp,type="l",col=site.cols[4],lwd=2)
mtext("I",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(3,10),xlab="",ylab="",axes=F)
mtext(expression('mean abs. humidity ('*g/m^3*')'),2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(cc.temp.rh.daily.means$date,cc.temp.rh.daily.means$mean.abs.hum,type="l",col=site.cols[1],lwd=2)
points(bt.temp.rh.daily.means$date,bt.temp.rh.daily.means$mean.abs.hum,type="l",col=site.cols[2],lwd=2)
points(gm.temp.rh.daily.means$date,gm.temp.rh.daily.means$mean.abs.hum,type="l",col=site.cols[3],lwd=2)
points(hm.temp.rh.daily.means$date,hm.temp.rh.daily.means$mean.abs.hum,type="l",col=site.cols[4],lwd=2)
mtext("J",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(0,27),xlab="",ylab="",axes=F)
mtext("mean daily rainfall (mm)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(cc.weath.daily.means$date,cc.weath.daily.means$mean.daily.rain,type="l",col=site.cols[1],lwd=2)
points(bt.weath.daily.means$date,bt.weath.daily.means$mean.daily.rain,type="l",col=site.cols[2],lwd=2)
points(gm.weath.daily.means$date,gm.weath.daily.means$mean.daily.rain,type="l",col=site.cols[3],lwd=2)
points(hm.weath.daily.means$date,hm.weath.daily.means$mean.daily.rain,type="l",col=site.cols[4],lwd=2)
mtext("K",side=3,adj=1,cex=2,line=.5)

par(mar=c(4,0,3,3),xpd=F)
plot(0,0,type="n",axes=F,xlab="",ylab="",bty="n")
legend("center",legend = c("HM","GM","BT","CC"),lwd=4,col=rev(site.cols),bty = "n",cex=1.75,seg.len = 4)

mean.temp.data.2020<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/Tair.csv"),header = F)
max.temp.data.2020<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/tasmax.csv"),header = F)
min.temp.data.2020<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/tasmin.csv"),header = F)
rh.data.2020<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/relHumid.csv"),header = F)
rainfall.data.2020<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/rainfall.csv"),header = F)
mean.temp.data.2045<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/Tair.csv"),header = F)
max.temp.data.2045<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/tasmax.csv"),header = F)
min.temp.data.2045<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/tasmin.csv"),header = F)
rh.data.2045<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/relHumid.csv"),header = F)
rainfall.data.2045<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/rainfall.csv"),header = F)
mean.temp.data.2070<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/Tair.csv"),header = F)
max.temp.data.2070<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/tasmax.csv"),header = F)
min.temp.data.2070<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/tasmin.csv"),header = F)
rh.data.2070<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/relHumid.csv"),header = F)
rainfall.data.2070<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/rainfall.csv"),header = F)

colnames(mean.temp.data.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.temp.data.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.temp.data.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(max.temp.data.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(max.temp.data.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(max.temp.data.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(min.temp.data.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(min.temp.data.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(min.temp.data.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(rh.data.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(rh.data.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(rh.data.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(rainfall.data.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(rainfall.data.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(rainfall.data.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")


weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2)])

par(mar=c(4,4,3,3),xpd=F)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(5,19),xlab="",ylab="",axes=F)
mtext("mean temperature (°C)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.temp.data.2020$cesm1.cam5.1.rcp45,type="l",col=weather.colors[2],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.temp.data.2020$cesm1.cam5.1.rcp85,type="l",col=weather.colors[3],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.temp.data.2045$cesm1.cam5.1.rcp45,type="l",col=weather.colors[4],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.temp.data.2045$cesm1.cam5.1.rcp85,type="l",col=weather.colors[5],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.temp.data.2070$cesm1.cam5.1.rcp45,type="l",col=weather.colors[6],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.temp.data.2070$cesm1.cam5.1.rcp85,type="l",col=weather.colors[7],lwd=2,lty=1)
mtext("L",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(11,36),xlab="",ylab="",axes=F)
mtext("max. temperature (°C)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),max.temp.data.2020$cesm1.cam5.1.rcp45,type="l",col=weather.colors[2],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),max.temp.data.2020$cesm1.cam5.1.rcp85,type="l",col=weather.colors[3],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),max.temp.data.2045$cesm1.cam5.1.rcp45,type="l",col=weather.colors[4],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),max.temp.data.2045$cesm1.cam5.1.rcp85,type="l",col=weather.colors[5],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),max.temp.data.2070$cesm1.cam5.1.rcp45,type="l",col=weather.colors[6],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),max.temp.data.2070$cesm1.cam5.1.rcp85,type="l",col=weather.colors[7],lwd=2,lty=1)
mtext("M",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(-6,12),xlab="",ylab="",axes=F)
mtext("min. temperature (°C)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),min.temp.data.2020$cesm1.cam5.1.rcp45,type="l",col=weather.colors[2],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),min.temp.data.2020$cesm1.cam5.1.rcp85,type="l",col=weather.colors[3],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),min.temp.data.2045$cesm1.cam5.1.rcp45,type="l",col=weather.colors[4],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),min.temp.data.2045$cesm1.cam5.1.rcp85,type="l",col=weather.colors[5],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),min.temp.data.2070$cesm1.cam5.1.rcp45,type="l",col=weather.colors[6],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),min.temp.data.2070$cesm1.cam5.1.rcp85,type="l",col=weather.colors[7],lwd=2,lty=1)
mtext("N",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(2,10),xlab="",ylab="",axes=F)
mtext(expression('mean abs. humidity ('*g/m^3*')'),2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
mean.abs.hums.2020.rcp45<-13.24732*exp((17.67*mean.temp.data.2020$cesm1.cam5.1.rcp45)/(mean.temp.data.2020$cesm1.cam5.1.rcp45+243.5))*rh.data.2020$cesm1.cam5.1.rcp45/(273.15+mean.temp.data.2020$cesm1.cam5.1.rcp45)
mean.abs.hums.2020.rcp85<-13.24732*exp((17.67*mean.temp.data.2020$cesm1.cam5.1.rcp85)/(mean.temp.data.2020$cesm1.cam5.1.rcp85+243.5))*rh.data.2020$cesm1.cam5.1.rcp85/(273.15+mean.temp.data.2020$cesm1.cam5.1.rcp85)
mean.abs.hums.2045.rcp45<-13.24732*exp((17.67*mean.temp.data.2045$cesm1.cam5.1.rcp45)/(mean.temp.data.2045$cesm1.cam5.1.rcp45+243.5))*rh.data.2045$cesm1.cam5.1.rcp45/(273.15+mean.temp.data.2045$cesm1.cam5.1.rcp45)
mean.abs.hums.2045.rcp85<-13.24732*exp((17.67*mean.temp.data.2045$cesm1.cam5.1.rcp85)/(mean.temp.data.2045$cesm1.cam5.1.rcp85+243.5))*rh.data.2045$cesm1.cam5.1.rcp85/(273.15+mean.temp.data.2045$cesm1.cam5.1.rcp85)
mean.abs.hums.2070.rcp45<-13.24732*exp((17.67*mean.temp.data.2070$cesm1.cam5.1.rcp45)/(mean.temp.data.2070$cesm1.cam5.1.rcp45+243.5))*rh.data.2070$cesm1.cam5.1.rcp45/(273.15+mean.temp.data.2070$cesm1.cam5.1.rcp45)
mean.abs.hums.2070.rcp85<-13.24732*exp((17.67*mean.temp.data.2070$cesm1.cam5.1.rcp85)/(mean.temp.data.2070$cesm1.cam5.1.rcp85+243.5))*rh.data.2070$cesm1.cam5.1.rcp85/(273.15+mean.temp.data.2070$cesm1.cam5.1.rcp85)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.abs.hums.2020.rcp45,type="l",col=weather.colors[2],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.abs.hums.2020.rcp85,type="l",col=weather.colors[3],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.abs.hums.2045.rcp45,type="l",col=weather.colors[4],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.abs.hums.2045.rcp85,type="l",col=weather.colors[5],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.abs.hums.2070.rcp45,type="l",col=weather.colors[6],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),mean.abs.hums.2070.rcp85,type="l",col=weather.colors[7],lwd=2,lty=1)
mtext("O",side=3,adj=1,cex=2,line=.5)

plot(0,0,type="n",xlim=c(as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC")),ylim=c(0,27),xlab="",ylab="",axes=F)
mtext("mean daily rainfall (mm)",2,line = 2.25,cex=1)
axis.POSIXct(1,seq(as.POSIXct("2020-06-12 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),by="1 day"))
axis(2)
box()
grid()
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),rainfall.data.2020$cesm1.cam5.1.rcp45,type="l",col=weather.colors[2],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),rainfall.data.2020$cesm1.cam5.1.rcp85,type="l",col=weather.colors[3],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),rainfall.data.2045$cesm1.cam5.1.rcp45,type="l",col=weather.colors[4],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),rainfall.data.2045$cesm1.cam5.1.rcp85,type="l",col=weather.colors[5],lwd=2,lty=1)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),rainfall.data.2070$cesm1.cam5.1.rcp45,type="l",col=weather.colors[6],lwd=2,lty=3)
points(as.POSIXct(paste0(mean.temp.data.2020$year,"-",mean.temp.data.2020$month,"-",mean.temp.data.2020$day," 00:00:00"),tz="utc"),rainfall.data.2070$cesm1.cam5.1.rcp85,type="l",col=weather.colors[7],lwd=2,lty=1)
mtext("P",side=3,adj=1,cex=2,line=.5)

par(mar=c(3.5,0,2.5,3),xpd=F)
plot(0,0,type="n",axes=F,xlab="",ylab="",bty="n")
legend("center",
       legend=c("2020 RCP4.5", "2020 RCP8.5", "2045 RCP4.5","2024 RCP8.5","2070 RCP4.5","2070 RCP8.5"),
       lwd=4,
       cex=1.75,
       seg.len = 4,
       lty=c(3,1,3,1,3,1,3),
       col=weather.colors[-1],
       bty="n"
)

# export at dimesions 1582x842