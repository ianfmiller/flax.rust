library(mgcv)
library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/transmission data set building.R")
transmission.data<-transmission.data[which(transmission.data$status==0),]
transmission.data$site<-as.factor(transmission.data$site)
for(i in 1:nrow(transmission.data))
{
  if(is.na(transmission.data[i,"tag"])) {transmission.data[i,"tag"]<-paste0(transmission.data[i,"site"],"X",transmission.data[i,"X"]+transmission.data[i,"x"],"Y",transmission.data[i,"Y"]+transmission.data[i,"y"])}
}
transmission.data$log.10.spore.deposition.per.day<-log10(transmission.data$tot.spore.deposition/transmission.data$time)

transmission.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")

site.cols<-viridis_pal(alpha=.5)(20)[c(20,15,6,1)]

layout(matrix(c(rep(16,10),1,1,1,1,2,3,3,3,3,4,4,4,4,rep(17,10),rep(10,10),1,1,1,1,2,3,3,3,3,4,4,4,4,rep(11,10),rep(10,10),12,5,5,5,5,6,6,6,6,7,7,7,7,rep(11,10),rep(10,10),12,5,5,5,5,6,6,6,6,7,7,7,7,rep(11,10),rep(10,10),13,14,14,8,8,8,8,9,9,9,9,15,15,rep(11,10),rep(18,10),13,14,14,8,8,8,8,9,9,9,9,15,15,rep(19,10)),6,33,byrow=T))

# B-I
par(mar=c(4,4,3,1))
plot(transmission.model,select = 1,scheme = 2,zlim=c(-12,12),too.far=1,hcolors=rev(heat.colors(101)),contour.col="black",cex.lab=1.5,cex.axis=1.5,xlab="",ylab="",main="",rug=F)
mtext(expression(log[10]*'(predicted spore deposition per day)'),1,line = 2.25,cex=.9)
mtext("plant height (cm)",2,line=2.25,cex=1)
points(transmission.data$log.10.spore.deposition.per.day,transmission.data$height.cm,pch=16,cex=.1)

par(mar=c(4,0,2.5,2.5))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-12-.12,12+.12),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-12,12,length.out=101)[i]
  rect(0,ii-.12,1,ii+.12,col=rev(heat.colors(101))[i],border = NA)
}
rect(0,-12-.12,1,12+.12)
mtext("B",cex=2)
axis(4,cex.axis=1,tck=-.5,padj=-1)

par(mar=c(4,4,3,1))
plot(transmission.model,select = 2,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("mean temperature (°C)",1,line = 2.25,cex=1)
mtext("s(mean temperature)",2,line=2.25,cex=1)
grid()
mtext("C",adj=1,cex=2)
plot(transmission.model,select = 3,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("max. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(max. temperature)",2,line=2.25,cex=1)
grid()
mtext("D",adj=1,cex=2)
plot(transmission.model,select = 4,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("min. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(min. temperature)",2,line=2.25,cex=1)
grid()
mtext("E",adj=1,cex=2)
plot(transmission.model,select = 5,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext(expression('mean abs. humidity ('*g/m^3*')'),1,line = 2.25,cex=1)
mtext("s(mean abs. humidity)",2,line=2.25,cex=1)
grid()
mtext("F",adj=1,cex=2)
plot(transmission.model,select = 6,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("total rainfall (mm)",1,line = 2.25,cex=1)
mtext("s(total rainfall)",2,line=2.25,cex=1)
grid()
mtext("G",adj=1,cex=2)
plot(transmission.model,select = 7,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(tag)",2,line=2.25,cex=1)
grid()
mtext("H",adj=1,cex=2)
plot(transmission.model,select = 8,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",col=site.cols[c(2,1,3,4)],cex=1.5/par()$cex,pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(site)",2,line=2.25,cex=1)
grid()
mtext("I",adj=1,cex=2)
legend("topleft",legend=c("CC","BT","GM","HM"),pch=16,col=site.cols,cex=1.25,bty="n",pt.cex = 1.5/par()$cex)

# A
set.seed(98579835)
transmission.data.shuffled= transmission.data[sample(1:nrow(transmission.data)), ]
site.indicies<-c(2,1,3,4)[as.numeric(transmission.data.shuffled$site)]
par(mar=c(5,6,5,2))
plot(jitter(log10(transmission.data.shuffled$tot.spore.deposition)),jitter(transmission.data.shuffled$status.next),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,cex.lab=2,col=site.cols[site.indicies],pch=16,cex=transmission.data.shuffled$time/2,panel.first = grid())
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
mtext("A",side=3,adj=1,cex=2)
legend("topleft",legend=c("CC","BT","GM","HM","2 days","4 days","6 days"),col=c(site.cols,"grey","grey","grey"),pt.cex=c(3,3,3,3,2/2,4/2,6/2),pch=16,cex=1.75,bty="n")

# J
dummy.log.10.tot.spore.deposition.per.day<-seq(-6.5,1,.01)

site<-"GM"

date0<-as.POSIXct("2020-07-07 00:00:00",tz="UTC")
date1<-as.POSIXct("2020-07-14 00:00:00",tz="UTC")
sim.dates<-date0

weath.data.vec<-c("observed","2020","2020","2045","2045","2070","2070")
weath.data.scenario.vec<-c(NA,"rcp45","rcp85","rcp45","rcp85","rcp45","rcp85")
weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2,1,1)])

plot.height<-25
plot(0,0,type="n",xlab=expression('predicted '*log[10]*'(spore deposition per day)'),ylab="odds of infection",cex.lab=2,cex.axis=2,ylim=c(0,.9),xlim=c(-6.5,1),panel.first = grid())
legend("topleft",
       legend=c("observed weather","2020 RCP4.5", "2020 RCP8.5", "2045 RCP4.5","2024 RCP8.5","2070 RCP4.5","2070 RCP8.5"),
       cex=1.75,
       lwd=4,
       seg.len = 2,
       lty=c(1,3,1,3,1,3,1,3),
       col=weather.colors,
       bty="n"
)

mtext("J",adj=1,cex=2)

for(j in 1:7)
{
  weath.data<-weath.data.vec[j]
  weath.data.scenario<-weath.data.scenario.vec[j]
  
  xcords<-rep(NA,length(sim.dates)) #time values
  inf.intens.cords<-rep(NA,length(sim.dates)) #infection intensity values
  height.cords<-rep(NA,length(sim.dates)) #height values
  
  if(!(weath.data=="observed"))
  {
    mean.temp.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/Tair.csv"),header = F)
    max.temp.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/tasmax.csv"),header = F)
    min.temp.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/tasmin.csv"),header = F)
    rh.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/relHumid.csv"),header = F)
    rainfall.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/rainfall.csv"),header = F)
    colnames(mean.temp.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
    colnames(max.temp.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
    colnames(min.temp.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
    colnames(rh.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
    colnames(rainfall.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
    
    mean.temp.data.sub.1<-mean.temp.data[which(format(as.Date(paste0(mean.temp.data$year,"-",mean.temp.data$month,"-",mean.temp.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),] #compare only month and day
    mean.temp.data.sub.2<-mean.temp.data.sub.1[which(format(as.Date(paste0(mean.temp.data.sub.1$year,"-",mean.temp.data.sub.1$month,"-",mean.temp.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
    
    max.temp.data.sub.1<-max.temp.data[which(format(as.Date(paste0(max.temp.data$year,"-",max.temp.data$month,"-",max.temp.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
    max.temp.data.sub.2<-max.temp.data.sub.1[which(format(as.Date(paste0(max.temp.data.sub.1$year,"-",max.temp.data.sub.1$month,"-",max.temp.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
    
    min.temp.data.sub.1<-min.temp.data[which(format(as.Date(paste0(min.temp.data$year,"-",min.temp.data$month,"-",min.temp.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
    min.temp.data.sub.2<-min.temp.data.sub.1[which(format(as.Date(paste0(min.temp.data.sub.1$year,"-",min.temp.data.sub.1$month,"-",min.temp.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
    
    rh.data.sub.1<-rh.data[which(format(as.Date(paste0(rh.data$year,"-",rh.data$month,"-",rh.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
    rh.data.sub.2<-rh.data.sub.1[which(format(as.Date(paste0(rh.data.sub.1$year,"-",rh.data.sub.1$month,"-",rh.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
    
    rainfall.data.sub.1<-rainfall.data[which(format(as.Date(paste0(rainfall.data$year,"-",rainfall.data$month,"-",rainfall.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
    rainfall.data.sub.2<-rainfall.data.sub.1[which(format(as.Date(paste0(rainfall.data.sub.1$year,"-",rainfall.data.sub.1$month,"-",rainfall.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
    
    new.mean.temp<-mean(mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
    new.max.temp<-max(max.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
    new.min.temp<-min(min.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
    abs.hums<-13.24732*exp((17.67*mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])/(mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)]+243.5))*rh.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)]/(273.15+mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
    new.mean.abs.hum<-mean(abs.hums)
    new.mean.daily.rain<-mean(rainfall.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
    
  }
  
  if(weath.data=="observed")
  {
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
  }
  
  new.data<-transmission.data[1:length(dummy.log.10.tot.spore.deposition.per.day),]
  new.data$log.10.spore.deposition.per.day<-dummy.log.10.tot.spore.deposition.per.day
  new.data$mean.temp<-new.mean.temp
  new.data$max.temp<-new.max.temp
  new.data$min.temp<-new.min.temp
  new.data$mean.abs.hum<-new.mean.abs.hum
  new.data$mean.daily.rain<-new.mean.temp
  new.data$site<-site
  new.data$height.cm<-plot.height
  
  points(new.data$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data,type = "response"),type="l",col=weather.colors[j],lwd=4,lty=c(1,3,1,3,1,3,1,3,1)[j])
}

# export at 1564x719

