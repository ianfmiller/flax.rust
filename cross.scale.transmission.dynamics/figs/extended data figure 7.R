# load data and packages
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

par(mfrow=c(3,3),mar=c(5,5,5,3))

# A
site.cols<-viridis_pal(alpha=.5)(20)[c(20,15,6,1)]
set.seed(98579835)
transmission.data.shuffled= transmission.data[sample(1:nrow(transmission.data)), ]
site.indicies<-c(2,1,3,4)[as.numeric(transmission.data.shuffled$site)]
plot(transmission.data.shuffled$log.10.spore.deposition.per.day,jitter(transmission.data.shuffled$time),xlab=expression(log[10]*' spore deposition per day'),ylab="days",col=site.cols[site.indicies],pch=16,panel.first = grid(),cex.lab=2,cex.axis=2,main="predicted data",cex.main=2)
legend("right",legend=c("CC","BT","GM","HM"),col=c(site.cols),pt.cex=2,pch=16,cex=1.5,bty="n")
mtext("A",side=3,adj=1,cex=1.5)

#B
dummy.log.10.tot.spore.deposition.per.day<-seq(-6.5,1,.01)

site<-"GM"
date0<-as.POSIXct("2020-07-07 00:00:00",tz="UTC")
date1<-as.POSIXct("2020-07-14 00:00:00",tz="UTC")
sim.dates<-date0

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

new.data<-transmission.data[1:length(dummy.log.10.tot.spore.deposition.per.day),]
new.data$log.10.spore.deposition.per.day<-dummy.log.10.tot.spore.deposition.per.day
new.data$time<-7
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
new.data.3$height.cm<-20
new.data.4<-new.data
new.data.4$height.cm<-30
new.data.5<-new.data
new.data.5$height.cm<-40
new.data.6<-new.data
new.data.6$height.cm<-50

exclude.vec<-c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(mean.daily.rain)")

colors<-colorRampPalette(c('palegreen1',"palegreen4"))(6)
plot(0,0,type="n",xlab=expression('predicted '*log[10]*' spore deposition per day'),ylab="odds of infection",cex.lab=2,cex.axis=2,ylim=c(0,1),xlim=c(-6.5,1),panel.first = grid(),main="plant height effect",cex.main=2)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.1,type = "response",exclude = exclude.vec),type="l",col=colors[1],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.2,type = "response",exclude = exclude.vec),type="l",col=colors[2],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.3,type = "response",exclude = exclude.vec),type="l",col=colors[3],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.4,type = "response",exclude = exclude.vec),type="l",col=colors[4],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.5,type = "response",exclude = exclude.vec),type="l",col=colors[5],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.6,type = "response",exclude = exclude.vec),type="l",col=colors[6],lwd=4,lty=1)

legend("topleft",legend=c("50 cm","40 cm","30 cm","20 cm","10 cm","5 cm"),col=rev(colors),lty=1,bty="n",cex=1.25,lwd=4)
mtext("B",side=3,adj=1,cex=1.5)

#C
dummy.log.10.tot.spore.deposition.per.day<-seq(-6.5,1,.01)

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

new.data<-transmission.data[1:length(dummy.log.10.tot.spore.deposition.per.day),]
new.data$log.10.spore.deposition.per.day<-dummy.log.10.tot.spore.deposition.per.day
new.data$height.cm<-25
new.data$mean.temp<-new.mean.temp
new.data$max.temp<-new.max.temp
new.data$min.temp<-new.min.temp
new.data$mean.abs.hum<-new.mean.abs.hum
new.data$mean.daily.rain<-new.mean.temp
new.data$site<-site

new.data.1<-new.data
new.data.1$time<-1
new.data.2<-new.data
new.data.2$time<-2
new.data.3<-new.data
new.data.3$time<-3
new.data.4<-new.data
new.data.4$time<-4
new.data.5<-new.data
new.data.5$time<-5
new.data.6<-new.data
new.data.6$time<-6
new.data.7<-new.data
new.data.7$time<-7

exclude.vec<-c("s(site)","s(tag)","s(mean.temp)","s(max.temp)","s(min.temp)","s(mean.abs.hum)","s(mean.daily.rain)")

colors<-colorRampPalette(c('skyblue1',"skyblue4"))(7)
plot(0,0,type="n",xlab=expression('predicted '*log[10]*' spore deposition per day'),ylab="odds of infection",cex.lab=2,cex.axis=2,ylim=c(0,.35),xlim=c(-6.5,1),panel.first = grid(),main="time effect",cex.main=2)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.1,type = "response",exclude = exclude.vec),type="l",col=colors[1],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.2,type = "response",exclude = exclude.vec),type="l",col=colors[2],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.3,type = "response",exclude = exclude.vec),type="l",col=colors[3],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.4,type = "response",exclude = exclude.vec),type="l",col=colors[4],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.5,type = "response",exclude = exclude.vec),type="l",col=colors[5],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.6,type = "response",exclude = exclude.vec),type="l",col=colors[6],lwd=4,lty=1)
points(new.data.1$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data.7,type = "response",exclude = exclude.vec),type="l",col=colors[7],lwd=4,lty=1)

legend("topleft",legend=c("7 days","6 days","5 days","4 days","3 days","2 days","1 day"),col=rev(colors),lty=1,bty="n",cex=1.25,lwd=4)
mtext("C",side=3,adj=1,cex=1.5)

# D-I
dummy.log.10.tot.spore.deposition.per.day<-seq(-6.5,1,.01)

date0<-as.POSIXct("2020-07-07 00:00:00",tz="UTC")
date1<-as.POSIXct("2020-07-14 00:00:00",tz="UTC")
sim.dates<-date0

weath.data.vec<-c("observed","2020","2020","2045","2045","2070","2070")
weath.data.scenario.vec<-c(NA,"rcp45","rcp85","rcp45","rcp85","rcp45","rcp85")
weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2,1,1)])

for(i in 1:6)
{
  plot.height<-c(5,10,20,30,40,50)[i]
  plot(0,0,type="n",xlab=expression('predicted '*log[10]*' spore deposition per day'),ylab="odds of infection",cex.lab=2,cex.axis=2,ylim=c(0,1),xlim=c(-6.5,1),main = paste0("plant height = ",plot.height," cm"),cex.main=2,panel.first = grid())
  if(i==1)
  {
    legend("topleft",
           legend=c("observed weather","2020 RCP4.5", "2020 RCP8.5", "2045 RCP4.5","2024 RCP8.5","2070 RCP4.5","2070 RCP8.5"),
           lwd=4,
           seg.len = 4,
           lty=c(1,3,1,3,1,3,1,3),
           col=weather.colors,
           bty="n",
           cex=1.75
    )
  }
  mtext(c("D","E","F","G","H","I")[i],adj=1,cex=1.5)
  
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
    
    points(new.data$log.10.spore.deposition.per.day,predict(transmission.model,newdata = new.data,type = "response"),type="l",col=weather.colors[j],lwd=3,lty=c(1,3,1,3,1,3,1,3,1)[j])
  }
}

#export at 1359 x 880
