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

dummy.log.10.tot.spore.deposition.per.day<-seq(-6.5,1,.01)
site<-"GM"
date0<-as.POSIXct("2020-07-07 00:00:00",tz="UTC")
date1<-as.POSIXct("2020-07-14 00:00:00",tz="UTC")
weath.data.scenario.vec<-c("rcp45","rcp85")
weath.data.vec<-c("2020","2030","2040","2045","2050","2060","2070")
plot.height<-25

colors<-c("lightblue4","lightblue2")
par(mar=c(4,6,2,2))
plot(0,0,type="n",xlim=c(2020,2070),ylim=c(0,.8),xlab="",ylab=expression('odds of infection for '*log[10]*' spore deposition = 0'),cex.lab=2,cex.axis=2)
grid()
legend("topright",legend=c("RCP4.5","RCP8.5"),cex=2,lwd=4,lty=c(2,1),col=colors,bty="n",seg.len = 3)
ltys<-c(2,1)
for(i in 1:length(weath.data.scenario.vec))
{
  weath.data.scenario<-weath.data.scenario.vec[i]
  final.preds<-c()

  predictions<-c()
  for(j in 1:length(weath.data.vec))
  {
    weath.data<-weath.data.vec[j]
    
    set.seed(874627)
    
      ycords.new<-c()
        
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
        
        new.prediction<-predict(transmission.model,newdata = new.data[which(new.data$log.10.spore.deposition.per.day==0),],type = "response")
     
        predictions<-rbind(predictions,new.prediction)  


  }
  points(as.numeric(weath.data.vec),predictions,type="l",lty=ltys[i],lwd=4,col=colors[i])
}

#export at 1202x777
