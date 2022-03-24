library(mgcv)
library(MASS)
library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")
delta.plant.heights<-subset(delta.plant.heights,time<=8)
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

site<-"GM"
start.date<-c(as.POSIXct("2020-06-23 00:00:00",tz="UTC"),as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-06-24 00:00:00",tz="UTC"),as.POSIXct("2020-06-26 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
end.date<-c(as.POSIXct("2020-07-27 00:00:00",tz="UTC"),as.POSIXct("2020-07-29 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),as.POSIXct("2020-07-10 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
sim.dates<-seq.POSIXt(start.date,end.date,"7 day")
weath.data.scenario.vec<-c("rcp45","rcp85")
weath.data.vec<-c("2020","2030","2040","2045","2050","2060","2070")
start.height<-10

colors<-c("lightblue4","lightblue2")
par(mar=c(4,6,2,2))
plot(0,0,type="n",xlim=c(2020,2070),ylim=c(15,30),xlab="",ylab="predicted plant height (cm) on July 22",cex.lab=2,cex.axis=2)
grid()
legend("topright",legend=c("RCP4.5","RCP8.5"),cex=2,lwd=4,lty=c(2,1),col=colors,bty="n",seg.len = 4)
ltys<-c(2,1)
for(i in 1:length(weath.data.scenario.vec))
{
  weath.data.scenario<-weath.data.scenario.vec[i]
  final.preds<-c()
  final.preds.1<-c()
  final.preds.9<-c()
  
  for(j in 1:length(weath.data.vec))
  {
    weath.data<-weath.data.vec[j]
    
    xcords<-rep(NA,length(sim.dates)) #time values
    ycords<-rep(NA,length(sim.dates)) #area values
    
    set.seed(874627)
    for(k in 1:100)
    {
      height<-start.height
      xcords.new<-c(sim.dates[1])
      ycords.new<-c(height)
      
      for (l in 1:(length(sim.dates)-1))
      {
        date0<-sim.dates[l]
        date1<-sim.dates[l+1]
        
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
        
        beta <- coef(plant.growth.model) ## posterior mean of coefs
        Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
        n <-2
        mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
        pred.data<-data.frame("height"=height,"inf.intens"=0,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.mean.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain,tag="NA",site=site)
        Xp <- predict(plant.growth.model, newdata = pred.data, exclude=c("s(plant ID)"),type="lpmatrix")
        ilink <- family(plant.growth.model)$linkinv
        preds <- rep(NA,n)
        for (l in seq_len(n)) { 
          preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
        }
        height.change<-preds[1]
        height<-height+height.change*as.numeric(date1-date0)
        if(height<5) {height<-5}
        
        xcords.new<-c(xcords.new,date1)
        ycords.new<-c(ycords.new,height)
      }
      xcords<-rbind(xcords,xcords.new)
      ycords<-rbind(ycords,ycords.new)  
    }
    final.preds.1<-c(final.preds.1,quantile(ycords[-1,ncol(ycords)],.1))
    final.preds<-c(final.preds,mean(ycords[-1,ncol(ycords)]))
    final.preds.9<-c(final.preds.9,quantile(ycords[-1,ncol(ycords)],.9))
    
  }
  polygon(c(2020,2030,2040,2045,2050,2060,2070,2070,2060,2050,2045,2040,2030,2020),c(final.preds.9,rev(final.preds.1)),col=colors[i],density=50,angle=c(45,135)[i])
  points(as.numeric(weath.data.vec),final.preds,type="l",lty=ltys[i],lwd=4,col=colors[i])
}

#export at 1202x777

