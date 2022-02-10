library(mgcv)
library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")
delta.pustules<-subset(delta.pustules,time<=8)
pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")

site.cols<-viridis_pal(alpha=.5)(20)[c(20,15,6,1)]
weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2,1,1)])

layout(matrix(c(rep(9,10),rep(12,13),
            rep(9,10),14,1,1,1,1,2,2,2,2,3,3,3,3,
            rep(9,10),14,1,1,1,1,2,2,2,2,3,3,3,3,
            rep(9,10),14,1,1,1,1,2,2,2,2,3,3,3,3,
            rep(9,10),14,4,4,4,4,5,5,5,5,6,6,6,6,
            rep(10,10),14,4,4,4,4,5,5,5,5,6,6,6,6,
            rep(10,10),14,4,4,4,4,5,5,5,5,6,6,6,6,
            rep(10,10),14,11,11,7,7,7,7,8,8,8,8,13,13,
            rep(10,10),14,11,11,7,7,7,7,8,8,8,8,13,13,
            rep(10,10),14,11,11,7,7,7,7,8,8,8,8,13,13),
    10,23,byrow=T))
par(mar=c(4,4,1,1))

plot(pustule.model,select = 1,scale=0,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext(expression('pustule area ('*mm^2*')'),1,line = 2.25,cex=1)
mtext("s(pustule area)",2,line=2.25,cex=1)
grid()
mtext("B",adj=1,cex=1.5)
plot(pustule.model,select = 2,ylim=c(-.025,.025),shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("mean temperature (°C)",1,line = 2.25,cex=1)
mtext("s(mean temperature)",2,line=2.25,cex=1)
grid()
mtext("C",adj=1,cex=1.5)
plot(pustule.model,select = 3,ylim=c(-.025,.025),shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("max. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(max. temperature)",2,line=2.25,cex=1)
grid()
mtext("D",adj=1,cex=1.5)
plot(pustule.model,select = 4,ylim=c(-.025,.025),shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("min. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(min. temperature)",2,line=2.25,cex=1)
grid()
mtext("E",adj=1,cex=1.5)
plot(pustule.model,select = 5,ylim=c(-.025,.025),shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext(expression('mean abs. humidity ('*g/m^3*')'),1,line = 2.25,cex=1)
mtext("s(mean abs. humidity)",2,line=2.25,cex=1)
grid()
mtext("F",adj=1,cex=1.5)
plot(pustule.model,select = 6,ylim=c(-.025,.025),shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("total rainfall (mm)",1,line = 2.25,cex=1)
mtext("s(total rainfall)",2,line=2.25,cex=1)
grid()
mtext("G",adj=1,cex=1.5)
par(col=NA) #hack to get rid of qqline
plot(pustule.model,select = 7,shade=T,main="",cex.lab=1.25,cex.axis=1,ylab="",xlab="",col="black")
par(col="black")
box()
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(plant ID)",2,line=2.25,cex=1)
grid()
mtext("H",adj=1,cex=1.5)
par(col=NA) #hack to get rid of qqline
plot(pustule.model,select = 8,shade=T,main="",cex.lab=1.25,cex.axis=1,ylab="",xlab="",col=site.cols[c(2,1,3,4)],cex=1.5/par()$cex,pch=16)
par(col="black")
box()
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(site)",2,line=2.25,cex=1)
grid()
mtext("I",adj=1,cex=1.5)
legend("topleft",legend=c("CC","BT","GM","HM"),pch=16,col=site.cols,cex=1.25,bty="n",pt.cex = 1.5/par()$cex)

mtext("generalized additive model",outer=T,adj=19.5/23,cex=2,line=-2.5)
mtext(expression('change in pustule area ('*mm^2*' per day)'),outer=T,adj=20/23,cex=1.75,line=-5.5)

site.indicies<-c(2,1,3,4)[as.numeric(delta.pustules$site)]
par(mar=c(5,6,5,2))
plot(delta.pustules$area,delta.pustules$area.next,xlab = "",ylab=expression('next observed pustule area '*(mm^2)),cex.lab=2,cex.axis=2,col=site.cols[site.indicies],pch=16,cex=delta.pustules$time/2,panel.first = {abline(0,1,lty=2);grid()})
mtext(expression('observed pustule area '*(mm^2)),side=1,cex=2*par()$cex,line=3.5)
mtext("A",side=3,adj=1,cex=1.5)
legend("bottomright",legend=c("CC","BT","GM","HM","2 days","4 days","6 days"),col=c(site.cols,"grey","grey","grey"),pt.cex=c(3,3,3,3,2/2,4/2,6/2),pch=16,cex=1.75,bty="n")
mtext("data",cex=2,line=1)

library("MASS")
library("viridis")
start.area<-0.1

site<-"GM"
start.date<-c(as.POSIXct("2020-06-23 00:00:00",tz="UTC"),as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-06-24 00:00:00",tz="UTC"),as.POSIXct("2020-06-26 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
end.date<-c(as.POSIXct("2020-07-27 00:00:00",tz="UTC"),as.POSIXct("2020-07-29 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),as.POSIXct("2020-07-10 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
sim.dates<-seq.POSIXt(start.date,end.date,"3 day")
weath.data.vec<-c("observed","2020","2020","2045","2045","2070","2070")
weath.data.scenario.vec<-c(NA,"rcp45","rcp85","rcp45","rcp85","rcp45","rcp85")

plot(0,0,xlim=c(start.date,end.date),ylim=c(.1,.85),type="n",xlab="date",ylab=expression('pustule area ('*mm^2*')'),cex.lab=2,axes=F)
grid()
mtext("J",side=3,adj=1,cex=1.5)
axis.POSIXct(1,sim.dates,cex.axis=2)
axis(2,cex.axis=2)
box()
mtext("projections",cex=2,line=1)

for(i in 2:7)
{
  weath.data<-weath.data.vec[i]
  weath.data.scenario<-weath.data.scenario.vec[i]
  
  xcords<-rep(NA,length(sim.dates)) #time values
  ycords<-rep(NA,length(sim.dates)) #area values
  
  set.seed(289988)
  
  for(j in 1:100)
  {
    area<-start.area
    xcords.new<-c(sim.dates[1])
    ycords.new<-c(area)
    
    for(k in 1:(length(sim.dates)-1))
    {
      date0<-sim.dates[k]
      date1<-sim.dates[k+1]
      
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
      beta <- coef(pustule.model) ## posterior mean of coefs
      Vb   <- vcov(pustule.model) ## posterior  cov of coefs
      n <-2
      mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
      pred.data<-data.frame("area"=area,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.mean.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain,tag="NA",site=site)
      Xp <- predict(pustule.model, newdata = pred.data, exclude=c("s(tag)"),type="lpmatrix")
      ilink <- family(pustule.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      area.change<-preds[1]
      area<-area+area.change*as.numeric(date1-date0)
      if(area<0) {area<-0}
      
      #reps<-reps+pred.window
      xcords.new<-c(xcords.new,date1)
      ycords.new<-c(ycords.new,area)
    }
    xcords<-rbind(xcords,xcords.new)
    ycords<-rbind(ycords,ycords.new)  
  }
  points(xcords[2,],colMeans(ycords[-1,]),col=weather.colors[i],type="l",lwd=5,lty=c(1,3,1,3,1,3,1,3,1)[i])
}
legend("topleft",
       legend=c("2020 RCP4.5", "2020 RCP8.5", "2045 RCP4.5","2045 RCP8.5","2070 RCP4.5","2070 RCP8.5"),
       cex=1.75,
       lwd=4,
       seg.len = 2,
       lty=c(1,3,1,3,1,3,1)[-1],
       col=weather.colors[-1],
       bty="n"
)

# export at 1291x812
