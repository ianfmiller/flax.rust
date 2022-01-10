source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")

setwd("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020")

par(mfrow=c(2,2))

mean.temp<-read.csv("Tair.csv",header = T)
dates<-as.Date(paste0(mean.temp$year,"-",mean.temp$month,"-",mean.temp$day))
plot(dates,mean.temp$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(0,20),xlab="date",ylab="mean temp")        
points(dates,mean.temp$cesm1.cam5.1.rcp85,type="l",col="red")

for(site in c("BT"))
{
  mean.temps<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      mean.temps<-rbind(mean.temps,data.frame(date,mean(sub.temp.rh.2$temp.c,na.rm = T)))
    }
  }
  points(as.Date(mean.temps[,1],origin="1970-01-01"),mean.temps[,2],type="l")
}

mean.rel.hum<-read.csv("relHumid.csv",header = T)
dates<-as.Date(paste0(mean.rel.hum$year,"-",mean.rel.hum$month,"-",mean.rel.hum$day))
plot(dates,mean.rel.hum$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(30,90),xlab="date",ylab="mean RH")        
points(dates,mean.rel.hum$cesm1.cam5.1.rcp85,type="l",col="red")

for(site in c("BT"))
{
  rel.hums<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      rel.hums<-rbind(rel.hums,data.frame(date,mean(sub.temp.rh.2$rh,na.rm = T)))
    }
  }
  points(as.Date(rel.hums[,1],origin="1970-01-01"),rel.hums[,2],type="l")
}


mean.mean.rel.hum<-read.csv("relHumid.csv",header = T)
mean.temp<-read.csv("Tair.csv",header = T)
dates<-as.Date(paste0(mean.rel.hum$year,"-",mean.rel.hum$month,"-",mean.rel.hum$day))
plot(dates,13.24732*exp((17.67*mean.temp$cesm1.cam5.1.rcp45)/(mean.temp$cesm1.cam5.1.rcp45+243.5))*mean.rel.hum$cesm1.cam5.1.rcp45/(273.15+mean.temp$cesm1.cam5.1.rcp45),type="l",col="orange",ylim=c(2,10),xlab="date",ylab="mean AH")        
points(dates,13.24732*exp((17.67*mean.temp$cesm1.cam5.1.rcp85)/(mean.temp$cesm1.cam5.1.rcp85+243.5))*mean.rel.hum$cesm1.cam5.1.rcp85/(273.15+mean.temp$cesm1.cam5.1.rcp85),type="l",col="red")

for(site in c("BT"))
{
  mean.abs.hums<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      abs.hums<-13.24732*exp((17.67*sub.temp.rh.2$temp.c)/(sub.temp.rh.2$temp.c+243.5))*sub.temp.rh.2$rh/(273.15+sub.temp.rh.2$temp.c)
      abs.hums<-6.112*exp((17.67*sub.temp.rh.2$temp.c/(sub.temp.rh.2$temp.c+243.5))*(sub.temp.rh.2$rh/100)*18.02)/((273.15+sub.temp.rh.2$temp.c)*100*0.08314)
      abs.hums<-6.112*exp(17.67*sub.temp.rh.2$temp.c/(sub.temp.rh.2$temp.c+243.5) )*sub.temp.rh.2$rh*2.1674/(273.15+sub.temp.rh.2$temp.c)
      mean.abs.hums<-rbind(mean.abs.hums,data.frame(date,mean(abs.hums,na.rm = T)))
    }
  }
  points(as.Date(mean.abs.hums[,1],origin="1970-01-01"),mean.abs.hums[,2],type="l")
}


mean.rain<-read.csv("rainfall.csv",header = T)
dates<-as.Date(paste0(mean.rain$year,"-",mean.rain$month,"-",mean.rain$day))
plot(dates,mean.rain$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(0,20),xlab="date",ylab="mean rainfall")        
points(dates,mean.rain$cesm1.cam5.1.rcp85,type="l",col="red")

for(site in c("BT"))
{
  mean.rains<-c()
  all.weath.1<-all.weath[which(all.weath$site==site),]
  for(date in dates)
  {
    if(date %in% as.Date(all.weath.1$date))
    {
      all.weath.2<-all.weath.1[which(as.Date(all.weath.1$date)==date),]
      mean.rains<-rbind(mean.rains,data.frame(date,mean(all.weath.2$rain,na.rm = T)*12*24))
    }
  }
  points(as.Date(mean.rains[,1],origin="1970-01-01"),mean.rains[,2],type="l")
}
