## Plant Growth
library(lubridate)
if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plant.heights.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS")))
{
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  
  # clean healthy plant data
  healthyraw <- read.csv("~/Documents/GitHub/flax.rust/data/healthyplants.csv")
  healthy <- healthyraw[!is.na(healthyraw$max.height),]
  tags <- unique(healthy$Tag)
  
  healthyindiv <- data.frame()
  for (i in 1:length(tags)) {
    tag <- tags[i]
    temp <- healthy[healthy$Tag == tag,]
    uniquedates <- unique(temp$Date)
    if(length(uniquedates) == 1) next
    for (j in 1:length(uniquedates)) {
      date <- uniquedates[j]
      temp2 <- temp[temp$Date == date,]
      temp2[1,]
      healthyindiv<- rbind(healthyindiv, temp2[1,])
    }
  }
  healthyindiv$Date <- mdy(healthyindiv$Date)
  healthyindiv$plant.inf.intens <- 0
  
  ## load diseased plant data
  plantsraw <- readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
  plants <- plantsraw[!is.na(plantsraw$max.height),]
  
  ## merge data
  subs <- c("Year", "Site", "Tag", "Date", "max.height", "plant.inf.intens")
  plant.heights <- rbind(healthyindiv[subs], plants[subs])
  colnames(plant.heights) <- c("year", "site", "tag", "date", "max.height", "plant.inf.intens")
  
  ## cut out incomplete records
  plant.heights<-plant.heights[-intersect(which(plant.heights$site=="HM"),which(plant.heights$date>"2020-07-10 00:00:00 UTC")),]

  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  sites<-c()
  start.vals<-c()
  end.vals<-c()
  inf.intens.vals<-c()
  days<-c()
  mean.temp<-c()
  max.temp<-c()
  min.temp<-c()
  mean.abs.hum<-c() #absolute humidity
  max.abs.hum<-c()
  min.abs.hum<-c()
  mean.vpd<-c() #vapor pressure deficit
  max.vpd<-c() 
  min.vpd<-c() 
  mean.wetness<-c()
  tot.rain<-c()
  mean.solar<-c()

  for (tag in unique(plant.heights$tag))
  {
    if(tag == "918") next # only one obs
    if(tag == "269") next # only one obs
    sub.plant.heights<-plant.heights[which(plant.heights$tag==tag),]
    
    for(i in 1:(dim(sub.plant.heights)[1]-1))
    {
      #pull reference data
      date0<-sub.plant.heights[i,"date"]
      date0<-as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")
      date1<-sub.plant.heights[i+1,"date"]
      date1<-as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")
      site<-sub.plant.heights[i,"site"]
      
      #subset temp data to relevant window
      temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
      temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #pull out relevant data
      temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #pull out relevant data
      temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #throw out NAs
      if(dim(temp.rh.sub)[1]==0) {next}
      temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
      
      #subset weather data to relevant window
      weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
      weath.sub<-subset(weath.sub,date<=date1) #pull out relevant data
      weath.sub<-subset(weath.sub,date>=date0) #pull out relevant data
      
      if(dim(weath.sub)[1]==0) {next}
      weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
      
      #calculate environmental variable metrics
      new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
      new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
      new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
      
      abs.hum<-6.112*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh*2.1674/(273.15+T)
      new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
      new.max.abs.hum<-max(abs.hum,na.rm=T)
      new.min.abs.hum<-min(abs.hum,na.rm=T)
      
      svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
      avps<- temp.rh.sub$rh / 100 * svps #actual vapor pressures 
      vpds<-avps-svps
      
      new.mean.vpd<-mean(vpds,na.rm=T)
      new.max.vpd<-max(vpds,na.rm=T)
      new.min.vpd<-min(vpds,na.rm=T)
      
      new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
      new.tot.rain<-sum(weath.sub$rain,na.rm=T)
      new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
      
      #pull out core predictors
      start.val<-sub.plant.heights[i,"max.height"]
      end.val<-sub.plant.heights[i+1,"max.height"]
      delta.days<-date1-date0
      
      #store values
      tags<-c(tags,tag)
      sites<-c(sites,site)
      
      start.vals<-c(start.vals,start.val)
      end.vals<-c(end.vals,end.val)
      inf.intens.vals<-c(inf.intens.vals,sub.plant.heights[i,"plant.inf.intens"])
      days<-c(days,delta.days)
      
      mean.temp<-c(mean.temp,new.mean.temp)
      max.temp<-c(max.temp,new.max.temp)
      min.temp<-c(min.temp,new.min.temp)
      mean.abs.hum<-c(mean.abs.hum,new.mean.abs.hum)
      max.abs.hum<-c(max.abs.hum,new.max.abs.hum)
      min.abs.hum<-c(min.abs.hum,new.min.abs.hum)
      mean.vpd<-c(mean.vpd,new.mean.vpd)
      max.vpd<-c(max.vpd,new.max.vpd)
      min.vpd<-c(min.vpd,new.min.vpd)
      mean.wetness<-c(mean.wetness,new.mean.wetness)
      tot.rain<-c(tot.rain,new.tot.rain)
      mean.solar<-c(mean.solar,new.mean.solar)
    }
  }
  delta.height<-data.frame(tag=factor(tags),site=factor(sites), height=start.vals, height.next=end.vals,inf.intens=inf.intens.vals,time=days,
                           mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,
                           mean.abs.hum=mean.abs.hum,max.abs.hum=max.abs.hum,min.abs.hum=min.abs.hum,
                           mean.vpd=mean.vpd,max.vpd=max.vpd,min.vpd=min.vpd,
                           mean.wetness=mean.wetness,tot.rain=tot.rain,mean.solar=mean.solar)
  
  saveRDS(plant.heights,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS")
  saveRDS(delta.height,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.height.RDS")
}
plant.heights<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS")
delta.height<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.height.RDS")
