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
  temp.days<-c()
  temp.days.16.22<-c()
  temp.days.7.30<-c()
  dew.point.days<-c()
  temp.dew.point.days<-c()
  temp.16.22.dew.point.days<-c()
  temp.7.30.dew.point.days<-c()
  wetness.days<-c()
  temp.wetness.days<-c()
  temp.16.22.wetness.days<-c()
  temp.7.30.wetness.days<-c()
  tot.rains<-c()
  solar.days<-c()
  wind.speed.days<-c()
  gust.speed.days<-c()

  for (tag in unique(plant.heights$tag))
  {
    if(tag == "918") next # only one obs
    if(tag == "269") next # only one obs
    sub.plant.heights<-plant.heights[which(plant.heights$tag==tag),]
    
    for(i in 1:(dim(sub.plant.heights)[1]-1))
    {
      #pull reference data
      date0<-sub.plant.heights[i,"date"]
      date1<-sub.plant.heights[i+1,"date"]
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
      new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
      new.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #time (in days) during which temp between 16 and 22 celsius
      new.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #time (in days) during which temp between 7 and 30 celsius
      new.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T) #Dew point days
      
      #calculate weather metrics
      new.wetness.days<-sum(weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      new.tot.rain<-sum(weath.sub$rain,na.rm=T)
      new.solar.days<-sum(weath.sub$solar.radiation*weath.sub$interval.length,na.rm = T)
      new.wind.speed.days<-sum(weath.sub$wind.speed*weath.sub$interval.length,na.rm = T)
      new.gust.speed.days<-sum(weath.sub$wind.direction*weath.sub$interval.length,na.rm = T)
      
      #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
      new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
      new.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
      new.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)
      new.temp.wetness.days<-sum(weath.sub$temp*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      new.temp.16.22.wetness.days<-sum(weath.sub$temp.16.22*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      new.temp.7.30.wetness.days<-sum(weath.sub$temp.7.30*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      
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
      
      temp.days.16.22<-c(temp.days.16.22,new.temp.days.16.22)
      temp.days.7.30<-c(temp.days.7.30,new.temp.days.7.30)
      temp.days<-c(temp.days,new.temp.days)
      dew.point.days<-c(dew.point.days,new.dew.point.days)
      temp.dew.point.days<-c(temp.dew.point.days,new.temp.dew.point.days)
      temp.16.22.dew.point.days<-c(temp.16.22.dew.point.days,new.temp.16.22.dew.point.days)
      temp.7.30.dew.point.days<-c(temp.7.30.dew.point.days,new.temp.7.30.dew.point.days)
      
      wetness.days<-c(wetness.days,new.wetness.days)
      temp.wetness.days<-c(temp.wetness.days,new.temp.wetness.days)
      temp.16.22.wetness.days<-c(temp.16.22.wetness.days,new.temp.16.22.wetness.days)
      temp.7.30.wetness.days<-c(temp.7.30.wetness.days,new.temp.7.30.wetness.days)
      tot.rains<-c(tot.rains,new.tot.rain)
      solar.days<-c(solar.days,new.solar.days)
      wind.speed.days<-c(wind.speed.days,new.wind.speed.days)
      gust.speed.days<-c(gust.speed.days,new.gust.speed.days)
    }
  }
  delta.height<-data.frame(tag=factor(tags),site=factor(sites), height=start.vals, height.next=end.vals,inf.intens=inf.intens.vals,time=days,
                             temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,
                             dew.point.days=dew.point.days,temp.dew.point.days=temp.dew.point.days,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,
                             wetness.days=wetness.days,temp.wetness.days=temp.wetness.days,temp.16.22.wetness.days=temp.16.22.wetness.days,temp.7.30.wetness.days=temp.7.30.wetness.days,
                             tot.rain=tot.rains,solar.days=solar.days,wind.speed.days=wind.speed.days,gust.speed.days=gust.speed.days)
  
  saveRDS(plant.heights,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS")
  saveRDS(delta.height,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.height.RDS")
}
plant.heights<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS")
delta.height<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.height.RDS")
