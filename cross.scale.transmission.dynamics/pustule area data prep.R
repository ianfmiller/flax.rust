if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.pustules.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pustules.RDS")))
{
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
    
  # clean data
  
  ## load data
  pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv")
  
  ## add in area
  area<-pi*(pustules$max.diam..mm./4)*(pustules$min.diam..mm./4)
  pustules$date<-as.POSIXct(pustules$date,tryFormats = "%m/%d/%y %H:%M",tz="UTC")
  pustules$area<-area
  
  ## clean data
  pustules<-pustules[which(pustules$pustule.ID.confidence=="Yes"),]
  pustules<-pustules[which(pustules$area>0),]
  
  ## cut out incomplete records
  pustules<-pustules[-intersect(which(pustules$site=="HM"),which(pustules$date>"2020-07-10 00:00:00 UTC")),]
  
  ## make new data object for change in pustule size
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  sites<-c()
  stem.iters<-c()
  leaf.iters<-c()
  pustule.nums<-c()
  start.vals<-c()
  end.vals<-c()
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
  measurer.ids<-c()
  
  
  for (tag in unique(pustules$tag))
  {
    sub.pustules1<-pustules[which(pustules$tag==tag),]
    
    for (color in unique(sub.pustules1$color))
    {
      sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
      
      for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
      {
        sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
        
        for(pustule.number in unique(sub.pustules3$pustule.number))
        {
          sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
          sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
          if(dim(sub.pustules4)[1]>=2)
          {
            if(length(unique(as.Date(sub.pustules4$date)))<dim(sub.pustules4)[1]) #get rid of duplicate data resulting from two pics of same pustule on same date
            {
              for(date in unique(as.Date(sub.pustules4$date)))
              {
                date.indicies<-which(as.Date(sub.pustules4$date)==date)
                if(length(date.indicies)>1) {sub.pustules4<-sub.pustules4[-date.indicies[2:length(date.indicies)],]}
              }
            }
            for(i in 1:(dim(sub.pustules4)[1]-1))
            {
              #pull reference data
              date0<-sub.pustules4[i,"date"]
              date1<-sub.pustules4[i+1,"date"]
              site<-sub.pustules4[i,"site"]
              
              #subset temp data to relevant window
              temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
              temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #pull out relevant data
              temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #pull out relevant data
              temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #throw out NAs
              temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
              
              #subset weather data to relevant window
              weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
              weath.sub<-subset(weath.sub,date<=date1) #pull out relevant data
              weath.sub<-subset(weath.sub,date>=date0) #pull out relevant data
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
              start.val<-sub.pustules4[i,"area"]
              end.val<-sub.pustules4[i+1,"area"]
              delta.days<-date1-date0
              measurer.id<-sub.pustules4[i,"who.entered"] #same person did all measurements for each pustule
              
              #store values
              tags<-c(tags,tag)
              sites<-c(sites,site)
              stem.iters<-c(stem.iters,color)
              leaf.iters<-c(leaf.iters,leaf.iteration)
              pustule.nums<-c(pustule.nums,pustule.number)
              
              start.vals<-c(start.vals,start.val)
              end.vals<-c(end.vals,end.val)
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
              
              measurer.ids<-c(measurer.ids,measurer.id)
              
            } 
          }
        }
      }
    }
  }
  
  delta.pustules<-data.frame(tag=factor(tags),site=factor(sites),stem.iter=stem.iters,leaf.iter=leaf.iters,pustule.num=pustule.nums,area=start.vals,area.next=end.vals,time=days,
                             temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,
                             dew.point.days=dew.point.days,temp.dew.point.days=temp.dew.point.days,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,
                             wetness.days=wetness.days,temp.wetness.days=temp.wetness.days,temp.16.22.wetness.days=temp.16.22.wetness.days,temp.7.30.wetness.days=temp.7.30.wetness.days,
                             tot.rain=tot.rains,solar.days=solar.days,wind.speed.days=wind.speed.days,gust.speed.days=gust.speed.days,who.measured=measurer.ids)
  
  saveRDS(pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pustules.RDS")
  saveRDS(delta.pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.pustules.RDS")
}

pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pustules.RDS")
delta.pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.pustules.RDS")


