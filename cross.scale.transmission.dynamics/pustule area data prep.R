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
  mean.temp<-c()
  max.temp<-c()
  min.temp<-c()
  mean.dew.point<-c()
  max.dew.point<-c()
  min.dew.point<-c()
  mean.vpd<-c() #vapor pressure deficit
  mean.wetness<-c()
  tot.rain<-c()
  mean.solar<-c()
  
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
              new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
              new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
              new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
              
              new.mean.dew.point<-mean(temp.rh.sub$dew.pt.c,na.rm = T) #mean dew point
              new.max.dew.point<-max(temp.rh.sub$dew.pt.c,na.rm = T) #mean dew point
              new.min.dew.point<-min(temp.rh.sub$dew.pt.c,na.rm = T) #mean dew point
              
              svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
              avps<- temp.rh.sub$rh.c / 100 * svps #actual vapor pressures 
              vpds<-avps-svps
              new.mean.vpd<-mean(vpds,na.rm=T)
              
              #calculate weather metrics
              new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
              new.tot.rain<-sum(weath.sub$rain,na.rm=T)
              new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
              
              #pull out core predictors
              start.val<-sub.pustules4[i,"area"]
              end.val<-sub.pustules4[i+1,"area"]
              delta.days<-date1-date0

              #store values
              tags<-c(tags,tag)
              sites<-c(sites,site)
              stem.iters<-c(stem.iters,color)
              leaf.iters<-c(leaf.iters,leaf.iteration)
              pustule.nums<-c(pustule.nums,pustule.number)
              
              start.vals<-c(start.vals,start.val)
              end.vals<-c(end.vals,end.val)
              days<-c(days,delta.days)
              
              mean.temp<-c(mean.temp,new.mean.temp)
              max.temp<-c(max.temp,new.max.temp)
              min.temp<-c(min.temp,new.min.temp)
              mean.dew.point<-c(mean.dew.point,new.mean.dew.point)
              max.dew.point<-c(max.dew.point,new.max.dew.point)
              min.dew.point<-c(min.dew.point,new.min.dew.point)
              mean.vpd<-c(mean.vpd,new.mean.vpd)
              mean.wetness<-c(mean.wetness,new.mean.wetness)
              tot.rain<-c(tot.rain,new.tot.rain)
              mean.solar<-c(mean.solar,new.mean.solar)
              
            } 
          }
        }
      }
    }
  }
  
  delta.pustules<-data.frame(tag=factor(tags),site=factor(sites),stem.iter=stem.iters,leaf.iter=leaf.iters,pustule.num=pustule.nums,area=start.vals,area.next=end.vals,time=days,
                             mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,mean.dew.point=mean.dew.point,max.dew.point=max.dew.point,min.dew.point=min.dew.point,mean.vpd=mean.vpd,
                             mean.wetness=mean.wetness,tot.rain=tot.rain,mean.solar=mean.solar)
  
  saveRDS(pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pustules.RDS")
  saveRDS(delta.pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.pustules.RDS")
}

pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pustules.RDS")
delta.pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.pustules.RDS")


