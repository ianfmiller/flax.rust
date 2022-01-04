if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS")))
{
  
  # load data
  pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv")
  
  # pull out relevant columns
  n.pustules<-pustules[,c("who.entered","date","picture","site","tag","color","leaf.iteration","N.pustules","N.pustule.count.confidence")]
  
  # format time
  n.pustules$date<-as.POSIXct(n.pustules$date,tryFormats = "%m/%d/%y %H:%M",tz="UTC")
  
  # cut out incomplete records
  n.pustules<-n.pustules[-intersect(which(n.pustules$site=="HM"),which(n.pustules$date>"2020-07-10 00:00:00 UTC")),]
  
  # subset to data of suitable quality
  n.pustules<-subset(n.pustules,N.pustule.count.confidence=="Yes")
  
  # remove NA values
  n.pustules<-n.pustules[-which(is.na(n.pustules$N.pustules)),]
  
  # subset to unique data (duplicates exist because a separate row was entered for each measurment of an individual pustule)
  n.pustules<-unique(n.pustules)
  
  
  # prep enviro data
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  
  # make new data object for change in pustule size
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  sites<-c()
  stem.iters<-c()
  leaf.iters<-c()
  start.vals<-c()
  end.vals<-c()
  days<-c()
  mean.temp<-c()
  max.temp<-c()
  min.temp<-c()
  mean.abs.hum<-c()
  max.abs.hum<-c()
  min.abs.hum<-c()
  mean.daily.rain<-c()
  mean.solar<-c()
  mean.wetness<-c()
  mean.windspeed<-c()
  mean.soil.moisture<-c()

  
  for (tag in unique(n.pustules$tag))
  {
    sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
    
    for (color in unique(sub.n.pustules1$color))
    {
      sub.n.pustules2<-sub.n.pustules1[which(sub.n.pustules1$color==color),]
      
      for(leaf.iteration in unique(sub.n.pustules2$leaf.iteration)) 
      {
        sub.n.pustules3<-sub.n.pustules2[which(sub.n.pustules2$leaf.iteration==leaf.iteration),]
        
        sub.n.pustules3<-sub.n.pustules3[order(sub.n.pustules3$date),]
        if(dim(sub.n.pustules3)[1]>=2)
        {
          if(length(unique(as.Date(sub.n.pustules3$date)))<dim(sub.n.pustules3)[1]) #get rid of duplicate data resulting from two pics of same pustule on same date
          {
            for(date in unique(as.Date(sub.n.pustules3$date)))
            {
              date.indicies<-which(as.Date(sub.n.pustules3$date)==date)
              if(length(date.indicies)>1) {sub.n.pustules3<-sub.n.pustules3[-date.indicies[2:length(date.indicies)],]}
            }
          }
          
          for(i in 1:(dim(sub.n.pustules3)[1]-1))
          {
            #pull reference data
            date0<-sub.n.pustules3[i,"date"]
            date1<-sub.n.pustules3[i+1,"date"]
            site<-sub.n.pustules3[i,"site"]
            
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
            
            abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
            new.mean.abs.hum<-mean(abs.hum,na.rm=T)
            new.max.abs.hum<-max(abs.hum,na.rm=T)
            new.min.abs.hum<-min(abs.hum,na.rm=T)
            
            new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(12*24)
            new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
            new.mean.wetness<-mean(weath.sub$wetness,na.rm=T)
            new.mean.windspeed<-mean(weath.sub$wind.speed,na.rm=T)
            new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
            
            #pull out core predictors
            start.val<-sub.n.pustules3[i,"N.pustules"]
            end.val<-sub.n.pustules3[i+1,"N.pustules"]
            delta.days<-as.numeric(date1-date0)
            measurer.id<-sub.n.pustules3[i,"who.entered"] #same person did all measurements for each pustule
            
            #store values
            tags<-c(tags,tag)
            sites<-c(sites,site)
            stem.iters<-c(stem.iters,color)
            leaf.iters<-c(leaf.iters,leaf.iteration)
  
            start.vals<-c(start.vals,start.val)
            end.vals<-c(end.vals,end.val)
            days<-c(days,delta.days)
            
            mean.temp<-c(mean.temp,new.mean.temp)
            max.temp<-c(max.temp,new.max.temp)
            min.temp<-c(min.temp,new.min.temp)
            mean.abs.hum<-c(mean.abs.hum,new.mean.abs.hum)
            max.abs.hum<-c(max.abs.hum,new.max.abs.hum)
            min.abs.hum<-c(min.abs.hum,new.min.abs.hum)

            mean.daily.rain<-c(mean.daily.rain,new.mean.daily.rain)
            mean.solar<-c(mean.solar,new.mean.solar)
            mean.wetness<-c(mean.wetness,new.mean.wetness)
            mean.windspeed<-c(mean.windspeed,new.mean.windspeed)
            mean.soil.moisture<-c(mean.soil.moisture,new.mean.soil.moisture)
            
          } 
        }
      }
    }
  }
  
  delta.n.pustules<-data.frame(tag=factor(tags),site=factor(sites),stem.iter=stem.iters,leaf.iter=leaf.iters,n.pustules=start.vals,n.pustules.next=end.vals,time=days,                             
                             mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,
                             mean.abs.hum=mean.abs.hum,max.abs.hum=max.abs.hum,min.abs.hum=min.abs.hum,
                             mean.daily.rain=mean.daily.rain,mean.solar=mean.solar,mean.wetness=mean.wetness,mean.windspeed=mean.windspeed,mean.soil.moisture=mean.soil.moisture)
  
  saveRDS(n.pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS")
  saveRDS(delta.n.pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS")
}

n.pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS")
delta.n.pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS")


