if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS")))
{
  # load model and data from pustule area analysis
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")
  delta.pustules<-subset(delta.pustules,time<=7)
  pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
  
  # prep enviro data
  source("prep.enviro.data.R")
  
  # subset 'pustules' data set to get unqiue number of pustules data
  
  ## pull out relevant columns
  n.pustules<-pustules[,c("who.entered","date","picture","site","tag","color","leaf.iteration","N.pustules","N.pustule.count.confidence")]
  
  ## subset to data of suitable quality
  n.pustules<-subset(n.pustules,N.pustule.count.confidence=="Yes")
  
  ## subset to unique data (duplicates exist because a separate row was entered for each measurment of an individual pustule)
  n.pustules<-unique(n.pustules)
  
  ## make new data object for change in pustule size
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  stem.iters<-c()
  leaf.iters<-c()
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
  pred.pustule.diam.growths<-c()
  
  
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
            start.val<-sub.n.pustules3[i,"N.pustules"]
            end.val<-sub.n.pustules3[i+1,"N.pustules"]
            delta.days<-as.numeric(date1-date0)
            measurer.id<-sub.n.pustules3[i,"who.entered"] #same person did all measurements for each pustule
            
            #predict pustule growth from pustule growth model and enviro conditions
            #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
            pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
            obs.time<-delta.days
            pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days/delta.days,"tot.rain"=new.tot.rain/delta.days,"gust.speed.days"=new.gust.speed.days/delta.days)
            pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
            
            #store values
            tags<-c(tags,tag)
            stem.iters<-c(stem.iters,color)
            leaf.iters<-c(leaf.iters,leaf.iteration)
  
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
            
            pred.pustule.diam.growths<-c(pred.pustule.diam.growths,pred.pustule.diam.growth)
          } 
        }
      }
    }
  }
  
  delta.n.pustules<-data.frame(tag=factor(tags),stem.iter=stem.iters,leaf.iter=leaf.iters,n.pustules=start.vals,n.pustules.next=end.vals,time=days,
                             temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,
                             dew.point.days=dew.point.days,temp.dew.point.days=temp.dew.point.days,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,
                             wetness.days=wetness.days,temp.wetness.days=temp.wetness.days,temp.16.22.wetness.days=temp.16.22.wetness.days,temp.7.30.wetness.days=temp.7.30.wetness.days,
                             tot.rain=tot.rains,solar.days=solar.days,wind.speed.days=wind.speed.days,gust.speed.days=gust.speed.days,who.measured=measurer.ids,pred.pustule.diam.growth=pred.pustule.diam.growths)
  
  saveRDS(n.pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS")
  saveRDS(delta.n.pustules,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS")
}

n.pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/n.pustules.RDS")
delta.n.pustules<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.n.pustules.RDS")


