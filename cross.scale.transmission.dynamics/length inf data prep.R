if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.length.inf.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/length.inf.RDS")))
{
  # load model and data from pustule area analysis
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")
  delta.pustules<-subset(delta.pustules,time<=7)
  pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
  
  # load model and data from n pustules analysis
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/n pustules data prep.R")
  delta.n.pustules<-subset(delta.pustules,time<=7)
  n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
  
  
  # prep enviro data
  source("prep.enviro.data.R")
  
  # prep data
  
  ## load data
  within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
  
  ## pull out relevant columns
  length.inf<-within.host[,c("Year","Site","Tag","Date","picture","stem.index","stem.height","percent.tissue.infected","length.tissue.infected")]
  
  ## subset to unique data (duplicates exist because a separate row was entered for each measurment of an individual pustule)
  length.inf<-length.inf[-which(length.inf$stem.index %in% c("","other")),]
  length.inf<-unique(length.inf)
  
  ## subset to 2020
  length.inf<-length.inf[which(length.inf$Year=="2020"),]
  
  ## finish cleaning
  length.inf$stem.height<-as.numeric(length.inf$stem.height)
  length.inf$length.tissue.infected<-as.numeric(length.inf$length.tissue.infected)
  ## make new data object for change in pustule size
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  stem.iters<-c()
  start.lengths<-c()
  end.lengths<-c()
  start.length.infs<-c()
  end.length.infs<-c()
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
  pred.pustule.diam.growths<-c()
  pred.pustule.num.increases<-c()
  
  
  for (tag in unique(length.inf$Tag))
  {
    sub.length.inf.1<-length.inf[which(length.inf$Tag==tag),]
    
    for (stem.index in unique(sub.length.inf.1$stem.index))
    {
      sub.length.inf.2<-sub.length.inf.1[which(sub.length.inf.1$stem.index==stem.index),]
      sub.length.inf.2<-sub.length.inf.2[-which(is.na(sub.length.inf.2$length.tissue.infected)),]
      
      if(dim(sub.length.inf.2)[1]>=2)
      {
        ### get rid of duplicate data resulting from two pics of same pustule on same date
        if(length(unique(as.Date(sub.length.inf.2$Date,tryFormats = c("%m/%d/%Y"))))<dim(sub.length.inf.2)[1])
        {
          for(date in unique(as.Date(sub.length.inf.2$Date,tryFormats = c("%m/%d/%Y"))))
          {
            date.indicies<-which(as.Date(sub.length.inf.2$Date,tryFormats = c("%m/%d/%Y"))==date)
            if(length(date.indicies)>1) {sub.length.inf.2<-sub.length.inf.2[-date.indicies[2:length(date.indicies)],]}
          }
        }
        
        for(i in 1:(dim(sub.length.inf.2)[1]-1))
        {
          ### pull reference data
          #### use pictures to get date times for plant, average
          date0.pic<-sub.length.inf.2[i,"picture"]
          date0<-pustules[which(pustules$picture==date0.pic)[1],"date"]
          date1.pic<-sub.length.inf.2[i+1,"picture"]
          date1<-pustules[which(pustules$picture==date1.pic)[1],"date"]
          #### if unable to line up pic with date time, use mean date time of all measurments at site on given day
          if(is.na(date0))
          {
            simple.date<-sub.length.inf.2[i,"Date"]
            sub.pustules<-pustules[which(pustules$site==site),]
            sub.pustules<-sub.pustules[which(as.Date(sub.pustules$date)==as.Date(simple.date,tryFormats = "%m/%d/%y")),]
            alt.date.0.pics<-unique(sub.pustules$picture)
            alt.dates<-pustules[which(pustules$picture %in% alt.date.0.pics),"date"]
            date0<-mean(unique(alt.dates))
          }
          if(is.na(date1))
          {
            simple.date<-sub.length.inf.2[i+1,"Date"]
            sub.pustules<-pustules[which(pustules$site==site),]
            sub.pustules<-sub.pustules[which(as.Date(sub.pustules$date)==as.Date(simple.date,tryFormats = "%m/%d/%y")),]
            alt.date.1.pics<-unique(sub.pustules$picture)
            alt.dates<-pustules[which(pustules$picture %in% alt.date.1.pics),"date"]
            date1<-mean(unique(alt.dates))
          }
          site<-sub.length.inf.2[i,"Site"]
          
          ### subset temp data to relevant window
          temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
          temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
          temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
          temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
          temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
          temp.rh.sub
          #### subset weather data to relevant window
          weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
          weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
          weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
          weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
          
          #### calculate environmental variable metrics
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
          start.length.inf<-sub.length.inf.2[i,"length.tissue.infected"]
          end.length.inf<-sub.length.inf.2[i+1,"length.tissue.infected"]
          start.length<-sub.length.inf.2[i,"stem.height"]
          end.length<-sub.length.inf.2[i+1,"stem.height"]
          delta.days<-as.numeric(date1-date0)
          
          #predict pustule growth from pustule growth model and enviro conditions
          #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
          pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
          obs.time<-delta.days
          pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days/delta.days)
          pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
          
          #predict change in number of pustules from enviro conditions
          #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
          n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
          obs.time<-delta.days
          n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days)
          pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
          
          #store values
          tags<-c(tags,tag)
          stem.iters<-c(stem.iters,stem.index)
          
          start.length.infs<-c(start.length.infs,start.length.inf)
          end.length.infs<-c(end.length.infs,end.length.inf)
          start.lengths<-c(start.lengths,start.length)
          end.lengths<-c(end.lengths,end.length)
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
          
          pred.pustule.diam.growths<-c(pred.pustule.diam.growths,pred.pustule.diam.growth)
          pred.pustule.num.increases<-c(pred.pustule.num.increases,pred.pustule.num.increase)
        } 
      }
    }    
  }
  
  delta.length.inf<-data.frame(tag=factor(tags),stem.iter=stem.iters,time=days,start.length=as.numeric(start.lengths),end.length=as.numeric(end.lengths),start.length.inf=as.numeric(start.length.infs),end.length.inf=as.numeric(end.length.infs),
                               temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,
                               dew.point.days=dew.point.days,temp.dew.point.days=temp.dew.point.days,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,
                               wetness.days=wetness.days,temp.wetness.days=temp.wetness.days,temp.16.22.wetness.days=temp.16.22.wetness.days,temp.7.30.wetness.days=temp.7.30.wetness.days,
                               tot.rain=tot.rains,solar.days=solar.days,wind.speed.days=wind.speed.days,gust.speed.days=gust.speed.days,pred.pustule.diam.growth=pred.pustule.diam.growths)
  
  saveRDS(length.inf,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/length.inf.RDS")
  saveRDS(delta.length.inf,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.length.inf.RDS")
}

length.inf<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/length.inf.RDS")
delta.length.inf<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.length.inf.RDS")


