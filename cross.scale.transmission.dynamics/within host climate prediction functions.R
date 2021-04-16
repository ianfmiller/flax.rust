source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
plant.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")

sites<-c("CC","BT","GM","HM")

pull.temp.mean.quantile.data<-function()
{
  daily.means<-data.frame("site"=character(),"date"=character(),"mean.temp"=numeric())
  for(site in sites)
  {
    sub.temp<-all.temp.rh[which(all.temp.rh$site==site),]
    dates<-unique(as.Date(sub.temp$date.time))[-c(1,length(unique(as.Date(sub.temp$date.time))))] # cut out incomplete days
    for (date in dates)
    {
      sub.sub.temp<-sub.temp[which(as.Date(sub.temp$date.time)==date),]
      daily.means<-rbind(daily.means,data.frame("site"=site,"date"=as.Date(date,origin = "1970-01-01"),"mean.temp"=mean(sub.sub.temp$temp.c,na.rm=T)))
    }
  }
  daily.means<-daily.means[order(daily.means$mean.temp),]
  daily.means
}

get.pred.data.temp.mean.quantile.pustule.model<-function(day.set,dummy.data,temp.addition=0)
{
  dim<-length(dummy.data)
  daily.means<-pull.temp.mean.quantile.data()
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

  date<-daily.means[day.set,"date"]
  date0<-as.POSIXct(strptime(paste0(date," 00:00:01"),"%Y-%m-%d %H:%M:%S"),tz="UTC")
  date1<-date0+60*60*24
  site<-daily.means[day.set,"site"]
  
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
  
  #pull out core predictors
  delta.days<-as.numeric(date1-date0)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  
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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.days.16.22"=new.temp.days.16.22/delta.days,"dew.point.days"=new.dew.point.days/delta.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days/delta.days,"temp.wetness.days"=new.temp.wetness.days/delta.days,"tot.rain"=new.tot.rain/delta.days)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"dew.point.days"=new.dew.point.days/delta.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days/delta.days,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  

  out.data<-data.frame(area=dummy.data,time=rep(delta.days,times=dim),
                             temp.days=rep(new.temp.days,times=dim),temp.days.16.22=rep(new.temp.days.16.22,times=dim),temp.days.7.30=rep(new.temp.days.7.30,times=dim),
                             dew.point.days=rep(new.dew.point.days,times=dim),temp.dew.point.days=rep(new.temp.dew.point.days,times=dim),temp.16.22.dew.point.days=rep(new.temp.16.22.dew.point.days,times=dim),temp.7.30.dew.point.days=rep(new.temp.7.30.dew.point.days,times=dim),
                             wetness.days=rep(new.wetness.days,times=dim),temp.wetness.days=rep(new.temp.wetness.days,times=dim),temp.16.22.wetness.days=rep(new.temp.16.22.wetness.days,times=dim),temp.7.30.wetness.days=rep(new.temp.7.30.wetness.days,times=dim),
                             tot.rain=rep(new.tot.rain,times=dim),solar.days=rep(new.solar.days,times=dim),wind.speed.days=rep(new.wind.speed.days,times=dim),gust.speed.days=rep(new.gust.speed.days),
                       pred.pustule.diam.growth=rep(pred.pustule.diam.growth,times=dim),pred.pustule.num.increase=rep(pred.pustule.num.increase,times=dim),pred.plant.inf.intens.increase=rep(pred.plant.inf.intens.increase,times=dim))
  out.data
}

get.pred.data.temp.mean.quantile.n.pustules.model<-function(day.set,dummy.data,temp.addition=0)
{
  dim<-length(dummy.data)
  daily.means<-pull.temp.mean.quantile.data()
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  date<-daily.means[day.set,"date"]
  date0<-as.POSIXct(strptime(paste0(date," 00:00:01"),"%Y-%m-%d %H:%M:%S"),tz="UTC")
  date1<-date0+60*60*24
  site<-daily.means[day.set,"site"]
  
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
  
  #pull out core predictors
  delta.days<-as.numeric(date1-date0)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  
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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.days.16.22"=new.temp.days.16.22/delta.days,"dew.point.days"=new.dew.point.days/delta.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days/delta.days,"temp.wetness.days"=new.temp.wetness.days/delta.days,"tot.rain"=new.tot.rain/delta.days)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"dew.point.days"=new.dew.point.days/delta.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days/delta.days,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  
  
  out.data<-data.frame(n.pustules=dummy.data,time=rep(delta.days,times=dim),
                       temp.days=rep(new.temp.days,times=dim),temp.days.16.22=rep(new.temp.days.16.22,times=dim),temp.days.7.30=rep(new.temp.days.7.30,times=dim),
                       dew.point.days=rep(new.dew.point.days,times=dim),temp.dew.point.days=rep(new.temp.dew.point.days,times=dim),temp.16.22.dew.point.days=rep(new.temp.16.22.dew.point.days,times=dim),temp.7.30.dew.point.days=rep(new.temp.7.30.dew.point.days,times=dim),
                       wetness.days=rep(new.wetness.days,times=dim),temp.wetness.days=rep(new.temp.wetness.days,times=dim),temp.16.22.wetness.days=rep(new.temp.16.22.wetness.days,times=dim),temp.7.30.wetness.days=rep(new.temp.7.30.wetness.days,times=dim),
                       tot.rain=rep(new.tot.rain,times=dim),solar.days=rep(new.solar.days,times=dim),wind.speed.days=rep(new.wind.speed.days,times=dim),gust.speed.days=rep(new.gust.speed.days),
                       pred.pustule.diam.growth=rep(pred.pustule.diam.growth,times=dim),pred.pustule.num.increase=rep(pred.pustule.num.increase,times=dim),pred.plant.inf.intens.increase=rep(pred.plant.inf.intens.increase,times=dim))
  out.data
}

get.pred.data.temp.mean.quantile.plants.model<-function(day.set,dummy.data,temp.addition=0)
{
  dim<-length(dummy.data)
  daily.means<-pull.temp.mean.quantile.data()
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  date<-daily.means[day.set,"date"]
  date0<-as.POSIXct(strptime(paste0(date," 00:00:01"),"%Y-%m-%d %H:%M:%S"),tz="UTC")
  date1<-date0+60*60*24
  site<-daily.means[day.set,"site"]
  
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
  
  #pull out core predictors
  delta.days<-as.numeric(date1-date0)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  
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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.days.16.22"=new.temp.days.16.22/delta.days,"dew.point.days"=new.dew.point.days/delta.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days/delta.days,"temp.wetness.days"=new.temp.wetness.days/delta.days,"tot.rain"=new.tot.rain/delta.days)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"dew.point.days"=new.dew.point.days/delta.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days/delta.days,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  
  
  out.data<-data.frame(plant.inf.intens=dummy.data,time=rep(delta.days,times=dim),
                       temp.days=rep(new.temp.days,times=dim),temp.days.16.22=rep(new.temp.days.16.22,times=dim),temp.days.7.30=rep(new.temp.days.7.30,times=dim),
                       dew.point.days=rep(new.dew.point.days,times=dim),temp.dew.point.days=rep(new.temp.dew.point.days,times=dim),temp.16.22.dew.point.days=rep(new.temp.16.22.dew.point.days,times=dim),temp.7.30.dew.point.days=rep(new.temp.7.30.dew.point.days,times=dim),
                       wetness.days=rep(new.wetness.days,times=dim),temp.wetness.days=rep(new.temp.wetness.days,times=dim),temp.16.22.wetness.days=rep(new.temp.16.22.wetness.days,times=dim),temp.7.30.wetness.days=rep(new.temp.7.30.wetness.days,times=dim),
                       tot.rain=rep(new.tot.rain,times=dim),solar.days=rep(new.solar.days,times=dim),wind.speed.days=rep(new.wind.speed.days,times=dim),gust.speed.days=rep(new.gust.speed.days),
                       pred.pustule.diam.growth=rep(pred.pustule.diam.growth,times=dim),pred.pustule.num.increase=rep(pred.pustule.num.increase,times=dim),pred.plant.inf.intens.increase=rep(pred.plant.inf.intens.increase,times=dim),site=rep(site,times=dim))
  out.data
}
