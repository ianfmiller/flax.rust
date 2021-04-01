# load data and model
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
plant.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")
n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")


# function for subsetting temp/rh
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

#function for predicting plant inf intens change
predict.plant.inf.intens<-function(plant.inf.intens.last,site,date0,date1)
{
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  ## subset weather data to relevant window
  weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
  weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
  weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
  weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
  
  # calculate environmental variable metrics
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
  
  delta.days<-as.numeric(date1-date0)
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
  
  # make forward prediction
  pred.data<-data.frame("plant.inf.intens"=plant.inf.intens.last,"dew.point.days"=new.dew.point.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days,"temp.wetness.days"=new.temp.wetness.days,"tot.rain"=new.tot.rain,"pred.pustule.num.increase"=pred.pustule.num.increase)
  change<-predict(plant.model,newdata=pred.data,re.form=NA)
  plant.inf.intens.next<-plant.inf.intens.last+change
  if(plant.inf.intens.next<.01) {plant.inf.intens.next<-.01}
  plant.inf.intens.next
}

#function for predicting plant inf intens change
predict.plant.inf.intens.last<-function(plant.inf.intens.next,site,date0,date1)
{
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  ## subset weather data to relevant window
  weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
  weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
  weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
  weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
  
  # calculate environmental variable metrics
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
  
  delta.days<-as.numeric(date1-date0)
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
  
  # make backwards prediction
  pred.func<-function(x)
  {
    plant.inf.intens.last.test<-x
    pred.data<-data.frame("plant.inf.intens"=plant.inf.intens.last.test,"dew.point.days"=new.dew.point.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days,"temp.wetness.days"=new.temp.wetness.days,"tot.rain"=new.tot.rain,"pred.pustule.num.increase"=pred.pustule.num.increase)
    change<-predict(plant.model,newdata=pred.data,re.form=NA)
    plant.inf.intens.next.pred<-plant.inf.intens.last.test+change
    abs(plant.inf.intens.next.pred-plant.inf.intens.next)
  }
  plant.inf.intens.last<-optim(c(plant.inf.intens.next),pred.func,method = "Brent",lower=0,upper=10e6)$par
  plant.inf.intens.last
}




