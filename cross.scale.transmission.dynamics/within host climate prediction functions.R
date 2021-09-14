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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"time"=delta.days,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"max.abs.hum"=new.max.abs.hum,"mean.solar"=new.mean.solar,"site"=site)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"time"=delta.days,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"tot.rain"=new.tot.rain,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"time"=delta.days,"min.vpd"=new.max.vpd,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  

  out.data<-data.frame(area=dummy.data,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                             mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                             mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),mean.vpd=rep(new.mean.vpd,times=dim),
                             max.vpd=rep(new.max.vpd,times=dim),min.vpd=rep(new.min.vpd,times=dim),mean.wetness=rep(new.mean.wetness,times=dim),tot.rain=rep(new.tot.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim),
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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"time"=delta.days,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"max.abs.hum"=new.max.abs.hum,"mean.solar"=new.mean.solar,"site"=site)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"time"=delta.days,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"tot.rain"=new.tot.rain,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"time"=delta.days,"min.vpd"=new.max.vpd,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  
  
  out.data<-data.frame(n.pustules=dummy.data,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),mean.vpd=rep(new.mean.vpd,times=dim),
                       max.vpd=rep(new.max.vpd,times=dim),min.vpd=rep(new.min.vpd,times=dim),mean.wetness=rep(new.mean.wetness,times=dim),tot.rain=rep(new.tot.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim),
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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"time"=delta.days,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"max.abs.hum"=new.max.abs.hum,"mean.solar"=new.mean.solar,"site"=site)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"time"=delta.days,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"tot.rain"=new.tot.rain,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"time"=delta.days,"min.vpd"=new.max.vpd,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  
  out.data<-data.frame(plant.inf.intens=dummy.data,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),mean.vpd=rep(new.mean.vpd,times=dim),
                       max.vpd=rep(new.max.vpd,times=dim),min.vpd=rep(new.min.vpd,times=dim),mean.wetness=rep(new.mean.wetness,times=dim),tot.rain=rep(new.tot.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim),
                       pred.pustule.diam.growth=rep(pred.pustule.diam.growth,times=dim),pred.pustule.num.increase=rep(pred.pustule.num.increase,times=dim),pred.plant.inf.intens.increase=rep(pred.plant.inf.intens.increase,times=dim))
  out.data
}

get.pred.data<-function(site,date0,date1,dummy.data,temp.addition=0)
{
  dim<-length(dummy.data)
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
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
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"time"=delta.days,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"max.abs.hum"=new.max.abs.hum,"mean.solar"=new.mean.solar,"site"=site)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"time"=delta.days,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"tot.rain"=new.tot.rain,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  #predict change in plant.inf.intensity from enviro conditions
  #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
  plant.model.new.plant.inf.intens<-.1
  obs.time<-delta.days
  plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"time"=delta.days,"min.vpd"=new.max.vpd,"site"=site)
  pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
  
  
  out.data<-data.frame(plant.inf.intens=dummy.data,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),mean.vpd=rep(new.mean.vpd,times=dim),
                       max.vpd=rep(new.max.vpd,times=dim),min.vpd=rep(new.min.vpd,times=dim),mean.wetness=rep(new.mean.wetness,times=dim),tot.rain=rep(new.tot.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim),
                       pred.pustule.diam.growth=rep(pred.pustule.diam.growth,times=dim),pred.pustule.num.increase=rep(pred.pustule.num.increase,times=dim),pred.plant.inf.intens.increase=rep(pred.plant.inf.intens.increase,times=dim))
  out.data
  
}
