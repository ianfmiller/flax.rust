source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")

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
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T)
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  

  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(60*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)

  out.data<-data.frame(area=dummy.data,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                             mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                             mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),
                             mean.wetness=rep(new.mean.wetness,times=dim),mean.daily.rain=rep(new.mean.daily.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim))
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
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  

  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(60*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  
  out.data<-data.frame(n.pustules=dummy.data,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),
                       mean.wetness=rep(new.mean.wetness,times=dim),mean.daily.rain=rep(new.mean.daily.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim))
  out.data
}

get.pred.data.temp.mean.quantile.infection.intensity.model<-function(day.set,dummy.data.infection.intensity,dummy.data.height,temp.addition=0)
{
  dim<-length(dummy.data.inf.intens)
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
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T) 
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  

  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(60*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  
  
  out.data<-data.frame(infection.intensity=dummy.data.inf.intens,max.height=dummy.data.height,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),
                       mean.wetness=rep(new.mean.wetness,times=dim),mean.daily.rain=rep(new.mean.daily.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim))
  out.data
}

get.pred.data.temp.mean.quantile.plant.growth.model<-function(day.set,dummy.data.height,temp.addition=0)
{
  dim<-length(dummy.data.height)
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
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T)
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  

  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(60*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  
  out.data<-data.frame(height=dummy.data.height,time=rep(delta.days,times=dim),site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),
                       #mean.vpd=rep(new.mean.vpd,times=dim),max.vpd=rep(new.max.vpd,times=dim),min.vpd=rep(new.min.vpd,times=dim),
                       mean.wetness=rep(new.mean.wetness,times=dim),mean.daily.rain=rep(new.mean.daily.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim),mean.soil.moisture=new.mean.soil.moisture)
  out.data
}

get.pred.data<-function(site,date0,date1,dummy.data,dummy.data.max.height=15,temp.addition=0)
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
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T)
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #adjust temp data
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition
  

  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(60*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  

  #area, n.pustules, and plant.inf.intens are all set to dummy data so that this function can be used interchangably for predictions at different scales without specifying the scale.
  out.data<-data.frame(area=dummy.data,n.pustules=dummy.data,plant.inf.intens=dummy.data,time=rep(delta.days,times=dim),max.height=dummy.data.max.height,site=rep(site,times=dim),
                       mean.temp=rep(new.mean.temp,times=dim),max.temp=rep(new.max.temp,times=dim),min.temp=rep(new.min.temp,times=dim),
                       mean.abs.hum=rep(new.mean.abs.hum,times=dim),max.abs.hum=rep(new.max.abs.hum,times=dim),min.abs.hum=rep(new.min.abs.hum,times=dim),
                       mean.wetness=rep(new.mean.wetness,times=dim),mean.daily.rain=rep(new.mean.daily.rain,times=dim),mean.solar=rep(new.mean.solar,times=dim),mean.soil.moisture=rep(new.mean.soil.moisture,times=dim))
  out.data
  
}
