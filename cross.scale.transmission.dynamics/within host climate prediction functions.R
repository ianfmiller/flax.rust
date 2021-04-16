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
  
  #pull out core predictors
  delta.days<-date1-date0

  out.data<-data.frame(area=dummy.data,time=rep(delta.days,times=dim),
                             temp.days=rep(new.temp.days,times=dim),temp.days.16.22=rep(new.temp.days.16.22,times=dim),temp.days.7.30=rep(new.temp.days.7.30,times=dim),
                             dew.point.days=rep(new.dew.point.days,times=dim),temp.dew.point.days=rep(new.temp.dew.point.days,times=dim),temp.16.22.dew.point.days=rep(new.temp.16.22.dew.point.days,times=dim),temp.7.30.dew.point.days=rep(new.temp.7.30.dew.point.days,times=dim),
                             wetness.days=rep(new.wetness.days,times=dim),temp.wetness.days=rep(new.temp.wetness.days,times=dim),temp.16.22.wetness.days=rep(new.temp.16.22.wetness.days,times=dim),temp.7.30.wetness.days=rep(new.temp.7.30.wetness.days,times=dim),
                             tot.rain=rep(new.tot.rain,times=dim),solar.days=rep(new.solar.days,times=dim),wind.speed.days=rep(new.wind.speed.days,times=dim),gust.speed.days=rep(new.gust.speed.days))
  out.data
}
