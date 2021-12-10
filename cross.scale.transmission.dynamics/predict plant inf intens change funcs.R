# load data and model
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
infection.intensity.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS")


# function for subsetting temp/rh
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

#function for predicting plant inf intens change
predict.inf.intens<-function(inf.intens.last,max.height.last,site,date0,date1)
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
  
  #calculate environmental variable metrics
  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T)
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("infection.intensity"=inf.intens.last,"max.height"=max.height.last,
                        "time"=delta.days,"site"=site,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.wetness"=new.mean.wetness,"mean.daily.rain"=new.mean.daily.rain,"mean.solar"=new.mean.solar,"mean.soil.moisture"=new.mean.soil.moisture)
  
  inf.intens.next<-inf.intens.last+delta.days*predict(infection.intensity.model,newdata = pred.data,type="response",exclude = c('s(site)','s(tag)'))
  if(inf.intens.next<.1) {inf.intens.next<-0}
  inf.intens.next
}

predict.inf.intens.boot<-function(inf.intens.last,max.height.last,site,date0,date1,temp.addition=0)
{
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition ### shift temperature
  
  ## subset weather data to relevant window
  weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
  weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
  weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
  weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
  
  # calculate environmental variable metrics
  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T)
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("infection.intensity"=inf.intens.last,"max.height"=max.height.last,
                        "time"=delta.days,"site"=site,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.wetness"=new.mean.wetness,"mean.daily.rain"=new.mean.daily.rain,"mean.solar"=new.mean.solar,"mean.soil.moisture"=new.mean.soil.moisture)
  
  Xp <- predict(infection.intensity.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
  beta <- coef(infection.intensity.model) ## posterior mean of coefs
  Vb   <- vcov(infection.intensity.model) ## posterior  cov of coefs
  n <- 2
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  preds <- rep(NA, n)
  ilink <- family(infection.intensity.model)$linkinv
  for (j in seq_len(n)) { 
    preds[j]   <- ilink(Xp %*% mrand[j, ])
  }
  inf.intens.next<-inf.intens.last+delta.days*preds[1]
  if(inf.intens.next<.1) {inf.intens.next<-0}
  inf.intens.next
}

#function for predicting plant inf intens change
predict.inf.intens.last<-function(inf.intens.next,max.height.last,site,date0,date1)
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
  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T)
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
  new.mean.daily.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  

  # make backwards prediction
  pred.func<-function(x)
  {
    inf.intens.last.test<-x
    pred.data<-data.frame("infection.intensity"=inf.intens.last.test,"max.height"=max.height.last,
                          "time"=delta.days,"site"=site,
                          "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                          "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                          "mean.wetness"=new.mean.wetness,"mean.daily.rain"=new.mean.daily.rain,"mean.solar"=new.mean.solar,"mean.soil.moisture"=new.mean.soil.moisture)
    inf.intens.next.pred<-inf.intens.last.test+delta.days*predict(infection.intensity.model,newdata=pred.data,exclude = c('s(site)','s(tag)'))
    abs(inf.intens.next.pred-inf.intens.next)
  }
  inf.intens.last<-optim(c(inf.intens.next),pred.func,method = "Brent",lower=0,upper=10e6)$par
  if(plant.inf.intens.last<0.1) {inf.intens.last<-.1}
  inf.intens.last
}




