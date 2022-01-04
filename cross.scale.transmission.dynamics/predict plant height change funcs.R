# load data and model
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")


# function for subsetting temp/rh
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

predict.plant.growth<-function(height.last,inf.intens.last,site,date0,date1,exclude.site=T)
{
  suppressWarnings(if(class(date0)=="Date") {date0<-as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")})
  suppressWarnings(if(class(date1)=="Date") {date1<-as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")})
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  #### subset weather data to relevant window
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
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(12*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)

  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("time"=delta.days,"height"=height.last,"inf.intens"=inf.intens.last,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.wetness"=new.mean.wetness,"mean.daily.rain"=new.mean.daily.rain,"mean.solar"=new.mean.solar,"mean.soil.moisture"=new.mean.soil.moisture,
                        "site"=site,"tag"="NA")
  
  if(exclude.site) {height.next<-height.last+(date1-date0)*predict(plant.growth.model,newdata=pred.data,exclude = c('s(site)','s(tag)'))} else {height.next<-height.last+predict(plant.growth.model,newdata=pred.data,exclude='s(tag)')}
  if(height.next<1) {height.next<-1}
  height.next
}

predict.plant.growth.boot<-function(height.last,inf.intens.last,site,date0,date1)
{
  if("Date" %in% class(date0)) {date0<-as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")}
  if("Date" %in% class(date1)) {date1<-as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")}
  
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  #### subset weather data to relevant window
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
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(12*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("time"=delta.days,"height"=height.last,"inf.intens"=inf.intens.last,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.wetness"=new.mean.wetness,"mean.daily.rain"=new.mean.daily.rain,"mean.solar"=new.mean.solar,"mean.soil.moisture"=new.mean.soil.moisture,
                        "site"=site,"tag"="NA")
  Xp <- predict(plant.growth.model, newdata = pred.data, exlude=c('s(site)','s(tag)'),type="lpmatrix")
  beta <- coef(plant.growth.model) ## posterior mean of coefs
  Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
  n <- 2
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  preds <- rep(NA, n)
  ilink <- family(plant.growth.model)$linkinv
  for (j in seq_len(n)) { 
    preds[j]   <- ilink(Xp %*% mrand[j, ])
  }
  height.next<-height.last+(date1-date0)*preds[1]
  if(height.next<1) {height.next<-1}
  height.next
}

predict.plant.growth.last<-function(height.next,inf.intens.last,site,date0,date1,exclude.site=T)
{
  if(class(date0)=="Date") {date0<-as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")}
  if(class(date1)=="Date") {date1<-as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")}
  
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  #### subset weather data to relevant window
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
  new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)/(12*24)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  new.mean.soil.moisture<-mean(weath.sub$soil.moisture,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  # make backwards prediction
  pred.func<-function(x)
  {
    plant.height.last.test<-x
    pred.data<-data.frame("time"=delta.days,"height"=plant.height.last.test,"inf.intens"=inf.intens.last, ## assume infection intensity did not change--for simplicity, and because tyring to simultaneously back-cast 
                           "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                           "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                          "mean.wetness"=new.mean.wetness,"mean.daily.rain"=new.mean.daily.rain,"mean.solar"=new.mean.solar,"mean.soil.moisture"=new.mean.soil.moisture,
                          "site"=site,"tag"="NA")
    if(exclude.site) {plant.height.next.pred<-plant.height.last.test+(date1-date0)*predict(plant.growth.model,newdata=pred.data,exclude = 's(site)')} else {plant.height.next.pred<-plant.height.last.test+(date1-date0)*predict(plant.growth.model,newdata=pred.data,exclude='s(tag)')}
    abs(plant.height.next.pred-height.next)
  }
  plant.height.last<-optim(c(height.next),pred.func,method = "Brent",lower=0,upper=10e6)$par
  if(plant.height.last<1) {plant.height.last<-1}
  plant.height.last
}

predict.plant.growth.and.infection.intensity.last<-function(height.next,inf.intens.next,site,date0,date1)
{
  pred.func.2<-function(z)
  {
    pred.height.last<-predict.plant.growth.last(height.next,z,site,date0,date1,exclude.site=T)
    pred.inf.intens.last<-predict.inf.intens.last(inf.intens.next,pred.height.last,site,date0,date1)
    abs(z-pred.inf.intens.last)
  }
  
  optim(c(height.next),pred.func.2,method = "Brent",lower=0,upper=10e6)$par
}



