# load data and model
library(MASS)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")


# function for subsetting temp/rh
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

#function for predicting plant height change
predict.plant.growth<-function(height.last,site,date0,date1,exclude.site=T)
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
  
  abs.hum<-6.112*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh*2.1674/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
  #avps<- temp.rh.sub$rh / 100 * svps #actual vapor pressures 
  #vpds<-avps-svps
  
  #new.mean.vpd<-mean(vpds,na.rm=T)
  #new.max.vpd<-max(vpds,na.rm=T)
  #new.min.vpd<-min(vpds,na.rm=T)
  
  new.tot.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)

  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("time"=delta.days,"height"=height.last,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                        "site"=site)
  
  if(exclude.site) {height.next<-height.last+predict(plant.growth.model,newdata=pred.data,exclude = 's(site)')} else {height.next<-height.last+predict(plant.growth.model,newdata=pred.data)}
  if(height.next<1) {height.next<-1}
  height.next
}

predict.plant.growth.boot<-function(height.last,site,date0,date1)
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
  
  #calculate environmental variable metrics
  new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
  new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
  new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
  
  abs.hum<-6.112*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh*2.1674/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
  #avps<- temp.rh.sub$rh / 100 * svps #actual vapor pressures 
  #vpds<-avps-svps
  
  #new.mean.vpd<-mean(vpds,na.rm=T)
  #new.max.vpd<-max(vpds,na.rm=T)
  #new.min.vpd<-min(vpds,na.rm=T)
  
  new.tot.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("time"=delta.days,"height"=height.last,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                        "site"=site)
  Xp <- predict(plant.growth.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
  beta <- coef(plant.growth.model) ## posterior mean of coefs
  Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
  n <- 2
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  preds <- rep(NA, n)
  ilink <- family(plant.growth.model)$linkinv
  for (j in seq_len(n)) { 
    preds[j]   <- ilink(Xp %*% mrand[j, ])
  }
  height.next<-height.last+preds[1]
  if(height.next<1) {height.next<-1}
  height.next
}

#function for predicting plant height change
predict.plant.growth.last<-function(height.next,site,date0,date1,exclude.site=T)
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
  
  abs.hum<-6.112*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh*2.1674/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
  #avps<- temp.rh.sub$rh / 100 * svps #actual vapor pressures 
  #vpds<-avps-svps
  
  #new.mean.vpd<-mean(vpds,na.rm=T)
  #new.max.vpd<-max(vpds,na.rm=T)
  #new.min.vpd<-min(vpds,na.rm=T)
  
  new.tot.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  # make backwards prediction
  pred.func<-function(x)
  {
    plant.height.last.test<-x
    pred.data<-data.frame("time"=delta.days,"height"=plant.height.last.test,
                           "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                           "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                           "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                           "site"=site)
    if(exclude.site) {plant.height.next.pred<-plant.height.last.test+predict(plant.growth.model,newdata=pred.data,exclude = 's(site)')} else {plant.height.next.pred<-plant.height.last.test+predict(plant.growth.model,newdata=pred.data)}
    abs(plant.height.next.pred-height.next)
  }
  plant.height.last<-optim(c(height.next),pred.func,method = "Brent",lower=0,upper=10e6)$par
  if(plant.height.last<1) {plant.height.last<-1}
  plant.height.last
}




