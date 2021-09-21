# load data and model
library(MASS)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
plant.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")
n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")


# function for subsetting temp/rh
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

#function for predicting plant inf intens change
predict.plant.inf.intens<-function(plant.inf.intens.last,max.height.last,site,date0,date1)
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
  
  new.tot.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,
                                      "time"=delta.days,"site"=site,
                                      "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                      "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                      "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-1 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,
                                         "time"=delta.days,"site"=site,
                                         "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                         "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                         "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                                         "pred.pustule.diam.growth"=pred.pustule.diam.growth)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  # make forward prediction
  pred.data<-data.frame("plant.inf.intens"=plant.inf.intens.last,"max.height"=max.height.last,
                        "time"=delta.days,"site"=site,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                        "pred.pustule.diam.growth"=pred.pustule.diam.growth,"pred.pustule.num.increase"=pred.pustule.num.increase)
  plant.inf.intens.next<-plant.inf.intens.last+predict(plant.model,newdata=pred.data,exclude = 's(site)')
  if(plant.inf.intens.next<.01) {plant.inf.intens.next<-.01}
  plant.inf.intens.next
}

predict.plant.inf.intens.boot<-function(plant.inf.intens.last,max.height.last,site,date0,date1)
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
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,
                                      "time"=delta.days,"site"=site,
                                      "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                      "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                      "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-1 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,
                                         "time"=delta.days,"site"=site,
                                         "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                         "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                         "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                                         "pred.pustule.diam.growth"=pred.pustule.diam.growth)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  # make forward prediction
  pred.data<-data.frame("plant.inf.intens"=plant.inf.intens.last,"max.height"=max.height.last,
                        "time"=delta.days,"site"=site,
                        "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                        "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                        "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                        "pred.pustule.diam.growth"=pred.pustule.diam.growth,"pred.pustule.num.increase"=pred.pustule.num.increase)
  
  Xp <- predict(plant.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
  beta <- coef(plant.model) ## posterior mean of coefs
  Vb   <- vcov(plant.model) ## posterior  cov of coefs
  n <- 2
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  preds <- rep(NA, n)
  ilink <- family(plant.model)$linkinv
  for (j in seq_len(n)) { 
    preds[j]   <- ilink(Xp %*% mrand[j, ])
  }
  plant.inf.intens.next<-preds[1]+plant.inf.intens.last
  if(plant.inf.intens.next<.01) {plant.inf.intens.next<-.01}
  plant.inf.intens.next
}

#function for predicting plant inf intens change
predict.plant.inf.intens.last<-function(plant.inf.intens.next,max.height.last,site,date0,date1)
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
  
  abs.hum<-6.112*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh*2.1674/(273.15+T)
  new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
  new.max.abs.hum<-max(abs.hum,na.rm=T)
  new.min.abs.hum<-min(abs.hum,na.rm=T)
  
  #svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
  #avps<- temp.rh.sub$rh / 100 * svps #actual vapor pressures 
  #vpds<-avps-svps
  
  new.mean.vpd<-mean(vpds,na.rm=T)
  new.max.vpd<-max(vpds,na.rm=T)
  new.min.vpd<-min(vpds,na.rm=T)
  
  new.tot.rain<-sum(weath.sub$rain,na.rm=T)
  new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
  
  delta.days<-as.numeric(date1-date0)
  
  #predict pustule growth from pustule growth model and enviro conditions
  #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
  pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
  obs.time<-delta.days
  pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,
                                      "time"=delta.days,"site"=site,
                                      "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                      "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                      "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain)
  pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,exclude = 's(site)')
  
  #predict change in number of pustules from enviro conditions
  #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
  n.pustules.model.new.n.pustules<-1 #included only for offset, picked 0 for ease of interpretability
  obs.time<-delta.days
  n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,
                                         "time"=delta.days,"site"=site,
                                         "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                         "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                         "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                                         "pred.pustule.diam.growth"=pred.pustule.diam.growth)
  pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,exclude = 's(site)')
  
  # make backwards prediction
  pred.func<-function(x)
  {
    plant.inf.intens.last.test<-x
    pred.data<-data.frame("plant.inf.intens"=plant.inf.intens.last.test,"max.height"=max.height.last,
                          "time"=delta.days,"site"=site,
                          "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                          "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                          "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain,
                          "pred.pustule.diam.growth"=pred.pustule.diam.growth,"pred.pustule.num.increase"=pred.pustule.num.increase)
    plant.inf.intens.next.pred<-plant.inf.intens.last.test+predict(plant.model,newdata=pred.data,exclude = 's(site)')
    abs(plant.inf.intens.next.pred-plant.inf.intens.next)
  }
  plant.inf.intens.last<-optim(c(plant.inf.intens.next),pred.func,method = "Brent",lower=0,upper=10e6)$par
  if(plant.inf.intens.last<0.1) {plant.inf.intens.last<-.1}
  plant.inf.intens.last
}




