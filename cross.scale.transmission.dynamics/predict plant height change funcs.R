# load data and model
library(MASS)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")


# function for subsetting temp/rh
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

#function for predicting plant height change
predict.plant.growth<-function(height.last,site,date0,date1,exclude.site=T)
{
  if(class(date0)=="Date") {date0<-as.POSIXct(date0)}
  if(class(date1)=="Date") {date1<-as.POSIXct(date1)}
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  # calculate environmental variable metrics
  new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
  new.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #time (in days) during which temp between 16 and 22 celsius
  new.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #time (in days) during which temp between 7 and 30 celsius
  new.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T) #Dew point days
  
  #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
  new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
  new.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
  new.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)

  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("height"=height.last,"dew.point.days"=new.dew.point.days,"temp.dew.point.days"=new.temp.dew.point.days,"site"=site)
  
  if(exclude.site) {height.next<-predict(plant.growth.model,newdata=pred.data,exclude = 's(site)')} else {height.next<-predict(plant.growth.model,newdata=pred.data)}
  if(height.next<1) {height.next<-1}
  height.next
}

predict.plant.growth.boot<-function(height.last,site,date0,date1)
{
  if(class(date0)=="Date") {date0<-as.POSIXct(date0)}
  if(class(date1)=="Date") {date1<-as.POSIXct(date1)}
  
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  # calculate environmental variable metrics
  new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
  new.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #time (in days) during which temp between 16 and 22 celsius
  new.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #time (in days) during which temp between 7 and 30 celsius
  new.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T) #Dew point days
  
  #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
  new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
  new.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
  new.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)

  delta.days<-as.numeric(date1-date0)
  
  # make forward prediction
  pred.data<-data.frame("height"=height.last,"dew.point.days"=new.dew.point.days,"temp.dew.point.days"=new.temp.dew.point.days,"site"=site)
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
  height.next<-preds[1]
  if(height.next<1) {height.next<-1}
  height.next
}

#function for predicting plant height change
predict.plant.growth.last<-function(height.next,site,date0,date1,exclude.site=T)
{
  if(class(date0)=="Date") {date0<-as.POSIXct(date0)}
  if(class(date1)=="Date") {date1<-as.POSIXct(date1)}
  
  # load weather data
  ## subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  # calculate environmental variable metrics
  new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
  new.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #time (in days) during which temp between 16 and 22 celsius
  new.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #time (in days) during which temp between 7 and 30 celsius
  new.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T) #Dew point days
  
  #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
  new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
  new.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
  new.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)

  delta.days<-as.numeric(date1-date0)
  
  # make backwards prediction
  pred.func<-function(x)
  {
    plant.height.last.test<-x
    pred.data<-data.frame("height"=plant.height.last.test,"dew.point.days"=new.dew.point.days,"temp.dew.point.days"=new.temp.dew.point.days,"site"=site)
    if(exclude.site) {plant.height.next.pred<-predict(plant.growth.model,newdata=pred.data,exclude = 's(site)')} else {plant.height.next.pred<-predict(plant.growth.model,newdata=pred.data)}
    abs(plant.height.next.pred-height.next)
  }
  plant.height.last<-optim(c(height.next),pred.func,method = "Brent",lower=0,upper=10e6)$par
  if(plant.height.last<1) {plant.height.last<-1}
  plant.height.last
}




