library(mgcv)


# load and prep data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/infection intensity data prep.R")
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
delta.infection.intensity<-subset(delta.infection.intensity,time<=10)
delta.plants.change<-cbind(delta.plants,"change"=ifelse(delta.plants$plant.inf.intens.next>=delta.plants$plant.inf.intens,1,0)) ## add binary variable indicating growth/stasis or shrinkage
delta.plants.plus<-delta.plants[which(delta.plants.change$change==1),] #subset data to increase
delta.plants.minus<-delta.plants[which(delta.plants.change$change==0),] #subset data to decrease

# visualize data
layout(matrix(c(1,2,4,4,3,3,3,3),2,4,byrow = T))
par(mar=c(5,5,3,0.5))
## histograms
hist(delta.infection.intensity$infection.intensity,main="",breaks=100,xlab="infection intensity",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.infection.intensity$infection.intensity.next-delta.infection.intensity$infection.intensity)/delta.infection.intensity$time,main="",breaks=100,xlab="change in infection intensity per day",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories

plot(c(min(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y")),max(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y"))),c(0,max(infection.intensity$infection.intensity)),type="n",xlab="date",ylab="infection intensity",main="infection intensity",cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
i<-0
plot.cols<-sample(rainbow(101))

for (tag in unique(infection.intensity$Tag))
{
  sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
  points(sub.infection.intensity.1$Date,sub.infection.intensity.1$infection.intensity,col=plot.cols[i],type="l",lwd=2)
  i<-i+1
}

#plot(c(min(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y")),max(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y"))),c(-1,max(log10(infection.intensity$infection.intensity))),type="n",xlab="date",ylab="log 10 plant infection intensity",main="log 10 plant infection intensity",cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
#i<-0
#plot.cols<-sample(rainbow(101))
#
#for (tag in unique(infection.intensity$Tag))
#{
#  sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
#  points(sub.infection.intensity.1$Date,log10(sub.infection.intensity.1$infection.intensity),col=plot.cols[i],type="l",lwd=2)
#  i<-i+1
#}

## plot change
par(mar=c(5,5,3,0.5))
plot(delta.infection.intensity$infection.intensity,delta.infection.intensity$infection.intensity.next,col="black",xlab = 'infection intensity',ylab='next obs. infection intensity',cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3.25,cex=2)
abline(0,1,lty=2)



## plot change
### Indicates need to constrain outcomes, as the change in plant infection intensity is bounded by current infection intensiy.
### To circumvent the constraint problem, we model log10(inf intens at t+1) rather than inf intens at t+1 - inf intens at t
par(mfrow=c(1,1))
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,xlab = "plant infection intensity",ylab="change in plant infection intensity")
abline(0,-1,lty=2)
par(fig=c(.6,.9,.5,.85),new=T,mar=c(0,0,0,0))
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,xlim=c(0,1),ylim=c(-1,10))
abline(0,-1,lty=2)
par(new=F,mar=c(5,5,5,5))

## plot transformed data to model
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next),xlab=expression(log[10]*' plant infection intensity'),ylab=expression(log[10]*' next observed plant infection intensity'),cex.lab=2)
abline(0,1)
## Indicates that most changes are small increases or decreases, but there are many large increases for low infection intensity (around -1), and several large decreases for medium infection intensity (around 0 to 2.5)
## Modeling the change in infection intensity as a single process resulted in a model that  1) often predicted log10(infection intensities) <  -1, 2) generally underpredicted infection intensity but 3) failed to capture large decreases
## As such, we first model whether the change in infection intensity is positive/0 or negative, and then model the growth/stasis or shrinkage in infection intensity as seperate processes.

# analyze data

## fit models--only if not already fit

if(any(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.change.model.RDS"),!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.growth.RDS"),!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.shrinkage.model.RDS")))
{
  ### model change
  change.mod0<-gam(change~
                   s(log10(plant.inf.intens),bs="cs",k=4)+
                   s(max.height,bs='cs',k=4)+
                   s(mean.temp,bs="cs",k=4)+
                   s(max.temp,bs="cs",k=4)+
                   s(min.temp,bs="cs",k=4)+
                   s(mean.abs.hum,bs="cs",k=4)+
                   s(max.abs.hum,bs="cs",k=4)+
                   s(min.abs.hum,bs="cs",k=4)+
                   s(mean.solar,bs="cs",k=4)+
                   s(tot.rain,bs="cs",k=4)+
                   s(site,bs="re",k=4),
                   data=delta.plants.change,
                   family = binomial())
  summary(change.mod0) #indicates that all but s(log10(plant.inf.intens)), s(max.temp), s(max.abs.hum) are insignificant, min.temp marginially significant
  
  change.mod1<-gam(change~
                   s(log10(plant.inf.intens),bs="cs",k=4)+
                   #s(max.height,bs='cs',k=4)+
                   #s(mean.temp,bs="cs",k=4)+
                   s(max.temp,bs="cs",k=4)+
                   s(min.temp,bs="cs",k=4)+
                   #s(mean.abs.hum,bs="cs",k=4)+
                   s(max.abs.hum,bs="cs",k=4)+
                   #s(min.abs.hum,bs="cs",k=4)+
                   #s(mean.solar,bs="cs",k=4)+
                   #(tot.rain,bs="cs",k=4)+
                   s(site,bs="re",k=4),
                   data=delta.plants.change,
                   family = binomial())
  summary(change.mod1) #Dropping marginially significant min.temp leads to model w/ lower AIC, so we retain it in the model. 
  
  saveRDS(change.mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.change.model.RDS")
  
  ### model growth
  
  growth.mod0<-gam(log10(plant.inf.intens.next)~
                     s(log10(plant.inf.intens),bs="cs",k=4)+
                     te(log10(plant.inf.intens),max.height,by=time,bs="cs",k=4)+
                     s(mean.temp,by=time,bs="cs",k=4)+
                     s(max.temp,by=time,bs="cs",k=4)+
                     s(min.temp,by=time,bs="cs",k=4)+
                     s(mean.abs.hum,by=time,bs="cs",k=4)+
                     s(max.abs.hum,by=time,bs="cs",k=4)+
                     s(min.abs.hum,by=time,bs="cs",k=4)+
                     s(mean.solar,by=time,bs="cs",k=4)+
                     s(tot.rain,bs="cs",k=4)+
                     s(site,bs="re",k=4),
                   data=delta.plants.plus)
  summary(growth.mod0) #indicates that all but s(log10(plant.inf.intens) and te(log10(plant.inf.intens),max.height) are insignificant 
  
  growth.mod1<-gam(log10(plant.inf.intens.next)~
                     s(log10(plant.inf.intens),bs="cs",k=4)+
                     te(log10(plant.inf.intens),max.height,by=time,bs="cs",k=4)+
                     #s(mean.temp,by=time,bs="cs",k=4)+
                     #s(max.temp,by=time,bs="cs",k=4)+
                     #s(min.temp,by=time,bs="cs",k=4)+
                     #s(mean.abs.hum,by=time,bs="cs",k=4)+
                     #s(max.abs.hum,by=time,bs="cs",k=4)+
                     #s(min.abs.hum,by=time,bs="cs",k=4)+
                     #s(mean.solar,by=time,bs="cs",k=4)+
                     #s(tot.rain,bs="cs",k=4)+
                     s(site,bs="re",k=4),
                   data=delta.plants.plus)
  summary(growth.mod1) #all predictors now significant
  
  saveRDS(growth.mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.growth.model.RDS")
  
  ## model shrinkage
  
  shrink.mod0<-gam(log10(plant.inf.intens.next)~
                     s(log10(plant.inf.intens),bs="cs",k=4)+
                     te(log10(plant.inf.intens),max.height,by=time,bs="cs",k=4)+
                     s(mean.temp,by=time,bs="cs",k=4)+
                     s(max.temp,by=time,bs="cs",k=4)+
                     s(min.temp,by=time,bs="cs",k=4)+
                     s(mean.abs.hum,by=time,bs="cs",k=4)+
                     s(max.abs.hum,by=time,bs="cs",k=4)+
                     s(min.abs.hum,by=time,bs="cs",k=4)+
                     s(mean.solar,by=time,bs="cs",k=4)+
                     s(tot.rain,bs="cs",k=4)+
                     s(site,bs="re",k=4),
                   data=delta.plants.minus)
  summary(shrink.mod0) #indicates that all but log10(plant.inf.intes), te(log10(plant.inf.intens,max.height)), max.temp, max.abs.hum insignificant, mean.temp marginally significant
  
  shrink.mod1<-gam(log10(plant.inf.intens.next)~
                     s(log10(plant.inf.intens),bs="cs",k=4)+
                     te(log10(plant.inf.intens),max.height,by=time,bs="cs",k=4)+
                     s(mean.temp,by=time,bs="cs",k=4)+
                     s(max.temp,by=time,bs="cs",k=4)+
                     #s(min.temp,by=time,bs="cs",k=4)+
                     #s(mean.abs.hum,by=time,bs="cs",k=4)+
                     s(max.abs.hum,by=time,bs="cs",k=4)+
                     #s(min.abs.hum,by=time,bs="cs",k=4)+
                     #s(mean.solar,by=time,bs="cs",k=4)+
                     #s(tot.rain,bs="cs",k=4)+
                     s(site,bs="re",k=4),
                   data=delta.plants.minus) 
  summary(shrink.mod1) #Dropping marginially significant mean.temp leads to model w/ lower AIC, so we retain it in the model. 
  
  saveRDS(shrink.mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.shrinkage.model.RDS")
  
}

## load best model

plants.change.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.change.model.RDS")
plants.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.growth.model.RDS")
plants.shrinkage.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.shrinkage.model.RDS")

## model checking
### overall results are acceptable
### change model
### plot smooths
#par(mfrow=c(2,3),mar=c(5,5,5,5))
#plot(plants.change.model,select = 1,scale=-1,xlab=expression(log[10]*'plant infection intensity'),ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2)
#plot(plants.change.model,select = 2,scale=-1,xlab="maximm temperature (C)",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")
#plot(plants.change.model,select = 3,scale=-1,xlab="minimum temperature (C)",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")
#plot(plants.change.model,select = 4,scale=-1,xlab="maximm absolute humidity",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")
#plot(plants.change.model,select = 5,scale=-1,xlab="site",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")

#par(mfrow=c(2,2))
#gam.check(plants.change.model) #indicates k vals are sufficient

### growth model
#par(mfrow=c(2,2),mar=c(5,6,5,5))
#plot(plants.growth.model,select = 1,scale=-1,xlab=expression(log[10]*' plant infection intensity'),ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2)
#plot(plants.growth.model,select=2,scheme=2,too.far=5,contour.col="black",labcex=1.5,xlab=expression(log[10]*' plant infection intensity'),cex.lab=2,main="")
#plot(plants.growth.model,select = 3,scale=-1,xlab="site",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")

#par(mfrow=c(1,1))
#plot(log10(delta.plants.plus$plant.inf.intens),log10(delta.plants.plus$plant.inf.intens.next))
#abline(0,1)

#par(mfrow=c(2,2))
#gam.check(plant.growth.model) # indicates k vals arre sufficient

### shrinkage model
#par(mfrow=c(2,3))
#plot(plants.shrinkage.model,select = 1,scale=0,xlab=expression(log[10]*' plant infection intensity'),ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2)
#plot(plants.shrinkage.model,select=2,scheme=2,too.far=5,contour.col="black",labcex=1.5,xlab=expression(log[10]*' plant infection intensity'),cex.lab=2,main="")
#plot(plants.shrinkage.model,select = 3,scale=0,xlab="mean temperature",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")
#plot(plants.shrinkage.model,select = 4,scale=0,xlab="maximum temperature",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")
#plot(plants.shrinkage.model,select = 5,scale=0,xlab="maximum absolute humidity",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")
#plot(plants.shrinkage.model,select = 6,scale=0,xlab="site",ylab="effect",shade=T,shade.col="palegreen",col="darkgreen",lwd=4,cex.lab=2,main="")


#par(mfrow=c(1,1))
#plot(log10(delta.plants.minus$plant.inf.intens),predict(plants.shrinkage.model)) #
#abline(0,1) 

#par(mfrow=c(2,2))
#gam.check(plants.shrinkage.model) # indicates k vals arre sufficient

## plot model pedictions against data

preds<-rep(NA,nrow(delta.plants))

for(i in 1:nrow(delta.plants))
{
  newdata=delta.plants[i,]
  odds<-predict(plants.change.model,newdata = newdata,type="response")
  if(odds >= .5)
  {
   preds[i]<-predict(plants.growth.model,newdata=newdata,type="response") 
  } else
  {
    preds[i]<-predict(plants.shrinkage.model,newdata=newdata,type="response") 
    
  }
}

plot(log10(delta.plants$plant.inf.intens.next),preds)


## visualize model predictions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
library("MASS")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(-1,5),ylim=c(-1,5),type="n",xlab="log 10 plant infection intensity",ylab="pred. change in plant infection intensity",cex.axis=1,cex.lab=1.2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in 10^seq(-1,5,.25))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plants.model(day,dummy.data.inf.intens=i,dummy.data.height=60)
    
    odds<-predict(plants.change.model,newdata = pred.data,type="response")
    n <- 1000
    n.growth<-round(odds*n)
    n.shrink<-1000-n.growth
    
    Xp.growth <- predict(plants.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta.growth <- coef(plants.growth.model) ## posterior mean of coefs
    Vb.growth  <- vcov(plants.growth.model) ## posterior  cov of coefs
    mrand.growth <- mvrnorm(n.growth, beta.growth, Vb.growth) ## simulate n rep coef vectors from posterior
    ilink <- family(plants.growth.model)$linkinv
    
    preds.growth<-rep(NA,n.growth)
    for (j in seq_len(n.growth)) { 
      preds.growth[j]   <- ilink(Xp.growth %*% mrand.growth[j, ])
    }
    
    Xp.shrink <- predict(plants.shrinkage.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta.shrink <- coef(plants.shrinkage.model) ## posterior mean of coefs
    Vb.shrink  <- vcov(plants.shrinkage.model) ## posterior  cov of coefs
    mrand.shrink <- mvrnorm(n.shrink, beta.shrink, Vb.shrink) ## simulate n rep coef vectors from posterior
    ilink <- family(plants.shrinkage.model)$linkinv
    
    preds.shrink<-rep(NA,n.shrink)
    for (j in seq_len(n.shrink)) { 
      preds.shrink[j]   <- ilink(Xp.shrink %*% mrand.shrink[j, ])
    }
    
    y<-c(preds.growth,preds.shrink)
    lower=rbind(lower,data.frame(x=log10(i),y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=log10(i),y=quantile(y,.95)))
    points(log10(i),mean(y),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
abline(0,1,lty=2)
legend("topleft",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")


plot(0,0,xlim=c(-1,5),ylim=c(-1,5),type="n",xlab="plant infection intensity",ylab="pred. next obs. plant infection intensity",cex.axis=1,cex.lab=1.2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in 10^seq(-1,5,.25))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plants.model(75,dummy.data.inf.intens=i,dummy.data.height=15,temp.addition = temp.addition)
    
    odds<-predict(plants.change.model,newdata = pred.data,type="response")
    n <- 1000
    n.growth<-round(odds*n)
    n.shrink<-1000-n.growth
    
    Xp.growth <- predict(plants.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta.growth <- coef(plants.growth.model) ## posterior mean of coefs
    Vb.growth  <- vcov(plants.growth.model) ## posterior  cov of coefs
    mrand.growth <- mvrnorm(n.growth, beta.growth, Vb.growth) ## simulate n rep coef vectors from posterior
    ilink <- family(plants.growth.model)$linkinv
    
    preds.growth<-rep(NA,n.growth)
    for (j in seq_len(n.growth)) { 
      preds.growth[j]   <- ilink(Xp.growth %*% mrand.growth[j, ])
    }
    
    Xp.shrink <- predict(plants.shrinkage.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta.shrink <- coef(plants.shrinkage.model) ## posterior mean of coefs
    Vb.shrink  <- vcov(plants.shrinkage.model) ## posterior  cov of coefs
    mrand.shrink <- mvrnorm(n.shrink, beta.shrink, Vb.shrink) ## simulate n rep coef vectors from posterior
    ilink <- family(plants.shrinkage.model)$linkinv
    
    preds.shrink<-rep(NA,n.shrink)
    for (j in seq_len(n.shrink)) { 
      preds.shrink[j]   <- ilink(Xp.shrink %*% mrand.shrink[j, ])
    }
    
    y<-c(preds.growth,preds.shrink)
    lower=rbind(lower,data.frame(x=log10(i),y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=log10(i),y=quantile(y,.95)))
    points(log10(i),mean(y),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
abline(0,1,lty=2)
legend("topleft",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

predict.plant.inf.trajectory<-function(site,temp.addition,plant.height,color,pred.window=2,plot=T,output=F)
{  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  weath.dat<-all.weath[which(all.weath$site==site),]
  temp.rh.dat<-all.temp.rh[which(all.temp.rh$site==site),]
  min.date<-max(min(unique(as.Date(weath.dat$date))),min(unique(as.Date(temp.rh.dat$date.time))))
  max.date<-min(max(unique(as.Date(weath.dat$date))),max(unique(as.Date(temp.rh.dat$date.time))))
  dates<-seq(min.date,max.date,pred.window)
  start.inf.intens<-1
  xcords<-rep(NA,length(dates)) #time values
  ycords<-rep(NA,length(dates)) #inf intensity values
  hcords<-rep(NA,length(dates)) #height values
  
  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.inf.intens
    h<-plant.height
    xcords.new<-c(1)
    ycords.new<-c(log10(i))
    hcords.new<-c(h)
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,i,dummy.data.max.height=h,temp.addition = temp.addition)
      pred.data.h<-pred.data
      colnames(pred.data.h)[which(colnames(pred.data.h)=="max.height")]<-"height"
      
      odds<-predict(plants.change.model,newdata = pred.data,type="response")
      draw<-runif(1)
      if(draw<=odds)
      {
        Xp.growth <- predict(plants.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
        beta.growth <- coef(plants.growth.model) ## posterior mean of coefs
        Vb.growth  <- vcov(plants.growth.model) ## posterior  cov of coefs
        n <-2
        mrand.growth <- mvrnorm(n, beta.growth, Vb.growth) ## simulate n rep coef vectors from posterior
        ilink <- family(plants.growth.model)$linkinv
        preds <- rep(NA,n)
        for (l in seq_len(n)) { 
          preds[l]   <- ilink(Xp.growth %*% mrand.growth[l, ])
        }
        y<-preds[1]
        i<-10^y
      } else
      {
        Xp.shrink <- predict(plants.shrinkage.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
        beta.shrink <- coef(plants.shrinkage.model) ## posterior mean of coefs
        Vb.shrink  <- vcov(plants.shrinkage.model) ## posterior  cov of coefs
        n <-2
        mrand.shrink <- mvrnorm(n, beta.shrink, Vb.shrink) ## simulate n rep coef vectors from posterior
        ilink <- family(plants.shrinkage.model)$linkinv
        preds <- rep(NA,n)
        for (l in seq_len(n)) { 
          preds[l]   <- ilink(Xp.shrink %*% mrand.shrink[l, ])
        }
        y<-preds[1]
        i<-10^y
      }
      
      beta.h <- coef(plant.growth.model) ## posterior mean of coefs
      Vb.h   <- vcov(plant.growth.model) ## posterior  cov of coefs
      n.h <-2
      mrand.h <- mvrnorm(n.h, beta.h, Vb.h) ## simulate n rep coef vectors from posterior
      Xp.h <- predict(plant.growth.model, newdata = pred.data.h, exlude="s(site)",type="lpmatrix")
      ilink <- family(plant.growth.model)$linkinv
      preds.h <- rep(NA,n.h)
      for (l in seq_len(n.h)) { 
        preds.h[l]   <- ilink(Xp.h %*% mrand.h[l, ])[1]
      }
      h<-h+preds.h[1]
      
      
      reps<-reps+pred.window
      xcords.new<-c(xcords.new,reps)
      ycords.new<-c(ycords.new,y)
      hcords.new<-c(hcords.new,h)
    }
    xcords<-rbind(xcords,xcords.new)
    ycords<-rbind(ycords,ycords.new)
    hcords<-rbind(hcords,hcords.new)
    print(j)
  }
  xcords<-xcords[-1,]
  ycords<-ycords[-1,]
  hcords<-hcords[-1,]
  
  if(plot)
  {
    for(k in 1:dim(xcords)[1])
    {
      points(xcords[k,],ycords[k,],type="l",col=color) 
    } 
  }
  if(output)
  {
    list(ycords,hcords)
  }
}

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

plot.purple<-t_col("purple",80)
plot.red<-t_col("red",80)
plot.orange<-t_col("orange",80)

plant.height<-15

### one day ahead projection
#layout(matrix(c(1,1,2),1,3,byrow = T))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,36),ylim=c(-1,5),ylab=expression(log[10]*'plant inf. intensity'),xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
dat1<-predict.plant.inf.trajectory("GM",0,plant.height=plant.height,pred.window=1,plot.orange,T,T) 
points(1:36,colMeans(dat1[[1]]),type="l",col="orange",lwd=4,lty=2)

dat2<-predict.plant.inf.trajectory("GM",1.8,plant.height=plant.height,pred.window=1,plot.red,T,T) 
points(1:36,colMeans(dat2[[1]]),type="l",col="red",lwd=4,lty=2)

dat3<-predict.plant.inf.trajectory("GM",3.7,plant.height=plant.height,pred.window=1,plot.purple,T,T) 
points(1:36,colMeans(dat3[[1]]),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",xlim=c(1,36),ylim=c(-1,5),ylab=expression(log[10]*'plant inf. intensity'),xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
polygon(c(1:36,36:1),c(apply(dat1[[1]],2,quantile,probs=.05),rev(apply(dat1[[1]],2,quantile,probs=.95))),col=plot.orange,density=100)
polygon(c(1:36,36:1),c(apply(dat2[[1]],2,quantile,probs=.05),rev(apply(dat2[[1]],2,quantile,probs=.95))),col=plot.red,density=100)
polygon(c(1:36,36:1),c(apply(dat3[[1]],2,quantile,probs=.05),rev(apply(dat3[[1]],2,quantile,probs=.95))),col=plot.purple,density=100)
points(1:36,colMeans(dat1[[1]]),type="l",col="orange",lwd=4,lty=2)
points(1:36,colMeans(dat2[[1]]),type="l",col="red",lwd=4,lty=2)
points(1:36,colMeans(dat3[[1]]),type="l",col="purple",lwd=4,lty=2)


