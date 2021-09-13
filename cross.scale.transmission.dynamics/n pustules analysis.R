library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/n pustules data prep.R")

delta.n.pustules<-subset(delta.n.pustules,time<=7)
delta.n.pustules.standardized<-delta.n.pustules
delta.n.pustules.standardized$time<-(delta.n.pustules$time-mean(delta.n.pustules$time,na.rm=T))/sd(delta.n.pustules$time)
delta.n.pustules.standardized$n.pustules<-(delta.n.pustules$n.pustules-mean(delta.n.pustules$n.pustules,na.rm=T))/sd(delta.n.pustules$n.pustules)
delta.n.pustules.standardized$mean.temp<-(delta.n.pustules$mean.temp-mean(delta.n.pustules$mean.temp,na.rm=T))/sd(delta.n.pustules$mean.temp)
delta.n.pustules.standardized$max.temp<-(delta.n.pustules$max.temp-mean(delta.n.pustules$max.temp,na.rm=T))/sd(delta.n.pustules$max.temp)
delta.n.pustules.standardized$min.temp<-(delta.n.pustules$min.temp-mean(delta.n.pustules$min.temp,na.rm=T))/sd(delta.n.pustules$min.temp)
delta.n.pustules.standardized$mean.abs.hum<-(delta.n.pustules$mean.abs.hum-mean(delta.n.pustules$mean.abs.hum,na.rm=T))/sd(delta.n.pustules$mean.abs.hum)
delta.n.pustules.standardized$max.abs.hum<-(delta.n.pustules$max.abs.hum-mean(delta.n.pustules$max.abs.hum,na.rm=T))/sd(delta.n.pustules$max.abs.hum)
delta.n.pustules.standardized$min.abs.hum<-(delta.n.pustules$min.abs.hum-mean(delta.n.pustules$min.abs.hum,na.rm=T))/sd(delta.n.pustules$min.abs.hum)
delta.n.pustules.standardized$mean.vpd<-(delta.n.pustules$mean.vpd-mean(delta.n.pustules$mean.vpd,na.rm=T))/sd(delta.n.pustules$mean.vpd)
delta.n.pustules.standardized$max.vpd<-(delta.n.pustules$max.vpd-mean(delta.n.pustules$max.vpd,na.rm=T))/sd(delta.n.pustules$max.vpd)
delta.n.pustules.standardized$min.vpd<-(delta.n.pustules$min.vpd-mean(delta.n.pustules$min.vpd,na.rm=T))/sd(delta.n.pustules$min.vpd)
delta.n.pustules.standardized$mean.wetness<-(delta.n.pustules$mean.wetness-mean(delta.n.pustules$mean.wetness,na.rm=T))/sd(delta.n.pustules$mean.wetness)
delta.n.pustules.standardized$tot.rain<-(delta.n.pustules$tot.rain-mean(delta.n.pustules$tot.rain,na.rm=T))/sd(delta.n.pustules$tot.rain)
delta.n.pustules.standardized$mean.solar<-(delta.n.pustules$mean.solar-mean(delta.n.pustules$mean.solar,na.rm=T))/sd(delta.n.pustules$mean.solar)



# visualize data

## histograms
par(mfrow=c(2,1))
hist(delta.n.pustules$n.pustules,main="n pustules",breaks=100,xlab="n pustules")
hist(delta.n.pustules$n.pustules.next-delta.n.pustules$n.pustules,main="change in n pustules",breaks=100,xlab="change in n pustules")

## plot trajectories
par(mfrow=c(1,1))
plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,max(n.pustules$N.pustules)),type="n",xlab="date",ylab="pustule area")
#plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,50),type="n",xlab="date",ylab="pustule area")

i<-0

plot.cols<-sample(rainbow(300))

for (tag in unique(n.pustules$tag))
{
  sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
  
  for (color in unique(sub.n.pustules1$color))
  {
    sub.n.pustules2<-sub.n.pustules1[which(sub.n.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.n.pustules2$leaf.iteration)) 
    {
      sub.n.pustules3<-sub.n.pustules2[which(sub.n.pustules2$leaf.iteration==leaf.iteration),]
      
      i<-i+1
      points(sub.n.pustules3$date,sub.n.pustules3$N.pustules,col=plot.cols[i],type="l",lwd=.5)
      
    }
  }
}

## one trajectory
par(mfrow=c(1,1),mar=c(6,6,6,6))
tags<-c(86,88,106,112,124)
i<-0
plot.cols<-sample(rainbow(36))
plot(c(min(n.pustules[which(n.pustules$tag %in% tags),]$date),max(n.pustules[which(n.pustules$tag %in% tags),]$date)),c(0,max(n.pustules[which(n.pustules$tag %in% tags),]$N.pustules)),type="n",xlab="date",ylab="N pustules",cex.lab=2,cex.axis=2)
sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
for (tag in tags)
{
  sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
  
  for (color in unique(sub.n.pustules1$color))
  {
    sub.n.pustules2<-sub.n.pustules1[which(sub.n.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.n.pustules2$leaf.iteration)) 
    {
      sub.n.pustules3<-sub.n.pustules2[which(sub.n.pustules2$leaf.iteration==leaf.iteration),]
      
      i<-i+1
      points(sub.n.pustules3$date,sub.n.pustules3$N.pustules,col=plot.cols[i],type="l",lwd=5)
      
    }
  }
}
#plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,50),type="n",xlab="date",ylab="pustule area")


## plot change
par(mfrow=c(1,1))
plot(delta.n.pustules$n.pustules,delta.n.pustules$n.pustules.next,col="black",xlab = "N pustules",ylab="next obs. N pustules",cex.lab=2,cex.axis=2)
abline(0,1,lty=2)
mtext(text="N = 650",cex=2)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/n.pustules.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("n.pustules.next ~ s(time,k=10) + s(n.pustules",predictors[x],'site,k=15,bs="re")'),collapse=",k=15) + s(")))

  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(gam(n.pustules.next~s(n.pustules,k=15)+s(time,k=15)+s(site,bs='re',k=15),data=delta.n.pustules)) #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-gam(model.set[[i]],data=delta.n.pustules))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }

  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
}

## load best model

n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")

## model checking
plot(n.pustules.model,scale=0,pages=1) #plot smooths
gam.check(n.pustules.model) #indicates more knots needed
#gam.check(more.knots.mod) #indicates that increasing number of knots doesn't solve the issue of unevenly distributed residuals. This indicates that the data is the root cause and original model is OK

## visualize model with standardized effects
standardized.var.model<-gam(best.model$formula,data=delta.n.pustules.standardized,)
plot(standardized.var.model,scale=0,pages=1) #plot standardized smooths

# predict climate change effect
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,26),ylim=c(-2,3),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(1,26,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.n.pustules.model(day,i)
    y<-bootMer(n.pustules.model, FUN=function(x)predict(x,newdata=pred.data, re.form=NA),nsim=100)$t
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=2,bty="n")

plot(0,0,xlim=c(0,26),ylim=c(-2,3),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(1,26,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.n.pustules.model(75,i,temp.addition = temp.addition)
    y<-bootMer(n.pustules.model, FUN=function(x)predict(x,newdata=pred.data, re.form=NA),nsim=100)$t
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),pch=16,cex=2,bty="n")
