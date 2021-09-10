library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")

delta.pustules<-subset(delta.pustules,time<=7)
delta.pustules.standardized<-delta.pustules
delta.pustules.standardized$time<-(delta.pustules$time-mean(delta.pustules$time,na.rm=T))/sd(delta.pustules$time)
delta.pustules.standardized$area<-(delta.pustules$area-mean(delta.pustules$area,na.rm=T))/sd(delta.pustules$area)
delta.pustules.standardized$mean.temp<-(delta.pustules$mean.temp-mean(delta.pustules$mean.temp,na.rm=T))/sd(delta.pustules$mean.temp)
delta.pustules.standardized$max.temp<-(delta.pustules$max.temp-mean(delta.pustules$max.temp,na.rm=T))/sd(delta.pustules$max.temp)
delta.pustules.standardized$min.temp<-(delta.pustules$min.temp-mean(delta.pustules$min.temp,na.rm=T))/sd(delta.pustules$min.temp)
delta.pustules.standardized$mean.abs.hum<-(delta.pustules$mean.abs.hum-mean(delta.pustules$mean.abs.hum,na.rm=T))/sd(delta.pustules$mean.abs.hum)
delta.pustules.standardized$max.abs.hum<-(delta.pustules$max.abs.hum-mean(delta.pustules$max.abs.hum,na.rm=T))/sd(delta.pustules$max.abs.hum)
delta.pustules.standardized$min.abs.hum<-(delta.pustules$min.abs.hum-mean(delta.pustules$min.abs.hum,na.rm=T))/sd(delta.pustules$min.abs.hum)
delta.pustules.standardized$mean.vpd<-(delta.pustules$mean.vpd-mean(delta.pustules$mean.vpd,na.rm=T))/sd(delta.pustules$mean.vpd)
delta.pustules.standardized$max.vpd<-(delta.pustules$max.vpd-mean(delta.pustules$max.vpd,na.rm=T))/sd(delta.pustules$max.vpd)
delta.pustules.standardized$min.vpd<-(delta.pustules$min.vpd-mean(delta.pustules$min.vpd,na.rm=T))/sd(delta.pustules$min.vpd)
delta.pustules.standardized$mean.wetness<-(delta.pustules$mean.wetness-mean(delta.pustules$mean.wetness,na.rm=T))/sd(delta.pustules$mean.wetness)
delta.pustules.standardized$tot.rain<-(delta.pustules$tot.rain-mean(delta.pustules$tot.rain,na.rm=T))/sd(delta.pustules$tot.rain)
delta.pustules.standardized$mean.solar<-(delta.pustules$mean.solar-mean(delta.pustules$mean.solar,na.rm=T))/sd(delta.pustules$mean.solar)



# visualize data

## histograms
par(mfrow=c(2,1))
hist(delta.pustules$area,main="pustule area",breaks=100,xlab="area")
hist(delta.pustules$area.next-delta.pustules$area,main="change in pustule area",breaks=100,xlab="change in area")

## plot trajectories
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="date",ylab="pustule area")
i<-0

plot.cols<-sample(rainbow(2042))

for (tag in unique(pustules$tag))
{
  sub.pustules1<-pustules[which(pustules$tag==tag),]
  
  for (color in unique(sub.pustules1$color))
  {
    sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
    {
      sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
      
      for(pustule.number in unique(sub.pustules3$pustule.number))
      {
        i<-i+1
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        points(sub.pustules4$date,sub.pustules4$area,col=plot.cols[i],type="l",lwd=.5)
      }
    }
  }
}

## plot change
plot(delta.pustules$area,delta.pustules$area.next,col="grey",xlab = "area",ylab="next obs. area")
abline(0,1)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ s(area",predictors[x],'site,bs="re")'),collapse=") + s(")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(gam(area.next~s(area)+s(site,bs="re"),data=delta.pustules.standardized)) #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-gam(model.set[[i]],data=delta.pustules.standardized,REML = F))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
}

## load best model

### LEFT OFF HERE

pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")

# predict climate change effect
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,.1),ylim=c(-.02,.02),type="n",xlab="area (cm)",ylab="pred. change in area (cm)",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(0,.1,.01))
  {
    pred.data<-get.pred.data.temp.mean.quantile.pustule.model(day,i)
    y<-bootMer(pustule.model, FUN=function(x)predict(x,newdata=pred.data, re.form=NA),nsim=100)$t
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=0)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=2,bty="n")

plot(0,0,xlim=c(0,.1),ylim=c(-.02,.02),type="n",xlab="area (cm)",ylab="pred. change in area (cm)",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(0,.1,.01))
  {
    pred.data<-get.pred.data.temp.mean.quantile.pustule.model(75,i,temp.addition = temp.addition)
    y<-bootMer(pustule.model, FUN=function(x)predict(x,newdata=pred.data, re.form=NA),nsim=100)$t
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=0)
}
legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),pch=16,cex=2,bty="n")
