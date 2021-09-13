library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.height<-subset(delta.height,time<=7)
delta.height.standardized<-delta.height
delta.height.standardized$time<-(delta.height$time-mean(delta.height$time,na.rm=T))/sd(delta.height$time)
delta.height.standardized$area<-(delta.height$area-mean(delta.height$area,na.rm=T))/sd(delta.height$area)
delta.height.standardized$mean.temp<-(delta.height$mean.temp-mean(delta.height$mean.temp,na.rm=T))/sd(delta.height$mean.temp)
delta.height.standardized$max.temp<-(delta.height$max.temp-mean(delta.height$max.temp,na.rm=T))/sd(delta.height$max.temp)
delta.height.standardized$min.temp<-(delta.height$min.temp-mean(delta.height$min.temp,na.rm=T))/sd(delta.height$min.temp)
delta.height.standardized$mean.abs.hum<-(delta.height$mean.abs.hum-mean(delta.height$mean.abs.hum,na.rm=T))/sd(delta.height$mean.abs.hum)
delta.height.standardized$max.abs.hum<-(delta.height$max.abs.hum-mean(delta.height$max.abs.hum,na.rm=T))/sd(delta.height$max.abs.hum)
delta.height.standardized$min.abs.hum<-(delta.height$min.abs.hum-mean(delta.height$min.abs.hum,na.rm=T))/sd(delta.height$min.abs.hum)
delta.height.standardized$mean.vpd<-(delta.height$mean.vpd-mean(delta.height$mean.vpd,na.rm=T))/sd(delta.height$mean.vpd)
delta.height.standardized$max.vpd<-(delta.height$max.vpd-mean(delta.height$max.vpd,na.rm=T))/sd(delta.height$max.vpd)
delta.height.standardized$min.vpd<-(delta.height$min.vpd-mean(delta.height$min.vpd,na.rm=T))/sd(delta.height$min.vpd)
delta.height.standardized$mean.wetness<-(delta.height$mean.wetness-mean(delta.height$mean.wetness,na.rm=T))/sd(delta.height$mean.wetness)
delta.height.standardized$tot.rain<-(delta.height$tot.rain-mean(delta.height$tot.rain,na.rm=T))/sd(delta.height$tot.rain)
delta.height.standardized$mean.solar<-(delta.height$mean.solar-mean(delta.height$mean.solar,na.rm=T))/sd(delta.height$mean.solar)



# visualize data
par(mfrow =c(1,1))
colors <- ifelse(delta.height$height.next > delta.height$height, "green", "red")
plot(height.next-height~height, data = delta.height, col = colors,xlab="height",ylab="change in heght")


## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant.growth.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("height.next ~ s(time, k=5) + s(height",predictors[x],'site,k=10,bs="re")'),collapse=",k=10) + s(")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(gam(height.next~s(time,k=5)+s(height,k=10)+s(site,bs="re",k=10),data=delta.height)) #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-gam(model.set[[i]],data=delta.height))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
plot(plant.growth.model,scale=0,pages=1) #plot smooths
#gam.check(pustule.model) #indicates no more knots needed

## visualize model with standardized effects
standardized.var.model<-gam(plant.growth.model$formula,data=delta.height.standardized,)
plot(standardized.var.model,scale=0,pages=1) #plot standardized smooths
