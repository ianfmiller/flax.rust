library(mgcv)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.plant.heights<-subset(delta.plant.heights,time<=8)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{
  set.seed(23094867)
  mod<-gam((height.next-height)/time~
             te(height,inf.intens)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(mean.daily.rain)+
             s(tag,bs="re")+
             s(site,bs="re"),
           select = T,
           method="REML",
           data=delta.plant.heights,
           control = list(nthreads=4))

  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
par(mfrow=c(2,2))
gam.check(plant.growth.model) #Indicates that basis dimension is sufficient
concurvity(plant.growth.model,full=F) #No obvious issues