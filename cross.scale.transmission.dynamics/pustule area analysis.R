library(mgcv)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")

delta.pustules<-subset(delta.pustules,time<=8)

# analyze data

## fit model--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS"))
{
  set.seed(34020968)
  mod<-gam((area.next-area)/time~
              s(area)+
              s(mean.temp)+
              s(max.temp)+
              s(min.temp)+
              s(mean.abs.hum)+
              s(mean.daily.rain)+
              s(tag,bs="re")+
              s(site,bs="re"),
            select = T,
            method="REML",
            data=delta.pustules,
            control = list(nthreads=4))
  
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
}

## load best model

pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")

## model checking
par(mfrow=c(2,2))
gam.check(pustule.model) #Indicates that k should be higher all smooths aside from area. Increasing k to significantly higher values is not possible due to terms having fewer unique covariate combinations than specified maximum degrees of freedom. 
concurvity(pustule.model,full=F) #No obvious issues.