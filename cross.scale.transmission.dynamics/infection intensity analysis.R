library(mgcv)

# load and prep data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/infection intensity data prep.R")
delta.infection.intensity<-subset(delta.infection.intensity,time<=8)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS"))
{
  set.seed(5708389)
  mod<-gam((infection.intensity.next-infection.intensity)/time~
             te(max.height,log10(infection.intensity),k=c(10,15))+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(mean.daily.rain)+
             s(tag,bs="re")+
             s(site,bs="re"),
           select = T,
           method="REML",
           data=delta.infection.intensity,
           control = list(nthreads=4))
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS")
  
}

## load best model

infection.intensity.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS")

## model checking 
par(mfrow=c(2,2))
dummy.delta.infection.intensity<-delta.infection.intensity
dummy.delta.infection.intensity$log.10.infection.intensity<-log10(delta.infection.intensity$infection.intensity)
set.seed(5708389)
dummy.infection.intensity.model<-gam((infection.intensity.next-infection.intensity)/time~
           te(max.height,log.10.infection.intensity,k=c(10,15))+
           s(mean.temp)+
           s(max.temp)+
           s(min.temp)+
           s(mean.abs.hum)+
           s(mean.daily.rain)+
           s(tag,bs="re")+
           s(site,bs="re"),
         select = T,
         method="REML",
         data=dummy.delta.infection.intensity,
         control = list(nthreads=4))
saveRDS(dummy.infection.intensity.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/dummy.infection.intensity.model.RDS")

gam.check(dummy.infection.intensity.model) #Indicates that k values are sufficient.
concurvity(infection.intensity.model,full=F) #no obvious issues
