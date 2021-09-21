library(mgcv)
library(lme4)
library(lmerTest)
library(progress)
library(gratia)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.height<-subset(delta.height,time<=7)

# visualize data
par(mfrow =c(1,1))
colors <- ifelse(delta.height$height.next > delta.height$height, "green", "red")
plot(height.next-height~height, data = delta.height, col = colors,xlab="height",ylab="change in heght")


## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{

  mod0<-gam(height.next~s(height,by=time,bs="cs",k=4)+
              s(inf.intens,by=time,bs="cs",k=4)+
              s(mean.temp,by=time,bs="cs",k=4)+
              s(max.temp,by=time,bs="cs",k=4)+
              s(min.temp,by=time,bs="cs",k=4)+
              s(mean.abs.hum,by=time,bs="cs",k=4)+
              s(max.abs.hum,by=time,bs="cs",k=4)+
              s(min.abs.hum,by=time,bs="cs",k=4)+
              s(mean.wetness,by=time,bs="cs",k=4)+
              s(tot.rain,by=time,bs="cs",k=4)+
              s(mean.solar,by=time,bs="cs",k=4)
            ,data=delta.height)
  summary(mod0) #indicates that max.temp, mean.abs.hum, max.abs.hum, min.abs.hum not significant
  
  mod1<-gam(height.next~s(height,by=time,bs="cs",k=4)+
                      #s(inf.intens,by=time,bs="cs",k=4)+
                      s(mean.temp,by=time,bs="cs",k=4)+
                      #s(max.temp,by=time,bs="cs",k=4)+
                      s(min.temp,by=time,bs="cs",k=4)+
                      #s(mean.abs.hum,by=time,bs="cs",k=4)+
                      #s(max.abs.hum,by=time,bs="cs",k=4)+
                      #s(min.abs.hum,by=time,bs="cs",k=4)+
                      s(mean.wetness,by=time,bs="cs",k=4)+
                      s(tot.rain,by=time,bs="cs",k=4)+
                      s(mean.solar,by=time,bs="cs",k=4)
                    ,data=delta.height)
  summary(mod1)
  saveRDS(mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
draw(plant.growth.model) #plot smooths
#gam.check(plant.growth.model) #indicates no more knots needed
