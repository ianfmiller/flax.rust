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
                      s(min.abs.hum,by=time,bs="cs",k=4)+
                      #s(tot.rain,by=time,bs="cs",k=4)+
                      s(mean.solar,by=time,bs="cs",k=4)
                    ,data=delta.height)
  summary(mod1) #indicates that min.temp is not significant
  
  mod2<-gam(height.next~s(height,by=time,bs="cs",k=4)+
              #s(inf.intens,by=time,bs="cs",k=4)+
              s(mean.temp,by=time,bs="cs",k=4)+
              #s(max.temp,by=time,bs="cs",k=4)+
              #s(min.temp,by=time,bs="cs",k=4)+
              #s(mean.abs.hum,by=time,bs="cs",k=4)+
              #s(max.abs.hum,by=time,bs="cs",k=4)+
              s(min.abs.hum,by=time,bs="cs",k=4)+
              #s(tot.rain,by=time,bs="cs",k=4)+
              s(mean.solar,by=time,bs="cs",k=4)
            ,data=delta.height)
  summary(mod2) #mean temp is marginally significant, comparison w/ model resulting from cutting marginally significant predictors (mean temp and then min.min.abs.hum) explains less deviance and has slightly lower AIC
  
  saveRDS(mod2,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
draw(plant.growth.model) #plot smooths
#gam.check(plant.growth.model) #indicates no more knots needed

# model vis
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
library("MASS")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,80),ylim=c(-50,50),type="n",xlab="log 10 plant infection intensity",ylab="pred. change in plant infection intensity",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(10,80,10))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plant.growth.model(day,dummy.data.height=i)
    Xp <- predict(plant.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta <- coef(plant.growth.model) ## posterior mean of coefs
    Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
    n <- 10000
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plant.growth.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)))
    points(i,mean(y),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")
