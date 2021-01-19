library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant data prep.R")
delta.plants<-subset(delta.plants,time<=10)

# visualize data

## histograms
par(mfrow=c(2,1))
hist(plants$plant.inf.intens,main="plant infection intensity",breaks=100,xlab="plant infection intensity")
hist(delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,main="change in plant infection intensity",breaks=100,xlab="change in plant infection intensity")

## plot trajectories
par(mfrow=c(1,1))
plot(c(min(as.Date(plants$Date,tryFormats = "%m/%d/%Y")),max(as.Date(plants$Date,tryFormats = "%m/%d/%Y"))),c(0,max(plants$plant.inf.intens)),type="n",xlab="date",ylab="plant infection intensity",main="plant infection intensity")
i<-0
plot.cols<-sample(rainbow(95))

for (tag in unique(plants$Tag))
{
  sub.plants.1<-plants[which(plants$Tag==tag),]
  points(sub.plants.1$Date,sub.plants.1$plant.inf.intens,col=plot.cols[i],type="l",lwd=.5)
  i<-i+1
}


## plot change
par(mfrow=c(2,1))
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next,col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity")
abline(0,1)
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,col="grey",xlab = "stem infection intensity",ylab="change in stem infection intensity")
abline(h=0)

## plot log10 change
### relationship doesn't appear obviously linear near 0--take polynomial fitting approach
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next),col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity")
abline(0,1)

## test different ordered polynomial fits
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next))
legend("bottomright",legend = c("linear","quadratic","cubic","quartic"),lty=1,col=c("blue","red","green","purple"))
test.mod.1<-lmer(log10(plant.inf.intens.next)~log10(plant.inf.intens)+(1|tag),data=delta.plants,REML = F)
abline(fixef(test.mod.1)[1],fixef(test.mod.1)[1],col="blue")

test.mod.2<-lmer(log10(plant.inf.intens.next)~poly(log10(plant.inf.intens),2,raw = T)+(1|tag),data=delta.plants,REML = F)
curve(fixef(test.mod.2)[1]+fixef(test.mod.2)[2]*x+fixef(test.mod.2)[3]*x^2,add=T,col="red")

test.mod.3<-lmer(log10(plant.inf.intens.next)~poly(log10(plant.inf.intens),3,raw = T)+(1|tag),data=delta.plants,REML = F)
curve(fixef(test.mod.3)[1]+fixef(test.mod.3)[2]*x+fixef(test.mod.3)[3]*x^2+fixef(test.mod.3)[4]*x^3,add=T,col="green")

test.mod.4<-lmer(log10(plant.inf.intens.next)~poly(log10(plant.inf.intens),4,raw = T)+(1|tag),data=delta.plants,REML = F)
curve(fixef(test.mod.4)[1]+fixef(test.mod.4)[2]*x+fixef(test.mod.4)[3]*x^2+fixef(test.mod.4)[4]*x^3+fixef(test.mod.4)[5]*x^4,add=T,col="purple")

AIC(test.mod.1,test.mod.2,test.mod.3,test.mod.4) ### quadratic and cubic give similar AICs, go with quadratic


# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("log10(plant.inf.intens.next) ~ poly(log10(plant.inf.intens),2,raw=T)",predictors[x],'(1|tag)'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(lmer(log10(plant.inf.intens.next)~poly(log10(plant.inf.intens),2,raw=T)+(1|tag),data=delta.plants,REML = F)) #cutoff to limit memory usage
  #AIC.benchmark<- 8025
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-lmer(model.set[[i]],data=delta.plants,REML = F))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")
}

## load best model

plants.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")

## visualize model
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next))

quant.dew.point.days<-quantile(delta.plants$dew.point.days,.5)
quant.temp.7.30.dew.point.days <-quantile(delta.plants$temp.7.30.dew.point.days ,.5)

curve.col<-"blue"
curve(fixef(plants.model)["(Intercept)"]+
        fixef(plants.model)["poly(log10(plant.inf.intens), 2, raw = T)1"]*x+
        fixef(plants.model)["poly(log10(plant.inf.intens), 2, raw = T)2"]*x^2+
        fixef(plants.model)["dew.point.days"]*quant.dew.point.days+
        fixef(plants.model)["temp.7.30.dew.point.days"]*quant.temp.7.30.dew.point.days
      ,add=T,col=curve.col)

