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
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next,col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity",xlim=c(0,1000),ylim=c(0,200))
abline(0,1)
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,col="grey",xlab = "plant infection intensity",ylab="change in plant infection intensity")
abline(h=0)

## plot log10 change
### relationship doesn't appear obviously linear near 0--take polynomial fitting approach
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next),col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity")
abline(0,1)

## test gam fit
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next),col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity")
test.mod<-gam(log10(plant.inf.intens.next)~s(log10(plant.inf.intens))+s(site,bs="re"),data = delta.plants)
points(log10(delta.plants$plant.inf.intens),predict(test.mod,exclude = 's(site)'),col="red")


# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("log10(plant.inf.intens.next) ~ s(log10(plant.inf.intens))",predictors[x],'s(site,bs="re")'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(gam(log10(plant.inf.intens.next)~s(log10(plant.inf.intens))+s(site,bs="re"),data = delta.plants)) #cutoff to limit memory usage
  #AIC.benchmark<- 8025
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-gam(model.set[[i]],data=delta.plants,REML = F))
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

new.data.plant.inf.intens<-10^seq(-1,5,.01)
new.data.dew.point.days<-rep(quantile(delta.plants$dew.point.days,.5),times=length(new.data.log10.plant.inf.intens))
new.data.temp.7.30.dew.point.days<-rep(quantile(delta.plants$temp.7.30.dew.point.days,.5),times=length(new.data.log10.plant.inf.intens))
new.data.pred.pustule.diam.growth<-rep(quantile(delta.plants$pred.pustule.diam.growth,.5),times=length(new.data.log10.plant.inf.intens))
new.data.site<-rep("BT",times=length(new.data.log10.plant.inf.intens))

new.data<-data.frame("plant.inf.intens"=new.data.plant.inf.intens,"dew.point.days"=new.data.dew.point.days,"temp.7.30.dew.point.days"=new.data.temp.7.30.dew.point.days,"pred.pustule.diam.growth"=new.data.pred.pustule.diam.growth,"site"=new.data.site)
points(log10(new.data.plant.inf.intens),predict(plants.model,newdata = new.data,exclude = 's(site)'),type="l",col="red")

