library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/stem data prep.R")
delta.stems<-subset(delta.stems,time<=10)
#delta.stems<-delta.stems[-which(delta.stems$end.stems>delta.stems$end.length),]
# visualize data

## histograms
par(mfrow=c(2,1))
hist(stems$stem.inf.intens,main="stem infection intensity",breaks=100,xlab="stem infection intensity")
hist(delta.stems$stem.inf.intens.next-delta.stems$stem.inf.intens,main="change in stem infection intensity",breaks=100,xlab="change in stem infection intensity")

## plot trajectories
par(mfrow=c(1,1))
plot(c(min(as.Date(stems$Date,tryFormats = "%m/%d/%Y")),max(as.Date(stems$Date,tryFormats = "%m/%d/%Y"))),c(0,2000),type="n",xlab="date",ylab="stem infection intensity",main="stem infection intensity")
i<-0
plot.cols<-sample(rainbow(238))

for (tag in unique(stems$Tag))
{
  sub.stems1<-stems[which(stems$Tag==tag),]
  
  for (stem.index in unique(sub.stems1$stem.index))
  {
    sub.stems2<-sub.stems1[which(sub.stems1$stem.index==stem.index),]
    if(any(is.na(sub.stems2$stem.inf.intens))) {sub.stems2<-sub.stems2[-which(is.na(sub.stems2$length.tissue.infected)),]}
    points(sub.stems2$Date,sub.stems2$stem.inf.intens,col=plot.cols[i],type="l",lwd=.5)
    i<-i+1
  }
}


## plot change
par(mfrow=c(2,1))
plot(delta.stems$stem.inf.intens,delta.stems$stem.inf.intens.next,col="grey",xlab = "stem infection intensity",ylab="next obs. stem infection intensity")
abline(0,1)
plot(delta.stems$stem.inf.intens,delta.stems$stem.inf.intens.next-delta.stems$stem.inf.intens,col="grey",xlab = "stem infection intensity",ylab="change in stem infection intensity")
abline(h=0)
# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/stems.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/stem.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("stem.inf.intens.next ~ offset(stem.inf.intens) + stem.inf.intens",predictors[x],'(1|tag)'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(lmer(stem.inf.intens.next~offset(stem.inf.intens)+stem.inf.intens+(1|tag),data=delta.stems,REML = F)) #cutoff to limit memory usage
  #AIC.benchmark<- 8025
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-lmer(model.set[[i]],data=delta.stems,REML = F))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/stems.model.RDS")
}

## load best model

stems.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/stems.model.RDS")

## visualize model
par(mfrow=c(1,1))
plot(delta.stems$stem.inf.intens,delta.stems$stem.inf.intens.next-delta.stems$stem.inf.intens)

quant.temp.days.16.22<-quantile(delta.stems$temp.days.16.22,.5)
quant.temp.16.22.dew.point.days<-quantile(delta.stems$temp.16.22.dew.point.days,.5)
quant.temp.7.30.wetness.days <-quantile(delta.stems$temp.16.22.wetness.days ,.5)
quant.tot.rain  <-quantile(delta.stems$tot.rain  ,.5)

curve.col<-"blue"
curve(fixef(stems.model)["(Intercept)"]+
        fixef(stems.model)["stem.inf.intens"]*x+
        fixef(stems.model)["temp.days.16.22"]*quant.temp.days.16.22+
        fixef(stems.model)["temp.16.22.dew.point.days"]*quant.temp.16.22.dew.point.days+
        fixef(stems.model)["temp.7.30.wetness.days"]*quant.temp.7.30.wetness.days+
        fixef(stems.model)["tot.rain"]*quant.tot.rain
        ,add=T,col=curve.col)

