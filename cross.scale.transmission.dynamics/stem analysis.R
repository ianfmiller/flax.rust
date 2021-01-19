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
par(mfrow=c(1,1))
plot(delta.stems$stem.inf.intens,delta.stems$stem.inf.intens.next,col="grey",xlab = "stem infection intensity",ylab="next obs. stem infection intensity")
abline(0,1)

## plot log10 change
### relationship doesn't appear obviously linear near 0--take polynomial fitting approach
par(mfrow=c(1,1))
plot(log10(delta.stems$stem.inf.intens),log10(delta.stems$stem.inf.intens.next),col="grey",xlab = "stem infection intensity",ylab="next obs. stem infection intensity")
abline(0,1)

## test different ordered polynomial fits
plot(log10(delta.stems$stem.inf.intens),log10(delta.stems$stem.inf.intens.next))
legend("bottomright",legend = c("linear","quadratic","cubic","quartic"),lty=1,col=c("blue","red","green","purple"))
test.mod.1<-lmer(log10(stem.inf.intens.next)~log10(stem.inf.intens)+(1|tag),data=delta.stems,REML = F)
abline(fixef(test.mod.1)[1],fixef(test.mod.1)[1],col="blue")

test.mod.2<-lmer(log10(stem.inf.intens.next)~poly(log10(stem.inf.intens),2,raw = T)+(1|tag),data=delta.stems,REML = F)
curve(fixef(test.mod.2)[1]+fixef(test.mod.2)[2]*x+fixef(test.mod.2)[3]*x^2,add=T,col="red")

test.mod.3<-lmer(log10(stem.inf.intens.next)~poly(log10(stem.inf.intens),3,raw = T)+(1|tag),data=delta.stems,REML = F)
curve(fixef(test.mod.3)[1]+fixef(test.mod.3)[2]*x+fixef(test.mod.3)[3]*x^2+fixef(test.mod.3)[4]*x^3,add=T,col="green")

test.mod.4<-lmer(log10(stem.inf.intens.next)~poly(log10(stem.inf.intens),4,raw = T)+(1|tag),data=delta.stems,REML = F)
curve(fixef(test.mod.4)[1]+fixef(test.mod.4)[2]*x+fixef(test.mod.4)[3]*x^2+fixef(test.mod.4)[4]*x^3+fixef(test.mod.4)[5]*x^4,add=T,col="purple")

AIC(test.mod.1,test.mod.2,test.mod.3,test.mod.4) ### quadratic and cubic give similar AICs, go with quadratic

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/stems.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/stem.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("log10(stem.inf.intens.next) ~ poly(log10(stem.inf.intens),2,raw=T)",predictors[x],'(1|tag)'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(lmer(log10(stem.inf.intens.next)~poly(log10(stem.inf.intens),2,raw = T)+(1|tag),data=delta.stems,REML = F)) #cutoff to limit memory usage
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
plot(log10(delta.stems$stem.inf.intens),log10(delta.stems$stem.inf.intens.next))

quant.temp.7.30.wetness.days <-quantile(delta.stems$temp.7.30.wetness.days,.5)
quant.tot.rain<-quantile(delta.stems$tot.rain,.5)

curve.col<-"blue"
curve(fixef(stems.model)["(Intercept)"]+
        fixef(stems.model)["poly(log10(stem.inf.intens), 2, raw = T)1"]*x+
        fixef(stems.model)["poly(log10(stem.inf.intens), 2, raw = T)2"]*x^2+
        fixef(stems.model)["temp.7.30.wetness.days"]*quant.temp.7.30.wetness.days+
        fixef(stems.model)["tot.rain"]*quant.tot.rain
        ,add=T,col=curve.col)

