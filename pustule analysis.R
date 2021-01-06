library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

source("~/Documents/GitHub/flax.rust/pustule analysis data prep.R")
# visualize data

## histograms
par(mfrow=c(2,1))
hist(delta.pustules$area,main="pustule area",breaks=100,xlab="area")
hist(delta.pustules$area.next-delta.pustules$area,main="change in pustule area",breaks=100,xlab="change in area")

## plot trajectories
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="date",ylab="pustule area")
i<-0

plot.cols<-sample(rainbow(2042))

for (tag in unique(pustules$tag))
{
  sub.pustules1<-pustules[which(pustules$tag==tag),]
  
  for (color in unique(sub.pustules1$color))
  {
    sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
    {
      sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
      
      for(pustule.number in unique(sub.pustules3$pustule.number))
      {
        i<-i+1
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        points(sub.pustules4$date,sub.pustules4$area,col=plot.cols[i],type="l",lwd=.5)
      }
    }
  }
}

## plot change
plot(delta.pustules$area,delta.pustules$area.next,col="grey",xlab = "area",ylab="next obs. area")
abline(0,1)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/data/models/pustule.model.RDS"))
{
  ### construct all combinations of predictors
  source("model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x],"(1|tag) + (1|who.measured)"),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- -17420 #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-lmer(model.set[[i]],data=delta.pustules))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  candidate.models<-unname(which(delta.AICs<4))
  model.set[candidate.models[order(AICs[candidate.models])]] #models to consider--offset(diam.last) not shown
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/data/models/pustule.model.RDS")
}

## load best model

pustule.model<-readRDS("~/Documents/GitHub/flax.rust/data/models/pustule.model.RDS")

## visualize model

par(mfrow=c(1,1))
plot(delta.pustules$area,delta.pustules$area.next-delta.pustules$area)

quant.temp.7.30.dew.point.days<-quantile(delta.pustules$temp.7.30.dew.point.days,.5)
quant.tot.rain<-quantile(delta.pustules$tot.rain,.5)
quant.gust.speed.days<-quantile(delta.pustules$gust.speed.days,.5)

curve.col<-"blue"
curve(fixef(pustule.model)["(Intercept)"]+
        fixef(pustule.model)["area"]*x+
        fixef(pustule.model)["temp.7.30.dew.point.days"]*quant.temp.7.30.dew.point.days+
        fixef(pustule.model)["tot.rain"]*quant.tot.rain+
        fixef(pustule.model)["gust.speed.days"]*quant.gust.speed.days
      ,add=T,col=curve.col)

