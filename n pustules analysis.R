library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/n pustules data prep.R")

delta.n.pustules<-subset(delta.n.pustules,time<=7)

# visualize data

## histograms
par(mfrow=c(2,1))
hist(delta.n.pustules$n.pustules,main="n pustules",breaks=100,xlab="n pustules")
hist(delta.n.pustules$n.pustules.next-delta.n.pustules$n.pustules,main="change in n pustules",breaks=100,xlab="change in n pustules")

## plot trajectories
par(mfrow=c(1,1))
plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,max(n.pustules$N.pustules)),type="n",xlab="date",ylab="pustule area")
#plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,50),type="n",xlab="date",ylab="pustule area")

i<-0

plot.cols<-sample(rainbow(300))

for (tag in unique(n.pustules$tag))
{
  sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
  
  for (color in unique(sub.n.pustules1$color))
  {
    sub.n.pustules2<-sub.n.pustules1[which(sub.n.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.n.pustules2$leaf.iteration)) 
    {
      sub.n.pustules3<-sub.n.pustules2[which(sub.n.pustules2$leaf.iteration==leaf.iteration),]
      
      i<-i+1
      points(sub.n.pustules3$date,sub.n.pustules3$N.pustules,col=plot.cols[i],type="l",lwd=.5)
      
    }
  }
}

## plot change
par(mfrow=c(1,2))
plot(delta.n.pustules$n.pustules,delta.n.pustules$n.pustules.next,col="grey",xlab = "area",ylab="next obs. area")
abline(0,1)
plot(delta.n.pustules$n.pustules,delta.n.pustules$n.pustules.next,col="grey",xlab = "area",ylab="next obs. area",xlim=c(0,20),ylim=c(0,20))
abline(0,1)


# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/data/models/n.pustule.model.RDS"))
{
  ### construct all combinations of predictors
  source("n.pustules.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("n.pustules.next ~ offset(n.pustules)",predictors[x],'(1|tag)'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(lmer(n.pustules.next~offset(n.pustules)+n.pustules+(1|tag),data=delta.n.pustules)) #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-lmer(model.set[[i]],data=delta.n.pustules))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/data/models/n.pustules.model.RDS")
}

## load best model

n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/data/models/n.pustules.model.RDS")

## visualize model

par(mfrow=c(1,1))
plot(delta.n.pustules$time,delta.n.pustules$n.pustules.next-delta.n.pustules$n.pustules)

quant.0.25.temp.days.16.22<-quantile(delta.n.pustules$temp.days.16.22,.25)
quant.0.5.temp.days.16.22<-quantile(delta.n.pustules$temp.days.16.22,.5)
quant.0.75.temp.days.16.22<-quantile(delta.n.pustules$temp.days.16.22,.75)

curve(fixef(n.pustules.model)["(Intercept)"]+fixef(n.pustules.model)["time"]*x+fixef(n.pustules.model)["temp.days.16.22"]*quant.0.25.temp.days.16.22,add=T,col="yellow")
curve(fixef(n.pustules.model)["(Intercept)"]+fixef(n.pustules.model)["time"]*x+fixef(n.pustules.model)["temp.days.16.22"]*quant.0.5.temp.days.16.22,add=T,col="orange")
curve(fixef(n.pustules.model)["(Intercept)"]+fixef(n.pustules.model)["time"]*x+fixef(n.pustules.model)["temp.days.16.22"]*quant.0.75.temp.days.16.22,add=T,col="red")

legend("topright",legend=c("25% quantile temp.days.16.22","50% quantile temp.days.16.22","75% quantile temp.days.16.22"),lty=1,col=c("yellow","orange","red"))
