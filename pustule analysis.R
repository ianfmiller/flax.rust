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

## fit models 

### construct all combinations of predictors
source("model.set.creation.R")

### create all sets of models
#model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x]),collapse=" + ")))
re.model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x],"(1|tag) + (1|who.measured)"),collapse=" + ")))
#gam.model.set <- apply(pred.mat, 1, function(x) as.formula( paste0(paste0(paste(c("area.next ~ offset(area",predictors[x]),collapse=") + s("),")")," + s(tag,bs='re')")))
  
#names(model.set)<-seq(1,length(model.set),1)
names(re.model.set)<-seq(1,length(re.model.set),1)
#names(gam.model.set)<-seq(1,length(gam.model.set),1)

## run to search for best lm model
#all.fit.models<-list()
#AIC.benchmark<-AIC(lm(area.next~offset(area)+area,data=delta.pustules))
#pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
#for (i in 1:length(model.set))
#{
#  new.mod<-lm(model.set[[i]],data=delta.pustules)
#  AIC.new.mod<-AIC(new.mod)
#  if(AIC.new.mod<=(AIC.benchmark+10)) {all.fit.models<-append(all.fit.models,list(new.mod))}
#  pb$tick()
#}

## run to search for best lmer model
#re.all.fit.models<-c()
#AIC.benchmark<-AIC(lmer(area.next~offset(area)+area +(1|tag),data=delta.pustules))
#AIC.benchmark<- -17400
#pb <- progress_bar$new(total = length(re.model.set),format = " fitting models [:bar] :percent eta: :eta")
#for (i in 1:length(re.model.set))
#{
#  suppressMessages(new.mod<-lmer(re.model.set[[i]],data=delta.pustules))
#  AIC.new.mod<-AIC(new.mod)
#  if(AIC.new.mod<=(AIC.benchmark)) {re.all.fit.models<-append(re.all.fit.models,list(new.mod))}
#  pb$tick()
#}

## run to search for best gam model
#gam.all.fit.models<-c()
#AIC.benchmark<-AIC(gam(area.next~offset(area)+s(area)+s(tag,bs="re"),data=delta.pustules,method="GCV.Cp"))
#AIC.benchmark<-(-17813)
#pb <- progress_bar$new(total = length(gam.model.set),format = " fitting models [:bar] :percent eta: :eta")
#for (i in 1:length(gam.model.set))
#{
#  suppressMessages(new.mod<-gam(gam.model.set[[i]],data=delta.pustules,method="GCV.Cp"))
#  AIC.new.mod<-AIC(new.mod)
#  if(AIC.new.mod<=(AIC.benchmark)) {gam.all.fit.models<-append(gam.all.fit.models,list(new.mod))}
#  pb$tick()
#}

## compare between models

### lms
#AICs<-unlist(lapply(all.fit.models,AIC))
#delta.AICs<-AICs-min(AICs)
#candidate.models<-unname(which(delta.AICs<4))
#model.set[candidate.models[order(AICs[candidate.models])]] #models to consider--offset(diam.last) not shown
#index<-1 #top two models have nearly identical AICs, #2 drops insignificant wind speed days predictor
#best.lm.model<-all.fit.models[[order(AICs)[index]]]
#best.lm.model<-lm(area.next~offset(area)+area+time+temp.days.7.30+dew.point.days+temp.7.30.dew.point.days+wetness.days+temp.wetness.days+wind.speed.days+gust.speed.days,data=delta.pustules)
#summary(best.lm.model)

#par(mfrow=c(1,1))
#plot(delta.pustules$area,delta.pustules$area.next-delta.pustules$area)
#quant.time<-quantile(delta.pustules$time,.5)
#quant.temp.days.7.30<-quantile(delta.pustules$temp.days.7.30,.5)
#quant.dew.point.days<-quantile(delta.pustules$dew.point.days,.5)
#quant.temp.7.30.dew.point.days<-quantile(delta.pustules$temp.7.30.dew.point.days,.5)
#quant.wetness.days<-quantile(delta.pustules$wetness.days,.5)
#quant.temp.wetness.days<-quantile(delta.pustules$temp.wetness.days,.5)
#quant.wind.speed.days<-quantile(delta.pustules$wind.speed.days,.5)
#quant.gust.speed.days<-quantile(delta.pustules$gust.speed.days,.5)

#curve.col<-"red"
#curve(best.lm.model$coefficients["(Intercept)"]+
#      best.lm.model$coefficients["area"]*x+
#      best.lm.model$coefficients["time"]*quant.time+
#      best.lm.model$coefficients["temp.days.7.30"]*quant.temp.days.7.30+
#      best.lm.model$coefficients["dew.point.days"]*quant.dew.point.days+
#      best.lm.model$coefficients["temp.7.30.dew.point.days"]*quant.temp.7.30.dew.point.days+
#      best.lm.model$coefficients["wetness.days"]*quant.wetness.days+
#      best.lm.model$coefficients["temp.wetness.days"]*quant.temp.wetness.days+
#      best.lm.model$coefficients["wind.speed.days"]*quant.wind.speed.days+
#     best.lm.model$coefficients["gust.speed.days"]*quant.gust.speed.days
#     ,add=T,col=curve.col)


### lmms
## if all models fit
#re.AICs<-unlist(lapply(re.all.fit.models,AIC))
#delta.re.AICs<-re.AICs-min(re.AICs)
#re.candidate.models<-unname(which(delta.re.AICs<4))
#re.model.set[re.candidate.models[order(re.AICs[re.candidate.models])]] #models to consider--offset(diam.last) not shown
#index<-1
#best.model<-re.all.fit.models[[order(re.AICs)[index]]]
if(!file.exists("~/Documents/GitHub/flax.rust/data/models/pustule.lmer.RDS")) {saveRDS(best.model,file="~/Documents/GitHub/flax.rust/data/models/pustule.lmer.RDS")}

##fitting just top model
best.lmer.model<-lmer(area.next~offset(area)+area+temp.7.30.dew.point.days+tot.rain+gust.speed.days+(1|tag),data=delta.pustules)
summary(best.lmer.model)

#par(mfrow=c(1,1))
#plot(delta.pustules$area,delta.pustules$area.next-delta.pustules$area)
quant.temp.7.30.dew.point.days<-quantile(delta.pustules$temp.7.30.dew.point.days,.5)
quant.tot.rain<-quantile(delta.pustules$tot.rain,.5)
quant.gust.speed.days<-quantile(delta.pustules$gust.speed.days,.5)

par(mfrow=c(1,1))
plot(delta.pustules$area,delta.pustules$area.next-delta.pustules$area)

curve.col<-"blue"
curve(fixef(best.lmer.model)["(Intercept)"]+
        fixef(best.lmer.model)["area"]*x+
        fixef(best.lmer.model)["temp.7.30.dew.point.days"]*quant.temp.7.30.dew.point.days+
        fixef(best.lmer.model)["tot.rain"]*quant.tot.rain+
        fixef(best.lmer.model)["gust.speed.days"]*quant.gust.speed.days
      ,add=T,col=curve.col)

### gams
## if all models fit
#gam.AICs<-unlist(lapply(gam.all.fit.models,AIC))
#delta.gam.AICs<-gam.AICs-min(gam.AICs)
#gam.candidate.models<-unname(which(delta.gam.AICs<4))
#gam.model.set[gam.candidate.models[order(delta.gam.AICs[gam.candidate.models])]] #models to consider--offset(diam.last) not shown
#index<-1
#best.model<-gam.all.fit.models[[order(gam.AICs)[index]]]

##fitting just top model
#best.gam.model<-gam(area.next~offset(area)+s(area)+s(temp.days)+s(temp.7.30.dew.point.days)+s(tot.rain)+s(solar.days)+s(wind.speed.days)+s(gust.speed.days)+s(tag,bs="re"),data=delta.pustules,method="GCV.Cp")
#summary(best.gam.model)

#pred.areas<-seq(0,.5,.001)
#new.preds<-c()
#quant.temp.days<-quantile(delta.pustules$temp.days,.5)
#quant.temp.7.30.dew.point.days<-quantile(delta.pustules$temp.7.30.dew.point.days,.5)
#quant.tot.rain<-quantile(delta.pustules$tot.rain,.5)
#quant.solar.days<-quantile(delta.pustules$solar.days,.5)
#quant.wind.speed.days<-quantile(delta.pustules$wind.speed.days,.5)
#quant.gust.speed.days<-quantile(delta.pustules$gust.speed.days,.5)
#for(i in 1:length(pred.areas))
#{
#  new.preds<-c(new.preds,predict(best.gam.model,
#                                 newdata = data.frame(area=pred.areas[i],temp.days=quant.temp.days,temp.7.30.dew.point.days=quant.temp.7.30.dew.point.days,tot.rain=quant.tot.rain,solar.days=quant.solar.days,wind.speed.days=quant.wind.speed.days,gust.speed.days=quant.gust.speed.days),   
#                                 type="response"))
#}
#points(pred.areas,new.preds-pred.areas,type="l",col="orange")

