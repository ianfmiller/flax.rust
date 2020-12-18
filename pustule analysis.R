library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

source("prep.enviro.data.R")

# clean data

## load data
pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv")

## add in area
area<-pi*(pustules$max.diam..mm./4)*(pustules$min.diam..mm./4)
pustules$date<-as.POSIXct(pustules$date,tryFormats = "%m/%d/%y %H:%M",tz="UTC")
pustules$area<-area

## clean data
pustules<-pustules[which(pustules$pustule.ID.confidence=="Yes"),]
pustules<-pustules[which(pustules$area>0),]

## cut out incomplete records
pustules<-pustules[-intersect(which(pustules$site=="HM"),which(pustules$date>"2020-07-10 00:00:00 UTC")),]

## make new data object for change in pustule size
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

tags<-c()
stem.iters<-c()
leaf.iters<-c()
pustule.nums<-c()
start.vals<-c()
end.vals<-c()
days<-c()
temp.days<-c()
temp.days.16.22<-c()
temp.days.7.30<-c()
dew.point.days<-c()
temp.dew.point.days<-c()
temp.16.22.dew.point.days<-c()
temp.7.30.dew.point.days<-c()
wetness.days<-c()
temp.wetness.days<-c()
temp.16.22.wetness.days<-c()
temp.7.30.wetness.days<-c()
tot.rains<-c()
solar.days<-c()
wind.speed.days<-c()
gust.speed.days<-c()
measurer.ids<-c()


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
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        if(dim(sub.pustules4)[1]>=2)
        {
          if(length(unique(sub.pustules4$date))<dim(sub.pustules4)[1]) #get rid of duplicate data resulting from two pics of same pustule on same date
          {
            for(date in unique(sub.pustules4$date))
            {
              date.indicies<-which(sub.pustules4$date==date)
              if(length(date.indicies)>1) {sub.pustules4<-sub.pustules4[-date.indicies[2:length(date.indicies)],]}
            }
          }
          for(i in 1:(dim(sub.pustules4)[1]-1))
          {
            #pull reference data
            date0<-sub.pustules4[i,"date"]
            date1<-sub.pustules4[i+1,"date"]
            site<-sub.pustules4[i,"site"]
            
            #subset temp data to relevant window
            temp.rh.sub<-subset(all.temp.rh,site==site) #pull out temp data for site
            temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #pull out relevant data
            temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #pull out relevant data
            temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #throw out NAs
            temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
            
            #subset weather data to relevant window
            weath.sub<-subset(all.weath,site==site)
            weath.sub<-subset(weath.sub,date<=date1) #pull out relevant data
            weath.sub<-subset(weath.sub,date>=date0) #pull out relevant data
            weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
            
            #calculate environmental variable metrics
            new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
            new.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #time (in days) during which temp between 16 and 22 celsius
            new.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #time (in days) during which temp between 7 and 30 celsius
            new.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T) #Dew point days

            #calculate weather metrics
            new.wetness.days<-sum(weath.sub$wetness*weath.sub$interval.length,na.rm = T)
            new.tot.rain<-sum(weath.sub$rain,na.rm=T)
            new.solar.days<-sum(weath.sub$solar.radiation*weath.sub$interval.length,na.rm = T)
            new.wind.speed.days<-sum(weath.sub$wind.speed*weath.sub$interval.length,na.rm = T)
            new.gust.speed.days<-sum(weath.sub$wind.direction*weath.sub$interval.length,na.rm = T)
            
            #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
            new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
            new.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
            new.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)
            new.temp.wetness.days<-sum(weath.sub$temp*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
            new.temp.16.22.wetness.days<-sum(weath.sub$temp.16.22*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
            new.temp.7.30.wetness.days<-sum(weath.sub$temp.7.30*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
            
            #pull out core predictors
            start.val<-sub.pustules4[i,"area"]
            end.val<-sub.pustules4[i+1,"area"]
            delta.days<-date1-date0
            measurer.id<-sub.pustules4[i,"who.entered"] #same person did all measurements for each pustule
            
            #store values
            tags<-c(tags,tag)
            stem.iters<-c(stem.iters,color)
            leaf.iters<-c(leaf.iters,leaf.iteration)
            pustule.nums<-c(pustule.nums,pustule.number)
            
            start.vals<-c(start.vals,start.val)
            end.vals<-c(end.vals,end.val)
            days<-c(days,delta.days)
            
            temp.days.16.22<-c(temp.days.16.22,new.temp.days.16.22)
            temp.days.7.30<-c(temp.days.7.30,new.temp.days.7.30)
            temp.days<-c(temp.days,new.temp.days)
            dew.point.days<-c(dew.point.days,new.dew.point.days)
            temp.dew.point.days<-c(temp.dew.point.days,new.temp.dew.point.days)
            temp.16.22.dew.point.days<-c(temp.16.22.dew.point.days,new.temp.16.22.dew.point.days)
            temp.7.30.dew.point.days<-c(temp.7.30.dew.point.days,new.temp.7.30.dew.point.days)
            
            wetness.days<-c(wetness.days,new.wetness.days)
            temp.wetness.days<-c(temp.wetness.days,new.temp.wetness.days)
            temp.16.22.wetness.days<-c(temp.16.22.wetness.days,new.temp.16.22.wetness.days)
            temp.7.30.wetness.days<-c(temp.7.30.wetness.days,new.temp.7.30.wetness.days)
            tot.rains<-c(tot.rains,new.tot.rain)
            solar.days<-c(solar.days,new.solar.days)
            wind.speed.days<-c(wind.speed.days,new.wind.speed.days)
            gust.speed.days<-c(gust.speed.days,new.gust.speed.days)
            
            measurer.ids<-c(measurer.ids,measurer.id)
            
          } 
        }
      }
    }
  }
}

delta.pustules<-data.frame(tag=factor(tags),stem.iter=stem.iters,leaf.iter=leaf.iters,pustule.num=pustule.nums,area=start.vals,area.next=end.vals,time=days,
                           temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,
                           dew.point.days=dew.point.days,temp.dew.point.days=temp.dew.point.days,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,
                           wetness.days=wetness.days,temp.wetness.days=temp.wetness.days,temp.16.22.wetness.days=temp.16.22.wetness.days,temp.7.30.wetness.days=temp.7.30.wetness.days,
                           tot.rain=tot.rains,solar.days=solar.days,wind.speed.days=wind.speed.days,gust.speed.days=gust.speed.days,who.measured=measurer.ids)

# visualize data

## histograms
par(mfrow=c(2,1))
hist(delta.pustules$area,main="pustule area",breaks=100,xlab="area")
hist(delta.pustules$area.next-delta.pustules$area,main="pustule area",breaks=100,xlab="change in area")

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
model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x]),collapse=" + ")))

re.model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x],"(1|tag) + (1|who.measured)"),collapse=" + ")))

#gam.model.set <- apply(pred.mat, 1, function(x) as.formula( paste0(paste0(paste(c("area.next ~ offset(area",predictors[x]),collapse=") + s("),")")," + s(tag,bs='re')")))
  
names(model.set)<-seq(1,length(model.set),1)
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
#AIC.benchmark<-AIC(lmer(area.next~offset(area)+area +(1|tag) + (1|who.measured),data=delta.pustules))
#AIC.benchmark<- -17555
#pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
#for (i in 1:length(re.model.set))
#{
#  suppressMessages(new.mod<-lmer(re.model.set[[i]],data=delta.pustules,REML=F))
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
best.lm.model<-lm(area.next~offset(area)+area+time+temp.days.7.30+dew.point.days+temp.7.30.dew.point.days+wetness.days+temp.wetness.days+wind.speed.days+gust.speed.days,data=delta.pustules)
summary(best.lm.model)

par(mfrow=c(1,1))
plot(delta.pustules$area,delta.pustules$area.next-delta.pustules$area)
quant.time<-quantile(delta.pustules$time,.5)
quant.temp.days.7.30<-quantile(delta.pustules$temp.days.7.30,.5)
quant.dew.point.days<-quantile(delta.pustules$dew.point.days,.5)
quant.temp.7.30.dew.point.days<-quantile(delta.pustules$temp.7.30.dew.point.days,.5)
quant.wetness.days<-quantile(delta.pustules$wetness.days,.5)
quant.temp.wetness.days<-quantile(delta.pustules$temp.wetness.days,.5)
quant.wind.speed.days<-quantile(delta.pustules$wind.speed.days,.5)
quant.gust.speed.days<-quantile(delta.pustules$gust.speed.days,.5)

curve.col<-"red"
curve(best.lm.model$coefficients["(Intercept)"]+
      best.lm.model$coefficients["area"]*x+
      best.lm.model$coefficients["time"]*quant.time+
      best.lm.model$coefficients["temp.days.7.30"]*quant.temp.days.7.30+
      best.lm.model$coefficients["dew.point.days"]*quant.dew.point.days+
      best.lm.model$coefficients["temp.7.30.dew.point.days"]*quant.temp.7.30.dew.point.days+
      best.lm.model$coefficients["wetness.days"]*quant.wetness.days+
      best.lm.model$coefficients["temp.wetness.days"]*quant.temp.wetness.days+
      best.lm.model$coefficients["wind.speed.days"]*quant.wind.speed.days+
      best.lm.model$coefficients["gust.speed.days"]*quant.gust.speed.days
      ,add=T,col=curve.col)


### lmms
## if all models fit
#re.AICs<-unlist(lapply(re.all.fit.models,AIC))
#delta.re.AICs<-re.AICs-min(re.AICs)
#re.candidate.models<-unname(which(delta.re.AICs<4))
#re.model.set[re.candidate.models[order(re.AICs[re.candidate.models])]] #models to consider--offset(diam.last) not shown
#index<-1
#best.model<-re.all.fit.models[[order(re.AICs)[index]]]

##fitting just top model
best.lmer.model<-lmer(area.next~offset(area)+area+temp.days.16.22+dew.point.days+temp.16.22.dew.point.days +wetness.days+temp.wetness.days+tot.rain+gust.speed.days+(1|tag)+(1|who.measured),data=delta.pustules)
summary(best.lmer.model)

#par(mfrow=c(1,1))
#plot(delta.pustules$area,delta.pustules$area.next-delta.pustules$area)
quant.temp.days.16.22<-quantile(delta.pustules$temp.days.16.22,.5)
quant.dew.point.days<-quantile(delta.pustules$dew.point.days,.5)
quant.temp.16.22.dew.point.days<-quantile(delta.pustules$temp.16.22.dew.point.days,.5)
quant.wetness.days<-quantile(delta.pustules$wetness.days,.5)
quant.temp.wetness.days<-quantile(delta.pustules$temp.wetness.days,.5)
quant.tot.rain<-quantile(delta.pustules$tot.rain,.5)
quant.gust.speed.days<-quantile(delta.pustules$gust.speed.days,.5)

curve.col<-"blue"
curve(fixef(best.lmer.model)["(Intercept)"]+
        fixef(best.lmer.model)["area"]*x+
        fixef(best.lmer.model)["temp.days.16.22"]*quant.temp.days.16.22+
        fixef(best.lmer.model)["dew.point.days"]*quant.dew.point.days+
        fixef(best.lmer.model)["temp.16.22.dew.point.days"]*quant.temp.16.22.dew.point.days+
        fixef(best.lmer.model)["wetness.days"]*quant.wetness.days+
        fixef(best.lmer.model)["temp.wetness.days"]*quant.temp.wetness.days+
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

