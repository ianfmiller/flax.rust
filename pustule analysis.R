library(mgcv)
library(lme4)

source("prep.enviro.data.R")

# clean data

## load data
pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv")

## add in area
area<-pi*(pustules$max.diam..mm./4)*(pustules$min.diam..mm./4)
pustules$date<-as.Date(pustules$date,tryFormats = "%m/%d/%y")
pustules$area<-area

## clean data
pustules<-pustules[which(pustules$pustule.ID.confidence=="Yes"),]
pustules<-pustules[which(pustules$area>0),]

## make new data object for change in pustule size
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

tags<-c()
stem.iters<-c()
leaf.iters<-c()
pustule.nums<-c()
start.vals<-c()
end.vals<-c()
days<-c()
temp.days.16.22<-c()
temp.days.7.30<-c()
temp.days<-c()
mean.temp<-c()
dew.point.days<-c()
mean.dew.point<-c()
temp.dew.point.days<-c()
temp.16.22.dew.point.days<-c()
temp.7.30.dew.point.days<-c()
mean.wetnesss<-c()
tot.rains<-c()
mean.solars<-c()
mean.wind.speeds<-c()
mean.gust.speeds<-c()



for (tag in unique(pustules$tag)) #916, 917, 920
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
            
            #calculate environmental variable metrics
            new.temp.days.16.22<-sum(temp.rh.sub.func(temp.rh.sub,16,22)$temp.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #temperature days for temp between 16 and 22 celsius
            new.temp.days.7.30<-sum(temp.rh.sub.func(temp.rh.sub,7,30)$temp.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #temperature days for temp between 7 and 30 celsius
            new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
            new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #temperature days
            new.dew.point.days<-mean(temp.rh.sub$dew.pt.c,na.rm = T) #Dew point days
            new.mean.dew.point<-mean(temp.rh.sub$dew.pt.c,na.rm = T) #temperature days

            #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
            new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
            new.temp.16.22.dew.point.days<-sum(temp.rh.sub.func(temp.rh.sub,16,22)$temp.c*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
            new.temp.7.30.dew.point.days<-sum(temp.rh.sub.func(temp.rh.sub,7,30)$temp.c*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)
            
            #calculate weather metrics
            new.mean.wetness<-mean(weath.sub$wetness)
            new.tot.rain<-sum(weath.sub$rain)
            new.mean.solar<-mean(weath.sub$solar.radiation)
            new.mean.wind.speed<-mean(weath.sub$wind.speed)
            new.mean.gust.speed<-mean(weath.sub$gust.speed)
            
            #pull out core predictors
            start.val<-sub.pustules4[i,"area"]
            end.val<-sub.pustules4[i+1,"area"]
            delta.days<-date1-date0
            
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
            mean.temp<-c(mean.temp,new.mean.temp)
            dew.point.days<-c(dew.point.days,new.dew.point.days)
            mean.dew.point<-c(mean.dew.point,new.mean.dew.point)
            temp.dew.point.days<-c(temp.dew.point.days,new.temp.dew.point.days)
            temp.16.22.dew.point.days<-c(temp.16.22.dew.point.days,new.temp.16.22.dew.point.days)
            temp.7.30.dew.point.days<-c(temp.7.30.dew.point.days,new.temp.7.30.dew.point.days)
            
            mean.wetnesss<-c(mean.wetnesss,new.mean.wetness)
            tot.rains<-c(tot.rains,new.tot.rain)
            mean.solars<-c(mean.solars,new.mean.solar)
            mean.wind.speeds<-c(mean.wind.speeds,new.mean.wind.speed)
            mean.gust.speeds<-c(mean.gust.speeds,new.mean.gust.speed)
            
            
          } 
        }
      }
    }
  }
}

delta.pustules<-data.frame(tag=factor(tags),stem.iter=stem.iters,leaf.iter=leaf.iters,pustule.num=pustule.nums,area=start.vals,area.next=end.vals,time=days,
                           temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,temp.days=temp.days,mean.temp=mean.temp,dew.point.days=dew.point.days,mean.dew.point=mean.dew.point,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,temp.dew.point.days=temp.dew.point.days,
                           mean.wetness=mean.wetnesss,tot.rain=tot.rains,mean.solar=mean.solars,mean.wind.speed=mean.wind.speeds,mean.gust.speed=mean.gust.speeds)
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
predictors<-c("area","time","mean.temp","mean.dew.point","mean.wetness","mean.solar",
              "mean.temp:time","temp.days",
              "mean.dew.point:time","dew.point.days",
              "mean.temp:mean.dew.point:time","temp.dew.point.days",
              "mean.wetness:time","mean.solar:time","mean.wind.speed","mean.wind.speed:time","tot.rain")

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

switch.vars<-list(
  list("mean.temp:time",c("temp.days")),
  list("temp.days",c("mean.temp:time")),
  list("mean.dew.point:time",c("dew.point.days")),
  list("dew.point.days",c("mean.dew.point:time")),
  list("mean.temp:mean.dew.point:time",c("temp.dew.point.days")),
  list("temp.dew.point.days",c("mean.temp:mean.dew.point:time"))
)

force.include.vars<-list(
  list("mean.temp:time",c("mean.temp","time")),
  list("mean.dew.point:time",c("mean.dew.point","time")),
  list("mean.temp:mean.dew.point:time",c("mean.temp","mean.dew.point","time")),
  list("mean.wetness:time",c("mean.wetness","time")),
  list("mean.solar:time",c("mean.solar","time")),
  list("mean.wind.speed:time",c("mean.wind.speed","time"))
)

for (i in 1:length(switch.vars))
{
  on.var<-unlist(switch.vars[[i]][1])
  off.vars<-unlist(switch.vars[[i]][2])
  off.rows<-which(pred.mat[,on.var]==T)
  pred.mat[off.rows,off.vars]<-F
}

for (i in 1:length(force.include.vars))
{
  on.var<-unlist(force.include.vars[[i]][1])
  force.on.vars<-unlist(force.include.vars[[i]][2])
  on.rows<-which(pred.mat[,on.var]==T)
  pred.mat[on.rows,force.on.vars]<-T
}

pred.mat<-unique(pred.mat)

### create all sets of models
model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x]),collapse=" + ")))
re.model.set <- apply(pred.mat, 1, function(x) as.formula( paste(c("area.next ~ offset(area)",predictors[x],"(1|tag)"),collapse=" + ")))
#gam.model.set <- apply(pred.mat, 1, function(x) as.formula( paste0(paste(c("area.next ~ offset(area",predictors[x]),collapse=") + s("),")")))
  
names(model.set)<-seq(1,length(model.set),1)
names(re.model.set)<-seq(1,length(re.model.set),1)
#names(gam.model.set)<-seq(1,length(gam.model.set),1)

all.fit.models<-lapply(model.set,function(x) lm(x,data=delta.pustules))
re.all.fit.models<-lapply(re.model.set,function(x) lmer(x,data=delta.pustules,REML=F))
#gam.all.fit.models<-lapply(gam.model.set,function(x) gam(x,data=delta.pustules,method = "REML"))

## compare between models

### lms
AICs<-unlist(lapply(all.fit.models,AIC))
delta.AICs<-AICs-min(AICs)
candidate.models<-unname(which(delta.AICs<4))
model.set[candidate.models[order(AICs[candidate.models])]] #models to consider--offset(diam.last) not shown
index<-1
best.model<-all.fit.models[[order(AICs)[index]]]
summary(best.model)


### lmms
re.AICs<-unlist(lapply(re.all.fit.models,AIC))
delta.re.AICs<-re.AICs-min(re.AICs)
re.candidate.models<-unname(which(delta.re.AICs<4))
re.model.set[re.candidate.models[order(re.AICs[re.candidate.models])]] #models to consider--offset(diam.last) not shown
index<-1
best.model<-re.all.fit.models[[order(re.AICs)[index]]]
summary(best.model)

### gams
#gam.AICs<-unlist(lapply(gam.all.fit.models,AIC))
#delta.gam.AICs<-gam.AICs-min(gam.AICs)
#gam.candidate.models<-unname(which(delta.gam.AICs<4))
#gam.model.set[gam.candidate.models[order(delta.gam.AICs[gam.candidate.models])]] #models to consider--offset(diam.last) not shown
#index<-1
#best.model<-gam.all.fit.models[[order(gam.AICs)[index]]]
#summary(best.model)

### lmm mumin approach w/ interactions
options(na.action="na.omit")
all.fit.lm<-lm(area.next~offset(area)+area+time+mean.temp+mean.temp:time+mean.dew.point+mean.dew.point:time+mean.wetness+mean.wetnesss:time+mean.solar+mean:solar,data=delta.pustules)
options(na.action="na.fail")
dd<-dredge(all.fit,fixed=~offset(area)+area)
best<-lmer(area.next~offset(area)+area+mean.dew.point+mean.temp+mean.wetnesss+mean.dew.point:mean.temp+(1|tag),data=delta.pustules)
summary(best)
