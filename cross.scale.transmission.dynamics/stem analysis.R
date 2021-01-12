library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/length inf data prep.R")
delta.length.inf<-subset(delta.length.inf,time<=7)
delta.length.inf<-delta.length.inf[-which(delta.length.inf$start.length.inf>delta.length.inf$start.length),]
#delta.length.inf<-delta.length.inf[-which(delta.length.inf$end.length.inf>delta.length.inf$end.length),]
# visualize data

## histograms
par(mfrow=c(2,1))
hist(as.numeric(length.inf$length.tissue.infected),main="length stem inf",breaks=100,xlab="length stem inf")
hist(delta.length.inf$end.length.inf-delta.length.inf$start.length.inf,main="change in length stem inf",breaks=100,xlab="change in length stem inf")

## plot trajectories
par(mfrow=c(2,1))

## plot length stem infected
plot(c(min(as.Date(length.inf$Date,tryFormats = "%m/%d/%Y")),max(as.Date(length.inf$Date,tryFormats = "%m/%d/%Y"))),c(0,max(length.inf$length.tissue.infected,na.rm = T)),type="n",xlab="date",ylab="length tissue infected",main="length tissue infected")
i<-0
plot.cols<-sample(rainbow(256))

for (tag in unique(length.inf$Tag))
{
  sub.length.inf1<-length.inf[which(length.inf$Tag==tag),]
  
  for (stem.index in unique(sub.length.inf1$stem.index))
  {
    sub.length.inf2<-sub.length.inf1[which(sub.length.inf1$stem.index==stem.index),]
    sub.length.inf2<-sub.length.inf2[-which(is.na(sub.length.inf2$length.tissue.infected)),]
    points(sub.length.inf2$Date,sub.length.inf2$length.tissue.infected,col=plot.cols[i],type="l",lwd=.5)
    i<-i+1
  }
}

## plot percent stem infected
plot(c(min(as.Date(length.inf$Date,tryFormats = "%m/%d/%Y")),max(as.Date(length.inf$Date,tryFormats = "%m/%d/%Y"))),c(0,1),type="n",xlab="date",ylab="percent tissue infected",main="percent tissue infected")
i<-0
plot.cols<-sample(rainbow(256))

for (tag in unique(length.inf$Tag))
{
  sub.length.inf1<-length.inf[which(length.inf$Tag==tag),]
  
  for (stem.index in unique(sub.length.inf1$stem.index))
  {
    sub.length.inf2<-sub.length.inf1[which(sub.length.inf1$stem.index==stem.index),]
    sub.length.inf2<-sub.length.inf2[-which(is.na(sub.length.inf2$length.tissue.infected)),]
    points(sub.length.inf2$Date,sub.length.inf2$length.tissue.infected/sub.length.inf2$stem.height,col=plot.cols[i],type="l",lwd=.5)
    i<-i+1
  }
}

## plot change
par(mfrow=c(1,1))
plot(delta.length.inf$start.length.inf,delta.length.inf$end.length.inf,col="grey",xlab = "length stem inf",ylab="next obs. length stem inf")
abline(0,1)


# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/length.inf.model.RDS"))
{
  ### construct all combinations of predictors
  source("length.inf.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("end.length.inf ~ offset(start.length.inf)",predictors[x],'(1|tag)'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  #AIC.benchmark<- AIC(lmer(end.length.inf~offset(start.length.inf)+start.length.inf+(1|tag),data=delta.length.inf)) #cutoff to limit memory usage
  AIC.benchmark<-1050
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-lmer(model.set[[i]],data=delta.length.inf))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/length.inf.model.RDS")
}

## load best model

length.inf.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/length.inf.model.RDS")

## visualize model

