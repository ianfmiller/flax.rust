library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/building foi dataset.R")
foi.data<-foi.data[which(foi.data$status==0),]

# visualize data

## histograms
par(mfrow=c(1,1))
hist(foi.data$foi,main="all sites",breaks=100,xlab="force of infection")
par(mfrow=c(2,2))
hist(foi.data[which(foi.data$site=="CC"),"foi"],breaks=100,xlab="force of infection",xlim=c(0,max(foi.data$foi)),col="yellow",main="CC")
hist(foi.data[which(foi.data$site=="BT"),"foi"],breaks=100,xlab="force of infection",xlim=c(0,max(foi.data$foi)),col="orange",main="BT")
hist(foi.data[which(foi.data$site=="GM"),"foi"],breaks=100,xlab="force of infection",xlim=c(0,max(foi.data$foi)),col="red",main="GM")
hist(foi.data[which(foi.data$site=="HM"),"foi"],breaks=100,xlab="force of infection",xlim=c(0,max(foi.data$foi)),col="purple",main="HM")

## outcome ~ foi

layout(matrix(c(1,6,2,7,3,8,4,9,5,10),5,2,byrow = T))
par(mar=c(2,3,2,3))
plot(jitter(foi.data$foi),jitter(foi.data$status.next),xlab="foi",ylab="outcome",axes=F,main="foi vs outcome")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(foi.data[which(foi.data$site=="CC"),"foi"]),jitter(foi.data[which(foi.data$site=="CC"),"status.next"]),xlab="foi",ylab="outcome",axes=F,col="yellow",xlim=c(0,max(foi.data$foi)),main="CC")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(foi.data[which(foi.data$site=="BT"),"foi"]),jitter(foi.data[which(foi.data$site=="BT"),"status.next"]),xlab="foi",ylab="outcome",axes=F,col="orange",xlim=c(0,max(foi.data$foi)),main="BT")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(foi.data[which(foi.data$site=="GM"),"foi"]),jitter(foi.data[which(foi.data$site=="GM"),"status.next"]),xlab="foi",ylab="outcome",axes=F,col="red",xlim=c(0,max(foi.data$foi)),main="GM")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(foi.data[which(foi.data$site=="HM"),"foi"]),jitter(foi.data[which(foi.data$site=="HM"),"status.next"]),xlab="foi",ylab="outcome",axes=F,col="purple",xlim=c(0,max(foi.data$foi)),main="HM")
axis(2,at=c(0,1),labels = c("healthy","infected"))

plot(jitter(log10(foi.data$foi)+1e-10),jitter(foi.data$status.next),xlab="log10 foi",ylab="outcome",axes=F,main="log10 foi vs outcome")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(log10(foi.data[which(foi.data$site=="CC"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="CC"),"status.next"]),xlab="log10 foi",ylab="outcome",axes=F,col="yellow",xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="CC")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(log10(foi.data[which(foi.data$site=="BT"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="BT"),"status.next"]),xlab="log10 foi",ylab="outcome",axes=F,col="orange",xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="BT")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(log10(foi.data[which(foi.data$site=="GM"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="GM"),"status.next"]),xlab="log10 foi",ylab="outcome",axes=F,col="red",xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="GM")
axis(2,at=c(0,1),labels = c("healthy","infected"))
plot(jitter(log10(foi.data[which(foi.data$site=="HM"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="HM"),"status.next"]),xlab="log10 foi",ylab="outcome",axes=F,col="purple",xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="HM")
axis(2,at=c(0,1),labels = c("healthy","infected"))

par(mfrow=c(2,1),mar=c(6,6,6,6))
plot(jitter(foi.data$foi),jitter(foi.data$status.next),xlab="predicted spore deposition",ylab="",axes=F,cex.lab=1.5)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5)
axis(1,cex.axis=1.5)

plot(jitter(log10(foi.data$foi)),jitter(foi.data$status.next),xlab=expression('predicted '*log[10]*'spore deposition'),ylab="",axes=F,cex.lab=1.5)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5)
axis(1,cex.axis=1.5)


if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS"))
{
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/foi model set creation.R") ### load models
  
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("status.next ~ foi",predictors[x]),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(glm(status.next~foi,data = foi.data)) #cutoff to limit memory usage
  #AIC.benchmark<- 8025
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-glm(model.set[[i]],data=foi.data))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
}

best.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")

par(mfrow=c(1,1),mar=c(6,8,6,6))
plot(jitter(foi.data$foi),jitter(foi.data$status.next),xlab="predicted spore deposition",ylab="",axes=F,cex.lab=1.5)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5,line=3,tick = F)
axis(2,cex.axis=1.5,col="red",col.axis="red")
mtext("odds of infection",side=2,line=2.5,cex=1.5,col="red")
axis(1,cex.axis=1.5)

fois<-seq(0,.15,.0001)
new.data<-data.frame("foi"=fois,
                     "mean.temp.days.16.22"=rep(mean(foi.data$mean.temp.days.16.22),times=length(fois)),
                     "mean.dew.point.days"=rep(mean(foi.data$mean.dew.point.days),times=length(fois)),
                     "mean.temp.dew.point.days"=rep(mean(foi.data$mean.temp.dew.point.days ),times=length(fois)),
                     "mean.temp.7.30.wetness.days"=rep(mean(foi.data$mean.temp.7.30.wetness.days),times=length(fois)),
                     "mean.tot.rain"=rep(mean(foi.data$mean.tot.rain ),times=length(fois))
                     )

points(new.data$foi,predict(best.model,newdata = new.data,type = "response"),type="l",col="red",lwd=4,lty=1)

