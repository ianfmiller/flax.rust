library(mgcv)
library(lme4)
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

par(mfrow=c(2,1),mar=c(5,5,2,2))
plot(jitter(foi.data$foi),jitter(foi.data$status.next),xlab="predicted spore deposition",ylab="",axes=F,cex.lab=1.5)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5)
axis(1,cex.axis=1.5)

plot(jitter(log10(foi.data$foi)),jitter(foi.data$status.next),xlab=expression('predicted '*log[10]*'spore deposition'),ylab="",axes=F,cex.lab=1.5)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5)
axis(1,cex.axis=1.5)


if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS"))
{
  mod0<-glmer(status.next~
                log10(foi)+
                #height.cm+
                #mean.temp + 
                #max.temp + 
                #min.temp + 
                #mean.abs.hum + 
                #max.abs.hum + 
                #min.abs.hum + 
                #tot.rain + 
                #mean.solar + 
                (1|site),
              data = foi.data,
              family = binomial(link="logit")
              )
  summary(mod0)
  
  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/foi model set creation.R") ### load models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("status.next ~  log10(foi)*height.cm + (1|site)",predictors[x]),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(glmer(status.next~log10(foi)*height.cm + (1|site),data = foi.data,family = binomial(link="logit"),nAGQ = 0)) #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-glmer(model.set[[i]],data=foi.data,family=binomial,nAGQ=0))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  summary(best.model)
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")
}

foi.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")

par(mfrow=c(1,1),mar=c(6,8,6,6))
plot(log10(foi.data$foi),jitter(foi.data$status.next,factor=.5),xlab="predicted spore deposition",ylab="",axes=F,cex.lab=1.5,ylim=c(-.2,1.2),xlim=c(-8,1))
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5,line=3,tick = F)
axis(2,at=c(0,.25,.5,.75,1),cex.axis=1.5,col="red",col.axis="red")
mtext("odds of infection",side=2,line=2.5,cex=1.5,col="red")
axis(1,cex.axis=1.5)

fois<-10^seq(-8,0,.25)
new.data.1<-data.frame("foi"=fois,
                     "height.cm"=5,
                     "mean.temp"=mean(foi.data$mean.temp),
                     "max.temp"=mean(foi.data$max.temp),
                     "min.temp"=mean(foi.data$min.temp),
                     "mean.abs.hum"=mean(foi.data$mean.abs.hum),
                     "tot.rain"=mean(foi.data$tot.rain)
                     )
new.data.2<-data.frame("foi"=fois,
                       "height.cm"=10,
                       "mean.temp"=mean(foi.data$mean.temp),
                       "max.temp"=mean(foi.data$max.temp),
                       "min.temp"=mean(foi.data$min.temp),
                       "mean.abs.hum"=mean(foi.data$mean.abs.hum),
                       "tot.rain"=mean(foi.data$tot.rain)
                      )
new.data.3<-data.frame("foi"=fois,
                       "height.cm"=25,
                       "mean.temp"=mean(foi.data$mean.temp),
                       "max.temp"=mean(foi.data$max.temp),
                       "min.temp"=mean(foi.data$min.temp),
                       "mean.abs.hum"=mean(foi.data$mean.abs.hum),
                       "tot.rain"=mean(foi.data$tot.rain)
                      )
new.data.4<-data.frame("foi"=fois,
                       "height.cm"=50,
                       "mean.temp"=mean(foi.data$mean.temp),
                       "max.temp"=mean(foi.data$max.temp),
                       "min.temp"=mean(foi.data$min.temp),
                       "mean.abs.hum"=mean(foi.data$mean.abs.hum),
                       "tot.rain"=mean(foi.data$tot.rain)
                      )
points(log10(new.data.1$foi),predict(foi.model,newdata = new.data.1,type = "response",re.form=NA),type="l",col="red",lwd=2,lty=1)
points(log10(new.data.2$foi),predict(foi.model,newdata = new.data.2,type = "response",re.form=NA),type="l",col="red",lwd=2,lty=2)
points(log10(new.data.3$foi),predict(foi.model,newdata = new.data.3,type = "response",re.form=NA),type="l",col="red",lwd=2,lty=3)
points(log10(new.data.4$foi),predict(foi.model,newdata = new.data.4,type = "response",re.form=NA),type="l",col="red",lwd=2,lty=4)
legend("topright",legend=c("50cm","25cm","10cm","5cm"),col="red",lty=c(4,3,2,1),bty="n",cex=1.5,lwd=2)

par(mfrow=c(1,1),mar=c(6,8,6,6))
plot(0,0,type = "n",xlab="predicted spore deposition",ylab="relative odds of infection",axes=F,cex.lab=1.5,ylim=c(0,10),xlim = c(-8,1))
axis(1)
axis(2)
base.pred<-predict(foi.model,newdata = new.data.1,type = "response",re.form=NA)[1]
points(log10(new.data.1$foi),predict(foi.model,newdata = new.data.1,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=2,lty=1)
points(log10(new.data.2$foi),predict(foi.model,newdata = new.data.2,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=2,lty=2)
points(log10(new.data.3$foi),predict(foi.model,newdata = new.data.3,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=2,lty=3)
points(log10(new.data.4$foi),predict(foi.model,newdata = new.data.4,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=2,lty=4)
legend("topleft",legend=c("50cm","25cm","10cm","5cm"),col="red",lty=c(4,3,2,1),bty="n",cex=1.5,lwd=2)


         