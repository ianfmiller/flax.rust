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
  ### create all sets of models--Only use foi and first obs height due to low sample size for environmental vars within several sites
  mod00<-glm(status.next~foi,data=foi.data,family = "binomial")
  mod0<-glm(status.next~foi,data=foi.data[which(!is.na(foi.data$first.obs.height)),],family = "binomial") #included as check to make sure coefficient of foi is relatively insensitive to data subsetting
  mod1<-glm(status.next~foi+first.obs.height,data=foi.data,family = "binomial") #subset data to make same n observations between mod0,mod1,mod2
  mod2<-glm(status.next~foi*first.obs.height,data=foi.data,family = "binomial")
  AIC(mod0,mod1,mod2)
  best.model<-mod1 #mod1 has lowest AIC, interaction insignificant in mod2
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")
}

par(mfrow=c(1,1),mar=c(6,8,6,6))
plot(jitter(foi.data$foi),jitter(foi.data$status.next),xlab="predicted spore deposition",ylab="",axes=F,cex.lab=1.5)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=1.5,line=3,tick = F)
axis(2,cex.axis=1.5,col="red",col.axis="red")
mtext("odds of infection",side=2,line=2.5,cex=1.5,col="red")
axis(1,cex.axis=1.5)

dummy.fois<-seq(0,.14,.0001)
new.data1<-data.frame("foi"=dummy.fois,"first.obs.height"=rep(5,length(dummy.fois)))
points(new.data1$foi,predict(best.model,newdata = new.data1,type = "response"),type="l",col="red",lwd=4,lty=1)

new.data1<-data.frame("foi"=dummy.fois,"first.obs.height"=rep(10,length(dummy.fois)))
points(new.data1$foi,predict(best.model,newdata = new.data1,type = "response"),type="l",col="red",lwd=4,lty=2)

new.data1<-data.frame("foi"=dummy.fois,"first.obs.height"=rep(25,length(dummy.fois)))
points(new.data1$foi,predict(best.model,newdata = new.data1,type = "response"),type="l",col="red",lwd=4,lty=3)

new.data1<-data.frame("foi"=dummy.fois,"first.obs.height"=rep(50,length(dummy.fois)))
points(new.data1$foi,predict(best.model,newdata = new.data1,type = "response"),type="l",col="red",lwd=4,lty=4)

legend("topright",legend=c("height = 5cm","height = 10cm","height = 25cm","height = 50cm"),col="red",lwd=4,lty=c(1,2,3,4),cex=2,bty="n")
