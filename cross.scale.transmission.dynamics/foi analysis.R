library(mgcv)
library(viridis)
# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/building foi dataset.R")
foi.data<-foi.data[which(foi.data$status==0),]
foi.data$site<-as.factor(foi.data$site)
for(i in 1:nrow(foi.data))
{
  if(is.na(foi.data[i,"tag"])) {foi.data[i,"tag"]<-paste0(foi.data[i,"site"],"X",foi.data[i,"X"]+foi.data[i,"x"],"Y",foi.data[i,"Y"]+foi.data[i,"y"])}
}
foi.data$tag<-as.factor(foi.data$tag)

## visualize data

layout(matrix(c(1,1,2,3,4,5),3,2,byrow = T))
par(mar=c(6,6,6,6))
plot(jitter(log10(foi.data$foi)),jitter(foi.data$status.next),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",main="all",axes=F,cex.lab=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(foi.data[which(foi.data$site=="CC"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="CC"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[2],xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="CC",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(foi.data[which(foi.data$site=="BT"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="BT"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[3],xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="BT",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(foi.data[which(foi.data$site=="GM"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="GM"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[4],xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="GM",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()
plot(jitter(log10(foi.data[which(foi.data$site=="HM"),"foi"]+1e-10)),jitter(foi.data[which(foi.data$site=="HM"),"status.next"]),xlab=expression('predicted '*log[10]*' spore deposition'),ylab="outcome",axes=F,col=magma(6)[5],xlim=c(min(log10(foi.data$foi+1e-10)),max(log10(foi.data$foi+1e-10))),main="HM",cex.lab=2,cex.main=2)
axis(2,at=c(0,1),labels = c("healthy","infected"),cex.axis=2)
axis(1,cex.axis=2)
box()


if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS"))
{
  mod<-gam(status.next~
             s(log10(foi),height.cm)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(max.abs.hum)+
             s(min.abs.hum)+
             s(mean.solar)+
             s(mean.daily.rain)+
             s(mean.wetness)+
             +offset(log(time))+
             s(tag,bs="re")+
             s(site,bs="re")+
           family=binomial(link="cloglog"),
           select = T,
           method="REML",
           data=foi.data,
           control = list(nthreads=4))
  
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")
}

foi.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")

par(mar=c(6,8,6,0))
layout(matrix(c(1,1,1,1,2),1,5))
plot(log10(foi.data$foi),jitter(foi.data$status.next*.1,factor=.5),xlab=expression('predicted '*log[10]*' spore depositin'),ylab="",axes=F,cex.lab=2,ylim=c(-.025,0.125),xlim=c(-6,3))
axis(2,at=c(0,.1),labels = c("healthy","infected"),cex.axis=2,line=3,tick = F)
axis(2,at=c(0,.05,.1),cex.axis=2,col="red",col.axis="red")
mtext("odds of infection",side=2,line=2.5,cex=2,col="red")
axis(1,cex.axis=1.5)

fois<-10^seq(-6,2.3,.1)
new.data.1<-data.frame("foi"=fois,
                     "height.cm"=5,
                     "site"="CC"
                     )
new.data.2<-data.frame("foi"=fois,
                       "height.cm"=10,
                       "site"="CC"
                      )
new.data.3<-data.frame("foi"=fois,
                       "height.cm"=25,
                       "site"="CC"
                      )
new.data.4<-data.frame("foi"=fois,
                       "height.cm"=50,
                       "site"="CC"
                      )
points(log10(new.data.1$foi),predict(foi.model,newdata = new.data.1,type = "response",exclude = "s(site)"),type="l",col="red",lwd=4,lty=1)
points(log10(new.data.2$foi),predict(foi.model,newdata = new.data.2,type = "response",exclude = "s(site)"),type="l",col="red",lwd=4,lty=2)
points(log10(new.data.3$foi),predict(foi.model,newdata = new.data.3,type = "response",exclude = "s(site)"),type="l",col="red",lwd=4,lty=3)
points(log10(new.data.4$foi),predict(foi.model,newdata = new.data.4,type = "response",exclude = "s(site)"),type="l",col="red",lwd=4,lty=4)
par(mar=c(0,0,0,0))
plot(0,0,type="n",axes=F,xlab = "",ylab="")
legend("left",legend=c("50cm","25cm","10cm","5cm"),col="red",lty=c(4,3,2,1),bty="n",cex=2,lwd=4)

par(mfrow=c(1,1),mar=c(6,8,6,6))
plot(0,0,type = "n",xlab=expression('predicted '*log[10]*' spore depositin'),ylab="relative odds of infection",axes=F,cex.lab=2,ylim=c(0,100),xlim = c(-6,3))
axis(1)
axis(2)
base.pred<-predict(foi.model,newdata = new.data.1,type = "response",re.form=NA)[1]
points(log10(new.data.1$foi),predict(foi.model,newdata = new.data.1,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=4,lty=1)
points(log10(new.data.2$foi),predict(foi.model,newdata = new.data.2,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=4,lty=2)
points(log10(new.data.3$foi),predict(foi.model,newdata = new.data.3,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=4,lty=3)
points(log10(new.data.4$foi),predict(foi.model,newdata = new.data.4,type = "response",re.form=NA)/base.pred,type="l",col="red",lwd=4,lty=4)
legend("topleft",legend=c("50cm","25cm","10cm","5cm"),col="red",lty=c(4,3,2,1),bty="n",cex=2,lwd=4)


         