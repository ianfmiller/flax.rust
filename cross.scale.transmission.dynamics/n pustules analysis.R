library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/n pustules data prep.R")

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

## one trajectory
par(mfrow=c(1,1),mar=c(6,6,6,6))
tags<-c(86,88,106,112,124)
i<-0
plot.cols<-sample(rainbow(36))
plot(c(min(n.pustules[which(n.pustules$tag %in% tags),]$date),max(n.pustules[which(n.pustules$tag %in% tags),]$date)),c(0,max(n.pustules[which(n.pustules$tag %in% tags),]$N.pustules)),type="n",xlab="date",ylab="N pustules",cex.lab=2,cex.axis=2)
sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
for (tag in tags)
{
  sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
  
  for (color in unique(sub.n.pustules1$color))
  {
    sub.n.pustules2<-sub.n.pustules1[which(sub.n.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.n.pustules2$leaf.iteration)) 
    {
      sub.n.pustules3<-sub.n.pustules2[which(sub.n.pustules2$leaf.iteration==leaf.iteration),]
      
      i<-i+1
      points(sub.n.pustules3$date,sub.n.pustules3$N.pustules,col=plot.cols[i],type="l",lwd=5)
      
    }
  }
}
#plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,50),type="n",xlab="date",ylab="pustule area")


## plot change
par(mfrow=c(1,1))
plot(delta.n.pustules$n.pustules,delta.n.pustules$n.pustules.next,col="black",xlab = "N pustules",ylab="next obs. N pustules",cex.lab=2,cex.axis=2)
abline(0,1,lty=2)
mtext(text="N = 650",cex=2)

# analyze data

## fit models--only if not already fit
### using k=15 after gam checks indicated that k=10 resulted in unevenly distributed residuals

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS"))
{
  mod0<-gam(n.pustules.next~s(n.pustules,bs="cs",k=4)+
              s(mean.temp,by=time,bs="cs",k=4)+
              s(max.temp,by=time,bs="cs",k=4)+
              s(min.temp,by=time,bs="cs",k=4)+
              s(mean.abs.hum,by=time,bs="cs",k=4)+
              s(max.abs.hum,by=time,bs="cs",k=4)+
              s(min.abs.hum,by=time,bs="cs",k=4)+
              s(mean.solar,by=time,bs="cs",k=4)+
              s(mean.wetness,by=time,bs="cs",k=4)+
              s(tot.rain,by=time,bs="cs",k=4)+
              s(site,bs="re",k=4),
            data=delta.n.pustules)
  
  summary(mod0) #indicates that max.abs.hum, mean.solar, and tot.rain are not significant predictors
  
  mod1<-gam(n.pustules.next~s(n.pustules,bs="cs",k=4)+
              #s(mean.temp,by=time,bs="cs",k=4)+
              s(max.temp,by=time,bs="cs",k=4)+
              #s(min.temp,by=time,bs="cs",k=4)+
              #s(mean.abs.hum,by=time,bs="cs",k=4)+
              #s(max.abs.hum,by=time,bs="cs",k=4)+
              #s(min.abs.hum,by=time,bs="cs",k=4)+
              #s(mean.solar,by=time,bs="cs",k=4)+
              s(mean.wetness,by=time,bs="cs",k=4)+
              #s(tot.rain,by=time,bs="cs",k=4)+
              s(site,bs="re",k=4),
            data=delta.n.pustules)
  summary(mod1) # All now significant. Stepwise removal yields the same result. There were two marginally significant (p<.1) terms. A model including these terms was ~+2 AIC points, indicating that the simpler model is likely the best.

  saveRDS(mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
}

## load best model

n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")

## model checking
plot(n.pustules.model,scale=0,pages=1) #plot smooths
gam.check(n.pustules.model) #indicates that the number of knots is sufficient

## better visualize model

par(mfrow=c(2,2))
plot(n.pustules.model,scale=0,select=1)
abline(0,1,col="red",lty=2)
vis.gam(n.pustules.model,view = c("max.temp","time"),n.grid=30,plot.type = "contour",zlim=c(0,11),color="topo",contour.col = "black")
vis.gam(n.pustules.model,view = c("mean.wetness","time"),n.grid=30,plot.type = "contour",zlim=c(0,11),color="topo",contour.col = "black")
plot(n.pustules.model,scale=0,select=4)

par(mfrow=c(2,2))
plot(n.pustules.model,scale=0,select=1)
vis.gam(n.pustules.model,view = c("max.temp","time"),n.grid=30,plot.type = "persp",zlim=c(0,11),se=1,theta=45,phi=15,ticktype="detailed")
vis.gam(n.pustules.model,view = c("mean.wetness","time"),n.grid=30,plot.type = "persp",zlim=c(0,11),se=1,theta=45,phi=15,ticktype="detailed")
plot(n.pustules.model,scale=0,select=4)

# predict climate change effect

## predict change in n.pustules by climate for different pustule sizes

### climate modeled as day matching ~ 50th/75th/90th quantile of mean temp

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,26),ylim=c(-4,2),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(1,26,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.n.pustules.model(day,i)
    
    Xp <- predict(n.pustules.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
    beta <- coef(n.pustules.model) ## posterior mean of coefs
    Vb   <- vcov(n.pustules.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(n.pustules.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

plot(0,0,xlim=c(0,26),ylim=c(-4,2),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(1,26,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.n.pustules.model(75,i,temp.addition = temp.addition)
    
    Xp <- predict(n.pustules.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
    beta <- coef(n.pustules.model) ## posterior mean of coefs
    Vb   <- vcov(n.pustules.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(n.pustules.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")
