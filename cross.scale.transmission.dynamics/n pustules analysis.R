library(mgcv)
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

## just a few trajectories
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
  mod0<-gam(n.pustules.next-n.pustules~0+s(n.pustules,by=time,bs="cs",k=4)+
              s(mean.temp,by=time,bs="cs",k=4)+
              s(max.temp,by=time,bs="cs",k=4)+
              s(min.temp,by=time,bs="cs",k=4)+
              s(mean.abs.hum,by=time,bs="cs",k=4)+
              s(max.abs.hum,by=time,bs="cs",k=4)+
              s(min.abs.hum,by=time,bs="cs",k=4)+
              s(mean.solar,by=time,bs="cs",k=4)+
              s(tot.rain,bs="cs",k=4)+
              s(site,bs="re",k=4),
            data=delta.n.pustules)
  
  summary(mod0)
  
  mod1<-gam(n.pustules.next-n.pustules~0+s(n.pustules,by=time,bs="cs",k=4),
              #s(mean.temp,by=time,bs="cs",k=4)+
              #s(max.temp,by=time,bs="cs",k=4)+
              #s(min.temp,by=time,bs="cs",k=4)+
              #s(mean.abs.hum,by=time,bs="cs",k=4)+
              #s(max.abs.hum,by=time,bs="cs",k=4)+
              #s(min.abs.hum,by=time,bs="cs",k=4)+
              #s(mean.solar,by=time,bs="cs",k=4)+
              #s(tot.rain,bs="cs",k=4)+
              #s(site,bs="re",k=4),
            data=delta.n.pustules)
  summary(mod1)    # All now significant. Mod1 fit via stepwise removal that resulted in 2 marginally signiificant terms. The AIC of this model is significantly less than that of the model without the marginally significant terms.
  saveRDS(mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
}

## load best model

n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")

## model checking
#plot(n.pustules.model,scale=0,pages=1) #plot smooths
#gam.check(n.pustules.model) #indicates that the number of knots is sufficient, except for mean temperature. Increasing knots doesn't resolve issue, which indicates the issue is with the data, and the model is OK.

# predict climate change effect

## predict change in n.pustules by climate for different pustule sizes

### climate modeled as day matching ~ 50th/75th/90th quantile of mean temp

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,26),ylim=c(-15,15),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
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

plot(0,0,xlim=c(0,26),ylim=c(-15,15),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
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

## predict size trajectory of n pustules across observation window
predict.n.pustules.trajectory<-function(site,temp.addition,color,pred.window=1,plot=T,output=F)
{  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  weath.dat<-all.weath[which(all.weath$site==site),]
  temp.rh.dat<-all.temp.rh[which(all.temp.rh$site==site),]
  min.date<-max(min(unique(as.Date(weath.dat$date))),min(unique(as.Date(temp.rh.dat$date.time))))
  max.date<-min(max(unique(as.Date(weath.dat$date))),max(unique(as.Date(temp.rh.dat$date.time))))
  dates<-seq(min.date,max.date,pred.window)
  start.n.pustules<-1
  xcords<-rep(NA,length(dates))
  ycords<-rep(NA,length(dates))
  
  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.n.pustules
    xcords.new<-c(1)
    ycords.new<-c(i)
    beta <- coef(n.pustules.model) ## posterior mean of coefs
    Vb   <- vcov(n.pustules.model) ## posterior  cov of coefs
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,i,temp.addition = temp.addition)
      Xp <- predict(n.pustules.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
      n <-2
      ilink <- family(n.pustules.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      y<-preds[1]
      if(y<0) {y<-0}
      reps<-reps+pred.window
      xcords.new<-c(xcords.new,reps)
      ycords.new<-c(ycords.new,y)
    }
    xcords<-rbind(xcords,xcords.new)
    ycords<-rbind(ycords,ycords.new)
    print(j)
  }
  xcords<-xcords[-1,]
  ycords<-ycords[-1,]
  
  if(plot)
  {
    for(k in 1:dim(xcords)[1])
    {
      points(xcords[k,],ycords[k,],type="l",col=color) 
    } 
  }
  if(output)
  {
    ycords
  }
}

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

plot.purple<-t_col("purple",80)
plot.red<-t_col("red",80)
plot.orange<-t_col("orange",80)

par(mar=c(6,6,2,2),mfrow=c(3,1))

### one day ahead projection
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,20),ylab='N.pustules',xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
dat<-predict.n.pustules.trajectory("GM",0,pred.window=1,plot.orange,T,T) 
points(1:36,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",1.8,pred.window=1,plot.red,T,T) 
points(1:36,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",3.7,pred.window=1,plot.purple,T,T) 
points(1:36,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)

### two days ahead projection
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,20),ylab='N. pustules',xlab="day",cex.lab=1.5,cex.axis=1.5,main="2 days ahead")
dat<-predict.n.pustules.trajectory("GM",0,pred.window=2,plot.orange,T,T) 
points(seq(1,35,2),colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",1.8,pred.window=2,plot.red,T,T) 
points(seq(1,35,2),colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",3.7,pred.window=2,plot.purple,T,T) 
points(seq(1,35,2),colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)

### seven days ahead projection
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,20),ylab='N pustules',xlab="week",cex.lab=1.5,cex.axis=1.5,main="1 week ahead")
dat<-predict.n.pustules.trajectory("GM",0,pred.window=7,plot.orange,T,T) 
points(seq(1,36,7),colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",1.8,pred.window=7,plot.red,T,T) 
points(seq(1,36,7),colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",3.7,pred.window=7,plot.purple,T,T) 
points(seq(1,36,7),colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)

