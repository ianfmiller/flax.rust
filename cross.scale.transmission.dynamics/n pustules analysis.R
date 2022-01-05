library(mgcv)
library(MASS)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/n pustules data prep.R")

delta.n.pustules<-subset(delta.n.pustules,time<=8)

# visualize data
layout(matrix(c(1,2,4,4,3,3,3,3),2,4,byrow = T))
par(mar=c(5,5,3,0.5))
## histograms
hist(delta.n.pustules$n.pustules,main="",breaks=100,xlab="N pustules",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.n.pustules$n.pustules.next-delta.n.pustules$n.pustules)/delta.n.pustules$time,main="",breaks=100,xlab="change in N pustules per day",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)
## plot trajectories
par(mar=c(2,5,0,0.5))
plot(c(min(n.pustules$date),max(n.pustules$date)),c(0,max(n.pustules$N.pustules)),type="n",xlab="date",ylab="N pustules")
mtext("D",side=3,adj=1,line=-3,cex=2)

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
#par(mfrow=c(1,1),mar=c(6,6,6,6))
#tags<-c(86,88,106,112,124)
#i<-0
#plot.cols<-sample(rainbow(28))
#plot(c(min(n.pustules[which(n.pustules$tag %in% tags),]$date),max(n.pustules[which(n.pustules$tag %in% tags),]$date)),c(0,max(n.pustules[which(n.pustules$tag %in% tags),]$N.pustules)),type="n",xlab="date",ylab="N pustules",cex.lab=2,cex.axis=2)
#sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
#for (tag in tags)
#{
#  sub.n.pustules1<-n.pustules[which(n.pustules$tag==tag),]
#  
#  for (color in unique(sub.n.pustules1$color))
#  {
#    sub.n.pustules2<-sub.n.pustules1[which(sub.n.pustules1$color==color),]
#    
#    for(leaf.iteration in unique(sub.n.pustules2$leaf.iteration)) 
#    {
#      sub.n.pustules3<-sub.n.pustules2[which(sub.n.pustules2$leaf.iteration==leaf.iteration),]
#      
#      i<-i+1
#      points(sub.n.pustules3$date,sub.n.pustules3$N.pustules,col=plot.cols[i],type="l",lwd=2)
#      
#    }
#  }
#}


## plot change
par(mar=c(5,5,3,0.5))
plot(delta.n.pustules$n.pustules,delta.n.pustules$n.pustules.next,col="black",xlab = 'N pustules',ylab='next obs. N pustules',cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3.25,cex=2)
abline(0,1,lty=2)

# analyze data

## fit models--only if not already fit
### using k=15 after gam checks indicated that k=10 resulted in unevenly distributed residuals

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS"))
{
  set.seed(23409687)
  mod<-gam((n.pustules.next-n.pustules)/time~0+
             s(n.pustules)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(max.abs.hum)+
             s(min.abs.hum)+
             s(mean.solar)+
             s(mean.daily.rain)+
             s(mean.wetness)+
             s(tag,bs="re")+
             s(site,bs="re"),
           select = T,
           method="REML",
           data=delta.n.pustules,
           control = list(nthreads=4))
  
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
  
}

## load best model

n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")

## model checking
par(mfrow=c(2,2))
gam.check(n.pustules.model) #Indicates that k is sufficient for all significant predictor terms
concurvity(n.pustules.model,full=F) #Concurvity is expected between all weather variables. The non-pessimistic estimations of concurvity (under $estimate and $observed) indicate that it is not a serious issue in the model fit.

## visualize model
par(mfrow=c(3,4),mar=c(4,4.5,3,4))
plot(n.pustules.model,select = 1,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="N pustules",ylab="s(N pustules)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 2,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean temperature (°C)",ylab="s(mean temperature)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 3,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="max. temperature (°C)",ylab="s(max. temperature)")
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 4,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="min. temperature (°C)",ylab="s(min. temperature)")
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 5,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean abs. humidity ('*g/m^3*')'),ylab="s(mean abs. humidity)")
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 6,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('max. abs. humidity ('*g/m^3*')'),ylab="s(max. abs. humidity)")
grid()
mtext("F",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 7,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('min abs. humidity ('*g/m^3*')'),ylab="s(min. abs. humidity)")
grid()
mtext("G",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean solar radiation ('*W/m^2*')'),ylab="s(mean solalr radiation")
grid()
mtext("H",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 9,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="total rainfall (mm)",ylab="s(total rainfall)")
grid()
mtext("I",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 10,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean leaf wetness (%)",ylab="s(mean leaf wetness)")
grid()
mtext("J",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 11,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(tag)")
grid()
mtext("K",adj=1,cex=1.25,font=2)
plot(n.pustules.model,select = 12,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(site)")
grid()
mtext("L",adj=1,cex=1.25,font=2)

# predict climate change effect

## predict change in n.pustules by climate for different pustule sizes

### climate modeled as day matching ~ 50th/75th/90th quantile of mean temp

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,206),ylim=c(-15,15),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(1,206,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.n.pustules.model(day,i)
    
    Xp <- predict(n.pustules.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
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
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)))
    points(i,mean(y),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

plot(0,0,xlim=c(0,206),ylim=c(-15,15),type="n",xlab="N pustules",ylab="pred. change in N pustules",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(1,206,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.n.pustules.model(75,i,temp.addition = temp.addition)
    
    Xp <- predict(n.pustules.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
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
    
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)))
    points(i,mean(y),col=colors[index],pch=16,cex=2)
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
  start.n.pustules<-10
  xcords<-rep(NA,length(dates))
  ycords<-rep(NA,length(dates))
  
  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.n.pustules
    xcords.new<-c(1)
    ycords.new<-c(i)
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,i,temp.addition = temp.addition)
      pred.data<-cbind(pred.data,"tag"=86) # add in tag to pred data. This term is ignored when simulating the model. It is included only to suppress warning messages that "tag" is missing in data
      beta <- coef(n.pustules.model) ## posterior mean of coefs
      Vb   <- vcov(n.pustules.model) ## posterior  cov of coefs
      n <-2
      mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
      Xp <- predict(n.pustules.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
      ilink <- family(n.pustules.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      y<-preds[1]
      i<-i+y
      if(i<1) {i<-1}
      reps<-reps+pred.window
      xcords.new<-c(xcords.new,reps)
      ycords.new<-c(ycords.new,i)
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


### one day ahead projection
par(mar=c(6,6,2,2),mfrow=c(2,2))

plot(0,0,type="n",xlim=c(1,36),ylim=c(0,40),ylab='N pustules',xlab="day",cex.lab=1.5,cex.axis=1.5,main="CC")
dat<-predict.n.pustules.trajectory("CC",0,pred.window=1,plot.orange,T,T) 
points(1:36,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("CC",1.8,pred.window=1,plot.red,T,T) 
points(1:36,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("CC",3.7,pred.window=1,plot.purple,T,T) 
points(1:36,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",xlim=c(1,41),ylim=c(0,40),ylab='N pustules',xlab="day",cex.lab=1.5,cex.axis=1.5,main="BT")
dat<-predict.n.pustules.trajectory("BT",0,pred.window=1,plot.orange,T,T) 
points(1:41,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("BT",1.8,pred.window=1,plot.red,T,T) 
points(1:41,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("BT",3.7,pred.window=1,plot.purple,T,T) 
points(1:41,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",xlim=c(1,41),ylim=c(0,40),ylab='N pustules',xlab="day",cex.lab=1.5,cex.axis=1.5,main="GM")
dat<-predict.n.pustules.trajectory("GM",0,pred.window=1,plot.orange,T,T) 
points(1:36,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",1.8,pred.window=1,plot.red,T,T) 
points(1:36,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("GM",3.7,pred.window=1,plot.purple,T,T) 
points(1:36,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",xlim=c(1,41),ylim=c(0,40),ylab='N pustules',xlab="day",cex.lab=1.5,cex.axis=1.5,main="HM")
dat<-predict.n.pustules.trajectory("HM",0,pred.window=1,plot.orange,T,T) 
points(1:16,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("HM",1.8,pred.window=1,plot.red,T,T) 
points(1:16,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.n.pustules.trajectory("HM",3.7,pred.window=1,plot.purple,T,T) 
points(1:16,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)

