set.seed(34020968)
library(mgcv)
library(lme4)
library(lmerTest)
library(progress)
library(MASS)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")

delta.pustules<-subset(delta.pustules,time<=7)

# visualize data
layout(matrix(c(1,2,4,4,3,3,3,3),2,4,byrow = T))
par(mar=c(5,5,3,0.5))
## histograms
hist(delta.pustules$area,main="",breaks=100,xlab="pustule area",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist(delta.pustules$area.next-delta.pustules$area,main="",breaks=100,xlab="change in pustule area",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories
par(mar=c(2,5,0,0.5))
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="",ylab="pustule area",cex.lab=2,cex.axis=2)
mtext("D",side=3,adj=1,line=-3,cex=2)
i<-0

plot.cols<-sample(rainbow(2007))

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

## just a few trajectories

#plot(c(min(pustules$date),max(pustules$date)),c(0,1),type="n",xlab="date",ylab="pustule area")
#i<-0

#plot.cols<-sample(rainbow(37))

#for (tag in unique(pustules$tag)[1:1])
#{
#  sub.pustules1<-pustules[which(pustules$tag==tag),]
#  
#  for (color in unique(sub.pustules1$color))
#  {
#    sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
#    
#    for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
#    {
#      sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
#      
#      for(pustule.number in unique(sub.pustules3$pustule.number))
#      {
#        i<-i+1
#        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
#        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
#        points(sub.pustules4$date,sub.pustules4$area,col=plot.cols[i],type="l",lwd=1)
#      }
#    }
#  }
#}

## plot change
par(mar=c(5,5,3,0.5))
plot(delta.pustules$area,delta.pustules$area.next,col="black",xlab = expression(area(mm^2)),ylab=expression('next obs. area '*(mm^2)),cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3.25,cex=2)
abline(0,1,lty=2)

# analyze data

## fit model--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS"))
{
  mod<-gam(area.next-area~0+
              te(area,time)+
              s(mean.temp)+
              s(max.temp)+
              s(min.temp)+
              s(mean.abs.hum)+
              s(max.abs.hum)+
              s(min.abs.hum)+
              s(mean.solar)+
              s(tot.rain)+
              s(mean.wetness)+
              s(tag,bs="re")+
              s(site,bs="re"),
            select = T,
            method="REML",
            data=delta.pustules,
            control = list(nthreads=4))
  
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
}

## load best model

pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")

## model checking
gam.check(pustule.model) #Indicates that k should be higher all smooths aside from the tensor. Increasing k to significantly higher values is not possible due to terms having fewer unique covariate combinations than specified maximum degrees of freedom. 
concurvity(pustule.model,full=F) #Concurvity is expected between all weather variables. The non-pessimistic estimations of concurvity (under $estimate and $observed) indicate that it is not a serious issue in the model fit.

## visualize model
layout(matrix(c(1,1,1,1,2,3,3,3,3,4,4,4,4,5,5,5,5,14,15,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,16,17,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,17),3,18,byrow = T))
par(mar=c(4,7,3,1.5))
options(warn=-1) ## suppress warnings due to passing levels to vis.gam
vis.gam(pustule.model,plot.type = "contour",type="response",labcex=.75,contour.col = "black",color="cm",zlim=c(-.35,.35),nCol = 100,too.far = .1,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('pustule area ('*mm^2*')'),ylab="time (days)")
points(delta.pustules$area,delta.pustules$time,pch=".")
par(mar=c(4,1,3,4))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-0.35-.0035,0.35+.035),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-.35,0.35,length.out=101)[i]
  rect(0,ii-.0035,1,ii+.0035,col=cm.colors(101)[i],border = NA)
}
rect(0,-.35-.0035,1,0.35+.0035)
mtext("A",cex=1.25,font=2)
mtext("te(pustule area, time)",side=4,line=.5)
axis(2,cex.axis=1.5)
par(mar=c(4,4.5,3,4))
plot(pustule.model,select = 2,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean temperature (°C)",ylab="s(mean temperature)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 3,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="max. temperature (°C)",ylab="s(max. temperature)")
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 4,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="min. temperature (°C)",ylab="s(min. temperature)")
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 5,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean abs. humidity ('*g/m^3*')'),ylab="s(mean abs. humidity)")
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 6,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('max. abs. humidity ('*g/m^3*')'),ylab="s(max. abs. humidity)")
grid()
mtext("F",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 7,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('min abs. humidity ('*g/m^3*')'),ylab="s(min. abs. humidity)")
grid()
mtext("G",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean solar radiation ('*W/m^2*')'),ylab="s(mean solalr radiation")
grid()
mtext("H",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 9,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="total rainfall (mm)",ylab="s(total rainfall)")
grid()
mtext("I",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 10,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean leaf wetness (%)",ylab="s(mean leaf wetness)")
grid()
mtext("J",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 11,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(tag)")
grid()
mtext("K",adj=1,cex=1.25,font=2)
plot(pustule.model,select = 12,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(site)")
grid()
mtext("L",adj=1,cex=1.25,font=2)


# predict climate change effect

## predict change in pustule area by climate for different pustule sizes

### climate modeled as day matching ~ 50th/75th/90th quantile of mean temp

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,1.5),ylim=c(-.5,.1),type="n",xlab="area (cm)",ylab="pred. change in area (mm^2)",cex.axis=1,cex.lab=1.2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(0,1.5,.1))
  {
    pred.data<-get.pred.data.temp.mean.quantile.pustule.model(day,i)
    
    Xp <- suppressWarnings(predict(pustule.model, newdata = pred.data, exlude='s(site)',type="lpmatrix"))
    beta <- coef(pustule.model) ## posterior mean of coefs
    Vb   <- vcov(pustule.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(pustule.model )$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds

    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)))
    points(i,mean(y),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=0)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

### climate modeled as ~50th quantile hotest day + temperature addition

plot(0,0,xlim=c(0,1.5),ylim=c(-.5,.1),type="n",xlab="area (cm)",ylab="pred. change in area (mm^2)",cex.axis=1,cex.lab=1.2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(0,1.5,.1))
  {
    pred.data<-get.pred.data.temp.mean.quantile.pustule.model(75,i,temp.addition = temp.addition)
    Xp <- predict(pustule.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
    beta <- coef(pustule.model) ## posterior mean of coefs
    Vb   <- vcov(pustule.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(pustule.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)))
    points(i,mean(y),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=0)
}
legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

## predict size trajectory of pustule across observation window
predict.pustule.trajectory<-function(site,temp.addition,color,start.date,end.date,pred.window=1,plot=T,output=F)
{  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  weath.dat<-all.weath[which(all.weath$site==site),]
  temp.rh.dat<-all.temp.rh[which(all.temp.rh$site==site),]
  min.date<-start.date
  max.date<-end.date
  dates<-seq(min.date,max.date,pred.window)
  start.area<-.1
  xcords<-rep(NA,length(dates))
  ycords<-rep(NA,length(dates))
  
  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.area
    xcords.new<-c(1)
    ycords.new<-c(i)
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,i,temp.addition = temp.addition)
      beta <- coef(pustule.model) ## posterior mean of coefs
      Vb   <- vcov(pustule.model) ## posterior  cov of coefs
      mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
      Xp <- predict(pustule.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
      n <-2
      ilink <- family(pustule.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      y<-preds[1]
      i<-i+y
      if(i<0) {i<-0}
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

plot.purple<-t_col("purple",90)
plot.red<-t_col("red",90)
plot.orange<-t_col("orange",90)

par(mar=c(6,6,2,2),mfrow=c(4,2))

### one day ahead projection
start.date<-as.Date("2020-06-23")
end.date<-as.Date("2020-06-27")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,0.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)

start.date<-as.Date("2020-06-28")
end.date<-as.Date("2020-07-02")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)

start.date<-as.Date("2020-07-03")
end.date<-as.Date("2020-07-07")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)

start.date<-as.Date("2020-07-08")
end.date<-as.Date("2020-07-12")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)

start.date<-as.Date("2020-07-13")
end.date<-as.Date("2020-07-17")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)


start.date<-as.Date("2020-07-18")
end.date<-as.Date("2020-07-22")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)

start.date<-as.Date("2020-07-23")
end.date<-as.Date("2020-07-27")
plot(0,0,type="n",xlim=c(1,5),ylim=c(0,.5),ylab='pustule area (mm^2)',xlab="day",cex.lab=1.75,cex.axis=1.75,main=paste0("start date: ",start.date),cex.main=1.5)
dat1<-predict.pustule.trajectory("GM",0,plot.orange,start.date,end.date,pred.window=1,T,T) 
dat2<-predict.pustule.trajectory("GM",1.8,plot.red,start.date,end.date,pred.window=1,T,T) 
dat3<-predict.pustule.trajectory("GM",3.7,plot.purple,start.date,end.date,pred.window=1,T,T) 
points(1:5,colMeans(dat1),type="l",col="orange",lwd=4,lty=2)
points(1:5,colMeans(dat2),type="l",col="red",lwd=4,lty=2)
points(1:5,colMeans(dat3),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",axes=F,bty="n")
legend("center",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=4,cex=3,bty="n")



