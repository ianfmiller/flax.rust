library(mgcv)
library(lme4)
library(lmerTest)
library(progress)
library(gratia)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.height<-subset(delta.height,time<=7)

# visualize data

layout(matrix(c(1,2,4,4,3,3,3,3),2,4,byrow = T))
par(mar=c(5,5,3,0.5))

## histograms
hist(delta.height$height,main="",breaks=100,xlab="plant height",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.height$height.next-delta.height$height)/delta.height$time,main="",breaks=100,xlab="change in plant height per day",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories
par(mar=c(2,5,0,0.5))
plot(c(min(plant.heights$date),max(plant.heights$date)),c(0,max(plant.heights$max.height)),type="n",xlab="",ylab="plant height",cex.lab=2,cex.axis=2)
mtext("D",side=3,adj=1,line=-3,cex=2)
i<-0

plot.cols<-sample(rainbow(180))

for (i in 1:length(unique(plant.heights$tag)))
{
  tag<-unique(plant.heights$tag)[i]
  sub.heights<-plant.heights[which(plant.heights$tag==tag),]
  sub.heights<-sub.heights[order(sub.heights$date),]
  points(sub.heights$date,sub.heights$max.height,col=plot.cols[i],type="l",lwd=.5)

}

## plot change
par(mar=c(5,5,3,0.5))
plot(delta.height$height,delta.height$height.next,col="black",xlab = "observed height",ylab="next observed height",cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3.25,cex=2)
abline(0,1,lty=2)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{
  mod<-gam((height.next-height)/time~0+
             te(height,inf.intens)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(max.abs.hum)+
             s(min.abs.hum)+
             s(mean.solar)+
             s(mean.daily.rain)+
             s(mean.wetness)+
             s(mean.soil.moisture)+
             s(tag,bs="re")+
             s(site,bs="re"),
           select = T,
           method="REML",
           data=delta.height,
           control = list(nthreads=4))
  summary(mod)

  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
par(mfrow=c(2,2))
gam.check(plant.growth.model) #indicates that basis dimension is sufficient
concurvity(plant.growth.model,full=F) #no obvious issues

## visualize model

layout(matrix(c(15,15,15,1,1,1,1,2,3,3,3,3,4,4,4,4,16,16,16,16,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14),3,20,byrow = T))

par(mar=c(4,7,3,1.5))
options(warn=-1) ## suppress warnings due to passing levels to vis.gam
vis.gam(plant.growth.model,view=c("height","inf.intens"),plot.type = "contour",type="response",labcex=.75,contour.col = "black",color="cm",zlim=c(-.75,.75),nCol = 100,main="",cex.lab=1.5,cex.axis=1.5,xlab="plant height (cm)",ylab="infection intensity")
points(delta.height$height,delta.height$inf.intens,pch=".")
par(mar=c(4,1,3,4))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-.75-0.0075,.75+0.0075),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-0.75,0.75,length.out=101)[i]
  rect(0,ii-0.0075,1,ii+0.0075,col=cm.colors(101)[i],border = NA)
}
rect(0,-0.75-0.0075,1,0.75+0.0075)
mtext("A",cex=1.25,font=2)
mtext("te(plant height, time)",side=4,line=.5)
axis(2,cex.axis=1.5)

par(mar=c(4,4.5,3,4))
plot(plant.growth.model,select = 2,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean temperature (°C)",ylab="s(mean temperature)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 3,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="max. temperature (°C)",ylab="s(max. temperature)")
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 4,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="min. temperature (°C)",ylab="s(min. temperature)")
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 5,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean abs. humidity ('*g/m^3*')'),ylab="s(mean abs. humidity)")
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 6,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('max. abs. humidity ('*g/m^3*')'),ylab="s(max. abs. humidity)")
grid()
mtext("F",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 7,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('min abs. humidity ('*g/m^3*')'),ylab="s(min. abs. humidity)")
grid()
mtext("G",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean solar radiation ('*W/m^2*')'),ylab="s(mean solalr radiation")
grid()
mtext("H",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 9,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="total rainfall (mm)",ylab="s(total rainfall)")
grid()
mtext("I",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 10,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean leaf wetness (%)",ylab="s(mean leaf wetness)")
grid()
mtext("J",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 11,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean soil moisture ('*m^3*' '*H[2]*0/m^3*' soil)'),ylab="s(mean soil moisture)")
grid()
mtext("K",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 12,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(tag)")
grid()
mtext("L",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 13,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(site)")
grid()
mtext("M",adj=1,cex=1.25,font=2)


# model vis
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
library("MASS")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(5,80),ylim=c(-5,5),type="n",xlab="plant height",ylab="pred. change in plant height per day",cex.axis=1,cex.lab=1.2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(5,80,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plant.growth.model(day,dummy.data.height=i)
    pred.data<-cbind(pred.data,inf.intens=0)
    Xp <- predict(plant.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta <- coef(plant.growth.model) ## posterior mean of coefs
    Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
    n <- 10000
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plant.growth.model)$linkinv
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

plot(0,0,xlim=c(5,80),ylim=c(-5,5),type="n",xlab="plant infection intensity",ylab="pred. change in plant infection intensity",cex.axis=1,cex.lab=1.2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(5,80,5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plant.growth.model(75,dummy.data.height=i,temp.addition = temp.addition)
    pred.data<-cbind(pred.data,inf.intens=0)
    Xp <- predict(plant.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta <- coef(plant.growth.model) ## posterior mean of coefs
    Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
    n <- 10000
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plant.growth.model)$linkinv
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


predict.plant.height.traj<-function(site,temp.addition,color,pred.window=1,plot=T,output=F)
{  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  weath.dat<-all.weath[which(all.weath$site==site),]
  temp.rh.dat<-all.temp.rh[which(all.temp.rh$site==site),]
  min.date<-max(min(unique(as.Date(weath.dat$date))),min(unique(as.Date(temp.rh.dat$date.time))))
  max.date<-min(max(unique(as.Date(weath.dat$date))),max(unique(as.Date(temp.rh.dat$date.time))))
  dates<-seq(min.date,max.date,pred.window)
  start.height<-20
  xcords<-rep(NA,length(dates)) #time values
  ycords<-rep(NA,length(dates)) #height values

  for(j in 1:10) #simulation iteration
  {
    reps<-1
    i<-start.height
    xcords.new<-c(1)
    ycords.new<-c(i)
    
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,dummy.data = NA,dummy.data.max.height=i,temp.addition = temp.addition)
      pred.data<-cbind(pred.data,inf.intens=0)
      colnames(pred.data)[which(colnames(pred.data)=="max.height")]<-"height"
      
      beta <- coef(plant.growth.model) ## posterior mean of coefs
      Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
      n <-2
      mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
      Xp <- predict(plant.growth.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
      ilink <- family(plant.growth.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      y<-preds[1]
      i<-i+y
      if(i<5) {i<-5}

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
    list(ycords)
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
plot(0,0,type="n",xlim=c(0,36),ylim=c(0,60),ylab='plant height (cm)',xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
dat1<-predict.plant.height.traj("GM",0,plot.orange,pred.window=1,T,T) 
points(1:36,colMeans(dat1[[1]]),type="l",col="orange",lwd=4,lty=2)

dat2<-predict.plant.height.traj("GM",1.8,plot.red,pred.window=1,T,T) 
points(1:36,colMeans(dat2[[1]]),type="l",col="red",lwd=4,lty=2)

dat3<-predict.plant.height.traj("GM",3.7,plot.purple,pred.window=1,T,T) 
points(1:36,colMeans(dat3[[1]]),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",xlim=c(1,36),ylim=c(0,60),ylab='plant inf. intens.',xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
points(1:36,colMeans(dat1[[1]]),type="l",col="orange",lwd=4,lty=2)
points(1:36,colMeans(dat2[[1]]),type="l",col="red",lwd=4,lty=2)
points(1:36,colMeans(dat3[[1]]),type="l",col="purple",lwd=4,lty=2)
