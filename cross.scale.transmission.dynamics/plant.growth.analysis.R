library(mgcv)
library(lme4)
library(lmerTest)
library(progress)
library(gratia)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.height<-subset(delta.height,time<=7)

# visualize data
par(mfrow =c(1,1))
colors <- ifelse(delta.height$height.next > delta.height$height, "green", "red")
plot(height.next-height~height, data = delta.height, col = colors,xlab="height",ylab="change in heght")


## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{

  mod0<-gam(height.next-height~s(height,by=time,bs="cs",k=4)+
              s(inf.intens,by=time,bs="cs",k=4)+
              s(mean.temp,by=time,bs="cs",k=4)+
              s(max.temp,by=time,bs="cs",k=4)+
              s(min.temp,by=time,bs="cs",k=4)+
              s(mean.abs.hum,by=time,bs="cs",k=4)+
              s(max.abs.hum,by=time,bs="cs",k=4)+
              s(min.abs.hum,by=time,bs="cs",k=4)+
              s(tot.rain,by=time,bs="cs",k=4)+
              s(mean.solar,by=time,bs="cs",k=4)+
              s(mean.soil.moisture,by=time,bs="cs",k=4)+
              s(site,bs="re",k=4)
            ,data=delta.height)
  summary(mod0) #indicates that all but min.temp, max.abs.hum, and min.abs.hum are not significant. min.temp is marginally significant in mod0 but becomes fully significaant in mod1.
  
  mod1<-gam(height.next-height~0+s(height,by=time,bs="cs",k=4)+
              #s(inf.intens,by=time,bs="cs",k=4)+
              #s(mean.temp,by=time,bs="cs",k=4)+
              #s(max.temp,by=time,bs="cs",k=4)+
              #s(min.temp,by=time,bs="cs",k=4)+
              #s(mean.abs.hum,by=time,bs="cs",k=4)+
              s(max.abs.hum,by=time,bs="cs",k=4)+
              #s(min.abs.hum,by=time,bs="cs",k=4)+
              #s(tot.rain,by=time,bs="cs",k=4)+
              #s(mean.solar,by=time,bs="cs",k=4)+
              #s(mean.soil.moisture,by=time,bs="cs",k=4)+
              s(site,bs="re",k=4)
            ,data=delta.height)
  summary(mod1) #indicates that all predictors are now significant

  saveRDS(mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
draw(plant.growth.model) #plot smooths
#gam.check(plant.growth.model) #indicates no more knots needed

# model vis
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
library("MASS")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(5,80),ylim=c(-2,2),type="n",xlab="plant height",ylab="pred. change in plant height per day",cex.axis=1,cex.lab=1.2)
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

plot(0,0,xlim=c(5,80),ylim=c(-2,2),type="n",xlab="plant infection intensity",ylab="pred. change in plant infection intensity",cex.axis=1,cex.lab=1.2)
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

  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.height
    xcords.new<-c(1)
    ycords.new<-c(i)
    
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,dummy.data = NA,dummy.data.max.height=i,temp.addition = temp.addition) #includes some irrelevant and meaningless predictors
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
      #if(i<.1) {i<-.1}

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
