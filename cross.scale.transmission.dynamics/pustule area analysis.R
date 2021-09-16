library(mgcv)
library(lme4)
library(lmerTest)
library(progress)
library(MASS)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")

delta.pustules<-subset(delta.pustules,time<=7)

# visualize data

## histograms
par(mfrow=c(2,1))
hist(delta.pustules$area,main="pustule area",breaks=100,xlab="area")
hist(delta.pustules$area.next-delta.pustules$area,main="change in pustule area",breaks=100,xlab="change in area")

## plot trajectories
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="date",ylab="pustule area")
i<-0

plot.cols<-sample(rainbow(2042))

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

## plot change
plot(delta.pustules$area,delta.pustules$area.next,col="grey",xlab = "area",ylab="next obs. area")
abline(0,1)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS"))
{

  mod0<-gam(area.next~area+
              te(mean.temp,area,by=time,bs="cs",k=4)+
              te(max.temp,area,by=time,bs="cs",k=4)+
              te(min.temp,area,by=time,bs="cs",k=4)+
              te(mean.abs.hum,area,by=time,bs="cs",k=4)+
              te(max.abs.hum,area,by=time,bs="cs",k=4)+
              te(min.abs.hum,area,by=time,bs="cs",k=4)+
              te(mean.solar,area,by=time,bs="cs",k=4)+
              te(mean.wetness,area,by=time,bs="cs",k=4)+
              te(tot.rain,area,by=time,bs="cs",k=4)+
              te(site,bs="re",k=4),
            data=delta.pustules)

  summary(mod0) #indicates that max.abs.hum, mean.solar, and tot.rain are not significant predictors
  
  mod1<-gam(area.next~area+
              te(mean.temp,area,by=time,bs="cs",k=4)+
              te(max.temp,area,by=time,bs="cs",k=4)+
              te(min.temp,area,by=time,bs="cs",k=4)+
              te(mean.abs.hum,area,by=time,bs="cs",k=4)+
              #te(max.abs.hum,area,by=time,bs="cs",k=4)+
              te(min.abs.hum,area,by=time,bs="cs",k=4)+
              te(mean.solar,area,by=time,bs="cs",k=4)+
              te(mean.wetness,area,by=time,bs="cs",k=4)+
              te(tot.rain,area,by=time,bs="cs",k=4)+
              te(site,bs="re",k=4),
            data=delta.pustules)
  
  summary(mod1) #indicates that max.abs.hum, mean.solar, and tot.rain are not significant predictors
  
  mod2<-gam(area.next~area+
              te(mean.temp,area,by=time,bs="cs",k=4)+
              #te(max.temp,area,by=time,bs="cs",k=4)+
              te(min.temp,area,by=time,bs="cs",k=4)+
              te(mean.abs.hum,area,by=time,bs="cs",k=4)+
              #te(max.abs.hum,area,by=time,bs="cs",k=4)+
              te(min.abs.hum,area,by=time,bs="cs",k=4)+
              te(mean.solar,area,by=time,bs="cs",k=4)+
              te(mean.wetness,area,by=time,bs="cs",k=4)+
              te(tot.rain,area,by=time,bs="cs",k=4)+
              te(site,bs="re",k=4),
            data=delta.pustules)
  
  summary(mod2) #indicates that max.abs.hum, mean.solar, and tot.rain are not significant predictors
  
  saveRDS(mod2,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
}

## load best model

pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")

## model checking
draw(pustule.model,dist=.5) #plot smooths
#gam.check(pustule.model) #indicates more knots needed
#more.knots.mod<-gam(area.next ~ s(area, bs = "cs", k = 20) + s(mean.temp, by = time, bs = "cs", k = 20) + s(max.temp, by = time, bs = "cs", k = 20) + s(min.temp, by = time, bs = "cs", k = 20) + s(mean.abs.hum, by = time, bs = "cs", k = 20) + s(min.abs.hum, by = time, bs = "cs", k = 20) + s(mean.wetness, by = time, bs = "cs", k = 20) + s(site, bs = "re", k = 20),data=delta.pustules)
#gam.check(more.knots.mod) #indicates that increasing number of knots doesn't solve the issue of unevenly distributed residuals. This indicates that the data is the root cause and original model is OK

## better visualize model
par(mfrow=c(2,4))
plot(pustule.model,scale=0,select=1)
vis.gam(pustule.model,view = c("mean.temp","time"),n.grid=30,plot.type = "contour",zlim=c(-.02,.15),color="topo",contour.col = "black")
vis.gam(pustule.model,view = c("max.temp","time"),n.grid=30,plot.type = "contour",zlim=c(-.02,.15),color="topo",contour.col = "black")
vis.gam(pustule.model,view = c("min.temp","time"),n.grid=30,plot.type = "contour",zlim=c(-.02,.15),color="topo",contour.col = "black")
vis.gam(pustule.model,view = c("mean.abs.hum","time"),n.grid=30,plot.type = "contour",zlim=c(-.02,.15),color="topo",contour.col = "black")
vis.gam(pustule.model,view = c("min.abs.hum","time"),n.grid=30,plot.type = "contour",zlim=c(-.02,.15),color="topo",contour.col = "black")
vis.gam(pustule.model,view = c("mean.wetness","time"),n.grid=30,plot.type = "contour",zlim=c(-.02,.15),color="topo",contour.col = "black")
plot(pustule.model,scale=0,select=8)

par(mfrow=c(2,4))
plot(pustule.model,scale=0,select=1)
vis.gam(pustule.model,view = c("mean.temp","time"),n.grid=30,plot.type = "persp",zlim=c(-.02,.15),se=1,theta=45,phi=15,ticktype="detailed")
vis.gam(pustule.model,view = c("max.temp","time"),n.grid=30,plot.type = "persp",zlim=c(-.02,.15),se=1,theta=45,phi=15,ticktype="detailed")
vis.gam(pustule.model,view = c("min.temp","time"),n.grid=30,plot.type = "persp",zlim=c(-.02,.15),se=1,theta=45,phi=15,ticktype="detailed")
vis.gam(pustule.model,view = c("mean.abs.hum","time"),n.grid=30,plot.type = "persp",zlim=c(-.02,.15),se=1,theta=45,phi=15,ticktype="detailed")
vis.gam(pustule.model,view = c("min.abs.hum","time"),n.grid=30,plot.type = "persp",zlim=c(-.02,.15),se=1,theta=45,phi=15,ticktype="detailed")
vis.gam(pustule.model,view = c("mean.wetness","time"),n.grid=30,plot.type = "persp",zlim=c(-.02,.15),se=1,theta=45,phi=15,ticktype="detailed")
plot(pustule.model,scale=0,select=8)

# predict climate change effect

## predict change in pustule area by climate for different pustule sizes

### climate modeled as day matching ~ 50th/75th/90th quantile of mean temp

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,.1),ylim=c(-.02,.02),type="n",xlab="area (cm)",ylab="pred. change in area (cm)",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(0,.1,.01))
  {
    pred.data<-get.pred.data.temp.mean.quantile.pustule.model(day,i)
    
    Xp <- predict(pustule.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
    beta <- coef(pustule.model) ## posterior mean of coefs
    Vb   <- vcov(pustule.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plant.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds

    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=0)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

### climate modeled as ~50th quantile hotest day + temperature addition

plot(0,0,xlim=c(0,.1),ylim=c(-.02,.02),type="n",xlab="area (cm)",ylab="pred. change in area (cm)",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in seq(0,.1,.01))
  {
    pred.data<-get.pred.data.temp.mean.quantile.pustule.model(75,i,temp.addition = temp.addition)
    Xp <- predict(pustule.model, newdata = pred.data, exlude='s(site)',type="lpmatrix")
    beta <- coef(pustule.model) ## posterior mean of coefs
    Vb   <- vcov(pustule.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plant.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    lower=rbind(lower,data.frame(x=i,y=quantile(y,.05)-i))
    upper=rbind(upper,data.frame(x=i,y=quantile(y,.95)-i))
    points(i,mean(y)-i,col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=0)
}
legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),pch=16,cex=1,bty="n")

## predict size trajectory of pustule across observation window
predict.pustule.trajectory<-function(site,temp.addition,color,pred.window=1,plot=T,output=F)
{  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  weath.dat<-all.weath[which(all.weath$site==site),]
  temp.rh.dat<-all.temp.rh[which(all.temp.rh$site==site),]
  min.date<-max(min(unique(as.Date(weath.dat$date))),min(unique(as.Date(temp.rh.dat$date.time))))
  max.date<-min(max(unique(as.Date(weath.dat$date))),max(unique(as.Date(temp.rh.dat$date.time))))
  dates<-seq(min.date,max.date,pred.window)
  start.area<-.01
  xcords<-rep(NA,length(dates))
  ycords<-rep(NA,length(dates))
  
  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.area
    xcords.new<-c(1)
    ycords.new<-c(i)
    beta <- coef(pustule.model) ## posterior mean of coefs
    Vb   <- vcov(pustule.model) ## posterior  cov of coefs
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,i,temp.addition = temp.addition)
      Xp <- predict(pustule.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
      n <-2
      ilink <- family(pustule.model)$linkinv
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
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,.13),ylab='pustule area (cm)',xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
dat<-predict.pustule.trajectory("GM",0,pred.window=1,plot.orange,T,T) 
points(1:36,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.pustule.trajectory("GM",1.8,pred.window=1,plot.red,T,T) 
points(1:36,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.pustule.trajectory("GM",3.7,pred.window=1,plot.purple,T,T) 
points(1:36,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)

### two days ahead projection
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,.13),ylab='pustule area (cm)',xlab="day",cex.lab=1.5,cex.axis=1.5,main="2 days ahead")
dat<-predict.pustule.trajectory("GM",0,pred.window=2,plot.orange,T,T) 
points(seq(1,35,2),colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.pustule.trajectory("GM",1.8,pred.window=2,plot.red,T,T) 
points(seq(1,35,2),colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.pustule.trajectory("GM",3.7,pred.window=2,plot.purple,T,T) 
points(seq(1,35,2),colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)

### seven days ahead projection
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,.13),ylab='pustule area (cm)',xlab="week",cex.lab=1.5,cex.axis=1.5,main="1 week ahead")
dat<-predict.pustule.trajectory("GM",0,pred.window=7,plot.orange,T,T) 
points(seq(1,36,7),colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.pustule.trajectory("GM",1.8,pred.window=7,plot.red,T,T) 
points(seq(1,36,7),colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.pustule.trajectory("GM",3.7,pred.window=7,plot.purple,T,T) 
points(seq(1,36,7),colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)





