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

  mod0<-gam(area.next~s(area,bs="cs",k=4)+
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
            data=delta.pustules)

  summary(mod0) #indicates that max.abs.hum, mean.solar, and tot.rain are not significant predictors
  
  mod1<-gam(area.next~s(area,bs="cs",k=4)+
              s(mean.temp,by=time,bs="cs",k=4)+
              s(max.temp,by=time,bs="cs",k=4)+
              s(min.temp,by=time,bs="cs",k=4)+
              s(mean.abs.hum,by=time,bs="cs",k=4)+
              #s(max.abs.hum,by=time,bs="cs",k=4)+
              s(min.abs.hum,by=time,bs="cs",k=4)+
              #s(mean.solar,by=time,bs="cs",k=4)+
              s(mean.wetness,by=time,bs="cs",k=4)+
              #s(tot.rain,by=time,bs="cs",k=4)+
              s(site,bs="re",k=4),
            data=delta.pustules)
  summary(mod1) # All now significant. Stepwise removal yields the same result.
  
  saveRDS(mod1,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
}

## load best model

pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")

## model checking
plot(pustule.model,scale=0,pages=1) #plot smooths
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
