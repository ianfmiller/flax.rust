library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant data prep.R")
delta.plants<-subset(delta.plants,time<=10)

# visualize data

## histograms
par(mfrow=c(2,1))
hist(plants$plant.inf.intens,main="plant infection intensity",breaks=100,xlab="plant infection intensity")
hist(delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,main="change in plant infection intensity",breaks=100,xlab="change in plant infection intensity")

## plot trajectories
par(mfrow=c(1,1))
plot(c(min(as.Date(plants$Date,tryFormats = "%m/%d/%Y")),max(as.Date(plants$Date,tryFormats = "%m/%d/%Y"))),c(-1,max(log10(plants$plant.inf.intens))),type="n",xlab="date",ylab="plant infection intensity",main="plant infection intensity",cex.lab=2,cex.axis=2,cex.main=2)
i<-0
plot.cols<-sample(rainbow(101))

for (tag in unique(plants$Tag))
{
  sub.plants.1<-plants[which(plants$Tag==tag),]
  points(sub.plants.1$Date,log10(sub.plants.1$plant.inf.intens),col=plot.cols[i],type="l",lwd=2)
  i<-i+1
}


## plot change
par(mfrow=c(2,1))
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next,col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity",xlim=c(0,1000),ylim=c(0,200))
abline(0,1)
plot(delta.plants$plant.inf.intens,delta.plants$plant.inf.intens.next-delta.plants$plant.inf.intens,col="grey",xlab = "plant infection intensity",ylab="change in plant infection intensity")
abline(h=0)

## plot log10 change
### relationship doesn't appear obviously linear near 0--take polynomial fitting approach
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next),col="black",xlab = "plant infection intensity",ylab="next obs. plant infection intensity",cex.lab=2,cex.axis=2,cex.main=2,main="N = 327")
abline(0,1,lty=2)

## test gam fit
par(mfrow=c(1,1))
plot(log10(delta.plants$plant.inf.intens),log10(delta.plants$plant.inf.intens.next),col="grey",xlab = "plant infection intensity",ylab="next obs. plant infection intensity")
test.mod<-gam(log10(plant.inf.intens.next)~s(log10(plant.inf.intens))+s(site,bs="re"),data = delta.plants)
points(log10(delta.plants$plant.inf.intens),predict(test.mod,exclude = 's(site)'),col="red")


# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("log10(plant.inf.intens.next) ~ s(log10(plant.inf.intens))",predictors[x],'s(site,bs="re")'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(gam(log10(plant.inf.intens.next)~s(log10(plant.inf.intens))+s(site,bs="re"),data = delta.plants)) #cutoff to limit memory usage
  #AIC.benchmark<- 8025
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-gam(model.set[[i]],data=delta.plants,REML = F))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")
}

## load best model

plants.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")

# predict climate change effect
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
library("MASS")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(-1,5),ylim=c(-1,2),type="n",xlab="plant infection intensity",ylab="pred. change in plant infection intensity",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in 10^seq(-1,5,.5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plants.model(day,i)
    Xp <- predict(plants.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta <- coef(plants.model) ## posterior mean of coefs
    Vb   <- vcov(plants.model) ## posterior  cov of coefs
    n <- 10000
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plants.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    lower=rbind(lower,data.frame(x=log10(i),y=quantile(y,.05)-log10(i)))
    upper=rbind(upper,data.frame(x=log10(i),y=quantile(y,.95)-log10(i)))
    points(log10(i),mean(y)-log10(i),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=2,bty="n")

plot(0,0,xlim=c(-1,5),ylim=c(-1,2),type="n",xlab="plant infection intensity",ylab="pred. change in plant infection intensity",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in 10^seq(-1,5,.5))
  {
    pred.data<-get.pred.data.temp.mean.quantile.plants.model(75,i,temp.addition = temp.addition)
    Xp <- predict(plants.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta <- coef(plants.model) ## posterior mean of coefs
    Vb   <- vcov(plants.model) ## posterior  cov of coefs
    n <- 10000
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(plants.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])
    }
    y<-preds
    lower=rbind(lower,data.frame(x=log10(i),y=quantile(y,.05)-log10(i)))
    upper=rbind(upper,data.frame(x=log10(i),y=quantile(y,.95)-log10(i)))
    points(log10(i),mean(y)-log10(i),col=colors[index],pch=16,cex=2)
  }
  polygon<-rbind(lower,upper[dim(upper)[1]:1,])
  polygon(polygon$x,polygon$y,col=colors[index],density=25)
}
legend("topright",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),pch=16,cex=2,bty="n")



predict.plant.inf.trajectory<-function(maxreps,day,temp.addition,color,plot=T,output=F)
{  
  reps<-1
  i<-rep(10^-1,times=1000)
  xcords<-rep(1,length(i))
  ycords<-i
  while(reps<maxreps)
  {
    pred.data<-get.pred.data.temp.mean.quantile.plants.model(day,i,temp.addition = temp.addition)
    Xp <- predict(plants.model, newdata = pred.data, exlude="s(site)",type="lpmatrix")
    beta <- coef(plants.model) ## posterior mean of coefs
    Vb   <- vcov(plants.model) ## posterior  cov of coefs
    n <- length(i)
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA,n)
    ilink <- family(plants.model)$linkinv
    for (j in seq_len(n)) { 
      preds[j]   <- ilink(Xp %*% mrand[j, ])[1]
    }
    y<-preds
    i<-10^y
    reps<-reps+1
    xcords<-cbind(xcords,rep(reps,length(i)))
    ycords<-cbind(ycords,i)
  }
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

plot.purple<-t_col("purple",95)
plot.red<-t_col("red",95)
plot.orange<-t_col("orange",95)

par(mfrow=c(1,1),mar=c(6,6,6,6))
plot(0,0,type="n",xlim=c(1,6),ylim=c(0,10),ylab=expression(log[10]*' plant infection intensity'),xlab="week",cex.lab=1.5,cex.axis=1.5)

dat<-predict.plant.inf.trajectory(6,75,0,plot.orange,T,T) 
points(1:6,colMeans(dat),type="l",col="orange",lwd=3,lty=2)

dat<-predict.plant.inf.trajectory(6,113,0,plot.red,T,T) 
points(1:6,colMeans(dat),type="l",col="red",lwd=3,lty=2)

dat<-predict.plant.inf.trajectory(6,135,0,plot.purple,T,T) 
points(1:6,colMeans(dat),type="l",col="purple",lwd=3,lty=2)

legend("topright",legend = c("50% quantile hottest days","75% quantile hottest days","90% quantile hottest days"),col = c("orange","red","purple"),pch=16,cex=1.5)


plot(0,0,type="n",xlim=c(1,6),ylim=c(0,5),ylab=expression(log[10]*' plant infection intensity'),xlab="week",cex.lab=1.5,cex.axis=1.5)
dat<-predict.plant.inf.trajectory(6,75,0,plot.orange,T,T) 
points(1:6,colMeans(dat),type="l",col="orange",lwd=4,lty=2)

dat<-predict.plant.inf.trajectory(6,75,1.8,plot.red,T,T) 
points(1:6,colMeans(dat),type="l",col="red",lwd=4,lty=2)

dat<-predict.plant.inf.trajectory(6,75,3.7,plot.purple,T,T) 
points(1:6,colMeans(dat),type="l",col="purple",lwd=4,lty=2)

legend("topleft",legend = c("+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col = c("orange","red","purple"),lty=2,lwd=2,cex=1.5)


