library(mgcv)


# load and prep data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/infection intensity data prep.R")
delta.infection.intensity<-subset(delta.infection.intensity,time<=8)

# visualize data
layout(matrix(c(1,2,4,4,3,3,3,3),2,4,byrow = T))
par(mar=c(5,5,3,0.5))
## histograms
hist(delta.infection.intensity$infection.intensity,main="",breaks=100,xlab="infection intensity",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.infection.intensity$infection.intensity.next-delta.infection.intensity$infection.intensity)/delta.infection.intensity$time,main="",breaks=100,xlab="change in infection intensity per day",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories

plot(c(min(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y")),max(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y"))),c(0,max(infection.intensity$infection.intensity)),type="n",xlab="date",ylab="infection intensity",main="infection intensity",cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
i<-0
plot.cols<-sample(rainbow(101))

for (tag in unique(infection.intensity$Tag))
{
  sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
  points(sub.infection.intensity.1$Date,sub.infection.intensity.1$infection.intensity,col=plot.cols[i],type="l",lwd=2)
  i<-i+1
}

#plot(c(min(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y")),max(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%Y"))),c(-1,max(log10(infection.intensity$infection.intensity))),type="n",xlab="date",ylab="log 10 plant infection intensity",main="log 10 plant infection intensity",cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
#i<-0
#plot.cols<-sample(rainbow(101))
#
#for (tag in unique(infection.intensity$Tag))
#{
#  sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
#  points(sub.infection.intensity.1$Date,log10(sub.infection.intensity.1$infection.intensity),col=plot.cols[i],type="l",lwd=2)
#  i<-i+1
#}

## plot change
par(mar=c(5,5,3,0.5))
plot(delta.infection.intensity$infection.intensity,delta.infection.intensity$infection.intensity.next,col="black",xlab = 'infection intensity',ylab='next obs. infection intensity',cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3.25,cex=2)
abline(0,1,lty=2)


# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS"))
{
  set.seed(5708389)
  mod<-gam((infection.intensity.next-infection.intensity)/time~
             te(max.height,infection.intensity)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(mean.daily.rain)+
             s(tag,bs="re")+
             s(site,bs="re"),
           select = T,
           method="REML",
           data=delta.infection.intensity,
           control = list(nthreads=4))
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS")
  
}

## load best model

infection.intensity.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS")

## model checking 
par(mfrow=c(2,2))
gam.check(infection.intensity.model) #Indicates that k should be higher in some smooths. Increasing k to higher value exacerbates the problem for the tensor. As such, we leave k at default values
concurvity(infection.intensity.model,full=F) #no obvious issues

## visualize model
layout(matrix(c(1,1,1,1,2,3,3,3,3,4,4,4,4,10,11,5,5,5,5,6,6,6,6,7,7,7,7,12,13,13,13,8,8,8,8,9,9,9,9,14,14,14),3,14,byrow = T))

par(mar=c(4,7,3,1.5))
options(warn=-1) ## suppress warnings due to passing levels to vis.gam
vis.gam(infection.intensity.model,view=c("max.height","infection.intensity"),plot.type = "contour",type="response",labcex=.75,contour.col = "black",color="cm",zlim=c(-2500,2500),nCol = 100,main="",cex.lab=1.5,cex.axis=1.5,xlab="plant height (cm)",ylab="infection intensity")
points(delta.infection.intensity$max.height,delta.infection.intensity$infection.intensity,pch=".")
par(mar=c(4,1,3,4))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-2500-25,2500+25),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-2500,2500,length.out=101)[i]
  rect(0,ii-25,1,ii+25,col=cm.colors(101)[i],border = NA)
}
rect(0,-2500-25,1,2500+25)
mtext("A",cex=1.25,font=2)
mtext("te(plant height, infection intensity)",side=4,line=.5)
axis(2,cex.axis=1.5)

par(mar=c(4,4.5,3,4))
plot(infection.intensity.model,select = 2,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean temperature (°C)",ylab="s(mean temperature)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(infection.intensity.model,select = 3,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="max. temperature (°C)",ylab="s(max. temperature)")
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(infection.intensity.model,select = 4,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="min. temperature (°C)",ylab="s(min. temperature)")
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(infection.intensity.model,select = 5,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean abs. humidity ('*g/m^3*')'),ylab="s(mean abs. humidity)")
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(infection.intensity.model,select = 6,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean daily rainfall (mm)",ylab="s(mean daily rainfall)")
grid()
mtext("I",adj=1,cex=1.25,font=2)
plot(infection.intensity.model,select = 7,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(tag)")
grid()
mtext("K",adj=1,cex=1.25,font=2)
plot(infection.intensity.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(site)")
grid()
mtext("L",adj=1,cex=1.25,font=2)

## visualize model predictions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
par(mfrow=c(1,2),mar=c(6,6,6,6))
plot(0,0,xlim=c(0,5000),ylim=c(-1000,1000),type="n",xlab="infection intensity",ylab="pred. change in infection intensity",cex.axis=2,cex.lab=2)
day.indicies<-c(75,113,135)
colors<-c("orange","red","purple")
for(day in day.indicies)
{
  index<-which(day.indicies==day)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in c(1,10,100,500,1000,5000))
  {
    pred.data<-get.pred.data.temp.mean.quantile.infection.intensity.model(day,dummy.data.infection.intensity=i,dummy.data.height = 15)
    
    Xp <- predict(infection.intensity.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
    beta <- coef(infection.intensity.model) ## posterior mean of coefs
    Vb   <- vcov(infection.intensity.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(infection.intensity.model)$linkinv
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

plot(0,0,xlim=c(0,5000),ylim=c(-1000,1000),type="n",xlab="infection intensity",ylab="pred. change in infection intensity",cex.axis=2,cex.lab=2)
temp.additions<-c(0,1.8,3.7)
colors<-c("orange","red","purple")
for(temp.addition in temp.additions)
{
  index<-which(temp.additions==temp.addition)
  lower<-data.frame(x=numeric(),y=numeric())
  upper<-data.frame(x=numeric(),y=numeric())
  for(i in c(1,10,100,500,1000,5000))
  {
    pred.data<-get.pred.data.temp.mean.quantile.infection.intensity.model(75,dummy.data.infection.intensity=i,dummy.data.height = 15)
    
    Xp <- predict(infection.intensity.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
    beta <- coef(infection.intensity.model) ## posterior mean of coefs
    Vb   <- vcov(infection.intensity.model) ## posterior  cov of coefs
    n <- 100
    mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
    preds <- rep(NA, n)
    ilink <- family(infection.intensity.model)$linkinv
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


predict.infection.intensity.trajectory<-function(site,temp.addition,plant.height,color,pred.window=1,plot=T,output=F)
{  
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  weath.dat<-all.weath[which(all.weath$site==site),]
  temp.rh.dat<-all.temp.rh[which(all.temp.rh$site==site),]
  min.date<-max(min(unique(as.Date(weath.dat$date))),min(unique(as.Date(temp.rh.dat$date.time))))
  max.date<-min(max(unique(as.Date(weath.dat$date))),max(unique(as.Date(temp.rh.dat$date.time))))
  dates<-seq(min.date,max.date,pred.window)
  start.infection.intensity<-1
  xcords<-rep(NA,length(dates)) #time values
  ycords<-rep(NA,length(dates)) #inf intensity values
  hcords<-rep(NA,length(dates)) #height values

  
  for(j in 1:100) #simulation iteration
  {
    reps<-1
    i<-start.infection.intensity
    h<-plant.height
    xcords.new<-c(1)
    ycords.new<-c(i)
    hcords.new<-c(h)
    for(k in 1:(length(dates)-1)) #date index
    {
      date0<-as.POSIXct(dates[k])
      date1<-as.POSIXct(dates[k+1])
      pred.data<-get.pred.data(site,date0,date1,i,dummy.data.max.height=h,temp.addition = temp.addition)
      pred.data.h<-pred.data
      colnames(pred.data.h)[which(colnames(pred.data.h)=="max.height")]<-"height"
      colnames(pred.data.h)[which(colnames(pred.data.h)=="infection.intensity")]<-"inf.intens"
      
      Xp <- predict(infection.intensity.model, newdata = pred.data, exclude=c("s(site)","s(tag)"),type="lpmatrix")
      beta <- coef(infection.intensity.model) ## posterior mean of coefs
      Vb   <- vcov(infection.intensity.model) ## posterior  cov of coefs
      n <-2
      mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
      ilink <- family(infection.intensity.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      y<-preds[1]
      i<-i+y*pred.window
      if(i<.1) {i<-.1}
      
      beta.h <- coef(infection.intensity.model) ## posterior mean of coefs
      Vb.h   <- vcov(infection.intensity.model) ## posterior  cov of coefs
      n.h <-2
      mrand.h <- mvrnorm(n.h, beta.h, Vb.h) ## simulate n rep coef vectors from posterior
      Xp.h <- predict(infection.intensity.model, newdata = pred.data.h, exclude=c("s(site)","s(tag)"),type="lpmatrix")
      ilink <- family(infection.intensity.model)$linkinv
      preds.h <- rep(NA,n.h)
      for (l in seq_len(n.h)) { 
        preds.h[l]   <- ilink(Xp.h %*% mrand.h[l, ])[1]
      }
      h<-h+preds.h[1]
      
      
      reps<-reps+pred.window
      xcords.new<-c(xcords.new,reps)
      ycords.new<-c(ycords.new,i)
      hcords.new<-c(hcords.new,h)
    }
    xcords<-rbind(xcords,xcords.new)
    ycords<-rbind(ycords,ycords.new)
    hcords<-rbind(hcords,hcords.new)
    print(j)
  }
  xcords<-xcords[-1,]
  ycords<-ycords[-1,]
  hcords<-hcords[-1,]
  
  if(plot)
  {
    for(k in 1:dim(xcords)[1])
    {
      points(xcords[k,],ycords[k,],type="l",col=color) 
    } 
  }
  if(output)
  {
    list(ycords,hcords)
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

plant.height<-15

### one day ahead projection
#layout(matrix(c(1,1,2),1,3,byrow = T))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,36),ylim=c(0,10000),ylab=expression(log[10]*'plant inf. intensity'),xlab="day",cex.lab=1.5,cex.axis=1.5,main="1 day ahead")
dat1<-predict.infection.intensity.trajectory("GM",0,plant.height=plant.height,pred.window=7,plot.orange,T,T) 
points(seq(1,36,length.out=length(colMeans(dat1[[1]]))),colMeans(dat1[[1]]),type="l",col="orange",lwd=4,lty=2)

dat2<-predict.infection.intensity.trajectory("GM",1.8,plant.height=plant.height,pred.window=7,plot.red,T,T) 
points(seq(1,36,length.out=length(colMeans(dat1[[1]]))),colMeans(dat2[[1]]),type="l",col="red",lwd=4,lty=2)

dat3<-predict.infection.intensity.trajectory("GM",3.7,plant.height=plant.height,pred.window=7,plot.purple,T,T) 
points(seq(1,36,length.out=length(colMeans(dat1[[1]]))),colMeans(dat3[[1]]),type="l",col="purple",lwd=4,lty=2)

plot(0,0,type="n",xlim=c(1,36),ylim=c(0,5000),ylab=expression(log[10]*'plant inf. intensity'),xlab="day",cex.lab=1.5,cex.axis=1.5,main="7 days ahead")
polygon(c(seq(1,36,length.out=length(colMeans(dat1[[1]]))),seq(36,1,length.out=length(colMeans(dat1[[1]])))),c(apply(dat1[[1]],2,quantile,probs=.05),rev(apply(dat1[[1]],2,quantile,probs=.95))),col=plot.orange,density=100)
polygon(c(seq(1,36,length.out=length(colMeans(dat1[[1]]))),seq(36,1,length.out=length(colMeans(dat1[[1]])))),c(apply(dat2[[1]],2,quantile,probs=.05),rev(apply(dat2[[1]],2,quantile,probs=.95))),col=plot.red,density=100)
polygon(c(seq(1,36,length.out=length(colMeans(dat1[[1]]))),seq(36,1,length.out=length(colMeans(dat1[[1]])))),c(apply(dat3[[1]],2,quantile,probs=.05),rev(apply(dat3[[1]],2,quantile,probs=.95))),col=plot.purple,density=100)
points(seq(1,36,length.out=length(colMeans(dat1[[1]]))),colMeans(dat1[[1]]),type="l",col="orange",lwd=4,lty=2)
points(seq(1,36,length.out=length(colMeans(dat1[[1]]))),colMeans(dat2[[1]]),type="l",col="red",lwd=4,lty=2)
points(seq(1,36,length.out=length(colMeans(dat1[[1]]))),colMeans(dat3[[1]]),type="l",col="purple",lwd=4,lty=2)


