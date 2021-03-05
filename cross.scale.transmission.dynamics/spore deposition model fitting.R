# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions.R")
## load enviro data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
## load spore dep data
spore.deposition<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts.csv")
spore.deposition[which(spore.deposition$Distance..cm.==0),"Distance..cm."]<-1 #set distance for "0" spore traps to 1cm
## load demog data
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
demog<-demog[which(demog$year==2020),] #subset to 2020
## load within host data
plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")

# functions
predict.kernel.decay.plume<-function(q,k,alphay,c,xtarget,ytarget,wind.data)
{
  tot.dep<-0
  for (i in 1:(dim(wind.data)[1]-1))
  {
    delta.t<-wind.data[i+1,"date"]-wind.data[i,"date"]
    cords<-get.plume.xy(2*pi*correct.wind.degree(wind.data[i,"wind.direction"],site = wind.data[i,"site"])/360,0,0,xtarget,ytarget,plot=F)
    tot.dep<-c(tot.dep,decay.plume(q=q,k=k,s=wind.data[i,"wind.speed"],x=cords[1],y=cords[2],alphay=alphay,c=c))
  }
  sum(tot.dep,na.rm = T)
}

predict.kernel.tilted.plume<-function(q,H,alphay,alphaz,Ws,xtarget,ytarget,wind.data)
{
  tot.dep<-0
  for (i in 1:(dim(wind.data)[1]-1))
  {
    delta.t<-wind.data[i+1,"date"]-wind.data[i,"date"]
    cords<-get.plume.xy(2*pi*correct.wind.degree(wind.data[i,"wind.direction"],site = wind.data[i,"site"])/360,0,0,xtarget,ytarget,plot=F)
    tot.dep<-c(tot.dep,tilted.plume(q=q,H=H,s=wind.data[i,"wind.speed"],x=cords[1],y=cords[2],alphay=alphay,alphaz=alphaz,Ws=Ws))
  }
  sum(tot.dep,na.rm = T)
}

param.search.optim.decay.plume<-function(x,return.vec=F)
{
  cval=x[1]
  alphayval=x[2]
  k=x[3]
  
  val<-c()
  preds<-c()
  obs<-c()
  #for(tag in unique(spore.deposition$Tag))
  for (tag in c(86,88))
  {
    site<-demog[which(demog$tag==tag),"Site"]
    plantx<-demog[which(demog$tag==tag),"X"]+demog[which(demog$tag==tag),"x"]
    planty<-demog[which(demog$tag==tag),"Y"]+demog[which(demog$tag==tag),"y"]
    sub.1.spore.deposition<-spore.deposition[which(spore.deposition$Tag==tag),]
    
    for(date in unique(sub.1.spore.deposition$Date.collected))
    {
      sub.2.spore.deposition<-sub.1.spore.deposition[which(sub.1.spore.deposition$Date.collected==date),]
      deploy.date<-sub.2.spore.deposition[1,"Date.deployed"]
      q<-plants[intersect(which(plants$Tag==tag),which(plants$Date==as.Date(deploy.date,tryFormats = c("%m/%d/%y")))),"plant.inf.intens"]
      
      wind.data<-all.weath[which(all.weath$site==site),]
      wind.data<-wind.data[which(wind.data$date>(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),]
      wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC")+60*60*24*2)),]
      #wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),]
      
      for(j in 1:dim(sub.2.spore.deposition)[1])
      {
        xtarget<-0
        ytarget<-0
        if(sub.2.spore.deposition[j,"Direction"]=="U") {ytarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="R") {xtarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="D") {ytarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="L") {xtarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        
        new.pred<-predict.kernel.decay.plume(q=q,k=k,alphay=alphayval,c=cval,xtarget=xtarget,ytarget=ytarget,wind.data=wind.data)
        new.obs<-sub.2.spore.deposition[j,"Pustules"]/sub.2.spore.deposition[j,"X..squares.counted"]
        val<-c(val,(new.obs-new.pred)^2)
        preds<-c(preds,new.pred)
        obs<-c(obs,new.obs)
      }
    }
  }
  if(return.vec==T) {data.frame(pred=preds,obs=obs)} else {sum(val)}
}

param.search.optim.tilted.plume<-function(x,return.vec=F)
{
  alphayval=x[1]
  alphazval=x[2]
  Wsval=x[3]
  
  val<-c()
  preds<-c()
  obs<-c()
  #for(tag in unique(spore.deposition$Tag))
  for (tag in 88)
  {
    site<-demog[which(demog$tag==tag),"Site"]
    plantx<-demog[which(demog$tag==tag),"X"]+demog[which(demog$tag==tag),"x"]
    planty<-demog[which(demog$tag==tag),"Y"]+demog[which(demog$tag==tag),"y"]
    sub.1.spore.deposition<-spore.deposition[which(spore.deposition$Tag==tag),]
    
    #for(date in unique(sub.1.spore.deposition$Date.collected))
    for(date in "7/16/20")
    {
      sub.2.spore.deposition<-sub.1.spore.deposition[which(sub.1.spore.deposition$Date.collected==date),]
      deploy.date<-sub.2.spore.deposition[1,"Date.deployed"]
      q<-plants[intersect(which(plants$Tag==tag),which(plants$Date==as.Date(deploy.date,tryFormats = c("%m/%d/%y")))),"plant.inf.intens"]
      H<-plants[intersect(which(plants$Tag==tag),which(plants$Date==as.Date(deploy.date,tryFormats = c("%m/%d/%y")))),"max.height"]
      
      wind.data<-all.weath[which(all.weath$site==site),]
      wind.data<-wind.data[which(wind.data$date>(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),]
      wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),]
      
      for(j in 1:dim(sub.2.spore.deposition)[1])
      {
        xtarget<-0
        ytarget<-0
        if(sub.2.spore.deposition[j,"Direction"]=="U") {ytarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="R") {xtarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="D") {ytarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="L") {xtarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        
        new.pred<-predict.kernel.tilted.plume(q=q,H=H,alphay=alphayval,alphaz=alphazval,Ws=Wsval,xtarget=xtarget,ytarget=ytarget,wind.data=wind.data)
        new.obs<-sub.2.spore.deposition[j,"Pustules"]/sub.2.spore.deposition[j,"X..squares.counted"]
        val<-c(val,(new.obs-new.pred)^2)
        preds<-c(preds,new.pred)
        obs<-c(obs,new.obs)
      }
    }
  }
  if(return.vec==T) {data.frame(pred=preds,obs=obs)} else {sum(val)}
}


# optimize decay plume
opt1<-optim(par=c(1.858651e+00,4.903428e-01,9.161020e-07),fn=param.search.optim.decay.plume,control=list(trace=1))

# x<-c(4.283863e-01,7.454097e-01,2.814154e-07) #OPT1 output!!! value = 22847.04 for full period
# x<-c(1.858651e+00,4.903428e-01,9.161020e-07) #OPT1 output!!! value = 22815.71 for two days
# x<-c(2.624685e+00,5.127938e-01,1.784182e-06) #OPT1 output!!! value = 22832.76 for one day

test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
out<-mapply(decay.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=4224.733,k=opt$par[3],s=1,alphay=opt$par[2],c=opt$par[1]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

pred.mat<-param.search.optim.decay.plume(x,return.vec=T)
plot(pred.mat$obs,pred.mat$pred)
plot(pred.mat$obs,(pred.mat$obs-pred.mat$pred)^2)


# optimize tilted plume
opt2<-optim(par=c(.15,.15,.5),fn=param.search.optim.tilted.plume,control=list(trace=1))

test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
out<-mapply(tilted.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=266.4167,H=.2,s=1,alphay=opt2$par[1],alphaz=opt2$par[2],Ws=opt2$par[3]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)


###### visualize parameter search
library(parallel)

param.search.decay.plume.plot<-function(cval,alphayval,k)
{
  param.search.optim.decay.plume(c(cval,alphayval,k))
}

param.search.tilted.plume.plot<-function(alphayval,alphazval,Wsval)
{
  param.search.optim.tilted.plume(c(alphayval,alphazval,Wsval))
}

## experiment results:

test.mat<-expand.grid(cval=seq(1.6,2,.04),alphayval=seq(.4,.6,.02),k=9.161020e-07) 
out<-mcmapply(param.search.decay.plume.plot,  cval = test.mat[,1],alphayval=test.mat[,2],k=test.mat[,3],mc.cores = 6)
par(mfrow=c(2,3))
for(i in seq(0,1,.1))
{
  res.mat<-matrix(out[which(test.mat$k==i)],11,11)
  contour(seq(1.6,2,.04),seq(.4,.6,.02),res.mat,xlab="cval",ylab="alphayval",main=paste("k=",i),nlevels = 20)
}

library(parallel)
test.mat<-expand.grid(cval=seq(0,.3,.03),alphayval=seq(.0,.3,.03),k=seq(0,2,.2))
out<-mcmapply(param.search.tilted.plume.plot,  alphayval = test.mat[,1],alphazval=test.mat[,2],Wsval=test.mat[,3],mc.cores = 4)
par(mfrow=c(2,3))
for(i in seq(0,2,.2))
{
  res.mat<-matrix(out[which(test.mat$k==i)],11,11)
  contour(seq(0,.05,.005),seq(.0,.3,.03),res.mat,xlab="cval",ylab="alphayval",main=paste("k=",i))
}
