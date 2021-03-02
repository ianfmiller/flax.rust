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


predict.kernel<-function(q,k,alphay,c,xtarget,ytarget,wind.data)
{
  tot.dep<-0
  for (i in 1:(dim(wind.data)[1]-1))
  {
    delta.t<-wind.data[i+1,"date"]-wind.data[i,"date"]
    cords<-get.plume.xy(2*pi*wind.data[i,"wind.direction"]/360,0,0,xtarget,ytarget,plot=F)
    tot.dep<-c(tot.dep,decay.plume(q=q,k=k,s=wind.data[i,"wind.speed"],x=cords[1],y=cords[2],alphay=alphay,c=c))
  }
  sum(tot.dep,na.rm = T)
}

param.search.optim<-function(x)
{
  cval=x[1]
  alphayval=x[2]
  k=x[3]
  
  val<-c()
  
  for(tag in unique(spore.deposition$Tag))
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
      wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),]
      
      for(j in 1:dim(sub.2.spore.deposition)[1])
      {
        xtarget<-0
        ytarget<-0
        if(sub.2.spore.deposition[j,"Direction"]=="U") {ytarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="R") {xtarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="D") {ytarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="L") {xtarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        
        pred<-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=xtarget,ytarget=ytarget,wind.data=wind.data)
        obs<-sub.2.spore.deposition[j,"Pustules"]/sub.2.spore.deposition[j,"X..squares.counted"]
        val<-c(val,(obs-pred)^2)
      }
    }
  }
  sum(val)
}


opt<-optim(par=c(.15,.03,0.00176),fn=param.search.optim,control=list(trace=1))

x<-c(1.563252e-01,2.931680e-02,4.718382e-07) #OPT output!!!

test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
out<-mapply(decay.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=266.4167,k=opt$par[3],s=5,alphay=opt$par[2],c=opt$par[1]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

###### visualize parameter search
library(parallel)

param.search.plot<-function(cval,alphayval,k)
{
  param.search.optim(c(cval,alphayval,k))
}

#test.mat<-data.frame(cval=rep(seq(0,.05,.01),each=11),alphayval=rep(seq(.05,.3,.025),times=11))
test.mat<-expand.grid(cval=seq(0,.05,.005),alphayval=seq(.0,.3,.03),k=seq(8e-4,.002,.00024))
out<-mcmapply(param.search.plot,  cval = test.mat[,1],alphayval=test.mat[,2],k=test.mat[,3],mc.cores = 4)
par(mfrow=c(2,3))
for(i in seq(8e-4,.002,.00024))
{
  res.mat<-matrix(out[which(test.mat$k==i)],11,11)
  contour(seq(0,.05,.005),seq(.0,.3,.03),res.mat,xlab="cval",ylab="alphayval",main=paste("k=",i))
}
