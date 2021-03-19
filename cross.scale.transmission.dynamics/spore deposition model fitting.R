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

param.search.optim.decay.plume<-function(x,return.vec=F)
{
  cval=x[1]
  alphayval=x[2]
  k=x[3]
  
  val<-c()
  preds<-c()
  obs<-c()
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
      wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC")+60*60*24*1)),] ### fit to one day post spore trap deploy
      #wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC")+60*60*24*2)),] ### fit to two days post spore trap deploy
      #wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),] ### fit to full spore trap deploy period
      
      for(j in 1:dim(sub.2.spore.deposition)[1])
      {
        xtarget<-0
        ytarget<-0
        if(sub.2.spore.deposition[j,"Direction"]=="U") {ytarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="R") {xtarget<-as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="D") {ytarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        if(sub.2.spore.deposition[j,"Direction"]=="L") {xtarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance..cm."])/100}
        
        new.pred<-predict.kernel.decay.plume(q=q,k=k,alphay=alphayval,c=cval,xtarget=xtarget,ytarget=ytarget,wind.data=wind.data)
        new.obs<-sub.2.spore.deposition[j,"Spores"]/sub.2.spore.deposition[j,"X..squares.counted"]
        val<-c(val,(new.obs-new.pred)^2)
        preds<-c(preds,new.pred)
        obs<-c(obs,new.obs)
      }
    }
  }
  if(return.vec==T) {data.frame(pred=preds,obs=obs)} else {sum(val)}
}


# optimize decay plume

## optimize
opt<-optim(par=c(1.226695e-01,1.066038e-01,1.205203e-05),fn=param.search.optim.decay.plume,control=list(trace=1))

## results for test-run fitting model to just tag %in% c(86,88)
### x<-c(6.845809e-03,7.640827e-01,2.265886e-06) #OPT1 output value = 3678.804 for full period <-pancake like distribution, looks unrealistic
### x<-c(6.929344e-02,7.717773e-02,5.565447e-06) #OPT1 output value = 3562.495 for two days
### x<-c(7.332391e-02,7.595204e-02,1.146244e-05) #OPT1 output value = 3512.003 for one day

## results for model fitting to full data set
### sum squared obs = 15833.73
### total sum of squares = 15563.23 (total sum of squares = sum((spores/squares - mean(spores/squares))^2) )
### x<-c(4.718717e-01,2.258340e-02,6.428169e-06) #OPT1 output value = 14652.07 for one day
### x<-c(5.514282e+00,1.886189e-02,3.285846e-06) #OPT1 output value = 14675.07 for two days ALMOST NO DECAY IN X DIRECTION but similar output value for all values of cval (opt$par[1])
### x<-c(1.470363e-02,1.355715e-02,1.128306e-06) #OPT1 output value = 14871.98 for full period

## visualize kernel
test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
out<-mapply(decay.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=4224.733,k=opt$par[3],s=1,alphay=opt$par[2],c=opt$par[1]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

## visualize fit
pred.mat<-param.search.optim.decay.plume(opt$par,return.vec=T)
plot(pred.mat$obs,pred.mat$pred)
plot(pred.mat$obs,pred.mat$obs-pred.mat$pred)


# visualize parameter search
library(parallel)

param.search.decay.plume.plot<-function(cval,alphayval,k)
{
  param.search.optim.decay.plume(c(cval,alphayval,k))
}


cvalset<-seq(.4,1,.06)
alphayvalset<-seq(0,.2,.02)
kset=6.428169e-06
test.mat<-expand.grid(cval=cvalset,alphayval=alphayvalset,k=kset) 
out<-mcmapply(param.search.decay.plume.plot,  cval = test.mat[,1],alphayval=test.mat[,2],k=test.mat[,3],mc.cores = 6)
res.mat<-matrix(out[which(test.mat$k==kset)],length(cvalset),length(alphayvalset))
contour(cvalset,alphayvalset,res.mat,xlab="cval",ylab="alphayval",main=paste("k=",i),nlevels = 1000)



