# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions.R")
## load enviro data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
## load spore dep data
spore.deposition<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts.csv")
spore.deposition[which(spore.deposition$Distance.cm==0),"Distance.cm"]<-1 #set distance for "0" spore traps to 1cm
## load demog data
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
demog<-demog[which(demog$year==2020),] #subset to 2020

# optimize decay plume

## optimize

opt<-optim(par=c(4.772329e-01,2.261641e-02,6.285425e-06),fn=param.search.optim.decay.plume,control=list(trace=1))

## results for model fitting
### sum squared obs = 14764.22
### total sum of squares = 15563.23 (total sum of squares = sum((spores/squares - mean(spores/squares))^2) )
### x<-c(4.772329e-01,2.261641e-02,6.285425e-06) #OPT1 output value = 14637.1 for one day
### x<-c(5.543098e+00,1.907743e-02,3.231231e-06) #OPT1 output value = 14651.75 for two days ALMOST NO DECAY IN X DIRECTION but similar output value for all values of cval (opt$par[1])
### x<-c(1.177116e-02,1.185136e-02,1.263117e-06) #OPT1 output value = 14843.28 for full period

## visualize kernel
test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
out<-mapply(decay.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=4224.733,k=opt$par[3],s=1,alphay=opt$par[2],c=opt$par[1]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

## visualize fit
pred.mat<-param.search.optim.decay.plume(opt$par,return.out=T)
plot(pred.mat$obs,pred.mat$pred)
plot(pred.mat$obs,pred.mat$pred-pred.mat$obs)


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



