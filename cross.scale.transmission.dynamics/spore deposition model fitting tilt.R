# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
## load enviro data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
## load spore dep data
spore.deposition<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts tmp.csv")
spore.deposition[which(spore.deposition$Distance.cm==0),"Distance.cm"]<-1 #set distance for "0" spore traps to 1cm
## load demog data
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
demog<-demog[which(demog$year==2020),] #subset to 2020

# optimize decay plume

## optimize

opt<-optim(par=c(1.470378e-06, 8.371566e-02, 1.065001e+01, 1.111622e+02),fn=param.search.optim.tilted.plume,control=list(trace=1))

## results for model fitting
### sum squared obs = 14991.02
### total sum of squares = 14764.22 (total sum of squares = sum((spores/squares - mean(spores/squares))^2) )
### x<-c(1.470378e-06, 8.371566e-02, 1.065001e+01, 1.111622e+02) # output value = 14232.8 for one day
### x<-c(5.288992e-01 , -3.461082e-02 , 2.512734e-06) # output value = 14319.83 for two days ALMOST NO DECAY IN X DIRECTION but similar output value for all values of cval (opt$par[1])
### x<-c(1.177116e-02,1.185136e-02,1.263117e-06) # output value = 14843.28 for full period

## visualize kernel
test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
#opt<-list(par=c(5e-5,.1,2,100))
out<-mapply(tilted.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=4224.733,H=.5,s=1,k=opt$par[1],alphay=opt$par[2],alphaz=opt$par[3],Ws=opt$par[4]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

## visualize fit
pred.mat<-param.search.optim.tilted.plume(opt$par,return.out=T)
plot(pred.mat$obs,pred.mat$pred)
plot(pred.mat$obs,pred.mat$pred-pred.mat$obs)
plot(pred.mat$dist,pred.mat$pred-pred.mat$obs)

# visualize parameter search
library(parallel)

param.search.tilted.plume<-function(kval,alphayval,alphazval,Wsval)
{
  param.search.optim.tilted.plume(c(kval,alphayval,alphazval,Wsval))
}

kset=1.470378e-06
alphayvalset<-seq(0,.1,.01)
alphazvalset<-seq(5,15,1)
Wsset= 1.111622e+02
test.mat<-expand.grid(kval=kset,alphayval=alphayvalset,alphazval=alphazvalset,Wsval=Wsset) 
out<-mcmapply(param.search.tilted.plume,  kval = test.mat[,1],alphayval=test.mat[,2],alphazval=test.mat[,3],Wsval=test.mat[,4],mc.cores = 6)
res.mat<-matrix(out[which(test.mat$k==kset)],length(alphayvalset),length(alphazvalset))
contour(alphayvalset,alphazvalset,res.mat,xlab="alphayval",ylab="Wsval",nlevels = 10000)



