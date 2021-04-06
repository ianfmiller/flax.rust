# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
## load enviro data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
## load spore dep data
spore.deposition<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts.csv")
spore.deposition[which(spore.deposition$Distance.cm==0),"Distance.cm"]<-5 #set distance for "0" spore traps to 1cm
## load demog data
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
demog<-demog[which(demog$year==2020),] #subset to 2020

# optimize decay plume

## optimize

#opt<-optim(par=c(5.803369e-07 , 1.596314e-01 , 1.100707e+00),fn=param.search.optim.tilted.plume,control=list(trace=1))

## NEED TO REFIT ALL--FIX k or get rid of k--three params making liklihood ridge
## results for model fitting
### sum squared obs = 14991.02
### total sum of squares = 14764.22 (total sum of squares = sum((spores/squares - mean(spores/squares))^2) )
### x<-c(9.318341e-07,1.897060e-01,1.085956e+00) # output value = 12297 for one day
### x<-c(5.803369e-07 , 1.596314e-01 , 1.100707e+00) # output value = 12283.01 for two days
### x<-c(4.355849e-08,1.714472e-06,1.384859e+00) # output value = 14084.46 for full period

opt<-list(par=c(5.803369e-07 , 1.596314e-01 , 1.100707e+00))

## visualize kernel
test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
out<-mapply(tilted.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=4224.733,H=.18,s=3,k=opt$par[1],alphaz=opt$par[2],Ws=opt$par[3]))
res.mat<-matrix(out,201,201,byrow = T)
contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)
points(c(.05,.25,.5,1),c(0,0,0,0),col="red")
points(0,0,col="blue")
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

## visualize fit
pred.mat<-param.search.optim.tilted.plume(opt$par,return.out=T)
plot(pred.mat$obs,pred.mat$pred)
plot(pred.mat$obs,pred.mat$pred-pred.mat$obs)
plot(pred.mat$dist,pred.mat$pred-pred.mat$obs)

# visualize parameter search
library(parallel)

param.search.tilted.plume<-function(kval,alphazval,Wsval)
{
  param.search.optim.tilted.plume(c(kval,alphazval,Wsval))
}

kset=4.355849e-08
alphazvalset<-seq(1e-6,2e-6,length.out = 11)
Wsset= seq(1.35,1.4,length.out = 11)
test.mat<-expand.grid(kval=kset,alphazval=alphazvalset,Wsval=Wsset) 
out<-mcmapply(param.search.tilted.plume,  kval = test.mat[,1],alphazval=test.mat[,2],Wsval=test.mat[,3],mc.cores = 6)
res.mat<-matrix(out[which(test.mat$k==kset)],length(alphazvalset),length(Wsset))
contour(alphazvalset,Wsset,res.mat,xlab="alphazval",ylab="Wsval",nlevels = 100)



