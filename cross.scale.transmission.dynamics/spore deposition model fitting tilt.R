# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
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

opt<-optim(par=c(5e-5,.1,2,100),fn=param.search.optim.tilted.plume,control=list(trace=1,maxit=200))

## results for model fitting
### sum squared obs = 14991.02
### total sum of squares = 14764.22 (total sum of squares = sum((spores/squares - mean(spores/squares))^2) )
### x<-c(0.0834738679, 0.1025881527, 0.0000049515) #OPT1 output value = 14305.1 for one day
### x<-c(5.288992e-01 , -3.461082e-02 , 2.512734e-06) #OPT1 output value = 14319.83 for two days ALMOST NO DECAY IN X DIRECTION but similar output value for all values of cval (opt$par[1])
### x<-c(1.177116e-02,1.185136e-02,1.263117e-06) #OPT1 output value = 14843.28 for full period

## visualize kernel
test.mat<-data.frame(x=rep(seq(-1,1,.01),each=201),y=rep(seq(-1,1,.01),times=201))
#opt<-list(par=c(5e-5,.1,2,100))
out<-mapply(tilted.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=4224.733,H=.5,s=4,k=opt$par[1],alphay=opt$par[2],alphaz=opt$par[3],Ws=opt$par[4]))
res.mat<-matrix(out,201,201,byrow = T)
filled.contour(x=seq(-1,1,.01),y=seq(-1,1,.01),z=res.mat)

## visualize fit
pred.mat<-param.search.optim.tilted.plume(opt$par,return.out=T)
plot(pred.mat$obs,pred.mat$pred)
plot(pred.mat$obs,pred.mat$pred-pred.mat$obs)


# visualize parameter search
library(parallel)

param.search.tilted.plume<-function(kval,alphayval,alphazval,Wsval)
{
  param.search.optim.tilted.plume(c(kval,alphayval,alphazval,Wsval))
}

kset=2.253107e-05
alphayvalset<-.25
alphazvalset<-seq(0,20,2)
Wsset=seq(0,2000,200)
test.mat<-expand.grid(kval=kset,alphayval=alphayvalset,alphazval=alphazvalset,Wsval=Wsset) 
out<-mcmapply(param.search.tilted.plume,  kval = test.mat[,1],alphayval=test.mat[,2],alphazval=test.mat[,3],Wsval=test.mat[,4],mc.cores = 6)
res.mat<-matrix(out[which(test.mat$k==kset)],length(alphazvalset),length(Wsset))
filled.contour(alphazvalset,Wsset,res.mat,xlab="alphayval",ylab="Wsval",nlevels = 1000)



