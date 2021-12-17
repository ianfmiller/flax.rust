# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
## load enviro data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
## load spore dep data
spore.deposition<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts.csv")
spore.deposition[which(spore.deposition$Distance.cm==0),"Distance.cm"]<-5 #set distance for "0" spore traps to 5cm
## load demog data
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
demog<-demog[which(demog$year==2020),] #subset to 2020

# optimize decay plume

## optimize

### optimizaiton
opt<-optim(par=c(5.514312e-07,7.658422e-01,.5),fn=param.search.optim.tilted.plume,control=list(trace=1))


opt<-list(par=c(5.514312e-07,7.658422e-01))

## visualize kernel
test.mat<-data.frame(x=rep(seq(-2,2,.01),each=401),y=rep(seq(-2,2,.01),times=401))
out<-mapply(tilted.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(I=4224.733,H=.18,s=10,k=opt$par[1],Ws=opt$par[2],A=opt$par[3]))
res.mat<-matrix(out,401,401,byrow = T)
contour(x=seq(-2,2,.01),y=seq(-2,2,.01),z=res.mat)
points(c(.05,.25,.5,1),c(0,0,0,0),col="red")
points(0,0,col="blue")
filled.contour(x=seq(-2,2,.01),y=seq(-2,2,.01),z=log10(res.mat),zlim=c(-10,3),xlim=c(-2,2),ylim=c(-2,2),key.title = mtext(expression(log[10]*'spores per mm'^2),cex=1.5),cex.lab=1.5,xlab="X (meters)",ylab="Y (meters)",)

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

kset=7.165041e-08
alphazvalset<-seq(1e-6,3e-6,length.out = 11)
Wsset= seq(.7,.8,length.out = 11)
test.mat<-expand.grid(kval=kset,alphazval=alphazvalset,Wsval=Wsset) 
out<-mcmapply(param.search.tilted.plume,  kval = test.mat[,1],alphazval=test.mat[,2],Wsval=test.mat[,3],mc.cores = 6)
res.mat<-matrix(out[which(test.mat$k==kset)],length(alphazvalset),length(Wsset))
contour(alphazvalset,Wsset,res.mat,xlab="alphazval",ylab="Wsval",nlevels = 100)
points(opt$par[2],opt$par[3])


