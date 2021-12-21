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
### will return silightly different result due to presence of a likelihood ridge
opt<-optim(par=c(5.739690e-07 ,4.451030e-02,7.777373e-02),fn=param.search.optim.tilted.plume,control=list(trace=1))


opt<-list(par=c(5.739690e-07,4.451030e-02, 7.777373e-02))

## visualize kernel
test.mat<-data.frame(x=rep(seq(-2,2,.01),each=401),y=rep(seq(-2,2,.01),times=401))
out<-mapply(tilted.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(I=4224.733,H=.25,s=1,k=opt$par[1],Ws=opt$par[2],A=opt$par[3]))
res.mat<-matrix(out,401,401,byrow = T)
contour(x=seq(-2,2,.01),y=seq(-2,2,.01),z=res.mat)
points(c(.05,.25,.5,1),c(0,0,0,0),col="red")
points(0,0,col="blue")
filled.contour(x=seq(-2,2,.01),y=seq(-2,2,.01),z=res.mat,xlim=c(-2,2),ylim=c(-2,2),key.title = mtext(expression(log[10]*'spores per mm'^2),cex=1.5),cex.lab=1.5,xlab="X (meters)",ylab="Y (meters)",)

## visualize fit
pred.mat<-param.search.optim.tilted.plume(opt$par,return.out=T)
plot(pred.mat$obs,pred.mat$pred,asp=1)
plot(log10(pred.mat$obs),log10(pred.mat$pred))
plot(pred.mat$obs,pred.mat$pred-pred.mat$obs)
plot(pred.mat$dist,pred.mat$pred-pred.mat$obs)

## visually check optimality
param.search.tilted.plume<-function(kval,Wsval,Aval)
{
  param.search.optim.tilted.plume(c(kval,Wsval,Aval))
}

Aset=7.777373e-02
kset<-seq(5e-7,6e-7,length.out = 11)
Wsset= seq(4e-2,5e-02,length.out = 11)
test.mat<-expand.grid(kval=kset,Aval=Aset,Wsval=Wsset) 
out<-mcmapply(param.search.tilted.plume,  kval = test.mat[,1],A=test.mat[,2],Wsval=test.mat[,3],mc.cores = 6)
res.mat<-matrix(out[which(test.mat$k==kset)],length(kset),length(Wsset))
contour(kset,Wsset,res.mat,xlab="k",ylab="Ws",nlevels = 100)
points(opt$par[1],opt$par[2])


