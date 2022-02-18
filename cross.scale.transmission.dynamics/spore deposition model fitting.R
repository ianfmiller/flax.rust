# load data
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions.R")
## load enviro data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
## load spore dep data
spore.deposition<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts.csv")
spore.deposition[which(spore.deposition$Distance.cm==0),"Distance.cm"]<-5 #set distance for "0" spore traps to 5cm
## load demog data
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
demog<-demog[which(demog$year==2020),] #subset to 2020

# optimize decay plume

## will return silightly different result due to presence of a likelihood ridge
opt<-optim(par=c(5.739690e-07 ,4.451030e-02,7.777373e-02),fn=param.search.optim.tilted.plume,control=list(trace=1))

## fiitted values
opt<-list(par=c(5.739690e-07,4.451030e-02, 7.777373e-02))

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


