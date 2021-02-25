library(parallel)

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")

wind.data<-subset(all.weath,site=="BT")
wind.data<-subset(wind.data,date>as.POSIXct("2020-07-09 12:00:00"))
wind.data<-subset(wind.data,date<=as.POSIXct("2020-07-16 12:00:00"))


predict.kernel<-function(q,k,alphay,c,xtarget,ytarget)
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

param.search.plot<-function(cval,alphayval)
{
  k=.001
  q=266.4167
  (0-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=.25,ytarget=0))^2+
  (0-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=-.25))^2+
  (232-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0.05,ytarget=0)^2)+
  (336-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=-.05))^2+
  (81-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=.05))^2+
  (134-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=-.05,ytarget=0))^2+
  (7-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=.25))^2+
  (18-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=-.25,ytarget=0))^2
}

test.mat<-data.frame(cval=rep(seq(0,.05,.005),each=11),alphayval=rep(seq(.05,.3,.025),times=11))
out<-mcmapply(param.search.plot,  cval = test.mat[,1],alphayval=test.mat[,2],mc.cores = 4)
contour(seq(.05,.3,.025),seq(0,.05,.005),res.mat,xlab="alphayval",ylab="cval")


param.search.optim<-function(x)
{
  alphayval=x[1]
  cval=x[2]
  k=.001
  q=266.4167
  (0-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=.25,ytarget=0))^2+
    (0-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=-.25))^2+
    (232-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0.05,ytarget=0)^2)+
    (336-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=-.05))^2+
    (81-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=.05))^2+
    (134-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=-.05,ytarget=0))^2+
    (7-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=0,ytarget=.25))^2+
    (18-predict.kernel(q=q,k=k,alphay=alphayval,c=cval,xtarget=-.25,ytarget=0))^2
}

library(optimParallel)
cl<-makeCluster(detectCores()-1,type="FORK")
setDefaultCluster(cl=cl)
opt<-optimParallel(par=c(.175,.035),fn=param.search.optim)
stopCluster(cl)
points(opt$par[1],opt$par[2])
