# load spore deposition functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")

# load fitted tilted gaussian plume model parameters
opt<-list(par=c(5.739690e-07,4.451030e-02, 7.777373e-02))

## visualize fit
pred.mat<-param.search.optim.tilted.plume(opt$par,return.out=T)

SS.tot<-sum((pred.mat$obs-mean(pred.mat$obs))^2)
SS.resid<-sum(pred.mat$obs-pred.mat$pred)^2
R.squared<-1-SS.resid/SS.tot

par(mfrow=c(2,2),mar=c(5,5,5,5))
plot(pred.mat$obs,pred.mat$pred,xlab="observed spore deposition",ylab="predicted spore deposition",cex=2,cex.lab=2)
grid()
text(65,27,expression(R^2*' = 0.409'),pos=4,offset=0,cex=2)
segments(65,23,77,23,lwd=1,lty=2)
text(78,23,"y = x",pos=4,offset=0,cex=2)
abline(0,1,lty=2,lwd=1)
mtext("A",adj=1,cex=1.5,font=2,line=1)

plot(log10(pred.mat$obs),log10(pred.mat$pred),xlab=expression(log[10]*" observed spore deposition"),ylab=expression(log[10]*" predicted spore deposition"),cex=2,cex.lab=2)
grid()
abline(0,1,lty=2,lwd=1)
mtext("B",adj=1,cex=1.5,font=2,line=1)

plot(pred.mat$obs,pred.mat$pred-pred.mat$obs,xlab="observed spore deposition",ylab="residual",cex=2,cex.lab=2)
grid()
abline(h=0,lty=2,lwd=1)
mtext("C",adj=1,cex=1.5,font=2,line=1)

plot(pred.mat$dist,pred.mat$pred-pred.mat$obs,xlab="distance",ylab="residual",cex=2,cex.lab=2)
grid()
abline(h=0,lty=2,lwd=1)
mtext("D",adj=1,cex=1.5,font=2,line=1)