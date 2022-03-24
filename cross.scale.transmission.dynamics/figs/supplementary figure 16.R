library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc data set building.R")

weath.data.scenario.vec<-c("rcp45","rcp85")
weath.data.vec<-c("2020","2030","2040","2045","2050","2060","2070")
weather.lty<-c(2,1)

par(mar=c(4,6,2,2),mfrow=c(1,2))
colors<-c("lightblue4","lightblue2")

for(site in c("BT","GM"))
{
  pred.epi.obs<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.2020.rcp45.",site,".step.size.7.RDS"))
  ylab=paste0("predicted prevalence on July ",c("29","23")[which(c("BT","GM")==site)])
  plot(0,0,xlim=c(2020,2070),ylim=c(0,c(.9,.45)[which(c("BT","GM")==site)]),xlab="",ylab=ylab,type="n",cex.axis=2,cex.lab=2)
  grid()
  mtext(c("A","B")[which(c("BT","GM")==site)],adj=1,cex=2)
  mtext(site,cex=2,font=2)
  for(i in 1:length(weath.data.scenario.vec))
  {
    scenario<-weath.data.scenario.vec[i]
    year.yvals<-c()
    year.yvals.1<-c()
    year.yvals.9<-c()
    
    for(j in 1:length(weath.data.vec))
    {
      year<-weath.data.vec[j]
      
      sim.data<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.",year,".",scenario,".",site,".step.size.7.RDS"))
      
      all.yvals<-c()
      
      for(k in 1:length(sim.data))
      {
        xvals<-c()
        yvals<-c()
        pred.epi<-sim.data[[k]]
        for(l in 1:length(unique(sim.data[[1]]$date)))
        {
          date<-unique(pred.epi$date)[l]
          xvals<-c(xvals,date)
          yvals<-c(yvals,sum(pred.epi[which(pred.epi$date==date),"status"])/length(which(pred.epi$date==pred.epi$date[1])))
        }
        all.yvals<-c(all.yvals,yvals[length(yvals)])
      }
      year.yvals<-c(year.yvals,mean(all.yvals))
      year.yvals.1<-c(year.yvals.1,quantile(all.yvals,.1))
      year.yvals.9<-c(year.yvals.9,quantile(all.yvals,.9))
    }
    polygon(c(2020,2030,2040,2045,2050,2060,2070,2070,2060,2050,2045,2040,2030,2020),c(year.yvals.9,rev(year.yvals.1)),col=colors[i],density=50,angle=c(45,135)[i])
    points(as.numeric(weath.data.vec),year.yvals,type="l",lty=weather.lty[i],lwd=4,col=colors[i])
  }
}
legend("topright",legend=c("RCP4.5","RCP8.5"),cex=2,lwd=4,lty=c(2,1),col=colors,bty="n",seg.len = 3)

#export at 1396x687