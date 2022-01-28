library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc data set building.R")

weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2)])
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}
weather.colors.trans<-unlist(lapply(weather.colors, t_col))

weather.lty=c(1,2,1,2,1,2,1)
weath.scenarios<-c("obs","2020.rcp45","2020.rcp85","2045.rcp45","2045.rcp85","2070.rcp45","2070.rcp85")

par(mfrow=c(1,2),mar=c(5,5,5,5))

for(site in c("BT","GM"))
{
  pred.epi.obs<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.obs.",site,".step.size.7.RDS"))
  plot(unique(pred.epi.obs[[1]]$date),rep(0,times=length(unique(pred.epi.obs[[1]]$date))),ylim=c(0,.6),xlab="date",ylab="prevalence",type="n",cex.axis=2,cex.lab=2,main=site,cex.main=2)
  grid()
  mtext(c("A","B")[which(c("BT","GM")==site)],adj=1,cex=2,font=2)
  
  sub.epi<-corrected.epi[which(corrected.epi$Site==site),]
  sub.locs<-corrected.locs[which(corrected.locs$Site==site),]
  
  xvals<-c()
  yvals<-c()
  for(i in 1:length(unique(sub.epi$Date.First.Observed.Diseased)))
  {
    date<-unique(sub.epi$Date.First.Observed.Diseased)[i]
    xvals<-c(xvals,date)
    yvals<-c(yvals,nrow(sub.epi[which(sub.epi$Date.First.Observed.Diseased<=date),])/nrow(sub.locs))
  }
  points(xvals,yvals,type="l",col="blue",lwd=4)
  
  for(scenario in weath.scenarios)
  {
    sim.data<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.",scenario,".",site,".step.size.7.RDS"))
    
    for(k in 1:length(sim.data))
    {
      xvals<-c()
      yvals<-c()
      pred.epi<-sim.data[[k]]
      for(i in 1:length(unique(sim.data[[1]]$date)))
      {
        date<-unique(pred.epi$date)[i]
        xvals<-c(xvals,date)
        yvals<-c(yvals,sum(pred.epi[which(pred.epi$date==date),"status"])/length(which(pred.epi$date==pred.epi$date[1])))
      }
      points(xvals,yvals,type="l",col=weather.colors.trans[which(weath.scenarios==scenario)])
    }
    
    xvals<-c()
    yvals<-c()
    for(j in 1:length(unique(sim.data[[1]]$date)))
    {
      sub.prevs<-c()
      date<-unique(sim.data[[1]]$date)[j]
      for(k in 1:length(sim.data))
      {
        sub.dat<-sim.data[[k]]
        sub.dat<-sub.dat[which(sub.dat$date==date),]
        sub.prevs<-c(sub.prevs,sum(sub.dat$status)/dim(sub.dat)[1])
      }
      yvals<-c(yvals,mean(sub.prevs))
      xvals<-c(xvals,date)
    }
    points(xvals,yvals,type="l",col=weather.colors[which(weath.scenarios==scenario)],lty=weather.lty[which(weath.scenarios==scenario)],lwd=4)
  }
}

legend("topright",
       legend=c("observed","predicted","predicted 2020 RCP 4.5", "2020 RCP 8.5 predicted", "2045 RCP 4.5 predicted","2024 RCP 8.5 predicted","2070 RCP 4.5 predicted","2070 RCP 8.5 predicted"),
       lwd=4,
       cex=1.25,
       seg.len = 4,
       lty=c(1,1,3,1,3,1,3,1,3),
       col=c("blue",weather.colors),
       bty="n"
)

