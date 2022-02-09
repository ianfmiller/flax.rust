library(mgcv)
library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")
delta.plant.heights<-subset(delta.plant.heights,time<=8)
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

site.cols<-viridis_pal(alpha=.5)(20)[c(20,15,6,1)]
weather.colors<-c("black",viridis_pal(option = "C")(5)[c(4,4,3,3,2,2,1,1)])

layout(matrix(c(rep(10,10),rep(16,13),
                rep(10,10),1,1,1,1,2,3,3,3,3,4,4,4,4,
                rep(10,10),1,1,1,1,2,3,3,3,3,4,4,4,4,
                rep(10,10),1,1,1,1,2,3,3,3,3,4,4,4,4,
                rep(10,10),12,5,5,5,5,6,6,6,6,7,7,7,7,
                rep(11,10),12,5,5,5,5,6,6,6,6,7,7,7,7,
                rep(11,10),12,5,5,5,5,6,6,6,6,7,7,7,7,
                rep(11,10),13,14,14,8,8,8,8,9,9,9,9,15,15,
                rep(11,10),13,14,14,8,8,8,8,9,9,9,9,15,15,
                rep(11,10),13,14,14,8,8,8,8,9,9,9,9,15,15),
              10,23,byrow=T))
par(mar=c(4,4,1,1))

vis.gam(plant.growth.model,view=c("height","inf.intens"),plot.type = "contour",type="response",contour.col = "black",color="cm",zlim=c(-.75,.75),xlim=c(5,75),nCol = 100,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("plant height (cm)",1,line = 2.25,cex=1)
mtext("infection intensity",2,line=2.25,cex=1)
points(delta.plant.heights$height,delta.plant.heights$inf.intens,pch=".")

par(mar=c(4,0,1,2.5))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-.75-0.0075,.75+0.0075),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-0.75,0.75,length.out=101)[i]
  rect(0,ii-0.0075,1,ii+0.0075,col=cm.colors(101)[i],border = NA)
}
rect(0,-0.75-0.0075,1,0.75+0.0075)
mtext("B",cex=1.5)
axis(4,cex.axis=1,tck=-.5,padj=-1)

par(mar=c(4,4,1,1))
plot(plant.growth.model,select = 2,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("mean temperature (°C)",1,line = 2.25,cex=1)
mtext("s(mean temperature)",2,line=2.25,cex=1)
grid()
mtext("C",adj=1,cex=1.5)
plot(plant.growth.model,select = 3,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("max. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(max. temperature)",2,line=2.25,cex=1)
grid()
mtext("D",adj=1,cex=1.5)
plot(plant.growth.model,select = 4,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("min. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(min. temperature)",2,line=2.25,cex=1)
grid()
mtext("E",adj=1,cex=1.5)
plot(plant.growth.model,select = 5,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext(expression('mean abs. humidity ('*g/m^3*')'),1,line = 2.25,cex=1)
mtext("s(mean abs. humidity)",2,line=2.25,cex=1)
grid()
mtext("F",adj=1,cex=1.5)
plot(plant.growth.model,select = 6,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("mean daily rainfall (mm)",1,line = 2.25,cex=1)
mtext("s(total rainfall)",2,line=2.25,cex=1)
grid()
mtext("G",adj=1,cex=1.5)
plot(plant.growth.model,select = 7,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(plant ID)",2,line=2.25,cex=1)
grid()
mtext("H",adj=1,cex=1.5)
plot(plant.growth.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1,xlab="",ylab="",col=site.cols[c(2,1,3,4)],cex=1.5/par()$cex,pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(site)",2,line=2.25,cex=1)
grid()
mtext("I",adj=1,cex=1.5)
legend("topleft",legend=c("CC","BT","GM","HM"),pch=16,col=site.cols,cex=1.25,bty="n",pt.cex = 1.5/par()$cex)

mtext("generalized additive model",outer=T,adj=19.5/23,cex=2,line=-2.5)
mtext("change in plant height (cm per day)",outer=T,adj=20/23,cex=1.75,line=-5.5)

site.indicies<-c(2,1,3,4)[as.numeric(delta.plant.heights$site)]
par(mar=c(5,6,5,2))
plot(delta.plant.heights$height,delta.plant.heights$height.next,xlab = "observed height (cm)",ylab="next observed height (cm)",cex.lab=2,cex.axis=2,col=site.cols[site.indicies],pch=16,cex=delta.plant.heights$time/2,panel.first = {abline(0,1,lty=2);grid()})
mtext("A",side=3,adj=1,cex=1.5)
legend("topleft",legend=c("CC","BT","GM","HM","2 days","4 days","6 days"),col=c(site.cols,"grey","grey","grey"),pt.cex=c(3,3,3,3,2/2,4/2,6/2),pch=16,cex=1.75,bty="n")
mtext("data",cex=2,line=1)

library("MASS")
library("viridis")
start.height<-10

site<-"GM"
start.date<-c(as.POSIXct("2020-06-23 00:00:00",tz="UTC"),as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-06-24 00:00:00",tz="UTC"),as.POSIXct("2020-06-26 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
end.date<-c(as.POSIXct("2020-07-27 00:00:00",tz="UTC"),as.POSIXct("2020-07-29 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),as.POSIXct("2020-07-10 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
sim.dates<-seq.POSIXt(start.date,end.date,"3 day")
weath.data.vec<-c("observed","2020","2020","2045","2045","2070","2070")
weath.data.scenario.vec<-c(NA,"rcp45","rcp85","rcp45","rcp85","rcp45","rcp85")

plot(0,0,xlim=c(start.date,end.date),ylim=c(10,30),type="n",xlab="date",ylab="plant height (cm)",cex.lab=2,axes=F)
grid()
mtext("J",side=3,adj=1,cex=1.5)
axis.POSIXct(1,sim.dates,cex.axis=2)
axis(2,cex.axis=2)
box()
mtext("projections",cex=2,line=1)

for(i in 2:7)
{
  weath.data<-weath.data.vec[i]
  weath.data.scenario<-weath.data.scenario.vec[i]
  
  xcords<-rep(NA,length(sim.dates)) #time values
  ycords<-rep(NA,length(sim.dates)) #height values
  
  set.seed(874627)
  
  for(j in 1:100)
  {
    height<-start.height
    xcords.new<-c(sim.dates[1])
    ycords.new<-c(height)
    
    for(k in 1:(length(sim.dates)-1))
    {
      date0<-sim.dates[k]
      date1<-sim.dates[k+1]
      
      if(!(weath.data=="observed"))
      {
        mean.temp.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/Tair.csv"),header = F)
        max.temp.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/tasmax.csv"),header = F)
        min.temp.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/tasmin.csv"),header = F)
        rh.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/relHumid.csv"),header = F)
        rainfall.data<-read.csv(paste0("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/",weath.data,"/rainfall.csv"),header = F)
        colnames(mean.temp.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
        colnames(max.temp.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
        colnames(min.temp.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
        colnames(rh.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
        colnames(rainfall.data)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
        
        mean.temp.data.sub.1<-mean.temp.data[which(format(as.Date(paste0(mean.temp.data$year,"-",mean.temp.data$month,"-",mean.temp.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),] #compare only month and day
        mean.temp.data.sub.2<-mean.temp.data.sub.1[which(format(as.Date(paste0(mean.temp.data.sub.1$year,"-",mean.temp.data.sub.1$month,"-",mean.temp.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
        
        max.temp.data.sub.1<-max.temp.data[which(format(as.Date(paste0(max.temp.data$year,"-",max.temp.data$month,"-",max.temp.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
        max.temp.data.sub.2<-max.temp.data.sub.1[which(format(as.Date(paste0(max.temp.data.sub.1$year,"-",max.temp.data.sub.1$month,"-",max.temp.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
        
        min.temp.data.sub.1<-min.temp.data[which(format(as.Date(paste0(min.temp.data$year,"-",min.temp.data$month,"-",min.temp.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
        min.temp.data.sub.2<-min.temp.data.sub.1[which(format(as.Date(paste0(min.temp.data.sub.1$year,"-",min.temp.data.sub.1$month,"-",min.temp.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
        
        rh.data.sub.1<-rh.data[which(format(as.Date(paste0(rh.data$year,"-",rh.data$month,"-",rh.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
        rh.data.sub.2<-rh.data.sub.1[which(format(as.Date(paste0(rh.data.sub.1$year,"-",rh.data.sub.1$month,"-",rh.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
        
        rainfall.data.sub.1<-rainfall.data[which(format(as.Date(paste0(rainfall.data$year,"-",rainfall.data$month,"-",rainfall.data$day)),"%m-%d")>=format(as.Date(date0),"%m-%d")),]
        rainfall.data.sub.2<-rainfall.data.sub.1[which(format(as.Date(paste0(rainfall.data.sub.1$year,"-",rainfall.data.sub.1$month,"-",rainfall.data.sub.1$day)),"%m-%d")<format(as.Date(date1),"%m-%d")),]
        
        new.mean.temp<-mean(mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
        new.max.temp<-max(max.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
        new.min.temp<-min(min.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
        abs.hums<-13.24732*exp((17.67*mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])/(mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)]+243.5))*rh.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)]/(273.15+mean.temp.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
        new.mean.abs.hum<-mean(abs.hums)
        new.mean.daily.rain<-mean(rainfall.data.sub.2[,paste0("cesm1.cam5.1.",weath.data.scenario)])
        
      }
      
      if(weath.data=="observed")
      {
        source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
        
        weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
        temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #### pull out temp data for site
        
        #### subst temp rh data to relevant window
        temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
        temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
        temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
        
        #### subset weather data to relevant window
        weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
        weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
        
        #### calculate environmental variable metrics
        new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
        new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
        new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
        abs.hum<-13.24732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+temp.rh.sub$temp.c)
        new.mean.abs.hum<-mean(abs.hum,na.rm=T)
        new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)*(12*24)
      }
      beta <- coef(plant.growth.model) ## posterior mean of coefs
      Vb   <- vcov(plant.growth.model) ## posterior  cov of coefs
      n <-2
      mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
      pred.data<-data.frame("height"=height,"inf.intens"=0,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.mean.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain,tag="NA",site=site)
      Xp <- predict(plant.growth.model, newdata = pred.data, exclude=c("s(plant ID)"),type="lpmatrix")
      ilink <- family(plant.growth.model)$linkinv
      preds <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds[l]   <- ilink(Xp %*% mrand[l, ])[1]
      }
      height.change<-preds[1]
      height<-height+height.change*as.numeric(date1-date0)
      if(height<5) {height<-5}
      
      #reps<-reps+pred.window
      xcords.new<-c(xcords.new,date1)
      ycords.new<-c(ycords.new,height)
    }
    xcords<-rbind(xcords,xcords.new)
    ycords<-rbind(ycords,ycords.new)  
  }
  points(xcords[2,],colMeans(ycords[-1,]),col=weather.colors[i],type="l",lwd=5,lty=c(1,3,1,3,1,3,1,3,1)[i])
}
legend("topleft",
       legend=c("2020 RCP 4.5", "2020 RCP 8.5", "2045 RCP 4.5","2045 RCP 8.5","2070 RCP 4.5","2070 RCP 8.5"),
       cex=1.75,
       lwd=4,
       seg.len = 2,
       lty=c(3,1,3,1,3,1,3),
       col=weather.colors[-1],
       bty="n"
       )

# export at 1291x812

