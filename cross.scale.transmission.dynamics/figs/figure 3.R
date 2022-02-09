library(mgcv)
library(viridis)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/infection intensity data prep.R")
delta.infection.intensity<-subset(delta.infection.intensity,time<=8)
infection.intensity.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/infection.intensity.model.RDS")
dummy.infection.intensity.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/dummy.infection.intensity.model.RDS")
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

vis.gam(dummy.infection.intensity.model,view=c("max.height","log.10.infection.intensity"),plot.type = "contour",type="response",contour.col = "black",color="cm",zlim=c(-4000,4000),xlim=c(5,60),nCol = 100,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",n.grid=100)
mtext("plant height (cm)",1,line = 2.25,cex=1)
mtext(expression(log[10]*' infection intensity'),2,line=2.25,cex=1)
points(delta.infection.intensity$max.height,log10(delta.infection.intensity$infection.intensity),pch=".")

par(mar=c(4,0,1,2.5))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-4000-40.40404,4000+40.40404),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-4000,4000,length.out=101)[i]
  rect(0,ii-40.40404,1,ii+40.40404,col=cm.colors(101)[i],border = NA)
}
rect(0,-4000-40.40404,1,4000+40.40404)
mtext("B",cex=2)
axis(4,cex.axis=1,tck=-.5,padj=-1)

par(mar=c(4,4,1,1))
plot(infection.intensity.model,select = 2,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("mean temperature (°C)",1,line = 2.25,cex=1)
mtext("s(mean temperature)",2,line=2.25,cex=1)
grid()
mtext("C",adj=1,cex=1.5)
plot(infection.intensity.model,select = 3,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("max. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(max. temperature)",2,line=2.25,cex=1)
grid()
mtext("D",adj=1,cex=1.5)
plot(infection.intensity.model,select = 4,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("min. temperature (°C)",1,line = 2.25,cex=1)
mtext("s(min. temperature)",2,line=2.25,cex=1)
grid()
mtext("E",adj=1,cex=1.5)
plot(infection.intensity.model,select = 5,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext(expression('mean abs. humidity ('*g/m^3*')'),1,line = 2.25,cex=1)
mtext("s(mean abs. humidity)",2,line=2.25,cex=1)
grid()
mtext("F",adj=1,cex=1.5)
plot(infection.intensity.model,select = 6,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="")
mtext("total rainfall (mm)",1,line = 2.25,cex=1)
mtext("s(total rainfall)",2,line=2.25,cex=1)
grid()
mtext("G",adj=1,cex=1.5)
plot(infection.intensity.model,select = 7,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(plant ID)",2,line=2.25,cex=1)
grid()
mtext("H",adj=1,cex=1.5)
plot(infection.intensity.model,select = 8,shade=T,main="",cex.lab=1.25,cex.axis=1,xlab="",ylab="",col=site.cols[c(2,1,3,4)],cex=1.5/par()$cex,pch=16)
mtext("Gaussian quantiles",1,line = 2.25,cex=1)
mtext("s(site)",2,line=2.25,cex=1)
grid()
mtext("I",adj=1,cex=1.5)
legend("bottomright",legend=c("CC","BT","GM","HM"),pch=16,col=site.cols,cex=1.25,bty="n",pt.cex = 1.5/par()$cex)

mtext("generalized additive model",outer=T,adj=19.5/23,cex=2,line=-2.5)
mtext("change in infection intensity per day",outer=T,adj=20/23,cex=1.75,line=-5.5)

site.indicies<-c(2,1,3,4)[as.numeric(delta.infection.intensity$site)]
par(mar=c(5,6,5,2))
plot(delta.infection.intensity$infection.intensity,delta.infection.intensity$infection.intensity.next,xlab = "observed infection intensity",ylab="next observed infection intensity",cex.lab=2,cex.axis=2,col=site.cols[site.indicies],pch=16,cex=delta.infection.intensity$time/2,panel.first = abline(0,1,lty=2))
grid()
mtext("A",side=3,adj=1,cex=1.5)
legend("bottomright",legend=c("CC","BT","GM","HM","2 days","4 days","6 days"),col=c(site.cols,"grey","grey","grey"),pt.cex=c(3,3,3,3,2/2,4/2,6/2),pch=16,cex=1.75,bty="n")
mtext("data",cex=2,line=1)

library("MASS")
library("viridis")
start.height<-25
start.inf.intens<-1

site<-"GM"
start.date<-c(as.POSIXct("2020-06-23 00:00:00",tz="UTC"),as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-06-24 00:00:00",tz="UTC"),as.POSIXct("2020-06-26 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
end.date<-c(as.POSIXct("2020-07-27 00:00:00",tz="UTC"),as.POSIXct("2020-07-29 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),as.POSIXct("2020-07-10 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
sim.dates<-seq.POSIXt(start.date,end.date,"7 day")
weath.data.vec<-c("observed","2020","2020","2045","2045","2070","2070")
weath.data.scenario.vec<-c(NA,"rcp45","rcp85","rcp45","rcp85","rcp45","rcp85")

plot(0,0,xlim=c(start.date,end.date),ylim=c(0,1500),type="n",xlab="date",ylab="infection intensity",cex.lab=2,axes=F)
grid()
mtext("J",side=3,adj=1,cex=2)
axis.POSIXct(1,sim.dates,cex.axis=2)
axis(2,cex.axis=2)
box()
mtext("projections",cex=2,line=1)

for(i in 2:7)
{
  weath.data<-weath.data.vec[i]
  weath.data.scenario<-weath.data.scenario.vec[i]
  
  xcords<-rep(NA,length(sim.dates)) #time values
  inf.intens.cords<-rep(NA,length(sim.dates)) #infection intensity values
  height.cords<-rep(NA,length(sim.dates)) #height values
  
  set.seed(73452749)
  
  for(j in 1:100)
  {
    inf.intens<-start.inf.intens
    height<-start.height
    xcords.new<-c(sim.dates[1])
    inf.intens.cords.new<-c(inf.intens)
    height.cords.new<-c(height)
    
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
      
      beta.inf.intens <- coef(infection.intensity.model) ## posterior mean of coefs
      Vb.inf.intens  <- vcov(infection.intensity.model) ## posterior  cov of coefs
      n <-2
      mrand.inf.intens <- mvrnorm(n, beta.inf.intens, Vb.inf.intens) ## simulate n rep coef vectors from posterior
      pred.data.inf.intens<-data.frame("max.height"=height,"infection.intensity"=inf.intens,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.mean.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain,tag="NA",site=site)
      Xp.inf.intens <- predict(infection.intensity.model, newdata = pred.data.inf.intens, exclude=c("s(tag)"),type="lpmatrix")
      ilink.inf.intens <- family(infection.intensity.model)$linkinv
      preds.inf.intens <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds.inf.intens[l]   <- ilink.inf.intens(Xp.inf.intens %*% mrand.inf.intens[l, ])[1]
      }
      inf.intens.change<-preds.inf.intens[1]
      inf.intens<-inf.intens+inf.intens.change*as.numeric(date1-date0)
      if(inf.intens<0.1) {inf.intens<-0.1}
      
      beta.growth <- coef(plant.growth.model) ## posterior mean of coefs
      Vb.growth  <- vcov(plant.growth.model) ## posterior  cov of coefs
      n <-2
      mrand.growth <- mvrnorm(n, beta.growth, Vb.growth) ## simulate n rep coef vectors from posterior
      pred.data.growth<-data.frame("height"=height,"inf.intens"=inf.intens,"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.mean.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain,tag="NA",site=site)
      Xp.growth <- predict(plant.growth.model, newdata = pred.data.growth, exclude=c("s(tag)"),type="lpmatrix")
      ilink.growth <- family(plant.growth.model)$linkinv
      preds.growth <- rep(NA,n)
      for (l in seq_len(n)) { 
        preds.growth[l]   <- ilink.growth(Xp.growth %*% mrand.growth[l, ])[1]
      }
      height.change<-preds.growth[1]
      height<-height+height.change*as.numeric(date1-date0)
      if(height<5) {height<-5}
      
      #reps<-reps+pred.window
      xcords.new<-c(xcords.new,date1)
      inf.intens.cords.new<-c(inf.intens.cords.new,inf.intens)
      height.cords.new<-c(height.cords.new,height)
    }
    xcords<-rbind(xcords,xcords.new)
    inf.intens.cords<-rbind(inf.intens.cords,inf.intens.cords.new)  
    height.cords<-rbind(height.cords,height.cords.new)
  }
  points(xcords[2,],colMeans(inf.intens.cords[-1,]),col=weather.colors[i],type="l",lwd=5,lty=c(1,3,1,3,1,3,1,3,1)[i])
}
legend("topleft",
       legend=c("2020 RCP 4.5", "2020 RCP 8.5", "2045 RCP 4.5","2045 RCP 8.5","2070 RCP 4.5","2070 RCP 8.5"),
       cex=1.75,
       lwd=4,
       seg.len = 2,
       lty=c(1,3,1,3,1,3,1,3),
       col=weather.colors[-1],
       bty="n"
)

# export at 1291x812
