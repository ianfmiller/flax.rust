library(mgcv)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant growth data prep.R")

delta.plant.heights<-subset(delta.plant.heights,time<=8)

# visualize data

layout(matrix(c(1,2,4,4,3,3,3,3),2,4,byrow = T))
par(mar=c(5,5,3,0.5))

## histograms
hist(delta.plant.heights$height,main="",breaks=100,xlab="plant height",cex.lab=2,cex.axis=2,cex.main=2)
mtext("A",side=3,adj=1,line=-3,cex=2)
hist((delta.plant.heights$height.next-delta.plant.heights$height)/delta.plant.heights$time,main="",breaks=100,xlab="change in plant height per day",cex.lab=2,cex.axis=2,cex.main=2)
mtext("B",side=3,adj=1,line=-3,cex=2)

## plot trajectories
par(mar=c(2,5,0,0.5))
plot(c(min(plant.heights$date),max(plant.heights$date)),c(0,max(plant.heights$max.height)),type="n",xlab="",ylab="plant height",cex.lab=2,cex.axis=2)
mtext("D",side=3,adj=1,line=-3,cex=2)
i<-0

plot.cols<-sample(rainbow(180))

for (i in 1:length(unique(plant.heights$tag)))
{
  tag<-unique(plant.heights$tag)[i]
  sub.heights<-plant.heights[which(plant.heights$tag==tag),]
  sub.heights<-sub.heights[order(sub.heights$date),]
  points(sub.heights$date,sub.heights$max.height,col=plot.cols[i],type="l",lwd=.5)

}

## plot change
par(mar=c(5,5,3,0.5))
plot(delta.plant.heights$height,delta.plant.heights$height.next,col="black",xlab = "observed height",ylab="next observed height",cex.lab=2,cex.axis=2)
mtext("C",side=3,adj=1,line=-3.25,cex=2)
abline(0,1,lty=2)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{
  set.seed(23094867)
  mod<-gam((height.next-height)/time~
             te(height,inf.intens)+
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(mean.daily.rain)+
             s(tag,bs="re")+
             s(site,bs="re"),
           select = T,
           method="REML",
           data=delta.plant.heights,
           control = list(nthreads=4))

  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

## model checking
par(mfrow=c(2,2))
gam.check(plant.growth.model) #Indicates that basis dimension is sufficient
concurvity(plant.growth.model,full=F) #No obvious issues

## visualize model

layout(matrix(c(1,1,1,1,2,3,3,3,3,4,4,4,4,10,11,5,5,5,5,6,6,6,6,7,7,7,7,12,13,13,13,8,8,8,8,9,9,9,9,14,14,14),3,14,byrow = T))

par(mar=c(4,7,3,1.5))
options(warn=-1) ## suppress warnings due to passing levels to vis.gam
vis.gam(plant.growth.model,view=c("height","inf.intens"),plot.type = "contour",type="response",labcex=.75,contour.col = "black",color="cm",zlim=c(-.75,.75),xlim=c(5,75),nCol = 100,main="",cex.lab=1.5,cex.axis=1.5,xlab="plant height (cm)",ylab="infection intensity")
points(delta.plant.heights$height,delta.plant.heights$inf.intens,pch=".")
par(mar=c(4,1,3,4))
plot(0,0,type="n",xlim=c(0,1),ylim=c(-.75-0.0075,.75+0.0075),axes=F,xlab="",ylab="")
for(i in 1:101)
{
  ii<-seq(-0.75,0.75,length.out=101)[i]
  rect(0,ii-0.0075,1,ii+0.0075,col=cm.colors(101)[i],border = NA)
}
rect(0,-0.75-0.0075,1,0.75+0.0075)
mtext("A",cex=1.25,font=2)
#mtext("te(plant height, infection intensity)",side=4,line=.5)
axis(2,cex.axis=1.5)

par(mar=c(4,4.5,3,4))
plot(plant.growth.model,select = 2,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="mean temperature (°C)",ylab="s(mean temperature)")
grid()
mtext("B",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 3,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="max. temperature (°C)",ylab="s(max. temperature)")
grid()
mtext("C",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 4,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="min. temperature (°C)",ylab="s(min. temperature)")
grid()
mtext("D",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 5,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab=expression('mean abs. humidity ('*g/m^3*')'),ylab="s(mean abs. humidity)")
grid()
mtext("E",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 6,shade=T,main="",cex.lab=1.5,cex.axis=1.5,xlab="total rainfall (mm)",ylab="s(total rainfall)")
grid()
mtext("I",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 7,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(tag)")
grid()
mtext("L",adj=1,cex=1.25,font=2)
plot(plant.growth.model,select = 8,shade=T,main="",cex.lab=1.5,cex.axis=1.5,ylab="s(site)")
grid()
mtext("M",adj=1,cex=1.25,font=2)


# predict climate change effects
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

set.seed(874627)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/within host climate prediction functions.R")
library("MASS")
library("viridis")
par(mfrow=c(2,2))
start.height<-10

for(site in c("CC","BT","GM","GM"))
{
  start.date<-c(as.POSIXct("2020-06-23 00:00:00",tz="UTC"),as.POSIXct("2020-06-20 00:00:00",tz="UTC"),as.POSIXct("2020-06-24 00:00:00",tz="UTC"),as.POSIXct("2020-06-26 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
  end.date<-c(as.POSIXct("2020-07-27 00:00:00",tz="UTC"),as.POSIXct("2020-07-29 00:00:00",tz="UTC"),as.POSIXct("2020-07-28 00:00:00",tz="UTC"),as.POSIXct("2020-07-10 00:00:00",tz="UTC"))[which(c("CC","BT","GM","HM")==site)]
  sim.dates<-seq.POSIXt(start.date,end.date,"3 day")
  weath.data.vec<-c("observed","2020","2020","2045","2045","2070","2070")
  weath.data.scenario.vec<-c(NA,"rcp45","rcp85","rcp45","rcp85","rcp45","rcp85")
  #colors<-c("black",magma(100)[c(20,30,45)],magma(100)[c(70,80,90)])
  colors<-c("black","darkorange","darkorange3","red1","darkred","purple1","purple4")
  plot(1:7,col=colors,pch=15,cex=5)
  plot(0,0,xlim=c(start.date,end.date),ylim=c(8,30),type="n",xlab="date",ylab="plant height",cex.axis=1,cex.lab=1.2,axes=F,main=site)
  axis.POSIXct(1,sim.dates)
  axis(2)
  box()
  
  for(i in 1:7)
  {
    weath.data<-weath.data.vec[i]
    weath.data.scenario<-weath.data.scenario.vec[i]
    
    xcords<-rep(NA,length(sim.dates)) #time values
    ycords<-rep(NA,length(sim.dates)) #height values
    
    for(j in 1:10)
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
        Xp <- predict(plant.growth.model, newdata = pred.data, exclude=c("s(tag)"),type="lpmatrix")
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
    
    for(m in 2:dim(xcords)[1])
    {
      points(xcords[m,],ycords[m,],col=t_col(colors[i],50),type="l",lwd=.5) 
    } 
    points(xcords[2,],colMeans(ycords[-1,]),col=colors[i],type="l",lwd=2)
  }
}


