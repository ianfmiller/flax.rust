source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")


par(mfrow=c(2,3))

mean.temp.2020<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/Tair.csv",header = F)
mean.temp.2045<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/Tair.csv",header = F)
mean.temp.2070<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/Tair.csv",header = F)
colnames(mean.temp.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.temp.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.temp.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
dates<-as.Date(paste0(mean.temp.2020$year,"-",mean.temp.2020$month,"-",mean.temp.2020$day))
plot.dates<-seq.Date(as.Date("2020-06-12"),as.Date("2020-07-28"),1)
plot(dates,mean.temp.2020$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(5,19),xlim=c(min(plot.dates),max(plot.dates)),xlab="date",ylab="mean temp",lty=1,lwd=2)        
axis.Date(1,plot.dates)
points(dates,mean.temp.2045$cesm1.cam5.1.rcp45,type="l",col="orange",lty=2,lwd=2)
points(dates,mean.temp.2070$cesm1.cam5.1.rcp45,type="l",col="orange",lty=3,lwd=2)
points(dates,mean.temp.2020$cesm1.cam5.1.rcp85,type="l",col="red",lty=1,lwd=2)
points(dates,mean.temp.2045$cesm1.cam5.1.rcp85,type="l",col="red",lty=2,lwd=2)
points(dates,mean.temp.2070$cesm1.cam5.1.rcp85,type="l",col="red",lty=3,lwd=2)

for(site in c("BT"))
{
  mean.temps<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in plot.dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      mean.temps<-rbind(mean.temps,data.frame(date,mean(sub.temp.rh.2$temp.c,na.rm = T)))
    }
  }
  points(as.Date(mean.temps[,1],origin="1970-01-01"),mean.temps[,2],type="l",lwd=2)
}


max.temp.2020<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/tasmax.csv",header = F)
max.temp.2045<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/tasmax.csv",header = F)
max.temp.2070<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/tasmax.csv",header = F)
colnames(max.temp.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(max.temp.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(max.temp.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
plot(dates,max.temp.2020$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(10,35),xlim=c(min(plot.dates),max(plot.dates)),xlab="date",ylab="max temp",lty=1,lwd=2)        
axis.Date(1,plot.dates)
points(dates,max.temp.2045$cesm1.cam5.1.rcp45,type="l",col="orange",lty=2,lwd=2)
points(dates,max.temp.2070$cesm1.cam5.1.rcp45,type="l",col="orange",lty=3,lwd=2)
points(dates,max.temp.2020$cesm1.cam5.1.rcp85,type="l",col="red",lty=1,lwd=2)
points(dates,max.temp.2045$cesm1.cam5.1.rcp85,type="l",col="red",lty=2,lwd=2)
points(dates,max.temp.2070$cesm1.cam5.1.rcp85,type="l",col="red",lty=3,lwd=2)

for(site in c("BT"))
{
  max.temps<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in plot.dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      max.temps<-rbind(max.temps,data.frame(date,max(sub.temp.rh.2$temp.c,na.rm = T)))
    }
  }
  points(as.Date(max.temps[,1],origin="1970-01-01"),max.temps[,2],type="l",lwd=2)
}

min.temp.2020<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/tasmin.csv",header = F)
min.temp.2045<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/tasmin.csv",header = F)
min.temp.2070<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/tasmin.csv",header = F)
colnames(min.temp.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(min.temp.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(min.temp.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
plot(dates,min.temp.2020$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(-10,10),xlim=c(min(plot.dates),max(plot.dates)),xlab="date",ylab="min temp",lty=1,lwd=2)        
axis.Date(1,plot.dates)
points(dates,min.temp.2045$cesm1.cam5.1.rcp45,type="l",col="orange",lty=2,lwd=2)
points(dates,min.temp.2070$cesm1.cam5.1.rcp45,type="l",col="orange",lty=3,lwd=2)
points(dates,min.temp.2020$cesm1.cam5.1.rcp85,type="l",col="red",lty=1,lwd=2)
points(dates,min.temp.2045$cesm1.cam5.1.rcp85,type="l",col="red",lty=2,lwd=2)
points(dates,min.temp.2070$cesm1.cam5.1.rcp85,type="l",col="red",lty=3,lwd=2)

for(site in c("BT"))
{
  min.temps<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in plot.dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      min.temps<-rbind(min.temps,data.frame(date,min(sub.temp.rh.2$temp.c,na.rm = T)))
    }
  }
  points(as.Date(min.temps[,1],origin="1970-01-01"),min.temps[,2],type="l",lwd=2)
}


mean.rh.2020<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/relHumid.csv",header = F)
mean.rh.2045<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/relHumid.csv",header = F)
mean.rh.2070<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/relHumid.csv",header = F)
colnames(mean.rh.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.rh.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.rh.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
plot(dates,13.24732*exp((17.67*mean.temp.2020$cesm1.cam5.1.rcp45)/(mean.temp.2020$cesm1.cam5.1.rcp45+243.5))*mean.rh.2020$cesm1.cam5.1.rcp45/(273.15+mean.temp.2020$cesm1.cam5.1.rcp45),type="l",col="orange",xlim=c(min(plot.dates),max(plot.dates)),ylim=c(2,10),xlab="date",ylab="mean AH",lty=1,lwd=2)  
axis.Date(1,plot.dates)
points(dates,13.24732*exp((17.67*mean.temp.2045$cesm1.cam5.1.rcp45)/(mean.temp.2045$cesm1.cam5.1.rcp45+243.5))*mean.rh.2045$cesm1.cam5.1.rcp45/(273.15+mean.temp.2045$cesm1.cam5.1.rcp45),type="l",col="orange",lty=2,lwd=2)
points(dates,13.24732*exp((17.67*mean.temp.2070$cesm1.cam5.1.rcp45)/(mean.temp.2070$cesm1.cam5.1.rcp45+243.5))*mean.rh.2070$cesm1.cam5.1.rcp45/(273.15+mean.temp.2070$cesm1.cam5.1.rcp45),type="l",col="orange",lty=3,lwd=2)
points(dates,13.24732*exp((17.67*mean.temp.2020$cesm1.cam5.1.rcp85)/(mean.temp.2020$cesm1.cam5.1.rcp85+243.5))*mean.rh.2020$cesm1.cam5.1.rcp85/(273.15+mean.temp.2020$cesm1.cam5.1.rcp85),type="l",col="red",lty=2,lwd=2)
points(dates,13.24732*exp((17.67*mean.temp.2045$cesm1.cam5.1.rcp85)/(mean.temp.2045$cesm1.cam5.1.rcp85+243.5))*mean.rh.2045$cesm1.cam5.1.rcp85/(273.15+mean.temp.2045$cesm1.cam5.1.rcp85),type="l",col="red",lty=2,lwd=2)
points(dates,13.24732*exp((17.67*mean.temp.2070$cesm1.cam5.1.rcp85)/(mean.temp.2070$cesm1.cam5.1.rcp85+243.5))*mean.rh.2070$cesm1.cam5.1.rcp85/(273.15+mean.temp.2070$cesm1.cam5.1.rcp85),type="l",col="red",lty=3,lwd=2)


for(site in c("BT"))
{
  mean.abs.hums<-c()
  sub.temp.rh.1<-all.temp.rh[which(all.temp.rh$site==site),]
  for(date in plot.dates)
  {
    if(date %in% as.Date(sub.temp.rh.1$date.time))
    {
      sub.temp.rh.2<-sub.temp.rh.1[which(as.Date(sub.temp.rh.1$date)==date),]
      abs.hums<-13.24732*exp((17.67*sub.temp.rh.2$temp.c)/(sub.temp.rh.2$temp.c+243.5))*sub.temp.rh.2$rh/(273.15+sub.temp.rh.2$temp.c)
      abs.hums<-6.112*exp((17.67*sub.temp.rh.2$temp.c/(sub.temp.rh.2$temp.c+243.5))*(sub.temp.rh.2$rh/100)*18.02)/((273.15+sub.temp.rh.2$temp.c)*100*0.08314)
      abs.hums<-6.112*exp(17.67*sub.temp.rh.2$temp.c/(sub.temp.rh.2$temp.c+243.5) )*sub.temp.rh.2$rh*2.1674/(273.15+sub.temp.rh.2$temp.c)
      mean.abs.hums<-rbind(mean.abs.hums,data.frame(date,mean(abs.hums,na.rm = T)))
    }
  }
  points(as.Date(mean.abs.hums[,1],origin="1970-01-01"),mean.abs.hums[,2],type="l",lwd=2)
}

mean.rainfall.2020<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2020/rainfall.csv",header = F)
mean.rainfall.2045<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2045/rainfall.csv",header = F)
mean.rainfall.2070<-read.csv("~/Documents/GitHub/flax.rust/data/enviro/climate.projections/2070/rainfall.csv",header = F)
colnames(mean.rainfall.2020)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.rainfall.2045)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
colnames(mean.rainfall.2070)<-c("year","month","day","cesm1.cam5.1.rcp45","cesm1.cam5.1.rcp85")
plot(dates,mean.rainfall.2020$cesm1.cam5.1.rcp45,type="l",col="orange",ylim=c(0,27),xlim=c(min(plot.dates),max(plot.dates)),xlab="date",ylab="mean rainfall",lty=1,lwd=2)        
axis.Date(1,plot.dates)
points(dates,mean.rainfall.2045$cesm1.cam5.1.rcp45,type="l",col="orange",lty=2,lwd=2)
points(dates,mean.rainfall.2070$cesm1.cam5.1.rcp45,type="l",col="orange",lty=3,lwd=2)
points(dates,mean.rainfall.2020$cesm1.cam5.1.rcp85,type="l",col="red",lty=1,lwd=2)
points(dates,mean.rainfall.2045$cesm1.cam5.1.rcp85,type="l",col="red",lty=2,lwd=2)
points(dates,mean.rainfall.2070$cesm1.cam5.1.rcp85,type="l",col="red",lty=3,lwd=2)


for(site in c("BT"))
{
  mean.rains<-c()
  all.weath.1<-all.weath[which(all.weath$site==site),]
  for(date in plot.dates)
  {
    if(date %in% as.Date(all.weath.1$date))
    {
      all.weath.2<-all.weath.1[which(as.Date(all.weath.1$date)==date),]
      mean.rains<-rbind(mean.rains,data.frame(date,mean(all.weath.2$rain,na.rm = T)*12*24))
    }
  }
  points(as.Date(mean.rains[,1],origin="1970-01-01"),mean.rains[,2],type="l",lwd=2)
}

plot(0,0,type="n",bty="n",xlab="",ylab="",axes=F,xlim=c(0,1),ylim=c(0,1))
segments(0,seq(0,1,length.out = 7)[7],.25,seq(0,1,length.out = 7)[7],col="black",lwd=2,lty=1)
segments(0,seq(0,1,length.out = 7)[6],.25,seq(0,1,length.out = 7)[6],col="orange",lwd=2,lty=1)
segments(0,seq(0,1,length.out = 7)[5],.25,seq(0,1,length.out = 7)[5],col="orange",lwd=2,lty=2)
segments(0,seq(0,1,length.out = 7)[4],.25,seq(0,1,length.out = 7)[4],col="orange",lwd=2,lty=3)
segments(0,seq(0,1,length.out = 7)[3],.25,seq(0,1,length.out = 7)[3],col="red",lwd=2,lty=1)
segments(0,seq(0,1,length.out = 7)[2],.25,seq(0,1,length.out = 7)[2],col="red",lwd=2,lty=2)
segments(0,seq(0,1,length.out = 7)[1],.25,seq(0,1,length.out = 7)[1],col="red",lwd=2,lty=3)
text(.25,seq(0,1,length.out = 7)[7],pos=4,"data",cex=1.5)
text(.25,seq(0,1,length.out = 7)[6],pos=4,"RCPP4.5 2020",cex=1.5)
text(.25,seq(0,1,length.out = 7)[5],pos=4,"RCPP4.5 2045",cex=1.5)
text(.25,seq(0,1,length.out = 7)[4],pos=4,"RCPP4.5 2070",cex=1.5)
text(.25,seq(0,1,length.out = 7)[3],pos=4,"RCPP8.5 2020",cex=1.5)
text(.25,seq(0,1,length.out = 7)[2],pos=4,"RCPP8.5 2045",cex=1.5)
text(.25,seq(0,1,length.out = 7)[1],pos=4,"RCPP8.5 2070",cex=1.5)

