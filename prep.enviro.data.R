### prep enviro data
### only 2020 data
### notes: GM soil moisture hits negative values, and also has some pretty extreme jumps--could potentially indicate issues w/ sensor

library(lubridate)

make.transparent<-function(col,alpha)
{
  rgb<-col2rgb(col)
  r<-rgb[1]/255
  g<-rgb[2]/255
  b<-rgb[3]/255
  rgb(r,g,b,alpha)
}

# temp/humidity data

## load
hm<-read.csv("~/Documents/Github/flax.rust/data/enviro/high meadow temp RH 2020.csv",skip=2)
gm<-read.csv("~/Documents/Github/flax.rust/data/enviro/gothic mountain temp RH 2020.csv",skip=2)
bt<-read.csv("~/Documents/Github/flax.rust/data/enviro/bus turnaround temp RH 2020.csv",skip=2)
cc<-read.csv("~/Documents/Github/flax.rust/data/enviro/cement creek temp RH 2020.csv",skip=2)

## clean
colnames(hm)[1:4]<-c("date.time","temp.c","rh","dew.pt.c")
colnames(gm)[1:4]<-c("date.time","temp.c","rh","dew.pt.c")
colnames(bt)[1:4]<-c("date.time","temp.c","rh","dew.pt.c")
colnames(cc)[1:4]<-c("date.time","temp.c","rh","dew.pt.c")

hm[,1]<-ymd_hms(hm[,1])
gm[,1]<-ymd_hms(gm[,1])
bt[,1]<-ymd_hms(bt[,1])
cc[,1]<-ymd_hms(cc[,1])

## join
all.temp.rh<-rbind(hm[,1:4],gm[,1:4],bt[,1:4],cc[,1:4])

# hm weather data
## load data, clean
hm.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/High_Meadow.csv",skip=2)[,-1]
colnames(hm.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
hm.weath[,1]<-parse_date_time(hm.weath[,1],'%m/%d/%y %I:%M:%S %p')
hm.weath<-hm.weath[-which(hm.weath$wetness==-888.88),]

## join with temperature data for calculating variables capturing co-occurence of weather vars and temp
temps<-c()
temps.16.22<-c()
temps.7.30<-c()
for(i in 1:(dim(hm.weath)[1]))
{
  date0<-as.POSIXct(hm.weath[i,"date"])
  date1<-as.POSIXct(hm.weath[i+1,"date"])
  interval.temps<-hm[intersect(which(hm$date.time>=date0),which(hm$date.time<=date1)),"temp.c"]
  ifelse(all(all(interval.temps>=16),all(interval.temps<=22)),interval.16.22<-1,interval.16.22<-0)
  ifelse(all(all(interval.temps>=7),all(interval.temps<=30)),interval.7.30<-1,interval.7.30<-0)
  temps<-c(temps,mean(interval.temps))
  temps.16.22<-c(temps.16.22,interval.16.22)
  temps.7.30<-c(temps.7.30,interval.7.30)
}
hm.weath<-cbind(hm.weath,temp=temps,temp.16.22=temps.16.22,temp.7.30=temps.7.30)

# gm weather data
## load data, clean
gm.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/Gothic_Mountain.csv",skip=2)[,-1]
colnames(gm.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
gm.weath[,1]<-parse_date_time(gm.weath[,1],'%m/%d/%y %I:%M:%S %p')
gm.weath<-gm.weath[-which(gm.weath$wetness==-888.88),]

## join with temperature data for calculating variables capturing co-occurence of weather vars and temp
temps<-c()
temps.16.22<-c()
temps.7.30<-c()
for(i in 1:(dim(gm.weath)[1]))
{
  date0<-as.POSIXct(gm.weath[i,"date"])
  date1<-as.POSIXct(gm.weath[i+1,"date"])
  interval.temps<-gm[intersect(which(gm$date.time>=date0),which(gm$date.time<=date1)),"temp.c"]
  ifelse(all(all(interval.temps>=16),all(interval.temps<=22)),interval.16.22<-1,interval.16.22<-0)
  ifelse(all(all(interval.temps>=7),all(interval.temps<=30)),interval.7.30<-1,interval.7.30<-0)
  temps<-c(temps,mean(interval.temps))
  temps.16.22<-c(temps.16.22,interval.16.22)
  temps.7.30<-c(temps.7.30,interval.7.30)
}
gm.weath<-cbind(gm.weath,temp=temps,temp.16.22=temps.16.22,temp.7.30=temps.7.30)

# bt weather data
## load data, clean
bt.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/Bus_Turnaround.csv",skip=2)[,-1]
colnames(bt.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
bt.weath[,1]<-parse_date_time(bt.weath[,1],'%m/%d/%y %I:%M:%S %p')
bt.weath<-bt.weath[-which(bt.weath$wetness==-888.88),]

## join with temperature data for calculating variables capturing co-occurence of weather vars and temp
temps<-c()
temps.16.22<-c()
temps.7.30<-c()
for(i in 1:(dim(bt.weath)[1]))
{
  date0<-as.POSIXct(bt.weath[i,"date"])
  date1<-as.POSIXct(bt.weath[i+1,"date"])
  interval.temps<-bt[intersect(which(bt$date.time>=date0),which(bt$date.time<=date1)),"temp.c"]
  ifelse(all(all(interval.temps>=16),all(interval.temps<=22)),interval.16.22<-1,interval.16.22<-0)
  ifelse(all(all(interval.temps>=7),all(interval.temps<=30)),interval.7.30<-1,interval.7.30<-0)
  temps<-c(temps,mean(interval.temps))
  temps.16.22<-c(temps.16.22,interval.16.22)
  temps.7.30<-c(temps.7.30,interval.7.30)
}
bt.weath<-cbind(bt.weath,temp=temps,temp.16.22=temps.16.22,temp.7.30=temps.7.30)


# cc weather data
## load data, clean
cc.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/Cement_Creek.csv",skip=2)[,-1]
colnames(cc.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
cc.weath[,1]<-parse_date_time(cc.weath[,1],'%m/%d/%y %I:%M:%S %p')
cc.weath<-cc.weath[-which(cc.weath$wetness==-888.88),]

## join with temperature data for calculating variables capturing co-occurence of weather vars and temp
temps<-c()
temps.16.22<-c()
temps.7.30<-c()
for(i in 1:(dim(cc.weath)[1]))
{
  date0<-as.POSIXct(cc.weath[i,"date"])
  date1<-as.POSIXct(cc.weath[i+1,"date"])
  interval.temps<-cc[intersect(which(cc$date.time>=date0),which(cc$date.time<=date1)),"temp.c"]
  ifelse(all(all(interval.temps>=16),all(interval.temps<=22)),interval.16.22<-1,interval.16.22<-0)
  ifelse(all(all(interval.temps>=7),all(interval.temps<=30)),interval.7.30<-1,interval.7.30<-0)
  temps<-c(temps,mean(interval.temps))
  temps.16.22<-c(temps.16.22,interval.16.22)
  temps.7.30<-c(temps.7.30,interval.7.30)
}
cc.weath<-cbind(cc.weath,temp=temps,temp.16.22=temps.16.22,temp.7.30=temps.7.30)

#join all weather data
all.weath<-rbind(hm.weath,gm.weath,bt.weath,cc.weath)

