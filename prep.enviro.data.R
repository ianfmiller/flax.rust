### prep enviro data
### only 2020 data
### notes: GM soil moisture hits negative values, and also has some pretty extreme jumps--could potentially indicate issues w/ sensor

library(lubridate)

setwd("~/Documents/Github/flax.rust/data/enviro")

make.transparent<-function(col,alpha)
{
  rgb<-col2rgb(col)
  r<-rgb[1]/255
  g<-rgb[2]/255
  b<-rgb[3]/255
  rgb(r,g,b,alpha)
}

hm<-read.csv("high meadow temp RH 2020.csv",skip=2)
gm<-read.csv("gothic mountain temp RH 2020.csv",skip=2)
bt<-read.csv("bus turnaround temp RH 2020.csv",skip=2)
cc<-read.csv("cement creek temp RH 2020.csv",skip=2)

hm[,1]<-ymd_hms(hm[,1])
gm[,1]<-ymd_hms(gm[,1])
bt[,1]<-ymd_hms(bt[,1])
cc[,1]<-ymd_hms(cc[,1])
all.temp.rh<-rbind(hm[,1:4],gm[,1:4],bt[,1:4],cc[,1:4])


hm.weath<-read.csv("High_Meadow.csv",skip=2)[,-1]
colnames(hm.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
hm.weath[,1]<-parse_date_time(hm.weath[,1],'%m/%d/%y %I:%M:%S %p')
hm.weath<-hm.weath[-which(hm.weath$wetness==-888.88),]

gm.weath<-read.csv("Gothic_Mountain.csv",skip=2)[,-1]
colnames(gm.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
gm.weath[,1]<-parse_date_time(gm.weath[,1],'%m/%d/%y %I:%M:%S %p')
gm.weath<-gm.weath[-which(gm.weath$wetness==-888.88),]

bt.weath<-read.csv("Bus_Turnaround.csv",skip=2)[,-1]
colnames(bt.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
bt.weath[,1]<-parse_date_time(bt.weath[,1],'%m/%d/%y %I:%M:%S %p')
bt.weath<-bt.weath[-which(bt.weath$wetness==-888.88),]

cc.weath<-read.csv("Cement_Creek.csv",skip=2)[,-1]
colnames(cc.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
cc.weath[,1]<-parse_date_time(cc.weath[,1],'%m/%d/%y %I:%M:%S %p')
cc.weath<-cc.weath[-which(cc.weath$wetness==-888.88),]

all.weath<-rbind(hm.weath,gm.weath,bt.weath,cc.weath)