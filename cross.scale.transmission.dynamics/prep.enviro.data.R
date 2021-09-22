if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/all.temp.rh.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/all.weath.RDS")))
{
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
  
  hm<-cbind("site"="HM",hm)
  gm<-cbind("site"="GM",gm)
  bt<-cbind("site"="BT",bt)
  cc<-cbind("site"="CC",cc)
  
  
  ## join
  all.temp.rh<-rbind(hm[,1:5],gm[,1:5],bt[,1:5],cc[,1:5])
  
  # hm weather data
  ## load data, clean
  hm.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/High_Meadow.csv",skip=2)[,-1]
  colnames(hm.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
  hm.weath[,1]<-parse_date_time(hm.weath[,1],'%m/%d/%y %I:%M:%S %p')
  hm.weath<-hm.weath[-which(hm.weath$wetness==-888.88),]
  hm.weath<-cbind("site"="HM",hm.weath)
  
  # gm weather data
  ## load data, clean
  gm.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/Gothic_Mountain.csv",skip=2)[,-1]
  colnames(gm.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
  gm.weath[,1]<-parse_date_time(gm.weath[,1],'%m/%d/%y %I:%M:%S %p')
  gm.weath<-gm.weath[-which(gm.weath$wetness==-888.88),]
  gm.weath<-cbind("site"="GM",gm.weath)
  
  # bt weather data
  ## load data, clean
  bt.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/Bus_Turnaround.csv",skip=2)[,-1]
  colnames(bt.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
  bt.weath[,1]<-parse_date_time(bt.weath[,1],'%m/%d/%y %I:%M:%S %p')
  bt.weath<-bt.weath[-which(bt.weath$wetness==-888.88),]
  bt.weath<-cbind("site"="BT",bt.weath)
  
  # cc weather data
  ## load data, clean
  cc.weath<-read.csv("~/Documents/Github/flax.rust/data/enviro/Cement_Creek.csv",skip=2)[,-1]
  colnames(cc.weath)<-c("date","wetness","rain","solar.radiation","wind.speed","gust.speed","wind.direction","soil.moisture")
  cc.weath[,1]<-parse_date_time(cc.weath[,1],'%m/%d/%y %I:%M:%S %p')
  cc.weath<-cc.weath[-which(cc.weath$wetness==-888.88),]
  cc.weath<-cbind("site"="CC",cc.weath)

  #join all weather data
  all.weath<-rbind(hm.weath,gm.weath,bt.weath,cc.weath)
  
  saveRDS(all.temp.rh,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/all.temp.rh.RDS")
  saveRDS(all.weath,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/all.weath.RDS")
}
all.temp.rh<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/all.temp.rh.RDS")
all.weath<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/all.weath.RDS")
