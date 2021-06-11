# script to build dataset of observed/predicted plant heights on dates relevant for foi dataset building

## setup

### load models and funcitons
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant height change funcs.R")
### convenience params
data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20"),"BT"=c("2020-06-24","2020-07-01"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09"))
sites<-c("CC","BT","GM","HM")
### load data
corrected.epi<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.epi.RDS")
corrected.locs<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.locs.RDS")
within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
#### subset to relevant records
within.host<-within.host[which(within.host$Year==2020),]
within.host<-unique(within.host[,c("Year","Site","Tag","Date","max.height")])
healthyplants<-read.csv("~/Documents/GitHub/flax.rust/data/healthyplants.csv")
#### subset to relevant records
healthyplants<-unique(healthyplants[,c("Site","Date","Tag","max.height")])

## build data set 
corrected.plant.heights<-data.frame("Site"=factor(),"tag"=factor(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"Date"=character(),"height.cm"=numeric(),"height.type"=factor())
for (site in sites)
{
  ### subset data by site
  sub.loc.data<-corrected.locs[which(corrected.locs$Site==site),]
  sub.data.dates<-unname(unlist(data.dates[which(names(data.dates)==site)]))
  
  ### for each plant in the transect
  
  for(i in 1:nrow(sub.loc.data))
  {
    #### if plant is tagged
    if(!is.na(sub.loc.data[i,"tag"]))
    {
      tag<-sub.loc.data[i,"tag"]
      for(date in sub.data.dates)
      {
        date<-as.Date(date)
        ##### get indicies of any existing records
        healthyplants.index<-intersect(which(healthyplants$Tag==tag),which(as.Date(healthyplants$Date,tryFormats = "%m/%d/%Y")==as.Date(date)))
        within.host.index<-intersect(which(within.host$Tag==tag),which(as.Date(within.host$Date,tryFormats = "%m/%d/%Y")==as.Date(date)))
        ##### if plant height was recorded
        if(length(healthyplants.index)>0 | length(within.host.index)>0)
        {
          ###### if plant height was recorded in healthy plants
          if(length(healthyplants.index)>0)
          {
            corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=healthyplants[healthyplants.index,"max.height"],"height.type"="observed"))
          }
          ###### if plant height was recorded in within host
          if(length(within.host.index)>0)
          {
            corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=within.host[within.host.index,"max.height"],"height.type"="observed"))
          }
        }
        ##### if plant height was not recorded
        
        if(length(healthyplants.index)==0 & length(within.host.index)==0)
        {
          print(tag)
          ###### look for closest observation or existing forecast/hindcast
          obs.dates<-data.frame("Date"=character(),"diff"=numeric(),"source"=character(),"height"=numeric(),"type"=character())
          
          if(length(within.host[which(within.host$Tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(within.host[which(within.host$Tag==tag),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-within.host[intersect(which(as.Date(within.host[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(within.host$Tag==tag)),"max.height"]
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="within.host","height"=height,"type"="obs"))
          }
          
          if(length(healthyplants[which(healthyplants$Tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(healthyplants[which(healthyplants$Tag==tag),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-healthyplants[intersect(which(as.Date(healthyplants[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(healthyplants$Tag==tag)),"max.height"]
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="healthyplants","height"=height,"type"="obs"))
          }
          
          if(length(corrected.epi[which(corrected.epi$Tag==tag),"Date.First.Observed.Diseased"])>0)
          {
            record.dates<-as.Date(corrected.epi[which(corrected.epi$Tag==tag),"Date.First.Observed.Diseased"])
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-corrected.epi[intersect(which(as.Date(corrected.epi[,"Date.First.Observed.Diseased"],)==min.date),which(corrected.epi$Tag==tag)),"max.height"]
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.epi","height"=height,"type"="obs"))
          }
          
          if(length(corrected.locs[which(corrected.locs$tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(corrected.locs[which(corrected.locs$tag==tag),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-corrected.locs[intersect(which(as.Date(corrected.locs[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(corrected.locs$tag==tag)),"height.cm"]
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.locs","height"=height,"type"="obs"))
          }
          
          if(length(corrected.plant.heights[which(corrected.plant.heights$tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(corrected.plant.heights[which(corrected.plant.heights$tag==tag),"Date"])
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-corrected.plant.heights[intersect(which(as.Date(corrected.plant.heights[,"Date"],)==min.date),which(corrected.plant.heights$tag==tag)),"height.cm"]
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.plant.heights","height"=height,"type"="obs"))
          }
          
          ###### pinpoint observation to use for forecast/hindcast, with preference for forecast
          if(length(which(is.na(obs.dates$height)))>0) {obs.dates<-obs.dates[-which(is.na(obs.dates$height)),]}
          obs.dates.index<-which(abs(obs.dates$diff)==min(abs(obs.dates$diff)))
          if(length(obs.dates.index>1)) {obs.dates.index<-obs.dates.index[which.min(obs.dates$diff[obs.dates.index])]}
          
          ###### forecast
          if(obs.dates$diff[obs.dates.index]<0)
          {
            forecast.height<-predict.plant.growth(obs.dates[obs.dates.index,"height"],site,obs.dates[obs.dates.index,"Date"],date,exclude.site = F)
            corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=forecast.height,"height.type"="forecast"))
          }
          
          ###### hindcast
          if(obs.dates$diff[obs.dates.index]>0)
          {
            hindcast.height<-predict.plant.growth.last(obs.dates[obs.dates.index,"height"],site,date,obs.dates[obs.dates.index,"Date"],exclude.site = F)
            corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=hindcast.height,"height.type"="hindcast"))
            
          }
        }
      }
    }
    #### if plant is untagged
    if(is.na(sub.loc.data[i,"tag"]))
    {
      
    }
  }
}
### for plant in corrected.locs
### for date in data.dates at plant site
### if height data exists--store it
### if height data doesn't exist--hindcast/forecast