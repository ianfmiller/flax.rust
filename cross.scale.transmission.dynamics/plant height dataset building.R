# script to build dataset of observed/predicted plant heights on dates relevant for foi dataset building

## setup

### load models and funcitons
plant.growth.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant height change funcs.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/starting plant inf intens model.R")
### convenience params
data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20","2020-07-27"),"BT"=c("2020-06-24","2020-07-01","2020-07-08"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09","2020-07-15"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09","2020-07-15"))
sites<-c("CC","BT","GM","HM")
### load data
corrected.epi<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.epi.RDS")
corrected.locs<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.locs.RDS")

corrected.locs$Date[which(corrected.locs$Date=="6/15/20")]<-"6/15/2020"
corrected.locs$Date[which(corrected.locs$Date=="6/16/20")]<-"6/16/2020"
corrected.locs$Date[which(corrected.locs$Date=="6/17/20")]<-"6/17/2020"
corrected.locs$Date[which(corrected.locs$Date=="6/18/20")]<-"6/18/2020"
corrected.locs$Date[which(corrected.locs$Date=="6/22/20")]<-"6/22/2020"

diseased.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/diseased.focal.plants.RDS")
healthy.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/healthy.focal.plants.RDS")

 ## build data set 
corrected.plant.heights<-data.frame("Site"=factor(),"tag"=factor(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"Date"=character(),"height.cm"=numeric(),"height.type"=factor())
for (site in sites)
{
  ### subset data by site
  sub.loc.data<-corrected.locs[which(corrected.locs$Site==site),]
  sub.epi.data<-corrected.epi[which(corrected.epi$Site==site),]
  sub.data.dates<-unname(unlist(data.dates[which(names(data.dates)==site)]))
  
  ### create ID strings
  loc.ID.strings<-c()
  for(i in 1:dim(sub.loc.data)[1]) 
  {
    loc.ID.strings<-c(loc.ID.strings,paste0("Tag"=sub.loc.data[i,"tag"],"X=",sub.loc.data[i,"X"],"Y=",sub.loc.data[i,"Y"],"x=",sub.loc.data[i,"x"],"y=",sub.loc.data[i,"y"]))
  }
  
  epi.ID.strings<-c()
  for(i in 1:dim(sub.epi.data)[1]) 
  {
    epi.ID.strings<-c(epi.ID.strings,paste0("Tag"=sub.loc.data[i,"tag"],"X=",sub.epi.data[i,"X"],"Y=",sub.epi.data[i,"Y"],"x=",sub.epi.data[i,"x"],"y=",sub.epi.data[i,"y"]))
  }
  
  ### for each plant in the transect
  
  for(i in 1:nrow(sub.loc.data))
  {
    print(paste0(site," ",i))
    #### if plant is tagged
    if(!is.na(sub.loc.data[i,"tag"]))
    {
      tag<-sub.loc.data[i,"tag"]
      for(date in sub.data.dates)
      {
        date<-as.Date(date)
        ##### get indicies of any existing records
        healthy.focal.index<-intersect(which(healthy.focal.plants$Tag==tag),which(as.Date(healthy.focal.plants$Date,tryFormats = "%m/%d/%Y")==as.Date(date)))
        diseased.focal.index<-intersect(which(diseased.focal.plants$Tag==tag),which(as.Date(diseased.focal.plants$Date,tryFormats = "%m/%d/%Y")==as.Date(date)))
        sub.epi.index<-intersect(which(sub.epi.data$Tag==tag),which(sub.epi.data$Date.First.Observed.Diseased==as.Date(date)))
        sub.loc.index<-intersect(which(sub.loc.data$tag==tag),which(as.Date(sub.loc.data$Date,tryFormats = "%m/%d/%Y")==as.Date(date)))
        ##### if plant height was recorded
        if(length(healthy.focal.index)>0 || length(diseased.focal.index)>0 || length(sub.epi.index)>0 || length(sub.loc.index>0))
        {
          ###### if plant height was recorded in healthy plants
          if(length(healthy.focal.index)>0)
          {
            if(!is.na(healthy.focal.plants$max.height[healthy.focal.index]))
            {
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=healthy.focal.plants[healthy.focal.index,"max.height"],"height.type"="observed"))
            }
          }
          ###### if plant height was recorded in within host
          if(length(diseased.focal.index)>0)
          {
            if(!is.na(diseased.focal.plants$max.height[diseased.focal.index]))
            {
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=diseased.focal.plants[diseased.focal.index,"max.height"],"height.type"="observed"))
            }
          }
          
          ###### if plant height was recorded in epi data
          if(length(sub.epi.index)>0)
          {
            if(!is.na(sub.epi.data$max.height[sub.epi.index]))
            {
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=sub.epi.data[sub.epi.index,"max.height"],"height.type"="observed"))
            }
          }
          ###### if plant height was recorded in location data (includes demog data)
          if(length(sub.loc.index)>0)
          {
            if(!is.na(sub.loc.data$height.cm[sub.loc.index]))
            {
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=sub.loc.data[sub.loc.index,"height.cm"],"height.type"="observed"))
            }
          }
        }
        ##### if plant height was not recorded
        
        if(length(healthy.focal.index)==0 &&  length(diseased.focal.index)==0 && length(sub.epi.index)==0 && length(sub.loc.index)==0)
        {
          ###### look for closest observation or existing forecast/hindcast
          obs.dates<-data.frame("Date"=character(),"diff"=numeric(),"source"=character(),"height"=numeric(),"height.type"=character(),"infection.intensity"=numeric(),"infection.intensity.type"=character())
          
          if(length(diseased.focal.plants[which(diseased.focal.plants$Tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(diseased.focal.plants[which(diseased.focal.plants$Tag==tag),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-diseased.focal.plants[intersect(which(as.Date(diseased.focal.plants[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(diseased.focal.plants$Tag==tag)),"max.height"]
            inf.intens<-diseased.focal.plants[intersect(which(as.Date(diseased.focal.plants[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(diseased.focal.plants$Tag==tag)),"inf.intens"]
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="diseased.focal.plants","height"=height,"height.type"="obs","infection.intensity"=inf.intens,"infection.intensity.type"="obs"))
          }
          
          if(length(healthy.focal.plants[which(healthy.focal.plants$Tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(healthy.focal.plants[which(healthy.focal.plants$Tag==tag),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-healthy.focal.plants[intersect(which(as.Date(healthy.focal.plants[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(healthy.focal.plants$Tag==tag)),"max.height"]
            inf.intens<-0
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="healthy.focal.plants","height"=height,"height.type"="obs","infection.intensity"=inf.intens,"infection.intensity.type"="obs"))
          }
          
          if(length(sub.epi.data[which(sub.epi.data$Tag==tag),"Date.First.Observed.Diseased"])>0)
          {
            record.dates<-as.Date(sub.epi.data[which(sub.epi.data$Tag==tag),"Date.First.Observed.Diseased"])
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-sub.epi.data[intersect(which(as.Date(sub.epi.data[,"Date.First.Observed.Diseased"])==min.date),which(sub.epi.data$Tag==tag)),"max.height"]
            inf.intens<-mean.starting.inf.intens
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.epi","height"=height,"height.type"="obs","infection.intensity"=inf.intens,"infection.intensity.type"="mean.start"))
          }
          
          if(length(sub.loc.data[which(sub.loc.data$tag==tag),"Date"])>0)
          {
            record.dates<-as.Date(sub.loc.data[which(sub.loc.data$tag==tag),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-sub.loc.data[intersect(which(as.Date(sub.loc.data[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(sub.loc.data$tag==tag)),"height.cm"]
            inf.intens<-0
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.locs","height"=height,"height.type"="obs","infection.intensity"=inf.intens,"infection.intensity.type"="infered.zero"))
          }
          
          ###### pinpoint observation to use for forecast/hindcast, with preference for forecast
          if(length(which(is.na(obs.dates$height)))>0) {obs.dates<-obs.dates[-which(is.na(obs.dates$height)),]}
          if(site=="HM" & any(obs.dates$Date>max(sub.data.dates))) {obs.dates<-obs.dates[-which(obs.dates$Date>max(sub.data.dates)),]}
          obs.dates.index<-which(abs(obs.dates$diff)==min(abs(obs.dates$diff)))
          if(length(obs.dates.index)>1) {obs.dates.index<-obs.dates.index[which.min(obs.dates$diff[obs.dates.index])]} ######## preference for forecast is set here--note absence of abs() and use of min()
          
          ###### forecast
          if(obs.dates$diff[obs.dates.index]<0)
          {
            forecast.height<-predict.plant.growth(obs.dates[obs.dates.index,"height"],obs.dates[obs.dates.index,"infection.intensity"],site,obs.dates[obs.dates.index,"Date"],date,exclude.site = F)
            corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=forecast.height,"height.type"="forecast"))
          }
          
          ###### hindcast
          if(obs.dates$diff[obs.dates.index]>0)
          {
            hindcast.height<-predict.plant.growth.last(obs.dates[obs.dates.index,"height"],obs.dates[obs.dates.index,"infection.intensity"],site,date,obs.dates[obs.dates.index,"Date"],exclude.site = F)
            corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=tag,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=hindcast.height,"height.type"="hindcast"))
          }
        }
      }
    }
    #### if plant is untagged
    if(is.na(sub.loc.data[i,"tag"]))
    {
      ID.string<-paste0("Tag"=sub.loc.data[i,"tag"],"X=",sub.loc.data[i,"X"],"Y=",sub.loc.data[i,"Y"],"x=",sub.loc.data[i,"x"],"y=",sub.loc.data[i,"y"])
      
      for(date in sub.data.dates)
      {
        date<-as.Date(date)
        ##### get indicies of any existing records
        sub.epi.index<-intersect(which(epi.ID.strings==ID.string),which(sub.epi.data$Date.First.Observed.Diseased==as.Date(date)))
        sub.loc.index<-intersect(which(loc.ID.strings==ID.string),which(as.Date(sub.loc.data$Date,tryFormats = "%m/%d/%Y")==as.Date(date)))
        if(isTRUE((is.na(sub.epi.data[sub.epi.index,"max.height"])))) {sub.epi.index<-integer(0)} ###### check to make sure na values aren't accidentally counted
        if(isTRUE((is.na(sub.loc.data[sub.loc.index,"height.cm"])))) {sub.loc.index<-integer(0)} ###### check to make sure na values aren't accidentally counted
            
        ##### if plant height was recorded
        if(length(sub.epi.index)>0 || length(sub.loc.index)>0)
        {
          ###### if plant height was recorded in epi data
          if(length(sub.epi.index)>0 ) 
          {
            if(!is.na(sub.epi.data$max.height[sub.epi.index]))
            {
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=NA,"X"=sub.epi.data[sub.epi.index,"X"],"Y"=sub.epi.data[sub.epi.index,"Y"],"x"=sub.epi.data[sub.epi.index,"x"],"y"=sub.epi.data[sub.epi.index,"y"],"Date"=date,"height.cm"=sub.epi.data[sub.epi.index,"max.height"],"height.type"="observed"))
            }
          }
          ###### if plant height was recorded in location data (includes demog data)
          if(length(sub.loc.index)>0) 
          {
            if(!is.na(sub.loc.data$height.cm[sub.loc.index]))
            {
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=NA,"X"=sub.loc.data[sub.loc.index,"X"],"Y"=sub.loc.data[sub.loc.index,"Y"],"x"=sub.loc.data[sub.loc.index,"x"],"y"=sub.loc.data[sub.loc.index,"y"],"Date"=date,"height.cm"=sub.loc.data[sub.loc.index,"height.cm"],"height.type"="observed"))
            }
          }
        }
        ##### if plant height was not recorded
        if(length(sub.epi.index)==0 && length(sub.loc.index)==0)
        {
          ###### look for closest observation or existing forecast/hindcast
          obs.dates<-data.frame("Date"=character(),"diff"=numeric(),"source"=character(),"height"=numeric(),"height.type"=character(),"infection.intensity"=numeric(),"infection.intensity.type"=character())
          
          if(length(sub.epi.data[which(epi.ID.strings==ID.string),"Date.First.Observed.Diseased"])>0)
          {
            record.dates<-as.Date(sub.epi.data[which(epi.ID.strings==ID.string),"Date.First.Observed.Diseased"])
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-sub.epi.data[intersect(which(as.Date(sub.epi.data[,"Date.First.Observed.Diseased"])==min.date),which(epi.ID.strings==ID.string)),"max.height"]
            inf.intens<-mean.starting.inf.intens
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.epi","height"=height,"height.type"="obs","infection.intensity"=inf.intens,"infection.intensity.type"="mean.start"))
          }
          
          if(length(sub.loc.data[which(loc.ID.strings==ID.string),"Date"])>0)
          {
            record.dates<-as.Date(sub.loc.data[which(loc.ID.strings==ID.string),"Date"],tryFormats = "%m/%d/%Y")
            min.date<-record.dates[which.min(abs(difftime(record.dates,date)))]
            height<-sub.loc.data[intersect(which(as.Date(sub.loc.data[,"Date"],tryFormats = "%m/%d/%Y")==min.date),which(loc.ID.strings==ID.string)),"height.cm"]
            inf.intens<-0
            obs.dates<-rbind(obs.dates,data.frame("Date"=min.date,"diff"=as.numeric(difftime(min.date,date,units = "days")),"source"="corrected.locs","height"=height,"height.type"="obs","infection.intensity"=inf.intens,"infection.intensity.type"="infered.zero"))
          }
          
          ###### pinpoint observation to use for forecast/hindcast, with preference for forecast
          if(any(!is.na(obs.dates$height)))
          {
            if(length(which(is.na(obs.dates$height)))>0) {obs.dates<-obs.dates[-which(is.na(obs.dates$height)),]}
            if(site=="HM" & any(obs.dates$Date>max(sub.data.dates))) {obs.dates<-obs.dates[-which(obs.dates$Date>max(sub.data.dates)),]}
            obs.dates.index<-which(abs(obs.dates$diff)==min(abs(obs.dates$diff)))
            if(length(obs.dates.index)>1) {obs.dates.index<-obs.dates.index[which.min(obs.dates$diff[obs.dates.index])]}
            ###### forecast
            if(obs.dates$diff[obs.dates.index]<0)
            {
              forecast.height<-predict.plant.growth(obs.dates[obs.dates.index,"height"],obs.dates[obs.dates.index,"infection.intensity"],site,obs.dates[obs.dates.index,"Date"],date,exclude.site = F)
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=NA,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=forecast.height,"height.type"="forecast"))
            }
            ###### hindcast
            if(obs.dates$diff[obs.dates.index]>0)
            {
              hindcast.height<-predict.plant.growth.last(obs.dates[obs.dates.index,"height"],obs.dates[obs.dates.index,"infection.intensity"],site,date,obs.dates[obs.dates.index,"Date"],exclude.site = F)
              corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=NA,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=hindcast.height,"height.type"="hindcast"))
            }
          } else {corrected.plant.heights<-rbind(corrected.plant.heights,data.frame("Site"=site,"tag"=NA,"X"=sub.loc.data[i,"X"],"Y"=sub.loc.data[i,"Y"],"x"=sub.loc.data[i,"x"],"y"=sub.loc.data[i,"y"],"Date"=date,"height.cm"=NA,"height.type"="observed"))}

        }
      }
    }
    if(any(corrected.plant.heights$height.cm>100,na.rm = T)) {print(paste0("ERROR; index == ",i))}
  }
}

### trim any duplicates
corrected.plant.heights<-unique(corrected.plant.heights)

###set missing values to mean for site at date

for (i in which(is.na(corrected.plant.heights$height.cm)))
{
  site<-corrected.plant.heights[i,"Site"]
  date<-corrected.plant.heights[i,"Date"]
  sub.data<-corrected.plant.heights[intersect(which(corrected.plant.heights$Site==site),which(corrected.plant.heights$Date==date)),]
  corrected.plant.heights[i,"height.cm"]<-mean(na.omit(sub.data$height.cm))
}
saveRDS(corrected.plant.heights,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.plant.heights.RDS")
