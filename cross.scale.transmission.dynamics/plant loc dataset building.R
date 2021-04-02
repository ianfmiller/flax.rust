# load + prep data and functions
## load functions
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant inf intens change funcs.R")
## load data
epi<-read.csv("~/Documents/GitHub/flax.rust/data/Epidemiology.csv")
within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
plant.locs<-read.csv("~/Documents/GitHub/flax.rust/data/plant locs.csv")
demog<-read.csv("~/Documents/GitHub/flax.rust/data/Demography.csv")
## prep data
epi<-epi[which(epi$Year==2020),]
epi$Date.First.Observed.Diseased<-as.Date(epi$Date.First.Observed.Diseased,tryFormats = "%m/%d/%Y")
within.host<-within.host[which(within.host$Year==2020),]
within.host$Date<-as.Date(within.host$Date,tryFormats = "%m/%d/%Y")
plant.locs<-plant.locs[which(plant.locs$Year==2020),]
demog<-demog[which(demog$year==2020),]
demog<-demog[which(demog$status %in% c("H","D")),]

# Time periods for fitting glm of infection~foi to data:
#CC: 6/22(first data)->7/27(last prediction)
#BT: 6/24(first data)->7/8(last prediction)
#GM: 6/16(first data)->7/15(last prediction)
#BT: 6/25(first data)->7/28(last prediction)
## These periods were selected such that the plant infection intensity metric was recorded for every infected plant in the transect for the observation prior to the observed changes in infection status
## Exceptions wer emade for:
### Infected seedling, which were assumed to have 0.1 starting inf intensity
### Miissing plant inf intensity measurments that could be fore/hindcasted from existing observations using using plant inf intens model. THe window for fore/hindcasting was limited to +/- 1 week frpm existing data
### Tag #15 in "CC" which was assumed to have 0 inf intensity as indicated by an observation on 6/22

sites<-c("CC","BT","GM","HM")
plant.loc.survey.dates<-list("CC"=c(as.Date("06/15/2020",tryFormats=c("%m/%d/%Y")),as.Date("06/22/2020",tryFormats=c("%m/%d/%Y"))),"BT"=as.Date("06/17/2020",tryFormats=c("%m/%d/%Y")),"GM"=as.Date("06/16/2020",tryFormats=c("%m/%d/%Y")),"HM"=as.Date("06/18/2020",tryFormats=c("%m/%d/%Y")))
epi.obs.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20","2020-07-27"),"BT"=c("2020-06-24","2020-07-01","2020-07-08"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09","2020-07-15"),"HM"=c("2020-06-25","2020-07-15","2020-07-21","2020-07-22","2020-07-28"))    
data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20"),"BT"=c("2020-06-24","2020-07-01"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09"),"HM"=c("2020-06-25","2020-07-15","2020-07-21","2020-07-22"))

corrected.plant.locs<-c()
for(site in sites)
{
  ## get plant location data
  loc.data<-plant.locs[which(plant.locs$Site==site),]
  loc.data<-loc.data[which(as.Date(loc.data$Date,tryFormats=c("%m/%d/%Y")) %in% plant.loc.survey.dates[which(sites==site)][[1]]),c("Site","X","Y","x","y","tag")]
  
  ## switch record for newly tagged plants not in loc.data for record of closest untagged plant
  for(date in plant.loc.survey.dates[which(sites==site)][[1]])
  {
    sub.epi.data<-epi[which(epi$Site==site),]
    sub.epi.data<-sub.epi.data[which(as.Date(sub.epi.data$Date.First.Observed.Diseased)==date),]
    
    for(index in which(!(sub.epi.data$Tag %in% loc.data$tag)))
    {
      ### condition to exclude data in unsurveyed region of CC
      if(
        any(
            c(
              (!(site=="CC")),
              (sub.epi.data[index,"Y"] %in% c(0:8,18,19)),
              all(sub.epi.data[index,"Y"] %in% c(15,17),sub.epi.data[index,"x"] %in% 7:9)
            )
        )
      )
      {
        min.dist.func<-function(i) {dist(rbind(c(sub.epi.data[index,"X"]+sub.epi.data[index,"x"],sub.epi.data[index,"Y"]+sub.epi.data[index,"y"]),c(loc.data[i,"X"]+loc.data[i,"x"],loc.data[i,"Y"]+loc.data[i,"y"])))[1]}
        distances<-sapply(which(is.na(loc.data$tag)),min.dist.func)
        if(min(distances)<=.25) 
        {
          #print(paste0("site = ",site," dist = ",min(distances)))
          loc.data<-loc.data[-which(is.na(loc.data$tag))[which.min(distances)],]
          loc.data<-rbind(loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=sub.epi.data[index,"Tag"]))
        } else {
          loc.data<-rbind(loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=sub.epi.data[index,"Tag"]))
        }
      }
    }
  }
  
  ## switch record for newly diseased plants not in loc.data for record of closest untagged plant
  for(date in epi.obs.dates[which(sites==site)][[1]])
  {
    sub.epi.data<-epi[which(epi$Site==site),]
    sub.epi.data<-sub.epi.data[which(as.Date(sub.epi.data$Date.First.Observed.Diseased)==date),]
    
    ### newly taggd plants
    #tagset<-sub.epi.data$Tag[-which(is.na(sub.epi.data$Tag))]
    tagset<-na.omit(sub.epi.data$Tag)
    for(index in which(!(tagset %in% loc.data$tag)))
    {
      ### condition to exclude data in unsurveyed region of CC
      if(
        any(
          c(
            (!(site=="CC")),
            (sub.epi.data[index,"Y"] %in% c(0:8,18,19)),
            all(sub.epi.data[index,"Y"] %in% c(15,17),sub.epi.data[index,"x"] %in% 7:9)
          )
        )
      )
      {
        min.dist.func<-function(i) {dist(rbind(c(sub.epi.data[index,"X"]+sub.epi.data[index,"x"],sub.epi.data[index,"Y"]+sub.epi.data[index,"y"]),c(loc.data[i,"X"]+loc.data[i,"x"],loc.data[i,"Y"]+loc.data[i,"y"])))[1]}
        distances<-sapply(which(is.na(loc.data$tag)),min.dist.func)
        if(min(distances)<=.25) 
        {
          #print(paste0("site = ",site," dist = ",min(distances)))
          loc.data<-loc.data[-which(is.na(loc.data$tag))[which.min(distances)],]
          loc.data<-rbind(loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=sub.epi.data[index,"Tag"]))
        } else {
          loc.data<-rbind(loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=sub.epi.data[index,"Tag"]))
        }
      }
    }
    
    ### untagged plants
    NAset<-which(is.na(sub.epi.data$Tag))
    for(index in NAset)
    {
      ### condition to exclude data in unsurveyed region of CC
      if(
        any(
          c(
            (!(site=="CC")),
            (sub.epi.data[index,"Y"] %in% c(0:8,18,19)),
            all(sub.epi.data[index,"Y"] %in% c(15,17),sub.epi.data[index,"x"] %in% 7:9)
          )
        )
      )
      {
        min.dist.func<-function(i) {dist(rbind(c(sub.epi.data[index,"X"]+sub.epi.data[index,"x"],sub.epi.data[index,"Y"]+sub.epi.data[index,"y"]),c(loc.data[i,"X"]+loc.data[i,"x"],loc.data[i,"Y"]+loc.data[i,"y"])))[1]}
        distances<-sapply(which(is.na(loc.data$tag)),min.dist.func)
        if(min(distances)<=.25) 
        {
          #print(paste0("site = ",site," dist = ",min(distances)))
          loc.data<-loc.data[-which(is.na(loc.data$tag))[which.min(distances)],]
          loc.data<-rbind(loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=NA))
        } else {
          loc.data<-rbind(loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=NA))
        }
      }

    }
  }
  corrected.plant.locs<-rbind(corrected.plant.locs,loc.data)
}

# vis loc adjustment
par(mfrow=c(2,2))
for(site0 in sites)
{
  xx<-plant.locs[which(plant.locs$Site==site0),]
  yy<-corrected.plant.locs[which(corrected.plant.locs$Site==site0),]
  plot(xx$X+xx$x,xx$Y+xx$y,xlim=c(0,10),ylim=c(0,20),xlab="X",ylab="Y",main=site0,pch=16,col="black",cex=.9)
  points(yy$X+yy$x,yy$Y+yy$y,col="red",cex=)
}



