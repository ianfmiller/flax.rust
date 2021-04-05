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

# loop for compiling corrected plant locations
corrected.epi<-epi
corrected.locs<-c()

for(site in sites)
{
  ## get plant location data
  sub.loc.data<-plant.locs[which(plant.locs$Site==site),]
  sub.loc.data<-sub.loc.data[which(as.Date(sub.loc.data$Date,tryFormats=c("%m/%d/%Y")) %in% plant.loc.survey.dates[which(sites==site)][[1]]),c("Site","X","Y","x","y","tag")]
  sub.loc.data<-data.frame(sub.loc.data,"matched"=F)
  sub.loc.data[which(!(is.na(sub.loc.data$tag))),"matched"]<-T
  
  dates<-unique(corrected.epi[which(corrected.epi$Site==site),"Date.First.Observed.Diseased"])

  for(date in dates)
  {
    ## for newly diseased and tagged plants not in loc.data: edit location in epi data to location of closest untagged plant
    sub.epi.data<-epi[which(epi$Site==site),]
    sub.epi.data<-sub.epi.data[which(as.Date(sub.epi.data$Date.First.Observed.Diseased)==date),]
    
    for(tag in sub.epi.data$Tag[!(sub.epi.data$Tag %in% sub.loc.data$tag)])
    {
      index<-which(sub.epi.data$Tag==tag)
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
        min.dist.func<-function(i) {dist(rbind(c(sub.epi.data[index,"X"]+sub.epi.data[index,"x"],sub.epi.data[index,"Y"]+sub.epi.data[index,"y"]),c(sub.loc.data[i,"X"]+sub.loc.data[i,"x"],sub.loc.data[i,"Y"]+sub.loc.data[i,"y"])))[1]}
        distances<-sapply(which(sub.loc.data$matched==F),min.dist.func)
        #replace data if there's a close match
        if(min(distances)<=.25) 
        {
          corrected.epi[which(corrected.epi$Tag==tag),c("X","Y","x","y")]<-sub.loc.data[which(sub.loc.data$matched==F)[which.min(distances)],c("X",'Y',"x","y")]
          sub.loc.data[intersect(which(is.na(sub.loc.data$tag)),which(sub.loc.data$matched==F))[which.min(distances)],"matched"]<-T
        } 
        # add a new record if there's not a close match
        if(min(distances)>.25) 
        {
          sub.loc.data<-rbind(sub.loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=sub.epi.data[index,"Tag"],"matched"=T))
        }
      }
    }
    
    ## for newly diseased and tagged plants in loc.data: edit location in epi data to location of plant in loc data
    for(tag in sub.epi.data$Tag[sub.epi.data$Tag %in% sub.loc.data$tag])
    {
      corrected.epi.index<-which(corrected.epi$Tag==tag)
      sub.loc.index<-which(sub.loc.data$tag==tag)
      corrected.epi[corrected.epi.index,c("X","Y","x","y")]<-sub.loc.data[sub.loc.index,c("X","Y","x","y")]
    }
    
    ## for newly diseased and untagged plants in loc.data: edit location in epi data to location of nearest unmatched plant in loc data
    for(index in which(is.na(sub.epi.data$Tag)))
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
        min.dist.func<-function(i) {dist(rbind(c(sub.epi.data[index,"X"]+sub.epi.data[index,"x"],sub.epi.data[index,"Y"]+sub.epi.data[index,"y"]),c(sub.loc.data[i,"X"]+sub.loc.data[i,"x"],sub.loc.data[i,"Y"]+sub.loc.data[i,"y"])))[1]}
        distances<-sapply(which(sub.loc.data$matched==F),min.dist.func)
        #replace data if there's a close match
        if(min(distances)<=.25) 
        {
          corrected.epi.index<-intersect(which(corrected.epi$Site==site),which(corrected.epi$X==sub.epi.data[index,"X"]))
          corrected.epi.index<-intersect(corrected.epi.index,which(corrected.epi$Y==sub.epi.data[index,"Y"]))
          corrected.epi.index<-intersect(corrected.epi.index,which(corrected.epi$x==sub.epi.data[index,"x"]))
          corrected.epi.index<-intersect(corrected.epi.index,which(corrected.epi$y==sub.epi.data[index,"y"]))
          corrected.epi.index<-corrected.epi.index[1]
          corrected.epi[corrected.epi.index,c("X","Y","x","y")]<-sub.loc.data[which(sub.loc.data$matched==F)[which.min(distances)],c("X",'Y',"x","y")]
          sub.loc.data[intersect(which(is.na(sub.loc.data$tag)),which(sub.loc.data$matched==F))[which.min(distances)],"matched"]<-T
        } 
        # add a new record if there's not a close match
        if(min(distances)>.25) 
        {
          sub.loc.data<-rbind(sub.loc.data,data.frame("Site"=site,"X"=sub.epi.data[index,"X"],"Y"=sub.epi.data[index,"Y"],"x"=sub.epi.data[index,"x"],"y"=sub.epi.data[index,"y"],"tag"=sub.epi.data[index,"Tag"],"matched"=T))
        }
      }
    }
  }
  corrected.locs<-rbind(corrected.locs,sub.loc.data)
}
  
# vis loc agreement between corrected epi and loc data
layout(matrix(c(1,2,5,3,4,5),2,3,byrow = T))
for(site0 in sites)
{
  xx<-corrected.locs[which(corrected.locs$Site==site0),]
  yy<-corrected.epi[which(corrected.epi$Site==site0),]
  plot(xx$X+xx$x,xx$Y+xx$y,xlim=c(0,10),ylim=c(0,20),xlab="X",ylab="Y",main=site0,pch=16,col="black",cex=.9)
  points(yy$X+yy$x,yy$Y+yy$y,col="red",cex=1)
  # shade out unsurveyed region of CC
  if(site0=="CC")
  {
    rect(0,9,15,15,col="black",density = 25,border=NA)
    rect(0,15,7,16,col="black",density = 25,border=NA)
    rect(0,16,10,17,col="black",density = 25,border=NA)
    rect(0,17,7,18,col="black",density = 25,border=NA) 
  }
}
plot(0,0,type="n",axes = F,xlab="",ylab="")
legend("center",legend = c("location data + additions","corrected epi data location"),col=c("black","red"),pch = c(16,1),cex=.9,1)

# view individual plants added to loc data
layout(matrix(c(1,2,5,3,4,5),2,3,byrow = T))
for(site0 in sites)
{
  xx<-plant.locs[which(plant.locs$Site==site0),]
  yy<-corrected.locs[which(corrected.locs$Site==site0),]
  plot(xx$X+xx$x,xx$Y+xx$y,xlim=c(0,10),ylim=c(0,20),xlab="X",ylab="Y",main=site0,pch=16,col="black",cex=.9)
  points(yy$X+yy$x,yy$Y+yy$y,col="red",cex=1)
}
plot(0,0,type="n",axes = F,xlab="",ylab="")
legend("center",legend = c("original location data","location data + additions"),col=c("black","red"),pch = c(16,1),cex=.9,1)

