# load plant loc dataset
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")

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

# convenience pars
sites<-c("CC","BT","GM","HM")
epi.obs.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20","2020-07-27"),"BT"=c("2020-06-24","2020-07-01","2020-07-08"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09","2020-07-15"),"HM"=c("2020-06-25","2020-07-15","2020-07-21","2020-07-22","2020-07-28"))    
data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20"),"BT"=c("2020-06-24","2020-07-01"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09"),"HM"=c("2020-06-25","2020-07-15","2020-07-21","2020-07-22"))

# build empty data object for foi vs outcome data
foi.data<-data.frame("site"=character(),"date"=as.Date(character()),"status"=character(),"status.next"=character(),"tag"=character(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"foi"=numeric())

# fill data object
for(site in sites)
{
  print(site)
  ## subset data by site
  sub.loc.data<-corrected.locs[which(corrected.locs$Site==site),]
  sub.1.epi.data<-corrected.epi[which(corrected.epi$Site==site),]
  
  ## get ID strings for each plant in sub.loc.data
  ID.strings<-c()
  for(i in 1:dim(sub.loc.data)[1]) 
  {
    if(is.na(sub.loc.data[i,"tag"])) {ID.strings<-c(ID.strings,paste0("X=",sub.loc.data[i,"X"],"Y=",sub.loc.data[i,"Y"],"x=",sub.loc.data[i,"x"],"y=",sub.loc.data[i,"y"]))}
    if(!is.na(sub.loc.data[i,"tag"])) {ID.strings<-c(ID.strings,sub.loc.data[i,"tag"])}
  }
  
  ## for each date for which we will be able to calculate foi:
  for(date.index in (1:length(data.dates[which(sites==site)][[1]])))
  {
    date0<-epi.obs.dates[which(sites==site)][[1]][date.index] ### date of current observation
    date1<-epi.obs.dates[which(sites==site)][[1]][date.index+1] ### date of next observation
    
    sub.2.epi.data<-sub.1.epi.data[which(sub.1.epi.data$Date.First.Observed.Diseased<=date0),] ### epi data up until date0
    sub.3.epi.data<-sub.1.epi.data[which(sub.1.epi.data$Date.First.Observed.Diseased<=date1),] ### epi date up until date1
    
    print(dim(sub.3.epi.data)[1]-dim(sub.2.epi.data)[1])
    ## build set of ID strings for plants that have become infected by date0
    epi.ID.strings<-c()
    for(i in 1:dim(sub.2.epi.data)[1]) 
    {
      if(is.na(sub.2.epi.data[i,"Tag"])) {epi.ID.strings<-c(epi.ID.strings,paste0("X=",sub.2.epi.data[i,"X"],"Y=",sub.2.epi.data[i,"Y"],"x=",sub.2.epi.data[i,"x"],"y=",sub.2.epi.data[i,"y"]))}
      if(!is.na(sub.2.epi.data[i,"Tag"])) {epi.ID.strings<-c(epi.ID.strings,sub.2.epi.data[i,"Tag"])}
    }
    
    ## build set of ID strings for plants that have become infected by date1
    epi.ID.strings.next<-c()
    for(i in 1:dim(sub.3.epi.data)[1]) 
    {
      if(is.na(sub.3.epi.data[i,"Tag"])) {epi.ID.strings.next<-c(epi.ID.strings.next,paste0("X=",sub.3.epi.data[i,"X"],"Y=",sub.3.epi.data[i,"Y"],"x=",sub.3.epi.data[i,"x"],"y=",sub.3.epi.data[i,"y"]))}
      if(!is.na(sub.3.epi.data[i,"Tag"])) {epi.ID.strings.next<-c(epi.ID.strings.next,sub.3.epi.data[i,"Tag"])}
    }
    
    ## for each plant in the transect:
    for(index in 1:dim(sub.loc.data)[1])
    {
      ### get ID string of plant
      if(!is.na(sub.loc.data[index,"tag"])) {ID.string<-sub.loc.data[index,"tag"]} else {ID.string<-paste0("X=",sub.loc.data[index,"X"],"Y=",sub.loc.data[index,"Y"],"x=",sub.loc.data[index,"x"],"y=",sub.loc.data[index,"y"])}
      ### if that ID string is unique
      if(length(which(ID.strings==ID.string))==1)
      {
        status<-ifelse(ID.string %in% epi.ID.strings,1,0) #### count plant as currently infected if its ID string shows up in the set of plants that have been infected by date0
        new.status<-ifelse(ID.string %in% epi.ID.strings.next,1,0) #### count plant as becoming infected by next obs if its ID string shows up in the set of plants that have been infected by date1
        foi<-NA #### calculate foi experienced
        foi.data<-rbind(foi.data,data.frame("site"=site,"date"=date0,"status"=status,"status.next"=new.status,"tag"=sub.loc.data[index,"tag"],"X"=sub.loc.data[index,"X"],"Y"=sub.loc.data[index,"Y"],"x"=sub.loc.data[index,"x"],"y"=sub.loc.data[index,"y"],"foi"=foi)) #### add new entry to data object
      }
      ### if that ID string is not unique
      if(length(which(ID.strings==ID.string))>1)
      {
        n.repeats<-length(which(ID.strings==ID.string)) #### how many times ID string is repeated
        order<-which(which(ID.strings==ID.string)==index) #### figure out which repeat (e.g. 2nd) this plant is
        status<-ifelse(length(which(epi.ID.strings==ID.string))>=order,1,0) #### for plant that is the nth repeat, count as currently infected if at least n matching ID strings show up in epi.ID.strings
        new.status<-ifelse(length(which(epi.ID.strings.next==ID.string))>=order,1,0) #### for plant that is the nth repeat, count as becoming infected by next obs if at least n matching ID strings show up in epi.ID.strings.next
        foi.data<-rbind(foi.data,data.frame("site"=site,"date"=date0,"status"=status,"status.next"=new.status,"tag"=sub.loc.data[index,"tag"],"X"=sub.loc.data[index,"X"],"Y"=sub.loc.data[index,"Y"],"x"=sub.loc.data[index,"x"],"y"=sub.loc.data[index,"y"],"foi"=foi)) #### add new entry to data object
      }
    }
  }
}



