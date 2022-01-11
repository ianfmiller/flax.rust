if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/transmission.data.RDS")))
{
  set.seed(569087)
  # load data sets and functions
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc data set building.R")
  corrected.locs<-corrected.locs[-which(corrected.locs$tag %in% c(911, 912, 913, 943, 934, 935, 942, 944, 945, 984, 985, 986)),] # Trim out plants in unsurveyed region of CC. They will still be included as sources of transmission, but should not be analyzed as targets of transmission.
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant height data set building.R")
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant inf intens change funcs.R")
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
  
  pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
  n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
  diseased.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/diseased.focal.plants.RDS")
  corrected.heights<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.plant.heights.RDS")
  corrected.epi<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.epi.RDS")
  mean.starting.inf.intens<-1
  # Time periods for fitting glm of infection~transmission to data:
  #CC: 6/22(first data)->7/27(last prediction)
  #BT: 6/24(first data)->7/8(last prediction)
  #GM: 6/16(first data)->7/15(last prediction)
  #BT: 6/25(first data)->7/28(last prediction)
  ## These periods were selected such that the plant infection intensity metric was recorded for every infected plant in the transect for the observation prior to the observed changes in infection status
  ## Exceptions wer emade for:
  ### Infected seedling, which were assumed to have 0.1 starting inf intensity
  ### Miissing plant inf intensity measurments that could be fore/hindcasted from existing observations using using plant inf intens model. The window for fore/hindcasting was limited to +/- 1 week frpm existing data
  ### Tag #15 in "CC" which was assumed to have 0 inf intensity as indicated by an observation on 6/22
  
  # function for calculating transmission
  
  total.spore.deposition.func<-function(site,date0,date1,epi.data=corrected.epi,xcord,ycord)
  {
    ## get wind data from site and dates
    
    wind.data<-all.weath[which(all.weath$site==site),]
    wind.data<-wind.data[which(wind.data$date>date0),]
    wind.data<-wind.data[which(wind.data$date<=date1),]
    
    ## get data about sources of spores
    source.data<-epi.data[intersect(which(epi.data$Site==site),which(epi.data$Date.First.Observed.Diseased<=as.Date(date0))),]
    
    ## substitute in corrected heights as average of height at date0 and date1
    
    ### loop to get corrected heights of source plants
    corrected.heights.vec<-c()
    for (k in 1:nrow(source.data))
    {
      #### tagged source plants
      if(!is.na(source.data[k,"Tag"])) 
      {
        height0<-as.numeric(corrected.heights[intersect(which(corrected.heights$tag==source.data[k,"Tag"]),which(corrected.heights$Date==as.Date(date0))),"height.cm"])
        height1<-as.numeric(corrected.heights[intersect(which(corrected.heights$tag==source.data[k,"Tag"]),which(corrected.heights$Date==as.Date(date1))),"height.cm"])
        corrected.height<-mean(height0,height1)
      }
      #### untagged source plants
      if(is.na(source.data[k,"Tag"])) 
      {
        index0<-Reduce(intersect,list(which(corrected.heights$X==source.data[k,"X"]),which(corrected.heights$Y==source.data[k,"Y"]),which(corrected.heights$x==source.data[k,"x"]),which(corrected.heights$y==source.data[k,"y"]),which(corrected.heights$Date==as.Date(date0))))
        index1<-Reduce(intersect,list(which(corrected.heights$X==source.data[k,"X"]),which(corrected.heights$Y==source.data[k,"Y"]),which(corrected.heights$x==source.data[k,"x"]),which(corrected.heights$y==source.data[k,"y"]),which(corrected.heights$Date==as.Date(date1))))
        height0<-as.numeric(corrected.heights[index0,"height.cm"])
        height1<-as.numeric(corrected.heights[index1,"height.cm"])
        corrected.height<-mean(height0,height1)
      }
      #### make vec
      corrected.heights.vec<-c(corrected.heights.vec,corrected.height)
    }
    
    ### add in courrected hieght data
    source.data<-data.frame(source.data,"corrected.height"=corrected.heights.vec)
    
    ## loop to calculate summed transmission from all sources
    tot.dep<-c()
    for(i in 1:dim(source.data)[1])
    {
      ### get cords of source
      sourceX<-source.data[i,"X"]+source.data[i,"x"]
      sourceY<-source.data[i,"Y"]+source.data[i,"y"]
      ### get I val
      #### if the source plant is tagged
      if(!is.na(source.data[i,"Tag"]))
      {
        if(source.data[i,"Tag"]==15 | (source.data[i,"Tag"]==928)) {I<-0; half.height=10} #### set  I = 0 for tags 15, 928, half.height doesn't matter. 928 was marked as diseased, but no pustules were observed until 7/21/20, which is outside of date range for this data set
        else{
          #### set I to 0.1 if the plant is a seedling
          if(!source.data[i,"Tag"] %in% diseased.focal.plants$Tag & (source.data[i,"corrected.height"]<=5 | grepl("seedling",source.data[i,"notes"]))) {I<-.1; if(is.na(source.data[i,"corrected.height"])) {half.height<-2.5} else {half.height<-source.data[i,"corrected.height"]*.5}}
          #### extract I from data
          else {
            diseased.focal.plants.index<-intersect(which(diseased.focal.plants$Date==as.Date(date0)),which(diseased.focal.plants$Tag==source.data[i,"Tag"]))
            ##### if there's an actual measurment use that
            if(length(diseased.focal.plants.index)==1) {I<-diseased.focal.plants[diseased.focal.plants.index,"inf.intens"]; half.height<-source.data[i,"corrected.height"]*.5}
            ##### if not, forecast or hindcast from closest observation
            if(length(diseased.focal.plants.index)==0)
            {
              obs.dates<-as.Date(diseased.focal.plants[which(diseased.focal.plants$Tag==source.data[i,"Tag"]),"Date"]) ###### get dates that have plant inf intens observations

              ###### make sure obs dates are included in temperature data
              if(site=="CC" & any(obs.dates<=as.Date("2020-06-21"))) {obs.dates<-obs.dates[-which(obs.dates<=as.Date("2020-06-21"))]}
              if(site=="BT" & any(obs.dates<=as.Date("2020-06-19"))) {obs.dates<-obs.dates[-which(obs.dates<=as.Date("2020-06-19"))]}
              if(site=="GM" & any(obs.dates<=as.Date("2020-06-22"))) {obs.dates<-obs.dates[-which(obs.dates<=as.Date("2020-06-22"))]}
              if(site=="HM" & any(obs.dates<=as.Date("2020-06-24"))) {obs.dates<-obs.dates[-which(obs.dates<=as.Date("2020-06-24"))]}
              
              if(site=="CC" & any(obs.dates>=as.Date("2020-07-28"))) {obs.dates<-obs.dates[-which(obs.dates>=as.Date("2020-07-28"))]}
              if(site=="BT" & any(obs.dates>=as.Date("2020-07-30"))) {obs.dates<-obs.dates[-which(obs.dates>=as.Date("2020-07-30"))]}
              if(site=="GM" & any(obs.dates>=as.Date("2020-07-29"))) {obs.dates<-obs.dates[-which(obs.dates>=as.Date("2020-07-29"))]}
              if(site=="HM" & any(obs.dates>=as.Date("2020-07-11"))) {obs.dates<-obs.dates[-which(obs.dates>=as.Date("2020-07-11"))]}
              
              date.diffs<-obs.dates-as.Date(date0) ### calculate time differences
              closest.date.index<-which(abs(date.diffs)==min(abs(date.diffs))) ###### find closest observation
              if(length(closest.date.index)>1) {closest.date.index<-closest.date.index[which(date.diffs[closest.date.index]<0)]} ###### if there's a tie, default to future observation
              closest.date<-obs.dates[closest.date.index]
              diseased.focal.plants.index<-intersect(which(diseased.focal.plants$Date==closest.date),which(diseased.focal.plants$Tag==source.data[i,"Tag"]))
              
              model.inf.intens<-diseased.focal.plants[diseased.focal.plants.index,"inf.intens"]
              
              if (!length(model.inf.intens)>1 & source.data[i,"Tag"]==928) {I<-0; half.height<-0} #### plant #928 was suspected to be diseased but no pustuels were actually observed until 7/21
              if (length(model.inf.intens)<1 & source.data[i,"Tag"]!=928) {predict.inf.intens(mean.starting.inf.intens,source.data[i,"corrected.height"],site,as.POSIXct(paste0(source.data[i,"Date.First.Observed.Diseased"]," 12:00:00"),tz="UTC"),date0)}
              if (length(model.inf.intens)>1) {
                if(closest.date>as.Date(date0)) ###### hindcast
                {
                  I<-predict.inf.intens.last(inf.intens.next=model.inf.intens,max.height.last=source.data[i,"corrected.height"],site=site,date0=date0,date1=as.POSIXct(paste0(closest.date," 12:00:00"),tz="UTC"))
                  half.height<-diseased.focal.plants[diseased.focal.plants.index,"max.height"]*.5
                }
                if(closest.date<as.Date(date0)) ###### forecast
                {
                  I<-predict.inf.intens(inf.intens.last=model.inf.intens,max.height.last=source.data[i,"corrected.height"],site=site,date0=as.POSIXct(paste0(closest.date," 12:00:00"),tz="UTC"),date1=date1)
                  half.height<-diseased.focal.plants[diseased.focal.plants.index,"max.height"]*.5
                }
              }
            }
          } 
        }
      }
      #### if the source plant is not tagged
      if(is.na(source.data[i,"Tag"]))
      {
        #### set I to 0.1 if the plant is a seedling
        if(source.data[i,"corrected.height"]<=5 | source.data[i,"notes"]=="seedling")
        {
          I<-.1
          half.height<-2.5
        } else {print('error--missing plant inf intens data'); print(i)}
      }
      
      ### calculate transmission from source plant
      tot.dep<-c(tot.dep,predict.kernel.tilted.plume(I=I,H=half.height/100,k=5.739690e-07,Ws=4.451030e-02,A=7.777373e-02,xtarget=xcord-sourceX,ytarget=ycord-sourceY,wind.data=wind.data)) 
    }
    sum(tot.dep)
  }
  
  # convenience pars
  sites<-c("CC","BT","GM","HM")
  epi.obs.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20","2020-07-27"),"BT"=c("2020-06-24","2020-07-01","2020-07-08"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09","2020-07-15"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09","2020-07-15"))    
  data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20"),"BT"=c("2020-06-24","2020-07-01"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09"))
  
  # build empty data object for transmission vs outcome data
  
  transmission.data<-data.frame("site"=character(),"date"=as.Date(character()),"status"=character(),"status.next"=character(),"tag"=character(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"time"=numeric(),"height.cm"=numeric(),"tot.spore.deposition"=numeric(),
                       "mean.temp"=numeric(),"max.temp"=numeric(),"min.temp"=numeric(),"mean.abs.hum"=numeric(),"max.abs.hum"=numeric(),"min.abs.hum"=numeric(),
                       "mean.daily.rain"=numeric(),"mean.solar"=numeric(),"mean.wetness"=numeric())

  
  # fill data object
  for(site in sites)
  {
    ## subset data by site
    sub.loc.data<-corrected.locs[which(corrected.locs$Site==site),]
    sub.loc.data<-sub.loc.data[which(sub.loc.data$added==F),] ### restrict analysis to plants who's location was recorded at beginning of field season so as to not introduce bias towards infection
    sub.1.epi.data<-corrected.epi[which(corrected.epi$Site==site),]
    
    ## get ID strings for each plant in sub.loc.data
    ID.strings<-c()
    for(i in 1:dim(sub.loc.data)[1]) 
    {
      if(is.na(sub.loc.data[i,"tag"])) {ID.strings<-c(ID.strings,paste0("X=",sub.loc.data[i,"X"],"Y=",sub.loc.data[i,"Y"],"x=",sub.loc.data[i,"x"],"y=",sub.loc.data[i,"y"]))}
      if(!is.na(sub.loc.data[i,"tag"])) {ID.strings<-c(ID.strings,sub.loc.data[i,"tag"])}
    }
    
    ## for each date for which we will be able to calculate total spore deposition:
    for(date.index in (1:length(data.dates[which(sites==site)][[1]])))
    {
      date0<-epi.obs.dates[which(sites==site)][[1]][date.index] ### date of current observation
      date0<-as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")
      date1<-epi.obs.dates[which(sites==site)][[1]][date.index+1] ### date of next observation
      date1<-as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")
      
      delta.days<-as.numeric(as.Date(date1)-as.Date(date0))
      ## get weather/enviro predictors
      
      ### subset temp data to relevant window
      temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #pull out temp data for site
      temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
      temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
      temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
      temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
      
      #### subset weather data to relevant window
      weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
      weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
      weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
      weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
      
      #### calculate environmental variable metrics
      new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
      new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
      new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
      
      abs.hum<-13.24732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+temp.rh.sub$temp.c)
      new.mean.abs.hum<-mean(abs.hum,na.rm=T) 
      
      new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)*(12*24)

      sub.2.epi.data<-sub.1.epi.data[which(sub.1.epi.data$Date.First.Observed.Diseased<=as.Date(date0)),] ### epi data up until date0
      sub.3.epi.data<-sub.1.epi.data[which(sub.1.epi.data$Date.First.Observed.Diseased<=as.Date(date1)),] ### epi date up until date1
      
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
        if(!
          any(
            c(
              (!(site=="CC")),
              (sub.loc.data[index,"Y"] %in% c(0:8,18,19)),
              all(sub.loc.data[index,"Y"] %in% c(15,17),sub.loc.data[index,"X"] %in% 7:9)
            )
          )
        ) {print(paste0("site = ",site," date = ",as.Date(date0)," ",index,"/",dim(sub.loc.data)[1]," complete"));next}
        ### get ID string of plant
        if(!is.na(sub.loc.data[index,"tag"])) {ID.string<-sub.loc.data[index,"tag"]} else {ID.string<-paste0("X=",sub.loc.data[index,"X"],"Y=",sub.loc.data[index,"Y"],"x=",sub.loc.data[index,"x"],"y=",sub.loc.data[index,"y"])}
        ### if that ID string is unique
        if(length(which(ID.strings==ID.string))==1)
        {
          if(is.na(sub.loc.data[index,"tag"]))
          {
            target.index<-Reduce(intersect,list(which(corrected.heights$X==sub.loc.data[index,"X"]),which(corrected.heights$Y==sub.loc.data[index,"Y"]),which(corrected.heights$x==sub.loc.data[index,"x"]),which(corrected.heights$y==sub.loc.data[index,"y"]),which(corrected.heights$Date==as.Date(date0))))
            target.height<-as.numeric(corrected.heights[target.index,"height.cm"])
          }
          if(!is.na(sub.loc.data[index,"tag"]))
          {
            target.height<-as.numeric(corrected.heights[intersect(which(corrected.heights$tag==sub.loc.data[index,"tag"]),which(corrected.heights$Date==as.Date(date0))),"height.cm"])
          }
          status<-ifelse(ID.string %in% epi.ID.strings,1,0) #### count plant as currently infected if its ID string shows up in the set of plants that have been infected by date0
          new.status<-ifelse(ID.string %in% epi.ID.strings.next,1,0) #### count plant as becoming infected by next obs if its ID string shows up in the set of plants that have been infected by date1
          tot.spore.deposition<-total.spore.deposition.func(site=site,date0=date0,date1=date1,epi.data=corrected.epi,xcord=sub.loc.data[index,"X"]+sub.loc.data[index,"x"],ycord=sub.loc.data[index,"Y"]+sub.loc.data[index,"y"]) #### calculate spore deposition experienced
          transmission.data<-rbind(transmission.data,data.frame("site"=site,"date"=as.Date(date0),"status"=status,"status.next"=new.status,"tag"=sub.loc.data[index,"tag"],"X"=sub.loc.data[index,"X"],"Y"=sub.loc.data[index,"Y"],"x"=sub.loc.data[index,"x"],"y"=sub.loc.data[index,"y"],"time"=delta.days,"height.cm"=target.height,"tot.spore.deposition"=tot.spore.deposition,
                                              "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain)) #### add new entry to data object
        }
        ### if that ID string is not unique
        if(length(which(ID.strings==ID.string))>1)
        {
          if(is.na(sub.loc.data[index,"tag"]))
          {
            target.index<-Reduce(intersect,list(which(corrected.heights$X==sub.loc.data[index,"X"]),which(corrected.heights$Y==sub.loc.data[index,"Y"]),which(corrected.heights$x==sub.loc.data[index,"x"]),which(corrected.heights$y==sub.loc.data[index,"y"]),which(corrected.heights$Date==as.Date(date0))))
            targt.height<-as.numeric(corrected.heights[target.index,"height.cm"])
          }
          if(!is.na(sub.loc.data[index,"tag"]))
          {
            target.height<-as.numeric(corrected.heights[intersect(which(corrected.heights$tag==sub.loc.data[index,"tag"]),which(corrected.heights$Date==as.Date(date0))),"height.cm"])
          }
          n.repeats<-length(which(ID.strings==ID.string)) #### how many times ID string is repeated
          order<-which(which(ID.strings==ID.string)==index) #### figure out which repeat (e.g. 2nd) this plant is
          status<-ifelse(length(which(epi.ID.strings==ID.string))>=order,1,0) #### for plant that is the nth repeat, count as currently infected if at least n matching ID strings show up in epi.ID.strings
          new.status<-ifelse(length(which(epi.ID.strings.next==ID.string))>=order,1,0) #### for plant that is the nth repeat, count as becoming infected by next obs if at least n matching ID strings show up in epi.ID.strings.next
          tot.spore.deposition<-total.spore.deposition.func(site=site,date0=date0,date1=date1,epi.data=corrected.epi,xcord=sub.loc.data[index,"X"]+sub.loc.data[index,"x"],ycord=sub.loc.data[index,"Y"]+sub.loc.data[index,"y"]) #### calculate spore deposition experienced
          transmission.data<-rbind(transmission.data,data.frame("site"=site,"date"=as.Date(date0),"status"=status,"status.next"=new.status,"tag"=sub.loc.data[index,"tag"],"X"=sub.loc.data[index,"X"],"Y"=sub.loc.data[index,"Y"],"x"=sub.loc.data[index,"x"],"y"=sub.loc.data[index,"y"],"time"=delta.days,"height.cm"=targt.height,"tot.spore.deposition"=tot.spore.deposition,
                                              "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"mean.daily.rain"=new.mean.daily.rain)) #### add new entry to data object
        }
        print(paste0("site = ",site," date = ",as.Date(date0)," ",index,"/",dim(sub.loc.data)[1]," complete"))
      }
    }
  }

  
  saveRDS(transmission.data,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/transmission.data.RDS")
}

transmission.data<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/transmission.data.RDS")

