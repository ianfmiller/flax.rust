if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/foi.data.RDS")))
{
  # load datasets and functions
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")
  corrected.locs<-corrected.locs[-which(corrected.locs$tag %in% c(911, 912, 913, 943, 934, 935, 942, 944, 945, 984, 985, 986)),] # Trim out plants in unsurveyed region of CC. They will still be included as sources of transmission, but should not be analyzed as targets of transmission.
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant inf intens change funcs.R")
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
  pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
  n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
  plant.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plants.model.RDS")
  plant.inf.intens<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
  corrected.heights<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/corrected.plant.heights.RDS")
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
  
  # function for calculating foi
  foi.func<-function(site,date0,date1,epi.data=corrected.epi,xcord,ycord)
  {
    ## get wind data from site and dates
    wind.data<-all.weath[which(all.weath$site==site),]
    wind.data<-wind.data[which(wind.data$date>as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")),]
    wind.data<-wind.data[which(wind.data$date<=as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")),]
    
    ## get data about sources of spores
    source.data<-epi.data[intersect(which(epi.data$Site==site),which(epi.data$Date.First.Observed.Diseased<=date0)),]
    
    ## substitute in corrected heights as average of height at date0 and date1
    
    ### loop to get corrected heights of source plants
    corrected.heights.vec<-c()
    for (k in 1:nrow(source.data))
    {
      #### tagged source plants
      if(!is.na(source.data[k,"Tag"])) 
      {
        height0<-as.numeric(corrected.heights[intersect(which(corrected.heights$tag==source.data[k,"Tag"]),which(corrected.heights$Date==date0)),"height.cm"])
        height1<-as.numeric(corrected.heights[intersect(which(corrected.heights$tag==source.data[k,"Tag"]),which(corrected.heights$Date==date1)),"height.cm"])
        corrected.height<-mean(height0,height1)
      }
      #### untagged source plants
      if(is.na(source.data[k,"Tag"])) 
      {
        index0<-Reduce(intersect,list(which(corrected.heights$X==source.data[k,"X"]),which(corrected.heights$Y==source.data[k,"Y"]),which(corrected.heights$x==source.data[k,"x"]),which(corrected.heights$y==source.data[k,"y"]),which(corrected.heights$Date==date0)))
        index1<-Reduce(intersect,list(which(corrected.heights$X==source.data[k,"X"]),which(corrected.heights$Y==source.data[k,"Y"]),which(corrected.heights$x==source.data[k,"x"]),which(corrected.heights$y==source.data[k,"y"]),which(corrected.heights$Date==date1)))
        height0<-as.numeric(corrected.heights[index0,"height.cm"])
        height1<-as.numeric(corrected.heights[index1,"height.cm"])
        corrected.height<-mean(height0,height1)
      }
      #### make vec
      corrected.heights.vec<-c(corrected.heights.vec,corrected.height)
    }
    
    ### add in courrected hieght data
    source.data<-data.frame(source.data,"corrected.height"=corrected.heights.vec)
    
    ## loop to calculate summed foi from all sources
    tot.dep<-c()
    for(i in 1:dim(source.data)[1])
    {
      ### get cords of source
      sourceX<-source.data[i,"X"]+source.data[i,"x"]
      sourceY<-source.data[i,"Y"]+source.data[i,"y"]
      ### get q val
      #### if the source plant is tagged
      if(!is.na(source.data[i,"Tag"]))
      {
        if(source.data[i,"Tag"]==15) {q<-0; half.height=10} #### set q = 0 for tag 15, half.height doesn't matter
        else{
          #### set q to 0.1 if the plant is a seedling
          if(!source.data[i,"Tag"] %in% plant.inf.intens$Tag & source.data[i,"corrected.height"]<=5) {q<-.1; if(is.na(source.data[i,"corrected.height"])) {half.height<-2.5} else {half.height<-source.data[i,"corrected.height"]*.5}}
          #### extract q from data
          else {
            plant.inf.intens.index<-intersect(which(plant.inf.intens$Date==date0),which(plant.inf.intens$Tag==source.data[i,"Tag"]))
            ##### if there's an actual measurment use that
            if(length(plant.inf.intens.index)==1) {q<-plant.inf.intens[plant.inf.intens.index,"plant.inf.intens"]; half.height<-source.data[i,"corrected.height"]*.5}
            ##### if not, forecast or hindcast from closest observation
            if(length(plant.inf.intens.index)==0)
            {
              obs.dates<-as.Date(plant.inf.intens[which(plant.inf.intens$Tag==source.data[i,"Tag"]),"Date"]) ###### get dates that have plant inf intens observations
              if(as.Date(date0) %in% obs.dates) {obs.dates<-obs.dates[-which(obs.dates==as.Date(date0))]}
              date.diffs<-as.Date(date0)-obs.dates ### calculate time differences
              closest.date.index<-which(abs(date.diffs)==min(abs(date.diffs))) ###### find closest observation
              if(length(closest.date.index)>1) {closest.date.index<-closest.date.index[which(date.diffs[closest.date.index]<0)]} ###### if there's a tie, default to future observation
              closest.date<-obs.dates[closest.date.index]
              plant.inf.intens.index<-intersect(which(plant.inf.intens$Date==closest.date),which(plant.inf.intens$Tag==source.data[i,"Tag"]))
              if(closest.date>date0) ###### hindcast
              {
                q<-predict.plant.inf.intens.last(plant.inf.intens.next=plant.inf.intens[plant.inf.intens.index,"plant.inf.intens"],site=site,date0=as.POSIXct(paste0(date0," 12:00:00")),date1=as.POSIXct(paste0(closest.date," 12:00:00")))
                half.height<-plant.inf.intens[plant.inf.intens.index,"max.height"]*.5
              }
              if(closest.date<date0) ###### forecast
              {
                q<-predict.plant.inf.intens(plant.inf.intens.last=plant.inf.intens[plant.inf.intens.index,"plant.inf.intens"],site=site,date0=as.POSIXct(paste0(closest.date," 12:00:00")),date1=as.POSIXct(paste0(date0," 12:00:00")))
                half.height<-plant.inf.intens[plant.inf.intens.index,"max.height"]*.5
              }
            }
          } 
        }
      }
      #### if the source plant is not tagged
      if(is.na(source.data[i,"Tag"]))
      {
        #### set q to 0.1 if the plant is a seedling
        if(source.data[i,"corrected.height"]<=5)
        {
          q<-.1
          half.height<-2.5
        } else {print('error--missing plant inf intens data'); print(i)}
      }
      
      ### calculate foi from source plant
      tot.dep<-c(tot.dep,predict.kernel.tilted.plume(q=q,H=half.height,k=5.803369e-07,alphaz=1.596314e-01,Ws=1.100707e+00,xtarget=xcord-sourceX,ytarget=ycord-sourceY,wind.data=wind.data)) 
    }
    sum(tot.dep)
  }
  
  # convenience pars
  sites<-c("CC","BT","GM","HM")
  epi.obs.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20","2020-07-27"),"BT"=c("2020-06-24","2020-07-01","2020-07-08"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09","2020-07-15"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09","2020-07-15"))    
  data.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20"),"BT"=c("2020-06-24","2020-07-01"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09"))
  
  # build empty data object for foi vs outcome data
  foi.data<-data.frame("site"=character(),"date"=as.Date(character()),"status"=character(),"status.next"=character(),"tag"=character(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"first.obs.height"=numeric(),"foi"=numeric(),
                       "temp.days"=numeric(),"temp.days.16.22"=numeric(),"temp.days.7.20"=numeric(),"dew.point.days"=numeric(),"wetness.days"=numeric(),"tot.rain"=numeric(),"solar.days"=numeric(),"wind.speed.days"=numeric(),
                       "gust.speed.days"=numeric(),"temp.16.22.dew.point.days"=numeric(),"temp.7.30.dew.point.days"=numeric(),"temp.wetness.days"=numeric(),"temp.16.22.wetness.days"=numeric(),"temp.7.30.wetness.days"=numeric(),
                       "pred.pustule.diam.growth"=numeric(),"pred.pustule.num.increase"=numeric(),"pred.plant.inf.intens.increase"=numeric())
  
  # fill data object
  for(site in sites)
  {
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
      new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
      new.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #time (in days) during which temp between 16 and 22 celsius
      new.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #time (in days) during which temp between 7 and 30 celsius
      new.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T) #Dew point days
      
      #calculate weather metrics
      new.wetness.days<-sum(weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      new.tot.rain<-sum(weath.sub$rain,na.rm=T)
      new.solar.days<-sum(weath.sub$solar.radiation*weath.sub$interval.length,na.rm = T)
      new.wind.speed.days<-sum(weath.sub$wind.speed*weath.sub$interval.length,na.rm = T)
      new.gust.speed.days<-sum(weath.sub$wind.direction*weath.sub$interval.length,na.rm = T)
      
      #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
      new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
      new.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
      new.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)
      new.temp.wetness.days<-sum(weath.sub$temp*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      new.temp.16.22.wetness.days<-sum(weath.sub$temp.16.22*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      new.temp.7.30.wetness.days<-sum(weath.sub$temp.7.30*weath.sub$wetness*weath.sub$interval.length,na.rm = T)
      
      #predict pustule growth from pustule growth model and enviro conditions
      #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
      pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
      obs.time<-delta.days
      pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.days.16.22"=new.temp.days.16.22/delta.days,"dew.point.days"=new.dew.point.days/delta.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days/delta.days,"temp.wetness.days"=new.temp.wetness.days/delta.days,"tot.rain"=new.tot.rain/delta.days)
      pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
      
      #predict change in number of pustules from enviro conditions
      #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
      n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
      obs.time<-delta.days
      n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
      pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
      
      #predict change in plant.inf.intensity from enviro conditions
      #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
      plant.model.new.plant.inf.intens<-.1
      obs.time<-delta.days
      plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"dew.point.days"=new.dew.point.days/delta.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days/delta.days,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
      pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
      
      sub.2.epi.data<-sub.1.epi.data[which(sub.1.epi.data$Date.First.Observed.Diseased<=date0),] ### epi data up until date0
      sub.3.epi.data<-sub.1.epi.data[which(sub.1.epi.data$Date.First.Observed.Diseased<=date1),] ### epi date up until date1
      
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
          foi<-foi.func(site=site,date0=date0,date1=date1,epi.data=corrected.epi,xcord=sub.loc.data[index,"X"]+sub.loc.data[index,"x"],ycord=sub.loc.data[index,"Y"]+sub.loc.data[index,"y"]) #### calculate foi experienced
          foi.data<-rbind(foi.data,data.frame("site"=site,"date"=date0,"status"=status,"status.next"=new.status,"tag"=sub.loc.data[index,"tag"],"X"=sub.loc.data[index,"X"],"Y"=sub.loc.data[index,"Y"],"x"=sub.loc.data[index,"x"],"y"=sub.loc.data[index,"y"],"first.obs.height"=sub.loc.data[index,"height.cm"],"foi"=foi,
                                              "temp.days"=new.temp.days,"temp.days.16.22"=new.temp.days.16.22,"temp.days.7.30"=new.temp.days.7.30,"dew.point.days"=new.dew.point.days,"wetness.days"=new.wetness.days,"tot.rain"=new.tot.rain,"solar.days"=new.solar.days,"wind.speed.days"=new.wind.speed.days,
                                              "gust.speed.days"=new.gust.speed.days,"temp.dew.point.days"=new.temp.dew.point.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days,"temp.wetness.days"=new.temp.wetness.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days,"temp.7.30.wetness.days"=new.temp.7.30.wetness.days,
                                              "pred.pustule.diam.growth"=pred.pustule.diam.growth,"pred.pustule.num.increase"=pred.pustule.num.increase,"pred.plant.inf.intens.increase"=pred.plant.inf.intens.increase)) #### add new entry to data object
        }
        ### if that ID string is not unique
        if(length(which(ID.strings==ID.string))>1)
        {
          n.repeats<-length(which(ID.strings==ID.string)) #### how many times ID string is repeated
          order<-which(which(ID.strings==ID.string)==index) #### figure out which repeat (e.g. 2nd) this plant is
          status<-ifelse(length(which(epi.ID.strings==ID.string))>=order,1,0) #### for plant that is the nth repeat, count as currently infected if at least n matching ID strings show up in epi.ID.strings
          new.status<-ifelse(length(which(epi.ID.strings.next==ID.string))>=order,1,0) #### for plant that is the nth repeat, count as becoming infected by next obs if at least n matching ID strings show up in epi.ID.strings.next
          foi<-foi.func(site=site,date0=date0,date1=date1,epi.data=corrected.epi,xcord=sub.loc.data[index,"X"]+sub.loc.data[index,"x"],ycord=sub.loc.data[index,"Y"]+sub.loc.data[index,"y"]) #### calculate foi experienced
          foi.data<-rbind(foi.data,data.frame("site"=site,"date"=date0,"status"=status,"status.next"=new.status,"tag"=sub.loc.data[index,"tag"],"X"=sub.loc.data[index,"X"],"Y"=sub.loc.data[index,"Y"],"x"=sub.loc.data[index,"x"],"y"=sub.loc.data[index,"y"],"first.obs.height"=sub.loc.data[index,"height.cm"],"foi"=foi,
                                              "temp.days"=new.temp.days,"temp.days.16.22"=new.temp.days.16.22,"temp.days.7.30"=new.temp.days.7.30,"dew.point.days"=new.dew.point.days,"wetness.days"=new.wetness.days,"tot.rain"=new.tot.rain,"solar.days"=new.solar.days,"wind.speed.days"=new.wind.speed.days,
                                              "gust.speed.days"=new.gust.speed.days,"temp.dew.point.days"=new.temp.dew.point.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days,"temp.7.30.dew.point.days"=new.temp.7.30.dew.point.days,"temp.wetness.days"=new.temp.wetness.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days,"temp.7.30.wetness.days"=new.temp.7.30.wetness.days,
                                              "pred.pustule.diam.growth"=pred.pustule.diam.growth,"pred.pustule.num.increase"=pred.pustule.num.increase,"pred.plant.inf.intens.increase"=pred.plant.inf.intens.increase)) #### add new entry to data object
        }
        print(paste0("site = ",site," date = ",date0," ",index,"/",dim(sub.loc.data)[1]," complete"))
      }
    }
  }
  
  ## rescale predictors to be time independent
  foi.data<-data.frame(foi.data,"delta.days"=NA)
  
  for (i in 1:dim(foi.data)[1])
  {
    date0<-foi.data[i,"date"]
    site<-foi.data[i,"site"]
    site.date.set<-epi.obs.dates[which(sites==site)][[1]]
    date1<-site.date.set[which(site.date.set==date0)+1]
    delta.days<-as.numeric(as.Date(date1)-as.Date(date0))
    foi.data[i,"delta.days"]<-delta.days
  }
  
  foi.data[,c("temp.days","temp.days.16.22","temp.days.7.30","dew.point.days","wetness.days","tot.rain","solar.days","wind.speed.days","gust.speed.days","temp.dew.point.days","temp.16.22.dew.point.days","temp.7.30.dew.point.days","temp.wetness.days","temp.16.22.wetness.days","temp.7.30.wetness.days")]<-foi.data[,c("temp.days","temp.days.16.22","temp.days.7.30","dew.point.days","wetness.days","tot.rain","solar.days","wind.speed.days","gust.speed.days","temp.16.22.dew.point.days","temp.dew.point.days","temp.7.30.dew.point.days","temp.wetness.days","temp.16.22.wetness.days","temp.7.30.wetness.days")]/foi.data[,"delta.days"]
  colnames(foi.data)[which(colnames(foi.data) %in% c("temp.days","temp.days.16.22","temp.days.7.30","dew.point.days","wetness.days","tot.rain","solar.days","wind.speed.days","gust.speed.days","temp.dew.point.days","temp.16.22.dew.point.days","temp.7.30.dew.point.days","temp.wetness.days","temp.16.22.wetness.days","temp.7.30.wetness.days"))]<-paste0("mean.",c("temp.days","temp.days.16.22","temp.days.7.30","dew.point.days","wetness.days","tot.rain","solar.days","wind.speed.days","gust.speed.days","temp.dew.point.days","temp.16.22.dew.point.days","temp.7.30.dew.point.days","temp.wetness.days","temp.16.22.wetness.days","temp.7.30.wetness.days"))
  
  saveRDS(foi.data,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/foi.data.RDS")
}

foi.data<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/foi.data.RDS")
