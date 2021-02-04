if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")))
{
  # load model and data from pustule area analysis
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/pustule area data prep.R")
  delta.pustules<-subset(delta.pustules,time<=7)
  pustule.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/pustule.model.RDS")
  
  # load model and data from n pustules analysis
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/n pustules data prep.R")
  delta.n.pustules<-subset(delta.n.pustules,time<=7)
  n.pustules.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/n.pustules.model.RDS")
  
  # load model and data from stem analysis
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/stem data prep.R")
  delta.stems<-subset(delta.stems,time<=7)
  stems.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/stems.model.RDS")
  
  
  # prep enviro data
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  
  # prep data
  
  ## load data
  within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
  
  ## subset to 2020
  within.host<-within.host[which(within.host$Year==2020),]
  
  ## cut out incomplete records
  within.host<-within.host[-intersect(which(within.host$Site=="HM"),which(as.Date(within.host$Date,tryFormats = "%m/%d/%y")>as.Date("2020-07-10"))),]
  
  ## pull out relevant columns
  plants<-within.host[,c("Year","Site","Tag","Date","N.Stems","N.D.Stems","max.height","picture","stem.index","stem.height","percent.tissue.infected","length.tissue.infected","N.pustules.middle")]
  
  ## get rid of data with missing  or NA length.tissue.infected or N.pustules.middle or stem.height
  if(any(is.na(plants$length.tissue.infected))) {plants<-plants[-which(is.na(plants$length.tissue.infected)),]}
  if(any(is.na(plants$length.tissue.infected))) {plants<-plants[-which(plants$length.tissue.infected==""),]}
  if(any(is.na(plants$N.pustules.middle))) {plants<-plants[-which(is.na(plants$N.pustules.middle)),]}
  if(any(is.na(plants$N.pustules.middle))) {plants<-plants[-which(plants$N.pustules.middle==""),]}
  if(any(is.na(plants$stem.height))) {plants<-plants[-which(is.na(plants$stem.height)),]}
  if(any(is.na(plants$stem.height))) {plants<-plants[-which(plants$stem.height==""),]}
  
  ## summarize by plant for each date
  years<-c()
  sites<-c()
  dates<-c()
  tags<-c()
  n.stems<-c()
  n.d.stems<-c()
  max.heights<-c()
  reference.pictures<-c()
  plant.inf.intens<-c()
  
  for (tag in unique(plants$Tag))
  {
    sub.plants.1<-plants[which(plants$Tag==tag),]
    
    for (date in unique(sub.plants.1$Date))
    {
      sub.plants.2<-sub.plants.1[which(sub.plants.1$Date==date),]
      
      if(any(!is.na(sub.plants.2$N.Stems)) & any(!is.na(sub.plants.2$N.D.Stems)))
      {
        ## new values
        new.year<-sub.plants.2[1,"Year"]
        new.site<-sub.plants.2[1,"Site"]
        new.tag<-sub.plants.2[1,"Tag"]
        new.date<-date
        new.n.stems<-sub.plants.2[1,"N.Stems"]
        new.n.d.stems<-sub.plants.2[1,"N.D.Stems"]
        new.max.height<-sub.plants.2[1,"max.height"]
        new.reference.picture<-sub.plants.2[1,"picture"]
        
        ## calculate plant infection intensity
        new.plant.inf.intens<-new.n.d.stems*mean(as.numeric(sub.plants.2$length.tissue.infected)*as.numeric(sub.plants.2$N.pustules.middle))
        
        ## store new values
        years<-c(years,new.year)
        sites<-c(sites,new.site)
        dates<-c(dates,new.date)
        tags<-c(tags,new.tag)
        n.stems<-c(n.stems,new.n.stems)
        n.d.stems<-c(n.d.stems,new.n.d.stems)
        max.heights<-c(max.heights,new.max.height)
        reference.pictures<-c(reference.pictures,new.reference.picture)
        plant.inf.intens<-c(plant.inf.intens,new.plant.inf.intens)
        
      }
    }
  }
  
  plants<-data.frame("Year"=years,"Site"=sites,"Tag"=tags,"Date"=dates,"N.Stems"=n.stems,"N.D.Stems"=n.d.stems,"max.height"=max.heights,"picture"=reference.pictures,"plant.inf.intens"=plant.inf.intens)

  ## correct 0 intensities to .1, logic being that this is a measure of tot infection load, and it shouldn't be less than 1 pustule (coded as .1cm infected tissue, 1 pustule/leaf)
  plants$plant.inf.intens[which(plants$plant.inf.intens<.1)]<-.1
  
  
  ## finish cleaning
  plants$Date<-as.Date(plants$Date,tryFormats = "%m/%d/%y")

  ## make new data object for change in pustule size
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  n.stems<-c()
  n.d.stems<-c()
  max.heights<-c()
  start.plant.inf.intens<-c()
  end.plant.inf.intens<-c()
  days<-c()
  temp.days<-c()
  temp.days.16.22<-c()
  temp.days.7.30<-c()
  dew.point.days<-c()
  temp.dew.point.days<-c()
  temp.16.22.dew.point.days<-c()
  temp.7.30.dew.point.days<-c()
  wetness.days<-c()
  temp.wetness.days<-c()
  temp.16.22.wetness.days<-c()
  temp.7.30.wetness.days<-c()
  tot.rains<-c()
  solar.days<-c()
  wind.speed.days<-c()
  gust.speed.days<-c()
  pred.pustule.diam.growths<-c()
  pred.pustule.num.increases<-c()
  pred.stem.inf.intens.increases<-c()
  
  for (tag in unique(plants$Tag))
  {
    sub.plants.1<-plants[which(plants$Tag==tag),]
    sub.plants.1<-sub.plants.1[order(sub.plants.1$Date),]
    if(dim(sub.plants.1)[1]>=2)
    {
      for(i in 1:(dim(sub.plants.1)[1]-1))
      {
        ### pull reference data
        #### use pictures to get date times for plant, average
        date0.pic<-sub.plants.1[i,"picture"]
        date0<-pustules[which(pustules$picture==date0.pic)[1],"date"]
        date1.pic<-sub.plants.1[i+1,"picture"]
        date1<-pustules[which(pustules$picture==date1.pic)[1],"date"]
        #### if unable to line up pic with date time, use mean date time of all measurments at site on given day
        if(is.na(date0))
        {
          simple.date<-sub.plants.1[i,"Date"]
          sub.pustules<-pustules[which(pustules$site==site),]
          sub.pustules<-sub.pustules[which(as.Date(sub.pustules$date)==as.Date(simple.date,tryFormats = "%m/%d/%y")),]
          alt.date.0.pics<-unique(sub.pustules$picture)
          alt.dates<-pustules[which(pustules$picture %in% alt.date.0.pics),"date"]
          date0<-mean(unique(alt.dates))
        }
        if(is.na(date1))
        {
          simple.date<-sub.plants.1[i+1,"Date"]
          sub.pustules<-pustules[which(pustules$site==site),]
          sub.pustules<-sub.pustules[which(as.Date(sub.pustules$date)==as.Date(simple.date,tryFormats = "%m/%d/%y")),]
          alt.date.1.pics<-unique(sub.pustules$picture)
          alt.dates<-pustules[which(pustules$picture %in% alt.date.1.pics),"date"]
          date1<-mean(unique(alt.dates))
        }
        site<-sub.plants.1[i,"Site"]
        
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
        
        #pull out core predictors
        new.start.plant.inf.intens<-sub.plants.1[i,"plant.inf.intens"]
        new.end.plant.inf.intens<-sub.plants.1[i+1,"plant.inf.intens"]
        delta.days<-as.numeric(date1-date0)
        
        new.n.stems<-sub.plants.1[i,"N.Stems"]
        new.n.d.stems<-sub.plants.1[i,"N.D.Stems"]
        new.max.heights<-sub.plants.1[i,"max.height"]
        
        #predict pustule growth from pustule growth model and enviro conditions
        #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
        pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
        pustule.model.new.time<-1
        obs.time<-delta.days
        pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"time"=pustule.model.new.time,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.dew.point.days"=new.temp.16.22.dew.point.days/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days,"tot.rain"=new.tot.rain/delta.days)
        pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
        
        #predict change in number of pustules from enviro conditions
        #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
        n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
        obs.time<-delta.days
        n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=new.temp.days.16.22/delta.days,"temp.16.22.wetness.days"=new.temp.16.22.wetness.days/delta.days)
        pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
        
        #predict change in stem infection intensity from enviro conditions
        #stem.model.vars<-names(fixef(stems.model))[2:length(names(fixef(stems.model)))]
        stems.model.stem.inf.intens<-1  #predict change for lightly infected stem, arbitrarily pick .01
        obs.time<-delta.days
        stems.model.pred.data<-data.frame("stem.inf.intens"=stems.model.stem.inf.intens,"temp.7.30.wetness.days"=new.temp.7.30.wetness.days/delta.days,"tot.rain"=new.tot.rain/delta.days)
        pred.stem.inf.intens.increase<-predict(stems.model,newdata=stems.model.pred.data,re.form=~0)
        
        #store values
        tags<-c(tags,tag)
        start.plant.inf.intens<-c(start.plant.inf.intens,new.start.plant.inf.intens)
        end.plant.inf.intens<-c(end.plant.inf.intens,new.end.plant.inf.intens)
        n.stems<-c(n.stems,new.n.stems)
        n.d.stems<-c(n.d.stems,new.n.d.stems)
        max.heights<-c(max.heights,new.max.heights)
        days<-c(days,delta.days)
        
        temp.days.16.22<-c(temp.days.16.22,new.temp.days.16.22)
        temp.days.7.30<-c(temp.days.7.30,new.temp.days.7.30)
        temp.days<-c(temp.days,new.temp.days)
        dew.point.days<-c(dew.point.days,new.dew.point.days)
        temp.dew.point.days<-c(temp.dew.point.days,new.temp.dew.point.days)
        temp.16.22.dew.point.days<-c(temp.16.22.dew.point.days,new.temp.16.22.dew.point.days)
        temp.7.30.dew.point.days<-c(temp.7.30.dew.point.days,new.temp.7.30.dew.point.days)
        
        wetness.days<-c(wetness.days,new.wetness.days)
        temp.wetness.days<-c(temp.wetness.days,new.temp.wetness.days)
        temp.16.22.wetness.days<-c(temp.16.22.wetness.days,new.temp.16.22.wetness.days)
        temp.7.30.wetness.days<-c(temp.7.30.wetness.days,new.temp.7.30.wetness.days)
        tot.rains<-c(tot.rains,new.tot.rain)
        solar.days<-c(solar.days,new.solar.days)
        wind.speed.days<-c(wind.speed.days,new.wind.speed.days)
        gust.speed.days<-c(gust.speed.days,new.gust.speed.days)
        
        pred.pustule.diam.growths<-c(pred.pustule.diam.growths,pred.pustule.diam.growth)
        pred.pustule.num.increases<-c(pred.pustule.num.increases,pred.pustule.num.increase)
        pred.stem.inf.intens.increases<-c(pred.stem.inf.intens.increases,pred.stem.inf.intens.increase)
      } 
    }
  }
  
  delta.plants<-data.frame(tag=factor(tags),time=days,N.stems=n.stems,N.D.Stems=n.d.stems,max.height=max.heights,plant.inf.intens=start.plant.inf.intens,plant.inf.intens.next=end.plant.inf.intens,
                          temp.days=temp.days,temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,
                          dew.point.days=dew.point.days,temp.dew.point.days=temp.dew.point.days,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,
                          wetness.days=wetness.days,temp.wetness.days=temp.wetness.days,temp.16.22.wetness.days=temp.16.22.wetness.days,temp.7.30.wetness.days=temp.7.30.wetness.days,
                          tot.rain=tot.rains,solar.days=solar.days,wind.speed.days=wind.speed.days,gust.speed.days=gust.speed.days,pred.pustule.diam.growth=pred.pustule.diam.growths,pred.pustule.num.increase=pred.pustule.num.increases,pred.stem.inf.intens.increase=pred.stem.inf.intens.increases)
  
  saveRDS(plants,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
  saveRDS(delta.plants,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")
}

plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
delta.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")


