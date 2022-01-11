if(!(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.infection.intensity.RDS")) | !(file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/infection.intensity.RDS")))
{

  # prep enviro data
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  
  # generate infection intensity data
  
  ## load raw data
  within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
  
  ## subset to 2020
  within.host<-within.host[which(within.host$Year==2020),]
  
  ## pull out relevant columns
  focal.dis.plants<-within.host[,c("Year","Site","Tag","Date","N.Stems","N.D.Stems","max.height","picture","stem.index","stem.height","percent.tissue.infected","length.tissue.infected","N.pustules.middle")]
  
  ## get rid of data with missing or NA length.tissue.infected or N.pustules.middle or stem.height
  if(any(is.na(focal.dis.plants$length.tissue.infected))) {focal.dis.plants<-focal.dis.plants[-which(is.na(focal.dis.plants$length.tissue.infected)),]}
  if(any(is.na(focal.dis.plants$length.tissue.infected))) {focal.dis.plants<-focal.dis.plants[-which(focal.dis.plants$length.tissue.infected==""),]}
  if(any(is.na(focal.dis.plants$N.pustules.middle))) {focal.dis.plants<-focal.dis.plants[-which(is.na(focal.dis.plants$N.pustules.middle)),]}
  if(any(is.na(focal.dis.plants$N.pustules.middle))) {focal.dis.plants<-focal.dis.plants[-which(focal.dis.plants$N.pustules.middle==""),]}
  if(any(is.na(focal.dis.plants$stem.height))) {focal.dis.plants<-focal.dis.plants[-which(is.na(focal.dis.plants$stem.height)),]}
  if(any(is.na(focal.dis.plants$stem.height))) {focal.dis.plants<-focal.dis.plants[-which(focal.dis.plants$stem.height==""),]}
  
  ## summarize by plant for each date
  
  years<-c()
  sites<-c()
  dates<-c()
  tags<-c()
  n.stems<-c()
  n.d.stems<-c()
  max.heights<-c()
  reference.pictures<-c()
  infection.intensity<-c()
  
  pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv") ### load pustule data for date referencing
  
  for (tag in unique(focal.dis.plants$Tag))
  {
    sub.focal.dis.plants.1<-focal.dis.plants[which(focal.dis.plants$Tag==tag),]
    
    for (date in unique(sub.focal.dis.plants.1$Date))
    {
      sub.focal.dis.plants.2<-sub.focal.dis.plants.1[which(sub.focal.dis.plants.1$Date==date),]
      
      if(any(!is.na(sub.focal.dis.plants.2$N.Stems)) & any(!is.na(sub.focal.dis.plants.2$N.D.Stems)))
      {
        ### new values
        new.year<-sub.focal.dis.plants.2[1,"Year"]
        new.site<-sub.focal.dis.plants.2[1,"Site"]
        new.tag<-sub.focal.dis.plants.2[1,"Tag"]
        new.date<-date
        new.n.stems<-sub.focal.dis.plants.2[1,"N.Stems"]
        new.n.d.stems<-sub.focal.dis.plants.2[1,"N.D.Stems"]
        new.max.height<-sub.focal.dis.plants.2[1,"max.height"]
        new.reference.picture<-sub.focal.dis.plants.2[1,"picture"]
        
        ## #calculate plant infection intensity
        new.infection.intensity<-new.n.d.stems*mean(as.numeric(sub.focal.dis.plants.2$length.tissue.infected)*as.numeric(sub.focal.dis.plants.2$N.pustules.middle))
        
        ### store new values
        years<-c(years,new.year)
        sites<-c(sites,new.site)
        dates<-c(dates,new.date)
        tags<-c(tags,new.tag)
        n.stems<-c(n.stems,new.n.stems)
        n.d.stems<-c(n.d.stems,new.n.d.stems)
        max.heights<-c(max.heights,new.max.height)
        reference.pictures<-c(reference.pictures,new.reference.picture)
        infection.intensity<-c(infection.intensity,new.infection.intensity)
        
      }
    }
  }
  
  infection.intensity<-data.frame("Year"=years,"Site"=sites,"Tag"=tags,"Date"=dates,"N.Stems"=n.stems,"N.D.Stems"=n.d.stems,"max.height"=max.heights,"picture"=reference.pictures,"infection.intensity"=infection.intensity)
  
  ## correct 0 intensities to .1, logic being that this is a measure of tot infection load, and it shouldn't be less than 1 pustule (coded as .1cm infected tissue, 1 pustule/leaf)
  infection.intensity$infection.intensity[which(infection.intensity$infection.intensity<.1)]<-.1
  
  ## finish cleaning
  infection.intensity$Date<-as.Date(infection.intensity$Date,tryFormats = "%m/%d/%y")
  
  ## save data
  saveRDS(infection.intensity,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/infection.intensity.RDS")
  
  # reformat to delta infection intensity
  
  ## trim out records w/o temp/rh data
  infection.intensity<-infection.intensity[-intersect(which(infection.intensity$Site=="HM"),which(as.Date(infection.intensity$Date,tryFormats = "%m/%d/%y")>as.Date("2020-07-10"))),]
  
  ## make new data object for change
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  sites<-c()
  max.heights<-c()
  n.stems<-c()
  n.d.stems<-c()
  max.heights<-c()
  start.infection.intensity<-c()
  end.infection.intensity<-c()
  days<-c()
  mean.temp<-c()
  max.temp<-c()
  min.temp<-c()
  mean.abs.hum<-c()
  mean.daily.rain<-c()

  for (tag in unique(infection.intensity$Tag))
  {
    sub.infection.intensity.1<-infection.intensity[which(infection.intensity$Tag==tag),]
    sub.infection.intensity.1<-sub.infection.intensity.1[order(sub.infection.intensity.1$Date),]
    site<-sub.infection.intensity.1[1,"Site"]
    if(dim(sub.infection.intensity.1)[1]>=2)
    {
      for(i in 1:(dim(sub.infection.intensity.1)[1]-1))
      {
        ### pull reference data
        #### use pictures to get date times for plant, average
        date0.pic<-sub.infection.intensity.1[i,"picture"]
        date0<-pustules[which(pustules$picture==date0.pic)[1],"date"][1]
        date0<-as.POSIXct(date0,tryFormats = "%m/%d/%y %H:%M",tz="UTC")
        date1.pic<-sub.infection.intensity.1[i+1,"picture"]
        date1<-pustules[which(pustules$picture==date1.pic)[1],"date"][1]
        date1<-as.POSIXct(date1,tryFormats = "%m/%d/%y %H:%M",tz="UTC")
        
        #### if unable to line up pic with date time, use mean date time of all measurments at site on given day
        if(is.na(date0))
        {
          simple.date<-sub.infection.intensity.1[i,"Date"]
          sub.pustules<-pustules[which(pustules$site==site),]
          sub.pustules<-sub.pustules[which(as.Date(sub.pustules$date,tryFormats = "%m/%d/%y")==as.Date(simple.date,tryFormats = "%m/%d/%y")),]
          alt.date.0.pics<-unique(sub.pustules$picture)
          alt.dates<-pustules[which(pustules$picture %in% alt.date.0.pics),"date"]
          date0<-mean(as.POSIXct(alt.dates,tryFormats = "%m/%d/%y %H:%M",tz="UTC"))
        }
        if(is.na(date1))
        {
          simple.date<-sub.infection.intensity.1[i+1,"Date"]
          sub.pustules<-pustules[which(pustules$site==site),]
          sub.pustules<-sub.pustules[which(as.Date(sub.pustules$date,tryFormats = "%m/%d/%y")==as.Date(simple.date,tryFormats = "%m/%d/%y")),]
          alt.date.1.pics<-unique(sub.pustules$picture)
          alt.dates<-pustules[which(pustules$picture %in% alt.date.1.pics),"date"]
          date1<-mean(as.POSIXct(alt.dates,tryFormats = "%m/%d/%y %H:%M",tz="UTC"))
        }
        
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
        
        new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
        new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
        new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
        
        abs.hum<-13.24732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+temp.rh.sub$temp.c)
        new.mean.abs.hum<-mean(abs.hum,na.rm=T)
        
        new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)*(12*24)
        
        #pull out core predictors
        new.start.infection.intensity<-sub.infection.intensity.1[i,"infection.intensity"]
        new.end.infection.intensity<-sub.infection.intensity.1[i+1,"infection.intensity"]
        delta.days<-as.numeric(date1-date0)
        
        new.n.stems<-sub.infection.intensity.1[i,"N.Stems"]
        new.n.d.stems<-sub.infection.intensity.1[i,"N.D.Stems"]
        new.max.height<-sub.infection.intensity.1[i,"max.height"]
        
        #store values
        tags<-c(tags,tag)
        sites<-c(sites,site)
        
        start.infection.intensity<-c(start.infection.intensity,new.start.infection.intensity)
        end.infection.intensity<-c(end.infection.intensity,new.end.infection.intensity)
        n.stems<-c(n.stems,new.n.stems)
        n.d.stems<-c(n.d.stems,new.n.d.stems)
        max.heights<-c(max.heights,new.max.height)
        days<-c(days,delta.days)
        
        mean.temp<-c(mean.temp,new.mean.temp)
        max.temp<-c(max.temp,new.max.temp)
        min.temp<-c(min.temp,new.min.temp)
        mean.abs.hum<-c(mean.abs.hum,new.mean.abs.hum)

        mean.daily.rain<-c(mean.daily.rain,new.mean.daily.rain)

      } 
    }
  }
  
  delta.infection.intensity<-data.frame(tag=factor(tags),site=factor(sites),max.height=max.heights,time=days,N.stems=n.stems,N.D.Stems=n.d.stems,infection.intensity=start.infection.intensity,infection.intensity.next=end.infection.intensity,
                           mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,
                           mean.abs.hum=mean.abs.hum,mean.daily.rain=mean.daily.rain)
  
  saveRDS(delta.infection.intensity,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.infection.intensity.RDS")
}

infection.intensity<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/infection.intensity.RDS")
delta.infection.intensity<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.infection.intensity.RDS")


