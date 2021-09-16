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
  
  # prep enviro data
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
  
  # prep data
  
  ## load data
  within.host<-read.csv("~/Documents/GitHub/flax.rust/data/Withinhost.csv")
  
  ## subset to 2020
  within.host<-within.host[which(within.host$Year==2020),]
  
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
  
  ## save data
  saveRDS(plants,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
  
  ## trim out records w/o temp/rh data
  plants<-plants[-intersect(which(plants$Site=="HM"),which(as.Date(plants$Date,tryFormats = "%m/%d/%y")>as.Date("2020-07-10"))),]
  
  
  ## make new data object for change in pustule size
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  sites<-c()
  n.stems<-c()
  n.d.stems<-c()
  max.heights<-c()
  start.plant.inf.intens<-c()
  end.plant.inf.intens<-c()
  days<-c()
  mean.temp<-c()
  max.temp<-c()
  min.temp<-c()
  mean.abs.hum<-c() #absolute humidity
  max.abs.hum<-c()
  min.abs.hum<-c()
  #mean.vpd<-c() #vapor pressure deficit
  #max.vpd<-c() 
  #min.vpd<-c() 
  mean.wetness<-c()
  tot.rain<-c()
  mean.solar<-c()
  pred.pustule.diam.growths<-c()
  pred.pustule.num.increases<-c()
  
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
        
        new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
        new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
        new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
        
        abs.hum<-6.112*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh*2.1674/(273.15+T)
        new.mean.abs.hum<-mean(abs.hum,na.rm=T) #absolute humidity, see https://www.medrxiv.org/content/10.1101/2020.02.12.20022467v1.full.pdf
        new.max.abs.hum<-max(abs.hum,na.rm=T)
        new.min.abs.hum<-min(abs.hum,na.rm=T)
        
        #svps<- 0.6108 * exp(17.27 * temp.rh.sub$temp.c / (temp.rh.sub$temp.c + 237.3)) #saturation vapor pressures
        #avps<- temp.rh.sub$rh / 100 * svps #actual vapor pressures 
        #vpds<-avps-svps
        
        #new.mean.vpd<-mean(vpds,na.rm=T)
        #new.max.vpd<-max(vpds,na.rm=T)
        #new.min.vpd<-min(vpds,na.rm=T)
        
        new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
        new.tot.rain<-sum(weath.sub$rain,na.rm=T)
        new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
        
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
        pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,
                                            "time"=delta.days,"site"=site,
                                            "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                            "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                            "mean.wetness"=new.mean.wetness,
                                            "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain)
        pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
        
        #predict change in number of pustules from enviro conditions
        #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
        n.pustules.model.new.n.pustules<-1  #predict change for small N  pustules, arbitrarily pick 1
        obs.time<-delta.days
        n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,
                                               "time"=delta.days,"site"=site,
                                               "mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,
                                               "mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,
                                               "mean.wetness"=new.mean.wetness,
                                               "mean.solar"=new.mean.solar,"tot.rain"=new.tot.rain)
        pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
        
        #store values
        tags<-c(tags,tag)
        sites<-c(sites,site)
        start.plant.inf.intens<-c(start.plant.inf.intens,new.start.plant.inf.intens)
        end.plant.inf.intens<-c(end.plant.inf.intens,new.end.plant.inf.intens)
        n.stems<-c(n.stems,new.n.stems)
        n.d.stems<-c(n.d.stems,new.n.d.stems)
        max.heights<-c(max.heights,new.max.heights)
        days<-c(days,delta.days)
        
        mean.temp<-c(mean.temp,new.mean.temp)
        max.temp<-c(max.temp,new.max.temp)
        min.temp<-c(min.temp,new.min.temp)
        mean.abs.hum<-c(mean.abs.hum,new.mean.abs.hum)
        max.abs.hum<-c(max.abs.hum,new.max.abs.hum)
        min.abs.hum<-c(min.abs.hum,new.min.abs.hum)
        #mean.vpd<-c(mean.vpd,new.mean.vpd)
        #max.vpd<-c(max.vpd,new.max.vpd)
        #min.vpd<-c(min.vpd,new.min.vpd)
        mean.wetness<-c(mean.wetness,new.mean.wetness)
        tot.rain<-c(tot.rain,new.tot.rain)
        mean.solar<-c(mean.solar,new.mean.solar)
        
        pred.pustule.diam.growths<-c(pred.pustule.diam.growths,pred.pustule.diam.growth)
        pred.pustule.num.increases<-c(pred.pustule.num.increases,pred.pustule.num.increase)
      } 
    }
  }
  
  delta.plants<-data.frame(tag=factor(tags),site=factor(sites),time=days,N.stems=n.stems,N.D.Stems=n.d.stems,max.height=max.heights,plant.inf.intens=start.plant.inf.intens,plant.inf.intens.next=end.plant.inf.intens,
                           mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,
                           mean.abs.hum=mean.abs.hum,max.abs.hum=max.abs.hum,min.abs.hum=min.abs.hum,
                           #mean.vpd=mean.vpd,max.vpd=max.vpd,min.vpd=min.vpd,
                           mean.wetness=mean.wetness,tot.rain=tot.rain,mean.solar=mean.solar,
                           pred.pustule.diam.growth=pred.pustule.diam.growths,pred.pustule.num.increase=pred.pustule.num.increases)
  
  saveRDS(delta.plants,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")
}

plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
delta.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")


