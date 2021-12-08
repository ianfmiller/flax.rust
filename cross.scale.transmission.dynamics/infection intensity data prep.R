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
  
  # load plants data
  
  plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")

  ## trim out records w/o temp/rh data
  plants<-plants[-intersect(which(plants$Site=="HM"),which(as.Date(plants$Date,tryFormats = "%m/%d/%y")>as.Date("2020-07-10"))),]
  
  
  ## make new data object for change
  
  temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}
  
  tags<-c()
  sites<-c()
  max.heights<-c()
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
  tot.rain<-c()
  mean.solar<-c()

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
        
        new.tot.rain<-sum(weath.sub$rain,na.rm=T)
        new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
        
        #pull out core predictors
        new.start.plant.inf.intens<-sub.plants.1[i,"plant.inf.intens"]
        new.end.plant.inf.intens<-sub.plants.1[i+1,"plant.inf.intens"]
        delta.days<-as.numeric(date1-date0)
        
        new.n.stems<-sub.plants.1[i,"N.Stems"]
        new.n.d.stems<-sub.plants.1[i,"N.D.Stems"]
        new.max.height<-sub.plants.1[i,"max.height"]
        
        #store values
        tags<-c(tags,tag)
        sites<-c(sites,site)
        
        start.plant.inf.intens<-c(start.plant.inf.intens,new.start.plant.inf.intens)
        end.plant.inf.intens<-c(end.plant.inf.intens,new.end.plant.inf.intens)
        n.stems<-c(n.stems,new.n.stems)
        n.d.stems<-c(n.d.stems,new.n.d.stems)
        max.heights<-c(max.heights,new.max.height)
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
        tot.rain<-c(tot.rain,new.tot.rain)
        mean.solar<-c(mean.solar,new.mean.solar)
              } 
    }
  }
  
  delta.plants<-data.frame(tag=factor(tags),site=factor(sites),max.height=max.heights,time=days,N.stems=n.stems,N.D.Stems=n.d.stems,max.height=max.heights,plant.inf.intens=start.plant.inf.intens,plant.inf.intens.next=end.plant.inf.intens,
                           mean.temp=mean.temp,max.temp=max.temp,min.temp=min.temp,
                           mean.abs.hum=mean.abs.hum,max.abs.hum=max.abs.hum,min.abs.hum=min.abs.hum,
                           #mean.vpd=mean.vpd,max.vpd=max.vpd,min.vpd=min.vpd,
                           tot.rain=tot.rain,mean.solar=mean.solar)
  
  saveRDS(delta.plants,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")
}

plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
delta.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/delta.plants.RDS")


