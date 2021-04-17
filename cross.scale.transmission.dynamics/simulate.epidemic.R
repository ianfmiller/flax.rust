# load + prep data
library(lme4)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant inf intens change funcs.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
foi.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")
plant.inf.intens<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")

# simulate

## setup
site<-"GM"

### subset data
sub.locs<-corrected.locs[which(corrected.locs$Site==site),]
sub.epi<-corrected.epi[which(corrected.epi$Site==site),]
start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased==min(sub.epi$Date.First.Observed.Diseased)),]

### data frame to fill
pred.epi<-data.frame("site"=factor(),"tag"=factor(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"date"=character(),"tot.stems"=numeric(),"status"=numeric(),"max.height"=numeric(),"plant.inf.intens"=numeric())

## initialize

for(i in 1:dim(sub.locs)[1])
{
  date0<-min(start.epi$Date.First.Observed.Diseased) ### first observed date
  tag<-sub.locs[i,"tag"] ### tag
  if(sub.locs[i,"tag"] %in% start.epi$Tag) ### if plant is starting diseased
  {
    status<-1 ### set status to diseased
    if(length(intersect(which(plant.inf.intens$Date==date0),which(plant.inf.intens$Tag==tag)))>0) ### if plant.inf.intens was observed on that date
    {
      new.plant.inf.intens<-plant.inf.intens[intersect(which(plant.inf.intens$Date==date0),which(plant.inf.intens$Tag==tag)),"plant.inf.intens"] #### record plant.inf.intens as observed value
    } else ### if plant.inf.intens was not observed on that date
    {
     target.date.set<-plant.inf.intens[which(plant.inf.intens$Tag==tag),"Date"] #### get set of dates w/ observations
     target.date.set<-target.date.set[which(target.date.set>date0)] #### look at those after date0
     target.date<-min(target.date.set) #### get closest date w/ observation after date0
     target.date.plant.inf.intens<-plant.inf.intens[intersect(which(plant.inf.intens$Date==target.date),which(plant.inf.intens$Tag==tag)),"plant.inf.intens"] #### get observed plant.inf.intens on target date
     new.plant.inf.intens<-predict.plant.inf.intens.last(target.date.plant.inf.intens,site,as.POSIXct(paste0(date0," 12:00:00")),as.POSIXct(paste0(target.date," 12:00:00"))) #### hindcast plant.inf.intens for date0
    }
  } else ### if plant is starting healthy
  {
    status<-0 ### set status to healthy
    new.plant.inf.intens<-NA  ### set plant.inf.intens to NA
  }
  new.row<-data.frame("site"=site,"tag"=sub.locs[i,"tag"],"X"=sub.locs[i,"X"],"Y"=sub.locs[i,"Y"],"x"=sub.locs[i,"x"],"y"=sub.locs[i,"y"],"date"= date0,"status"=status,"max.height"=max(sub.locs[i,"height.cm"],5,na.rm = T),"plant.inf.intens"=new.plant.inf.intens) ### new data row; set max height to 5 for seedlings w/o max height recorded         
  pred.epi<-rbind(pred.epi,new.row) ### join data
}

## loopt to simulate epi process

for(date.index in 1:(length(unique(sub.epi$Date.First.Observed.Diseased))-1)) 
{
  date0<-unique(sub.epi$Date.First.Observed.Diseased)[date.index] ### last date
  date0<-as.POSIXct(paste0(date0," 12:00:00")) ### convert format
  date1<-unique(sub.epi$Date.First.Observed.Diseased)[date.index+1] ### next date
  date1<-as.POSIXct(paste0(date1," 12:00:00")) ### convert format
  
  ### load weather data
  #### subst temp rh data to relevant window
  temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #### pull out temp data for site
  temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
  temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
  temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
  
  #### subset weather data to relevant window
  weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
  weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
  weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
  weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
  
  ### get date from previous date
  last.epi<-pred.epi[which(pred.epi$date==as.Date(date0)),]
  ### get set of diseased plants from last date
  last.epi.dis.set<-last.epi[which(last.epi$status==1),]
  
  for (i in 1:dim(sub.locs)[1]) ### for each plant
  {
    ### if plant is healthy calulate odds of becoming infected. If infected, set p.inf.intens to 0.1
    if(last.epi[i,"status"]==0)
    {
      foi<-0
      xcord<-sub.locs[i,"X"]+sub.locs[i,"x"]
      ycord<-sub.locs[i,"Y"]+sub.locs[i,"y"]
      for(j in 1:dim(last.epi.dis.set)[1])
      {
        sourceX<-last.epi.dis.set[j,"X"]+last.epi.dis.set[j,"x"]
        sourceY<-last.epi.dis.set[j,"Y"]+last.epi.dis.set[j,"y"]
        half.height<-last.epi.dis.set[j,"max.height"]/2
        q<-last.epi.dis.set[j,"plant.inf.intens"]
        foi<-foi+predict.kernel.tilted.plume(q=q,H=half.height,k=5.803369e-07,alphaz=1.596314e-01,Ws=1.100707e+00,xtarget=xcord-sourceX,ytarget=ycord-sourceY,wind.data=weath.sub)
      }
      
      pred.inf.odds<-predict(foi.model,newdata = data.frame("foi"=foi,"first.obs.height"=last.epi[i,"max.height"]),type="response") ### predict odds of becoming infected
      draw<-runif(1) ### draw random number between 0 and 1
      
      if(draw<=pred.inf.odds) ### if draw <= odds make plant infected
      {
        new.status<-1
        new.plant.inf.intens<-0.1
      } else { ### else plant stays healthy
        new.status<-0
        new.plant.inf.intens<-NA
      }
    }
    
    ### if plant is diseased increase infection intensity DETERMINISTIC FOR NOW
    if(last.epi[i,"status"]==1)
    {
      new.status<-1
      
      new.plant.inf.intens<-predict.plant.inf.intens.boot(plant.inf.intens.last = last.epi[i,"plant.inf.intens"],site = site,date0 = date0,date1 = date1)
    }
    new.row<-data.frame("site"=site,"tag"=sub.locs[i,"tag"],"X"=sub.locs[i,"X"],"Y"=sub.locs[i,"Y"],"x"=sub.locs[i,"x"],"y"=sub.locs[i,"y"],"date"= as.Date(date1),"status"=new.status,"max.height"=last.epi[i,"max.height"],"plant.inf.intens"=new.plant.inf.intens) ### new data row
    pred.epi<-rbind(pred.epi,new.row) ### join data
    print(paste0("finished ",site," ",as.Date(date1)," tag ",i,"/",dim(sub.locs)[1]))
  }
}

plot(unique(sub.epi$Date.First.Observed.Diseased),rep(0,times=length(unique(sub.epi$Date.First.Observed.Diseased))),ylim=c(0,1),xlab="date",ylab="prev",type="n")
for(i in 2:7)
{
  date<-unique(sub.epi$Date.First.Observed.Diseased)[i]
  points(date,dim(sub.epi[which(sub.epi$Date.First.Observed.Diseased<=date),])[1]/dim(sub.locs)[1],col="black")
  points(date,sum(pred.epi[which(pred.epi$date==date),"status"])/dim(sub.locs)[1],col="red")
}

par(mfrow=c(3,3))

for(date in unique(pred.epi$date))
{
  dat<-pred.epi[which(pred.epi$date==date),]
  sub.dat<-dat[which(dat$status==1),]
  plot(dat$X+dat$x,dat$Y+dat$y,xlab="X",ylab="y",main=as.Date(date,origin="1970-01-01"),pch=16,cex=.9)
  points(sub.dat$X+sub.dat$x,sub.dat$Y+sub.dat$y,col="red",cex=2)
}

