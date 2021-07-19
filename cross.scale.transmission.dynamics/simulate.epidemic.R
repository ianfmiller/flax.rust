# load + prep data
library(lme4)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant inf intens change funcs.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant height change funcs.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/starting plant inf intens model.R")
foi.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/foi.model.RDS")
plant.inf.intens<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")



# simulate epi function

simulate.epi<-function(site,temp.addition,print.progress=T)
{
  ## setup
  
  ### subset data
  sub.locs<-corrected.locs[which(corrected.locs$Site==site),]
  sub.epi<-corrected.epi[which(corrected.epi$Site==site),]
  #start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased==min(sub.epi$Date.First.Observed.Diseased)),]
  start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased<=unique(sub.epi$Date.First.Observed.Diseased)[2]),] ## option for shorter sim
  
  ### data frame to fill
  pred.epi<-data.frame("site"=factor(),"tag"=factor(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"date"=character(),"tot.stems"=numeric(),"status"=numeric(),"max.height"=numeric(),"plant.inf.intens"=numeric())
  
  ## initialize
  
  for(i in 1:dim(sub.locs)[1])
  {
    date0<-max(start.epi$Date.First.Observed.Diseased) ### first observed date
    tag<-sub.locs[i,"tag"] ### tag
    if(sub.locs[i,"tag"] %in% start.epi$Tag) ### if plant is starting diseased
    {
      status<-1 ### set status to diseased
      if(length(intersect(which(plant.inf.intens$Date==date0),which(plant.inf.intens$Tag==tag)))>0) ### if plant.inf.intens was observed on that date
      {
        new.plant.inf.intens<-plant.inf.intens[intersect(which(plant.inf.intens$Date==date0),which(plant.inf.intens$Tag==tag)),"plant.inf.intens"] #### record plant.inf.intens as observed value
      } else ### if plant.inf.intens was not observed on that date
      {
        if(!tag %in% plant.inf.intens$Tag)
        {
          new.plant.inf.intens<-mean(plant.inf.intens$plant.inf.intens[intersect(which(plant.inf.intens$Date<=date0),which(plant.inf.intens$Site==site))])
        } else {
          target.date.set<-plant.inf.intens[which(plant.inf.intens$Tag==tag),"Date"] #### get set of dates w/ observations
          target.date.set<-target.date.set[which(target.date.set>date0)] #### look at those after date0
          target.date<-min(target.date.set) #### get closest date w/ observation after date0
          target.date.plant.inf.intens<-plant.inf.intens[intersect(which(plant.inf.intens$Date==target.date),which(plant.inf.intens$Tag==tag)),"plant.inf.intens"] #### get observed plant.inf.intens on target date
          new.plant.inf.intens<-predict.plant.inf.intens.last(target.date.plant.inf.intens,site,as.POSIXct(paste0(date0," 12:00:00")),as.POSIXct(paste0(target.date," 12:00:00"))) #### hindcast plant.inf.intens for date0 
        }
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
  
  for(date.index in 1:(length(unique(sub.epi$Date.First.Observed.Diseased))-2)) 
  {
    date0<-unique(sub.epi$Date.First.Observed.Diseased)[date.index+1] ### last date
    date0<-as.POSIXct(paste0(date0," 12:00:00")) ### convert format
    date1<-unique(sub.epi$Date.First.Observed.Diseased)[date.index+2] ### next date
    date1<-as.POSIXct(paste0(date1," 12:00:00")) ### convert format
    delta.days<-as.numeric(as.Date(date1)-as.Date(date0))
    
    ### load weather data
    #### subst temp rh data to relevant window
    temp.rh.sub<-all.temp.rh[which(all.temp.rh$site==site),] #### pull out temp data for site
    temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #### pull out relevant data
    temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #### pull out relevant data
    temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #### throw out NAs
    temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
    
    temp.rh.sub$temp.c<-temp.rh.sub$temp.c+temp.addition ### shift temperature
    #### subset weather data to relevant window
    weath.sub<-all.weath[which(all.weath$site==site),] #pull out weath data for site
    weath.sub<-subset(weath.sub,date<=date1) #### pull out relevant data
    weath.sub<-subset(weath.sub,date>=date0) #### pull out relevant data
    weath.sub<-cbind(weath.sub,interval.length=c(diff(as.numeric(weath.sub$date))/(60*60*24),NA))
    
    #### calculate environmental variable metrics
    mean.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T)/delta.days #temperature days
    mean.temp.days.16.22<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)/delta.days  #time (in days) during which temp between 16 and 22 celsius
    mean.temp.days.7.30<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)/delta.days  #time (in days) during which temp between 7 and 30 celsius
    mean.dew.point.days<-sum(temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)/delta.days  #Dew point days
    
    #calculate weather metrics
    mean.wetness.days<-sum(weath.sub$wetness*weath.sub$interval.length,na.rm = T)/delta.days 
    mean.tot.rain<-sum(weath.sub$rain,na.rm=T)/delta.days 
    mean.solar.days<-sum(weath.sub$solar.radiation*weath.sub$interval.length,na.rm = T)/delta.days 
    mean.wind.speed.days<-sum(weath.sub$wind.speed*weath.sub$interval.length,na.rm = T)/delta.days 
    mean.gust.speed.days<-sum(weath.sub$wind.direction*weath.sub$interval.length,na.rm = T)/delta.days 
    
    #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
    mean.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)/delta.days 
    mean.temp.16.22.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)/delta.days 
    mean.temp.7.30.dew.point.days<-sum(1*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)/delta.days 
    mean.temp.wetness.days<-sum(weath.sub$temp*weath.sub$wetness*weath.sub$interval.length,na.rm = T)/delta.days 
    mean.temp.16.22.wetness.days<-sum(weath.sub$temp.16.22*weath.sub$wetness*weath.sub$interval.length,na.rm = T)/delta.days 
    mean.temp.7.30.wetness.days<-sum(weath.sub$temp.7.30*weath.sub$wetness*weath.sub$interval.length,na.rm = T)/delta.days 
    
    #calculate weather metrics
    new.wetness.days<-sum(weath.sub$wetness*weath.sub$interval.length,na.rm = T)
    new.tot.rain<-sum(weath.sub$rain,na.rm=T)
    new.solar.days<-sum(weath.sub$solar.radiation*weath.sub$interval.length,na.rm = T)
    new.wind.speed.days<-sum(weath.sub$wind.speed*weath.sub$interval.length,na.rm = T)
    new.gust.speed.days<-sum(weath.sub$wind.direction*weath.sub$interval.length,na.rm = T)
    
    
    #predict pustule growth from pustule growth model and enviro conditions
    #pustule.model.vars<-names(fixef(pustule.model))[2:length(names(fixef(pustule.model)))]
    pustule.model.new.area<-.01 #predict change for small pustule, arbitrarily pick .01
    obs.time<-delta.days
    pustule.model.pred.data<-data.frame("area"=pustule.model.new.area,"temp.days.16.22"=mean.temp.days.16.22,"dew.point.days"=mean.dew.point.days,"temp.16.22.dew.point.days"=mean.temp.16.22.dew.point.days,"temp.wetness.days"=mean.temp.wetness.days,"tot.rain"=new.tot.rain/delta.days)
    pred.pustule.diam.growth<-predict(pustule.model,newdata=pustule.model.pred.data,re.form=~0)
    
    #predict change in number of pustules from enviro conditions
    #n.pustule.model.vars<-names(fixef(n.pustule.model))[2:length(names(fixef(n.pustule.model)))]
    n.pustules.model.new.n.pustules<-0 #included only for offset, picked 0 for ease of interpretability
    obs.time<-delta.days
    n.pustules.model.pred.data<-data.frame("n.pustules"=n.pustules.model.new.n.pustules,"temp.days.16.22"=mean.temp.days.16.22,"temp.16.22.wetness.days"=mean.temp.16.22.wetness.days)
    pred.pustule.num.increase<-predict(n.pustules.model,newdata=n.pustules.model.pred.data,re.form=~0)
    
    #predict change in plant.inf.intensity from enviro conditions
    #plant.model.vars<-names(fixef(plant.model))[2:length(names(fixef(plant.model)))]
    plant.model.new.plant.inf.intens<-.1
    obs.time<-delta.days
    plant.model.pred.data<-data.frame("plant.inf.intens"=plant.model.new.plant.inf.intens,"dew.point.days"=mean.dew.point.days,"temp.7.30.dew.point.days"=mean.temp.7.30.dew.point.days,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"site"=site)
    pred.plant.inf.intens.increase<-10^predict(plant.model,newdata=plant.model.pred.data,exclude = 's(site)')
    
    ### compare weather data coverage to make sure foi is not being underestimated
    time.diff.weather.data<-weath.sub$date[nrow(weath.sub)]-weath.sub$date[1]
    units(time.diff.weather.data)<-"days"
    time.diff.epi.data<-date1-date0
    units(time.diff.epi.data)<-"days"
    if(time.diff.epi.data>time.diff.weather.data) {foi.mod<-(as.numeric(time.diff.weather.data)/as.numeric(time.diff.epi.data))^-1} else{foi.mod<-1}
    
    ### get data from previous date
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
          foi<-foi+predict.kernel.tilted.plume(q=q,H=half.height,k=4.828517e-07,alphaz=1.687830e-06,Ws=9.299220e-01,xtarget=xcord-sourceX,ytarget=ycord-sourceY,wind.data=weath.sub)
        }
        foi<-foi*foi.mod ### correct for any gaps in weath data
        pred.inf.odds<-predict(foi.model,newdata = data.frame("foi"=foi,"height.cm"=last.epi[i,"max.height"],"mean.temp.days.16.22"=mean.temp.days.16.22,"mean.temp.days.7.30"=mean.temp.days.7.30,"mean.dew.point.days"=mean.dew.point.days,"mean.temp.7.30.dew.point.days"=mean.temp.7.30.dew.point.days,"mean.wetness.days"=mean.wetness.days,"mean.temp.wetness.days"=mean.temp.wetness.days,"pred.pustule.diam.growth"=pred.pustule.diam.growth,"pred.pustule.num.increase"=pred.pustule.num.increase),type="response") ### predict odds of becoming infected
        #pred.inf.odds<-predict(foi.model,newdata = data.frame("foi"=foi,"height.cm"=last.epi[i,"max.height"],"mean.temp.days.16.22"=mean.temp.days.16.22,"mean.temp.7.30.dew.point.days"=mean.temp.7.30.dew.point.days),type="response",re.form=NA) ### predict odds of becoming infected
        #pred.inf.odds<-predict(foi.model,newdata=data.frame("foi"=foi),type="response")
        draw<-runif(1) ### draw random number between 0 and 1
        
        if(draw<=pred.inf.odds) ### if draw <= odds make plant infected
        {
          new.status<-1
          new.plant.inf.intens<-starting.plant.inf.intens.mod(runif(1)) # starting infection intensity draw from distribution fitted to empirical data
        } else { ### else plant stays healthy
          new.status<-0
          new.plant.inf.intens<-NA
        }
      }
      
      ### if plant is diseased increase infection intensity 
      if(last.epi[i,"status"]==1)
      {
        new.status<-1
        
        new.plant.inf.intens<-predict.plant.inf.intens.boot(plant.inf.intens.last = last.epi[i,"plant.inf.intens"],site = site,date0 = date0,date1 = date1)
      }
      
      ### plant growth
      new.height<-predict.plant.growth.boot(last.epi[i,"max.height"],site,date0,date1)
      
      new.row<-data.frame("site"=site,"tag"=sub.locs[i,"tag"],"X"=sub.locs[i,"X"],"Y"=sub.locs[i,"Y"],"x"=sub.locs[i,"x"],"y"=sub.locs[i,"y"],"date"= as.Date(date1),"status"=new.status,"max.height"=new.height,"plant.inf.intens"=new.plant.inf.intens) ### new data row
      pred.epi<-rbind(pred.epi,new.row) ### join data
      if(print.progress) {print(paste0("finished ",site," ",as.Date(date1)," tag ",i,"/",dim(sub.locs)[1]))}
    }
  }
  pred.epi
}

simulate.epi("GM",0,print.progress = T)->pred.epi

par(mfrow=c(3,3))
for(date in unique(pred.epi$date))
{
  hist(pred.epi[which(pred.epi$date==date),"plant.inf.intens"],breaks=20)
}

library(parallel)
library(doParallel)
n.cores<-5
registerDoParallel(n.cores)
site<-"GM"
pred.epi.all.0<-foreach(k = 1:1, .multicombine = T) %dopar% simulate.epi(site,0,print.progress = T)
pred.epi.all.1.8<-foreach(k = 1:1, .multicombine = T) %dopar% simulate.epi(site,1.8,print.progress = F)
pred.epi.all.3.7<-foreach(k = 1:1, .multicombine = T) %dopar% simulate.epi(site,3.7,print.progress = F)

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

plot.purple<-t_col("purple",75)
plot.red<-t_col("red",75)
plot.orange<-t_col("orange",75)

sub.epi<-corrected.epi[which(corrected.epi$Site==site),]
sub.locs<-corrected.locs[which(corrected.locs$Site==site),]

par(mfrow=c(1,1))
par(mar=c(6,6,6,6))
plot(unique(pred.epi.all.0[[1]]$date),rep(0,times=length(unique(pred.epi.all.0[[1]]$date))),ylim=c(0,1),xlab="date",ylab="prev",type="n",cex.axis=2,cex.lab=2)
xvals<-c()
yvals<-c()
for(i in 1:9)
{
  date<-unique(sub.epi$Date.First.Observed.Diseased)[i]
  xvals<-c(xvals,date)
  yvals<-c(yvals,dim(sub.epi[which(sub.epi$Date.First.Observed.Diseased<=date),])[1]/dim(sub.locs)[1])
}
points(xvals,yvals,type="l",col="black",lwd=4)


for(k in 1:length(pred.epi.all.0))
{
  xvals<-c()
  yvals<-c()
  pred.epi<-pred.epi.all.0[[k]]
  for(i in 1:9)
  {
    date<-unique(pred.epi$date)[i]
    xvals<-c(xvals,date)
    yvals<-c(yvals,sum(pred.epi[which(pred.epi$date==date),"status"])/dim(sub.locs)[1])
  }
  points(xvals,yvals,type="l",col=plot.orange)
}

xvals<-c()
yvals<-c()
for(k in 1:length(unique(pred.epi.all.0[[1]]$date)))
{
  sub.prevs<-c()
  date<-unique(pred.epi.all.0[[1]]$date)[k]
  for(j in 1:length(pred.epi.all.0))
  {
   sub.dat<-pred.epi.all.0[[j]]
   sub.dat<-sub.dat[which(sub.dat$date==date),]
   sub.prevs<-c(sub.prevs,sum(sub.dat$status)/dim(sub.dat)[1])
  }
  yvals<-c(yvals,mean(sub.prevs))
  xvals<-c(xvals,date)
}
points(xvals,yvals,type="l",col="orange",lwd=4)


for(k in 1:length(pred.epi.all.1.8))
{
  xvals<-c()
  yvals<-c()
  pred.epi<-pred.epi.all.1.8[[k]]
  for(i in 1:9)
  {
    date<-unique(pred.epi$date)[i]
    xvals<-c(xvals,date)
    yvals<-c(yvals,sum(pred.epi[which(pred.epi$date==date),"status"])/dim(sub.locs)[1])
  }
  points(xvals,yvals,type="l",col=plot.red)
}

xvals<-c()
yvals<-c()
for(k in 1:length(unique(pred.epi.all.1.8[[1]]$date)))
{
  sub.prevs<-c()
  date<-unique(pred.epi.all.1.8[[1]]$date)[k]
  for(j in 1:length(pred.epi.all.0))
  {
    sub.dat<-pred.epi.all.1.8[[j]]
    sub.dat<-sub.dat[which(sub.dat$date==date),]
    sub.prevs<-c(sub.prevs,sum(sub.dat$status)/dim(sub.dat)[1])
  }
  yvals<-c(yvals,mean(sub.prevs))
  xvals<-c(xvals,date)
}
points(xvals,yvals,type="l",col="red",lwd=4)


for(k in 1:length(pred.epi.all.3.7))
{
  xvals<-c()
  yvals<-c()
  pred.epi<-pred.epi.all.3.7[[k]]
  for(i in 1:9)
  {
    date<-unique(pred.epi$date)[i]
    xvals<-c(xvals,date)
    yvals<-c(yvals,sum(pred.epi[which(pred.epi$date==date),"status"])/dim(sub.locs)[1])
  }
  points(xvals,yvals,type="l",col=plot.purple)
}

xvals<-c()
yvals<-c()
for(k in 1:length(unique(pred.epi.all.3.7[[1]]$date)))
{
  sub.prevs<-c()
  date<-unique(pred.epi.all.3.7[[1]]$date)[k]
  for(j in 1:length(pred.epi.all.3.7))
  {
    sub.dat<-pred.epi.all.3.7[[j]]
    sub.dat<-sub.dat[which(sub.dat$date==date),]
    sub.prevs<-c(sub.prevs,sum(sub.dat$status)/dim(sub.dat)[1])
  }
  yvals<-c(yvals,mean(sub.prevs))
  xvals<-c(xvals,date)
}
points(xvals,yvals,type="l",col="purple",lwd=4)
legend("bottomright",legend=c("data","+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col=c("black","orange","red","purple"),lwd=4,cex=2)




par(mfrow=c(2,4))

library(viridis)
#pred.epi<-pred.epi.all.3.7[[1]]
#layout(matrix(c(1,1,2,2,3,3,6,4,4,5,5,7),2,6,byrow = T))
#par(mar=c(6,6,6,6))
data<-pred.epi.all.0[[1]]
for(date in unique(data$date))
{
  dat<-data[which(data$date==date),]
  sub.dat<-dat[which(dat$status==1),]
  plot(dat$X+dat$x,dat$Y+dat$y,xlab="X",ylab="Y",main=as.Date(date,origin="1970-01-01"),pch=16,cex=.9,cex.lab=1.25,cex.axis=1.25,cex.main=1.5,asp=1)
  points(sub.dat$X+sub.dat$x,sub.dat$Y+sub.dat$y,col="red",cex=2,lwd=2)
}

