# load + prep data
library(mgcv)
library(MASS)
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc data set building.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant inf intens change funcs.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/spore deposition functions tilt.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/predict plant height change funcs.R")
transmission.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")
transmission.data<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/transmission.data.RDS")
heights<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plant.heights.RDS")
diseased.focal.plants<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/diseased.focal.plants.RDS")
rm(site)


# simulate epi function

simulate.epi<-function(site,temp.addition,seed,step.size=7,print.progress=T)
{
  ## setup
  set.seed(seed)
  ### subset data
  sub.locs<-corrected.locs[which(corrected.locs$Site==site),]
  sub.epi<-corrected.epi[which(corrected.epi$Site==site),]
  start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased==min(sub.epi$Date.First.Observed.Diseased)),]
  if(site=="GM") {start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased<=as.Date("2020-06-23")),]}
  if(site=="HM") {start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased<=as.Date("2020-06-25")),]}
  if(site=="BT") {start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased<=as.Date("2020-06-24")),]}
  if(site=="CC") {start.epi<-sub.epi[which(sub.epi$Date.First.Observed.Diseased<=as.Date("2020-06-22")),]}
  if(site=="CC") {sub.epi<-sub.epi[-which(sub.epi$Y %in% c(0:8,18,19)),]; sub.epi<-sub.epi[-which((sub.epi$Y %in% c(15:17)) & (sub.epi$X %in% 7:9)),]}
  
  ### data frame to fill
  pred.epi<-data.frame("site"=factor(),"tag"=factor(),"X"=numeric(),"Y"=numeric(),"x"=numeric(),"y"=numeric(),"date"=character(),"tot.stems"=numeric(),"status"=numeric(),"max.height"=numeric(),"inf.intens"=numeric())
  
  ## initialize
  
  for(i in 1:dim(sub.locs)[1])
  {
    date0<-max(start.epi$Date.First.Observed.Diseased) ### first observed date
    
    tag<-sub.locs[i,"tag"] ### tag
    if (as.Date(sub.locs[i,"Date"],tryFormats = "%m/%d/%y")==date0) ### simulate height on date0
    {
      new.max.height<-max(sub.locs[i,"height.cm"],5,na.rm = T)
    } else 
    {
      if((as.Date(sub.locs[i,"Date"],tryFormats = "%m/%d/%y") >= as.Date(date0) ) & is.na(sub.locs[i,"height.cm"])) {new.max.height<-5} else  #set max height to 5 for seedlings w/o max height recorded
        {
          new.max.height<-predict.plant.growth(height.last=max(sub.locs[i,"height.cm"],5,na.rm = T),inf.intens.last=0,site=site,date0=as.Date(sub.locs[i,"Date"],tryFormats = "%m/%d/%y"),date1=as.Date(date0)) #otherwise, predict
        }
    }
    
    if(sub.locs[i,"tag"] %in% start.epi$Tag) ### if plant is starting diseased
    {
      status<-1 ### set status to diseased
      if(length(intersect(which(diseased.focal.plants$Date==date0),which(diseased.focal.plants$Tag==tag)))>0) ### if inf.intens was observed on that date
      {
        new.inf.intens<-diseased.focal.plants[intersect(which(diseased.focal.plants$Date==date0),which(diseased.focal.plants$Tag==tag)),"inf.intens"] #### record inf.intens as observed value
      } else ### if inf.intens was not observed on that date
      {
        if(!tag %in% diseased.focal.plants$Tag)
        {
          new.inf.intens<-mean(diseased.focal.plants$inf.intens[intersect(which(diseased.focal.plants$Date<=date0),which(diseased.focal.plants$Site==site))])
        } else {
          target.date.set<-diseased.focal.plants[which(diseased.focal.plants$Tag==tag),"Date"] #### get set of dates w/ observations
          target.date.set<-target.date.set[which(target.date.set>date0)] #### look at those after date0
          target.date<-min(target.date.set) #### get closest date w/ observation after date0
          target.date.inf.intens<-diseased.focal.plants[intersect(which(diseased.focal.plants$Date==target.date),which(diseased.focal.plants$Tag==tag)),"inf.intens"] #### get observed inf.intens on target date
          new.inf.intens<-predict.inf.intens.last(inf.intens.next = target.date.inf.intens, max.height.last = new.max.height, site,as.POSIXct(paste0(date0," 12:00:00")),as.POSIXct(paste0(target.date," 12:00:00"))) #### hindcast inf.intens for date0, warnings generated because site/tag=NA not included in factors for original data, nothing to worry about 
        }
      }
    } else ### if plant is starting healthy
    {
      status<-0 ### set status to healthy
      new.inf.intens<-NA  ### set inf.intens to NA
    }
    new.row<-data.frame("site"=site,"tag"=sub.locs[i,"tag"],"X"=sub.locs[i,"X"],"Y"=sub.locs[i,"Y"],"x"=sub.locs[i,"x"],"y"=sub.locs[i,"y"],"date"= date0,"status"=status,"max.height"=new.max.height,"inf.intens"=new.inf.intens) ### new data row       
    pred.epi<-rbind(pred.epi,new.row) ### join data
  }
  
  ## loop to simulate epi process
  start.date<-min(pred.epi$date)
  end.date<-max(unique(sub.epi$Date.First.Observed.Diseased))
  sim.dates<-seq(start.date,end.date,step.size)
  if(site=="HM") {sim.dates<-sim.dates[which(sim.dates<as.Date("2020-07-10"))]}
  
  for(date.index in 1:(length(sim.dates)-1)) 
  {
    date0<-sim.dates[date.index] ### last date
    date0<-as.POSIXct(paste0(date0," 12:00:00")) ### convert format
    date1<-sim.dates[date.index+1] ### next date
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
    new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #mean temperature
    new.max.temp<-max(temp.rh.sub$temp.c,na.rm = T) #max temperature
    new.min.temp<-min(temp.rh.sub$temp.c,na.rm = T) #min temperature
    
    abs.hum<-0.1324732*exp((17.67*temp.rh.sub$temp.c)/(temp.rh.sub$temp.c+243.5))*temp.rh.sub$rh/(273.15+T)
    new.mean.abs.hum<-mean(abs.hum,na.rm=T)
    new.max.abs.hum<-max(abs.hum,na.rm=T)
    new.min.abs.hum<-min(abs.hum,na.rm=T)
    
    new.mean.wetness<-mean(weath.sub$wetness,na.rm = T)
    new.mean.daily.rain<-mean(weath.sub$rain,na.rm=T)*(12*24)
    new.mean.solar<-mean(weath.sub$solar.radiation,na.rm=T)
    
    ### compare weather data coverage to make sure spore deposition is not being underestimated
    time.diff.weather.data<-weath.sub$date[nrow(weath.sub)]-weath.sub$date[1]
    units(time.diff.weather.data)<-"days"
    time.diff.epi.data<-date1-date0
    units(time.diff.epi.data)<-"days"
    if(time.diff.epi.data-time.diff.weather.data) {transmission.mod<-(as.numeric(time.diff.weather.data)/as.numeric(time.diff.epi.data))^-1} else{transmission.mod<-1}
    
    ### get data from previous date
    last.epi<-pred.epi[which(pred.epi$date==as.Date(date0)),]
    ### get set of diseased plants from last date
    last.epi.dis.set<-last.epi[which(last.epi$status==1),]
    
    for (i in 1:dim(sub.locs)[1]) ### for each plant
    {
      ### if plant is healthy calulate odds of becoming infected. If infected, set p.inf.intens to 0.1
      if(last.epi[i,"status"]==0)
      {
        spore.deposition<-0
        xcord<-sub.locs[i,"X"]+sub.locs[i,"x"]
        ycord<-sub.locs[i,"Y"]+sub.locs[i,"y"]
        target.tag<-ifelse(is.na(last.epi[i,"tag"]),paste0(site,"X",xcord,"Y",ycord),last.epi[i,"tag"])
        for(j in 1:dim(last.epi.dis.set)[1])
        {
          sourceX<-last.epi.dis.set[j,"X"]+last.epi.dis.set[j,"x"]
          sourceY<-last.epi.dis.set[j,"Y"]+last.epi.dis.set[j,"y"]
          half.height<-last.epi.dis.set[j,"max.height"]/2
          I<-last.epi.dis.set[j,"inf.intens"]
          spore.deposition<-spore.deposition+predict.kernel.tilted.plume(I=I,H=half.height/100,k=5.739690e-07,Ws=4.451030e-02,A=7.777373e-02,xtarget=xcord-sourceX,ytarget=ycord-sourceY,wind.data=weath.sub)
        }
        spore.deposition<-spore.deposition*transmission.mod ### correct for any gaps in weath data
        if(spore.deposition==0) spore.deposition<-10e-10
        #newdata = transmission.data[1,c("site","tag","tot.spore.deposition","height.cm","mean.temp","max.temp","min.temp","mean.abs.hum","max.abs.hum","min.abs.hum","mean.solar","mean.daily.rain","mean.wetness","time")]
        #newdata[1,c("site","tag","tot.spore.deposition","height.cm","mean.temp","max.temp","min.temp","mean.abs.hum","max.abs.hum","min.abs.hum","mean.solar","mean.daily.rain","mean.wetness","time")]<-c(site,as.character(target.tag),spore.deposition,last.epi[i,"max.height"],new.mean.temp,new.max.temp,new.min.temp,new.mean.abs.hum,new.max.abs.hum,new.min.abs.hum,new.mean.solar,new.mean.daily.rain,new.mean.wetness,delta.days)
        #class(newdata[,3:14])<-"numeric"
        if(target.tag %in% unique(transmission.data$tag))
        {
          pred.inf.odds<-predict(transmission.model,
                                 newdata = data.frame("site"=site,"tag"=target.tag,"tot.spore.deposition"=spore.deposition,"height.cm"=last.epi[i,"max.height"],"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,"mean.solar"=new.mean.solar,"mean.daily.rain"=new.mean.daily.rain,"mean.wetness"=new.mean.wetness,"time"=delta.days),
                                 type="response",
                                 newdata.guaranteed=T,discrete = T) ### predict odds of becoming infected
        } else 
        {
          pred.inf.odds<-predict(transmission.model,
                                 newdata = data.frame("site"=site,"tag"=1,"tot.spore.deposition"=spore.deposition,"height.cm"=last.epi[i,"max.height"],"mean.temp"=new.mean.temp,"max.temp"=new.max.temp,"min.temp"=new.min.temp,"mean.abs.hum"=new.mean.abs.hum,"max.abs.hum"=new.max.abs.hum,"min.abs.hum"=new.min.abs.hum,"mean.solar"=new.mean.solar,"mean.daily.rain"=new.mean.daily.rain,"mean.wetness"=new.mean.wetness,"time"=delta.days),
                                 type="response",exclude = "s(tag)",
                                 newdata.guaranteed=T,discrete = T) ### predict odds of becoming infected
          
        }
        draw<-runif(1) ### draw random number between 0 and 1
        
        if(draw<=pred.inf.odds) ### if draw <= odds make plant infected
        {
          new.status<-1
          new.inf.intens<-1
        } else { ### else plant stays healthy
          new.status<-0
          new.inf.intens<-NA
        }
      }
      
      ### if plant is diseased increase infection intensity 
      if(last.epi[i,"status"]==1)
      {
        new.status<-1
        
        new.inf.intens<-predict.inf.intens.boot(inf.intens.last = last.epi[i,"inf.intens"], max.height.last = last.epi[i,"max.height"], site = site,date0 = date0,date1 = date1,temp.addition=temp.addition)
      }
      
      ### plant growth
      new.height<-predict.plant.growth.boot(last.epi[i,"max.height"],ifelse(is.na(last.epi[i,"inf.intens"]),0,last.epi[i,"max.height"]),site,date0,date1)
      
      new.row<-data.frame("site"=site,"tag"=sub.locs[i,"tag"],"X"=sub.locs[i,"X"],"Y"=sub.locs[i,"Y"],"x"=sub.locs[i,"x"],"y"=sub.locs[i,"y"],"date"= as.Date(date1),"status"=new.status,"max.height"=new.height,"inf.intens"=new.inf.intens) ### new data row
      pred.epi<-rbind(pred.epi,new.row) ### join data
      if(print.progress) {print(paste0("finished ",site," ",as.Date(date1)," tag ",i,"/",dim(sub.locs)[1]))}
    }
  }
  pred.epi
}

# run simulations
step.size<-7
site<-"GM"


if(any((!file.exists(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.0.site.",site,".step.size.",step.size,".RDS"))),
        (!file.exists(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.1.8.site.",site,".step.size.",step.size,".RDS"))),
         (!file.exists(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.3.7.site.",site,".step.size.",step.size,".RDS")))
)) {
  library(parallel)
  library(doParallel)
  library(foreach)
  library(doRNG)
  
  n.cores<-4
  registerDoParallel(n.cores)
  pred.epi.all.0<-foreach(k = 1:10, .multicombine = T, .options.RNG=2389572) %dorng% simulate.epi(site,0,seed=249685,step.size=step.size,print.progress = F)
  saveRDS(pred.epi.all.0,file=paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.0.site.",site,".step.size.",step.size,".RDS"))
  
  pred.epi.all.1.8<-foreach(k = 1:10, .multicombine = T, .options.RNG=2389572) %dorng% simulate.epi(site,1.8,seed=875345,step.size=step.size,print.progress = F)
  saveRDS(pred.epi.all.1.8,file=paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.1.8.site.",site,".step.size.",step.size,".RDS"))
  
  pred.epi.all.3.7<-foreach(k = 1:10, .multicombine = T, .options.RNG=2389572) %dorng% simulate.epi(site,3.7,seed=7465932,step.size=step.size,print.progress = F)
  saveRDS(pred.epi.all.3.7,file=paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.3.7.site.",site,".step.size.",step.size,".RDS"))
}

# plot simulations

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

plot.func<-function(site,step.size)
{
  
  pred.epi.all.0<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.0.site.",site,".step.size.",step.size,".RDS"))
  pred.epi.all.1.8<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.1.8.site.",site,".step.size.",step.size,".RDS"))
  pred.epi.all.3.7<-readRDS(paste0("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/pred.epi.all.3.7.site.",site,".step.size.",step.size,".RDS"))
  
  sub.epi<-corrected.epi[which(corrected.epi$Site==site),]
  sub.locs<-corrected.locs[which(corrected.locs$Site==site),]
  
  if(site=="GM")
  {
    plot(unique(pred.epi.all.0[[1]]$date),rep(0,times=length(unique(pred.epi.all.0[[1]]$date))),ylim=c(0,.4),xlab="date",ylab="prevalence",type="n",cex.axis=2,cex.lab=2)
  }
  
  if(site=="HM")
  {
    plot(unique(pred.epi.all.0[[1]]$date),rep(0,times=length(unique(pred.epi.all.0[[1]]$date))),ylim=c(0,.05),xlab="date",ylab="prevalence",type="n",cex.axis=2,cex.lab=2)
  }
  
  
  xvals<-c()
  yvals<-c()
  for(i in 1:length(unique(sub.epi$Date.First.Observed.Diseased)))
  {
    date<-unique(sub.epi$Date.First.Observed.Diseased)[i]
    xvals<-c(xvals,date)
    yvals<-c(yvals,nrow(sub.epi[which(sub.epi$Date.First.Observed.Diseased<=date),])/nrow(sub.locs))
  }
  points(xvals,yvals,type="l",col="black",lwd=4)
  
  
  for(k in 1:length(pred.epi.all.0))
  {
    xvals<-c()
    yvals<-c()
    pred.epi<-pred.epi.all.0[[k]]
    for(i in 1:length(unique(pred.epi.all.0[[1]]$date)))
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
    for(i in 1:length(unique(pred.epi.all.1.8[[1]]$date)))
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
    for(j in 1:length(pred.epi.all.1.8))
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
    for(i in 1:length(unique(pred.epi.all.3.7[[1]]$date)))
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
  legend("topleft",legend=c("data","+0 degrees C","+1.8 degrees C","+3.7 degrees C"),col=c("black","orange","red","purple"),lwd=4,cex=2,bty="n")
}

par(mfrow=c(1,2),mar=c(5,5,3,3))
plot.func("GM",2)
mtext("GM",cex = 2)
plot.func("HM",2)
mtext("HM",cex=2)



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

