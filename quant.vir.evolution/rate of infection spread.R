library(lubridate)
library(viridis)

setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")
within.host<-read.csv("Withinhost.csv")

#### within host data

### clean

within.host$Date<-as.Date(within.host$Date,tryFormats=c("%m/%d/%Y"))
column.index<-c(1:7,9,11,14,15,17)
within.host<-within.host[,column.index]
within.host[,"length.tissue.infected"]<-NA

## subset to relevant year

within.host<-subset(within.host,Year=="2019")
within.host<-droplevels(within.host)

## for stems with height recorded as '<5' set height to 5
within.host[which(within.host[,"stem.height"]=="<5"),"stem.height"]<-5
within.host<-transform(within.host,stem.height=as.numeric(as.character(stem.height)))

## for stems with height recorded as '<5' set height to 5
within.host[which(within.host[,"percent.tissue.infected"]=="<.05"),"percent.tissue.infected"]<-.05
within.host[which(within.host[,"percent.tissue.infected"]==""),"percent.tissue.infected"]<-NA
within.host<-transform(within.host,percent.tissue.infected=as.numeric(as.character(percent.tissue.infected)))

### restructure

tags<-c()
sites<-c()
dates<-c()
N.Stems<-c()
N.D.Stems<-c()
max.heights<-c()
tot.inf.metric.measures<-c()

for(tag in unique(within.host$Tag))
{
  sub.data<-subset(within.host,Tag==tag)
  
  for(date in unique(sub.data$Date)[order(unique(sub.data$Date))])
  {
    sub.data.week<-subset(sub.data,Date==date)
    
    tag.data<-as.character(sub.data.week[1,"Tag"])
    site<-as.character(sub.data.week[1,"Site"])
    date.data<-sub.data.week[1,"Date"]
    N.Stem<-sub.data.week[1,"N.Stems"]
    N.D.Stem<-sub.data.week[1,"N.D.Stems"]
    max.height<-sub.data.week[1,"max.height"]
    
    tags<-c(tags,tag.data)
    sites<-c(sites,site)
    dates<-c(dates,paste(date.data))
    N.Stems<-c(N.Stems,N.Stem)
    N.D.Stems<-c(N.D.Stems,N.D.Stem)
    max.heights<-c(max.heights,max.height)
    
    inf.metric.measured<-c()
    for(i in 1:dim(sub.data.week)[1])
    {
      if(is.na(sub.data.week[i,"length.tissue.infected"]) & !is.na(sub.data.week[i,"percent.tissue.infected"])) {sub.data.week[i,"length.tissue.infected"]<-sub.data.week[i,"percent.tissue.infected"]*sub.data.week[i,"stem.height"]}
      new.inf.metric<-sub.data.week[i,"length.tissue.infected"]*sub.data.week[i,"N.pustules.middle"]
      inf.metric.measured<-c(inf.metric.measured,new.inf.metric)
    }
    tot.inf.metric.measure<-c()
    if(!all(is.na(inf.metric.measured))) {tot.inf.metric.measure<-sub.data.week$N.D.Stems[1]*sum(inf.metric.measured,na.rm = T)/length(inf.metric.measured)} else{tot.inf.metric.measure<-NA}
    tot.inf.metric.measures<-c(tot.inf.metric.measures,tot.inf.metric.measure)    
  }
  
}

predictor.data<-data.frame(tag=as.character(tags),site=sites,date=as.Date(dates),N.Stems=N.Stems,N.D.Stems=N.D.Stems,p.D.stems=N.D.Stems/(N.D.Stems+N.Stems),max.height=max.heights,tot.inf.metric=tot.inf.metric.measures)

layout(matrix(c(1,1,5,5,9,9,2,2,6,6,10,10,3,4,7,8,11,12),3,6,byrow = T))
par(mar=c(5.1,4.1,1.1,2.1),oma=c(0,0,4,0))

### N.D.stems 
plot(0,0,xlim=c(min(predictor.data$date),max(predictor.data$date)),ylim=c(0,65),axes = F,main="all data",xlab="",ylab = "N D stems")
mtext("N stems infected",3,line=2,font=2)
axis(2)
axis.Date(1,predictor.data$date)
box()

for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  col<-viridis(length(unique(predictor.data$tag)))[which(unique(predictor.data$tag)==tag.index)]
  points(plot.data$date,plot.data$N.D.Stems,type="l",col=col)
}

plot(0,0,xlim=c(min(predictor.data$date),max(predictor.data$date)),ylim=c(0,65),axes = F,main="end points",xlab="",ylab = "N D stems")
axis(2)
axis.Date(1,predictor.data$date)
box()

for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  col<-viridis(length(unique(predictor.data$tag)))[which(unique(predictor.data$tag)==tag.index)]
  points(c(plot.data[1,"date"],plot.data[dim(plot.data)[1],"date"]),c(plot.data[1,"N.D.Stems"],plot.data[dim(plot.data)[1],"N.D.Stems"]),type="l",col=col)
}

slopes<-c()
for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  if(dim(plot.data)[1]>1)
  {
    slope<-coef(lm(N.D.Stems~date,data = plot.data))[2]
    slopes<-c(slopes,slope)
  } else{slopes<-c(slopes,NA)}
}
hist(na.omit(slopes),xlab="slope",main="fit to all data",breaks=10)
abline(v=0,col="red")

slopes<-c()
for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  if(dim(plot.data)[1]>1)
  {
    slope<-coef(lm(N.D.Stems~date,data = plot.data[c(min(which(!is.na(plot.data$N.D.Stems))),max(which(!is.na(plot.data$N.D.Stems)))),] ))[2]
    slopes<-c(slopes,slope)
  } else{slopes<-c(slopes,NA)}
}
hist(na.omit(slopes),xlab="slope",main="fit to endpoints",breaks=10)
abline(v=0,col="red")

### p.D.stems

plot(0,0,xlim=c(min(predictor.data$date),max(predictor.data$date)),ylim=c(0,.6),axes = F,main="all data",xlab="",ylab = "% stems infected")
mtext("% stems infected",3,line=2,font=2)
axis(2)
axis.Date(1,predictor.data$date)
box()

for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  col<-viridis(length(unique(predictor.data$tag)))[which(unique(predictor.data$tag)==tag.index)]
  points(plot.data$date,plot.data$p.D.stems,type="l",col=col)
}

plot(0,0,xlim=c(min(predictor.data$date),max(predictor.data$date)),ylim=c(0,.6),axes = F,main="end points",xlab="",ylab = "% stems infected")
axis(2)
axis.Date(1,predictor.data$date)
box()

for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  col<-viridis(length(unique(predictor.data$tag)))[which(unique(predictor.data$tag)==tag.index)]
  points(c(plot.data[1,"date"],plot.data[dim(plot.data)[1],"date"]),c(plot.data[1,"p.D.stems"],plot.data[dim(plot.data)[1],"p.D.stems"]),type="l",col=col)
}

slopes<-c()
for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  if(dim(plot.data)[1]>1)
  {
    slope<-coef(lm(p.D.stems~date,data = plot.data))[2]
    slopes<-c(slopes,slope)
  } else{slopes<-c(slopes,NA)}
}
hist(na.omit(slopes),xlab="slope",main="fit to all data",breaks=10)
abline(v=0,col="red")

slopes<-c()
for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  if(dim(plot.data)[1]>1)
  {
    slope<-coef(lm(p.D.stems~date,data = plot.data[c(min(which(!is.na(plot.data$p.D.stems))),max(which(!is.na(plot.data$p.D.stems)))),] ))[2]
    slopes<-c(slopes,slope)
  } else{slopes<-c(slopes,NA)}
}
hist(na.omit(slopes),xlab="slope",main="fit to endpoints",breaks=10)
abline(v=0,col="red")

### total infection load metric

plot(0,0,xlim=c(min(predictor.data$date),max(predictor.data$date)),ylim=c(0,32000),axes = F,main="all data",xlab="",ylab = "total infection load")
mtext("total infection load",3,line=2,font=2)
axis(2)
axis.Date(1,predictor.data$date)
box()

for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  col<-viridis(length(unique(predictor.data$tag)))[which(unique(predictor.data$tag)==tag.index)]
  points(plot.data$date,plot.data$tot.inf.metric,type="l",col=col)
}

plot(0,0,xlim=c(min(predictor.data$date),max(predictor.data$date)),ylim=c(0,32000),axes = F,main="end points",xlab="",ylab = "total infection load")
axis(2)
axis.Date(1,predictor.data$date)
box()

for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  col<-viridis(length(unique(predictor.data$tag)))[which(unique(predictor.data$tag)==tag.index)]
  points(c(plot.data[1,"date"],plot.data[dim(plot.data)[1],"date"]),c(plot.data[1,"tot.inf.metric"],plot.data[dim(plot.data)[1],"tot.inf.metric"]),type="l",col=col)
}

slopes<-c()
for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  if(dim(plot.data)[1]>1)
  {
    slope<-coef(lm(tot.inf.metric~date,data = plot.data))[2]
    slopes<-c(slopes,slope)
  } else{slopes<-c(slopes,NA)}
}
hist(na.omit(slopes),xlab="slope",main="fit to all data",breaks=10)
abline(v=0,col="red")

slopes<-c()
for(tag.index in unique(predictor.data$tag))
{
  plot.data<-subset(predictor.data,tag==tag.index)
  if(dim(plot.data)[1]>1)
  {
    slope<-coef(lm(tot.inf.metric~date,data = plot.data[c(min(which(!is.na(plot.data$tot.inf.metric))),max(which(!is.na(plot.data$tot.inf.metric)))),] ))[2]
    slopes<-c(slopes,slope)
  } else{slopes<-c(slopes,NA)}
}
hist(na.omit(slopes),xlab="slope",main="fit to endpoints",breaks=10)
abline(v=0,col="red")
