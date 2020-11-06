library(lubridate)

setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")
within.host<-read.csv("Withinhost.csv")


#### clean demog data

### pull relevant info out of demog
tag<-c()
status<-c()
death<-c()
year<-c()

for(index in unique(demog$tag))
{
  
  stat1<-NA
  stat2<-NA
  death1<-NA
  death2<-NA
 
  sub.data<-subset(demog,tag==index)
  if(all("2018" %in% sub.data$year,"2019" %in% sub.data$year))
  {
    stat1<-as.character(sub.data[which(sub.data$year=="2018"),"final.status"])
    if(sub.data[which(sub.data$year=="2019"),"status"] %in% c("H","D","X"))
    {
      death1<-ifelse(sub.data[which(sub.data$year=="2019"),"status"]=="X",1,0)
    } else {death1<-NA}
    
  }
  if(all("2019" %in% sub.data$year,"2020" %in% sub.data$year))
  {
    stat2<-as.character(sub.data[which(sub.data$year=="2019"),"final.status"])
    if(sub.data[which(sub.data$year=="2020"),"status"] %in% c("H","D","X"))
    {
      death2<-ifelse(sub.data[which(sub.data$year=="2020"),"status"]=="X",1,0)
    } else {death2<-NA}
  }
  
  if(!any(is.na(c(stat1,death1)))) {tag<-c(tag,index);status<-c(status,stat1);death<-c(death,death1);year<-c(year,"2018")}
  if(!any(is.na(c(stat2,death2)))) {tag<-c(tag,index);status<-c(status,stat2);death<-c(death,death2);year<-c(year,"2019")}            
  if(!(length(tag)==length(status))) {print(index)}
}

surv.data<-data.frame(tag=tag,year=year,status=status,death=death)
surv.data<-droplevels(surv.data)

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

### extract plant metrics and total infected tissue at last observation

tags<-c()
sites<-c()
N.Stems<-c()
N.D.Stems<-c()
max.heights<-c()
tot.inf.metric.measures<-c()

for(tag in unique(within.host$Tag))
{
  sub.data<-subset(within.host,Tag==tag)
  sub.data<-subset(sub.data,Date==max(unique(sub.data$Date)))
  
  tag<-as.character(sub.data[1,"Tag"])
  site<-as.character(sub.data[1,"Site"])
  N.Stem<-sub.data[1,"N.Stems"]
  N.D.Stem<-sub.data[1,"N.D.Stems"]
  max.height<-sub.data[1,"max.height"]
  
  tags<-c(tags,tag)
  sites<-c(sites,site)
  N.Stems<-c(N.Stems,N.Stem)
  N.D.Stems<-c(N.D.Stems,N.D.Stem)
  max.heights<-c(max.heights,max.height)
  
  inf.metric.measured<-c()
  for(i in 1:dim(sub.data)[1])
  {
    if(is.na(sub.data[i,"length.tissue.infected"]) & !is.na(sub.data[i,"percent.tissue.infected"])) {sub.data[i,"length.tissue.infected"]<-sub.data[i,"percent.tissue.infected"]*sub.data[i,"stem.height"]}
    new.inf.metric<-sub.data[i,"length.tissue.infected"]*sub.data[i,"N.pustules.middle"]
    inf.metric.measured<-c(inf.metric.measured,new.inf.metric)
  }
  tot.inf.metric.measure<-c()
  if(!all(is.na(inf.metric.measured))) {tot.inf.metric.measure<-sub.data$N.D.Stems[1]*sum(inf.metric.measured,na.rm = T)/length(inf.metric.measured)} else{tot.inf.metric.measure<-NA}
  tot.inf.metric.measures<-c(tot.inf.metric.measures,tot.inf.metric.measure)
}

predictor.data<-data.frame(tag=as.character(tags),site=sites,N.Stems=N.Stems,N.D.Stems=N.D.Stems,max.height=max.heights,tot.inf.metric=tot.inf.metric.measures)

#### analyze total tissue infected, including healhty plant data

### join data
sub.surv.data<-subset(surv.data,year==2019)

new.metrics<-c()
for(i in 1:dim(sub.surv.data)[1])
{
  if(sub.surv.data[i,"status"]=="H") {new.metric<-0}
  if(sub.surv.data[i,"status"]=="D") 
  {
    if(as.character(sub.surv.data[i,"tag"]) %in% as.character(predictor.data$tag))
    {
      new.metric<-predictor.data[which(as.character(predictor.data$tag)==as.character(sub.surv.data[i,"tag"])),"tot.inf.metric"]
    } else {new.metric<-NA}
  }
  new.metrics<-c(new.metrics,new.metric)
}

final.data1<-data.frame(sub.surv.data,"tot.inf.metric"=new.metrics)

### analysis
plot(final.data1$tot.inf.metric,final.data1$death,xlab="total infection metric",ylab="death")
summary(glm(death~tot.inf.metric,data=final.data1))

#### analyze by plant metrics, including healhty plant data

### join data

sub.surv.data<-subset(surv.data,year==2019)

death.metrics<-c()
for(i in 1:dim(predictor.data)[1])
{
  index<-which(as.character(sub.surv.data$tag)==as.character(predictor.data[i,"tag"]))
  death.metrics<-c(death.metrics,sub.surv.data[index,"death"])
  print(paste0("tag = ",as.character(predictor.data[i,"tag"])," death = ",sub.surv.data[index,"death"]))
}

final.data2<-data.frame(predictor.data,"death"=death.metrics)

