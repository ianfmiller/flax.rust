library(lubridate)

setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")
within.host<-read.csv("Withinhost.csv")


#### clean demog data

demog<-subset(demog,status %in% c("H","D","X"))
demog<-subset(demog,final.status %in% c("H","D","X"))


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
    death1<-ifelse(sub.data[which(sub.data$year=="2019"),"status"]=="X",1,0)
  }
  if(all("2019" %in% sub.data$year,"2020" %in% sub.data$year))
  {
    stat2<-as.character(sub.data[which(sub.data$year=="2019"),"final.status"])
    death2<-ifelse(sub.data[which(sub.data$year=="2020"),"status"]=="X",1,0)
  }
  
  if(!any(is.na(c(stat1,death1)))) {tag<-c(tag,index);status<-c(status,stat1);death<-c(death,death1);year<-c(year,"2018")}
  if(!any(is.na(c(stat2,death2)))) {tag<-c(tag,index);status<-c(status,stat2);death<-c(death,death2);year<-c(year,"2019")}            
  if(!(length(tag)==length(status))) {print(index)}
}

surv.data<-data.frame(tag=tag,year=year,status=status,death=death)
surv.data<-droplevels(surv.data)

#### within host data

### clean

within.host<-read.csv("Withinhost.csv")
within.host$Date<-as.Date(within.host$Date,tryFormats=c("%m/%d/%Y"))
column.index<-c(1:7,9,11,14,15,17)
within.host<-within.host[,column.index]
within.host[,"length.tissue.infected"]<-NA

## subset to relevant year

within.host<-subset(within.host,Year=="2019")
within.host<-droplevels(within.host)

## for stems with height recorded as '<5' set height to 5
within.host[which(within.host[,"stem.height"]=="<5"),"stem.height"]<-5
within.host<-transform(within.host,stem.height=as.numeric(stem.height))

## for stems with height recorded as '<5' set height to 5
within.host[which(within.host[,"percent.tissue.infected"]=="<.05"),"percent.tissue.infected"]<-.05
within.host[which(within.host[,"percent.tissue.infected"]==""),"percent.tissue.infected"]<-NA
within.host<-transform(within.host,percent.tissue.infected=as.numeric(as.character(percent.tissue.infected)))

### extract total infected tissue at last observation

tags<-c()

tot.inf.tissue.measures<-c()

for(tag in unique(within.host$Tag))
{
  sub.data<-subset(within.host,Tag==tag)
  sub.data<-subset(sub.data,Date==max(unique(sub.data$Date)))
  
  inf.tissue.measured<-c()
  for(i in 1:dim(sub.data)[1])
  {
    if(is.na(sub.data[i,"length.tissue.infected"]) & !is.na(sub.data[i,"percent.tissue.infected"])) {sub.data[i,"length.tissue.infected"]<-sub.data[i,"percent.tissue.infected"]*sub.data[i,"stem.height"]}
    new.inf.tissue.measured<-sub.data[i,"length.tissue.infected"]*sub.data[i,"N.pustules.middle"]
    inf.tissue.measured<-c(inf.tissue.measured,new.inf.tissue.measured)
  }
  if(!all(is.na(inf.tissue.measured))) {tot.inf.tissue.measure<-sub.data$N.D.Stems[1]*sum(inf.tissue.measured,na.rm = T)/length(inf.tissue.measured)} else{tot.inf.tissue.measure<-NA}
  tags<-c(tag,tags)
  tot.inf.tissue.measures<-c(tot.inf.tissue.measures,tot.inf.tissue.measure)
}

tot.inf.tissue.measures.frame<-data.frame(tag=tags,tot.inf.tissue=tot.inf.tissue.measures)

#### join survival data, tot inf.load
sub.surv.data<-subset(surv.data,year==2019)

new.metrics<-c()
for(i in 1:dim(sub.surv.data)[1])
{
  if(sub.surv.data[i,"status"]=="H") {new.metric<-0}
  if(sub.surv.data[i,"status"]=="D") 
  {
    if(as.character(sub.surv.data[i,"tag"]) %in% as.character(tot.inf.tissue.measures.frame$tag))
    {
      new.metric<-tot.inf.tissue.measures.frame[which(as.character(tot.inf.tissue.measures.frame$tag)==as.character(sub.surv.data[i,"tag"])),"tot.inf.tissue"]
    } else {new.metric<-NA}
  }
  new.metrics<-c(new.metrics,new.metric)
}

final.data<-data.frame(sub.surv.data,"tot.inf.tissue"=new.metrics)
plot(final.data$tot.inf.tissue,final.data$death,xlab="total infection load",ylab="death")
summary(glm(death~tot.inf.tissue,data=final.data))
