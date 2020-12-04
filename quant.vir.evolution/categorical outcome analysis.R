library(lubridate)
library(viridis)
library(MASS)

setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")
within.host<-read.csv("Withinhost.csv")


#### clean demog data

### pull relevant info out of demog
tag<-c()
site<-c()
status<-c()
nextstat<-c()
year<-c()

for(index in unique(demog$tag))
{
  
  stat1<-NA
  stat2<-NA
  nextstat1<-NA
  nextstat2<-NA
  site1<-NA
  site2<-NA
  
  sub.data<-subset(demog,tag==index)
  if(all("2018" %in% sub.data$year,"2019" %in% sub.data$year))
  {
    stat1<-as.character(sub.data[which(sub.data$year=="2018"),"final.status"])
    site1<-as.character(sub.data[which(sub.data$year=="2018"),"Site"])
    if(sub.data[which(sub.data$year=="2019"),"status"] %in% c("H","D","X"))
    {
      nextstat1<-sub.data[which(sub.data$year=="2019"),"status"]
    } else {nextstat1<-NA}
    
  }
  if(all("2019" %in% sub.data$year,"2020" %in% sub.data$year))
  {
    stat2<-as.character(sub.data[which(sub.data$year=="2019"),"final.status"])
    site2<-as.character(sub.data[which(sub.data$year=="2019"),"Site"])
    if(sub.data[which(sub.data$year=="2020"),"status"] %in% c("H","D","X"))
    {
      nextstat2<-sub.data[which(sub.data$year=="2020"),"status"]
    } else {nextstat2<-NA}
  }
  
  if(!any(is.na(c(stat1,nextstat1)))) {tag<-c(tag,index);status<-c(status,stat1);nextstat<-c(nextstat,nextstat1);year<-c(year,"2018");site<-c(site,site1)}
  if(!any(is.na(c(stat2,nextstat2)))) {tag<-c(tag,index);status<-c(status,stat2);nextstat<-c(nextstat,nextstat2);year<-c(year,"2019");site<-c(site,site2)}            
  if(!(length(tag)==length(status))) {print(index)}
}

surv.data<-data.frame(tag=tag,year=year,site=site,status=status,nextstat=nextstat)
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
dates<-c()
N.Stems<-c()
N.D.Stems<-c()
max.heights<-c()
tot.inf.metrics<-c()

for(tag in unique(within.host$Tag))
{
  sub.data<-subset(within.host,Tag==tag)
  
  tag.data<-as.character(sub.data[1,"Tag"])
  site<-as.character(sub.data[1,"Site"])
  date.data<-sub.data[1,"Date"]
  
  tags<-c(tags,tag.data)
  sites<-c(sites,site)
  dates<-c(dates,paste(date.data))
  
  tot.inf.metric.measures.week<-c()
  N.Stems.week<-c()
  N.D.Stems.week<-c()
  max.heights.week<-c()
  
  for(date in unique(sub.data$Date)[order(unique(sub.data$Date))])
  {
    sub.data.week<-subset(sub.data,Date==date)
    
    inf.metric.measured<-c()
    for(i in 1:dim(sub.data.week)[1])
    {
      if(is.na(sub.data.week[i,"length.tissue.infected"]) & !is.na(sub.data.week[i,"percent.tissue.infected"])) {sub.data.week[i,"length.tissue.infected"]<-sub.data.week[i,"percent.tissue.infected"]*sub.data.week[i,"stem.height"]}
      new.inf.metric<-sub.data.week[i,"length.tissue.infected"]*sub.data.week[i,"N.pustules.middle"]
      inf.metric.measured<-c(inf.metric.measured,new.inf.metric)
    }
    tot.inf.metric.measure.week<-c()
    if(!all(is.na(inf.metric.measured))) {tot.inf.metric.measure.week<-sub.data.week$N.D.Stems[1]*sum(inf.metric.measured,na.rm = T)/length(inf.metric.measured)} else{tot.inf.metric.measure.week<-NA}
    tot.inf.metric.measures.week<-c(tot.inf.metric.measures.week,tot.inf.metric.measure.week)    
    
    N.Stems.week<-c(N.Stems.week,sub.data.week[1,"N.Stems"])
    N.D.Stems.week<-c(N.D.Stems.week,sub.data.week[1,"N.D.Stems"])
    max.heights.week<-c(max.heights.week,sub.data.week[1,"max.height"])
    
  }
  
  if(any(!is.na(tot.inf.metric.measures.week))) {tot.inf.metrics<-c(tot.inf.metrics,tot.inf.metric.measures.week[max(which(!is.na(tot.inf.metric.measures.week)))])} else {tot.inf.metrics<-c(tot.inf.metrics,NA)}
  if(any(!is.na(N.Stems.week))) {N.Stems<-c(N.Stems,N.Stems.week[max(which(!is.na(N.D.Stems.week)))])} else {N.Stems<-c(N.Stems,NA)} #get N.stems and N.D.Stems from same time point
  if(any(!is.na(N.D.Stems.week))) {N.D.Stems<-c(N.D.Stems,N.D.Stems.week[max(which(!is.na(N.D.Stems.week)))])} else {N.D.Stems<-c(N.D.Stems,NA)} #get N.stems and N.D.Stems from same time point
  if(any(!is.na(max.heights.week))) {max.heights<-c(max.heights,max.heights.week[max(which(!is.na(max.heights.week)))])} else {max.heights<-c(max.heights,NA)}
}

predictor.data<-data.frame(tag=as.character(tags),site=sites,date=as.Date(dates),N.Stems=N.Stems,N.D.Stems=N.D.Stems,p.D.stems=N.D.Stems/N.Stems,max.height=max.heights,tot.inf.metric=tot.inf.metrics)

### join data

sub.surv.data<-subset(surv.data,year==2019)

nextstat.metrics<-c()
for(i in 1:dim(predictor.data)[1])
{
  index<-which(as.character(sub.surv.data$tag)==as.character(predictor.data[i,"tag"]))
  if(length(index)>0) {nextstat.metrics<-c(nextstat.metrics,sub.surv.data[index,"nextstat"])} else {nextstat.metrics<-c(nextstat.metrics,NA)}
}

nextstat.metrics<-factor(nextstat.metrics,ordered = T,levels=c("H","D","X"))

final.data<-data.frame(predictor.data,"nextstat"=nextstat.metrics)

#ordinal logistic regression
mod1<-polr(nextstat~p.D.stems*max.height*N.Stems,data=final.data)
mod2<-polr(nextstat~p.D.stems+max.height*N.Stems,data=final.data)
mod3<-polr(nextstat~p.D.stems*max.height+N.Stems,data=final.data)
mod4<-polr(nextstat~p.D.stems*N.Stems+max.height,data=final.data)
mod5<-polr(nextstat~p.D.stems+N.Stems+max.height,data=final.data)
mod6<-polr(nextstat~p.D.stems*N.Stems,data=final.data)
mod7<-polr(nextstat~p.D.stems+N.Stems,data=final.data)
mod8<-polr(nextstat~p.D.stems*max.height,data=final.data)
mod9<-polr(nextstat~p.D.stems+max.height,data=final.data)
mod10<-polr(nextstat~N.Stems*max.height,data=final.data)
mod11<-polr(nextstat~N.Stems+max.height,data=final.data)
mod12<-polr(nextstat~p.D.stems,data=final.data)
mod13<-polr(nextstat~N.Stems,data=final.data)
mod14<-polr(nextstat~max.height,data=final.data)
AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14)

##get p values
ctable <- coef(summary(mod9))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = p)
confint.default(mod1)
