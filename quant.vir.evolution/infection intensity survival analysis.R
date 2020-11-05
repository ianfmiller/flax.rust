setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")
within.host<-read.csv("Withinhost.csv")

demog<-subset(demog,status %in% c("H","D","X"))
demog<-subset(demog,final.status %in% c("H","D","X"))

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