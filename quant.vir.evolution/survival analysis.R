setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")

demog<-subset(demog,status %in% c("H","D","X"))
demog<-subset(demog,final.status %in% c("H","D","X"))

tag<-c()
status<-c()
surv<-c()

for(index in unique(demog$tag))
{
  
  stat1<-NA
  stat2<-NA
  surv1<-NA
  surv2<-NA
  sub.data<-subset(demog,tag==index)
  if(all("2018" %in% sub.data$year,"2019" %in% sub.data$year))
    {
      stat1<-as.character(sub.data[which(sub.data$year=="2018"),"final.status"])
      surv1<-ifelse(sub.data[which(sub.data$year=="2019"),"status"]=="X",0,1)
    }
  if(all("2019" %in% sub.data$year,"2020" %in% sub.data$year))
    {
    stat2<-as.character(sub.data[which(sub.data$year=="2019"),"final.status"])
    surv2<-ifelse(sub.data[which(sub.data$year=="2020"),"status"]=="X",0,1)
  }
  
  if(!any(is.na(c(stat1,surv1)))) {tag<-c(tag,index);status<-c(status,stat1);surv<-c(surv,surv1)}
  if(!any(is.na(c(stat2,surv2)))) {tag<-c(tag,index);status<-c(status,stat2);surv<-c(surv,surv2)}
  
  if(!(length(tag)==length(status))) {print(index)}
}

surv.data<-data.frame(tag=tag,status=status,surv=surv)
surv.data<-droplevels(surv.data)
plot(jitter(as.numeric(surv.data$status)),jitter(surv.data$surv),xlab="status",ylab="",axes=F)
axis(1,at=c(1,2),labels=c("diseased","healthy"))
axis(2,at=c(0,1),labels = c("died","survived"))

mod<-glm(surv~status,data=surv.data)
