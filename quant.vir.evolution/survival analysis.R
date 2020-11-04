setwd("~/Documents/GitHub/flax.rust/data")
demog<-read.csv("Demography.csv")

demog<-subset(demog,status %in% c("H","D","X"))
demog<-subset(demog,final.status %in% c("H","D","X"))

tag<-c()
status<-c()
surv<-c()
height<-c()
ninfl<-c()

for(index in unique(demog$tag))
{
  
  stat1<-NA
  stat2<-NA
  surv1<-NA
  surv2<-NA
  height1<-NA
  height2<-NA
  ninfl1<-NA
  ninfl2<-NA
  sub.data<-subset(demog,tag==index)
  if(all("2018" %in% sub.data$year,"2019" %in% sub.data$year))
    {
      stat1<-as.character(sub.data[which(sub.data$year=="2018"),"final.status"])
      surv1<-ifelse(sub.data[which(sub.data$year=="2019"),"status"]=="X",0,1)
      height1<-as.numeric(sub.data[which(sub.data$year=="2018"),"height.cm"])
      ninfl1<-as.numeric(sub.data[which(sub.data$year=="2018"),"num.infl"])
    }
  if(all("2019" %in% sub.data$year,"2020" %in% sub.data$year))
    {
    stat2<-as.character(sub.data[which(sub.data$year=="2019"),"final.status"])
    surv2<-ifelse(sub.data[which(sub.data$year=="2020"),"status"]=="X",0,1)
    height2<-as.numeric(sub.data[which(sub.data$year=="2019"),"height.cm"])
    ninfl2<-as.numeric(sub.data[which(sub.data$year=="2019"),"num.infl"])
  }
  
  if(!any(is.na(c(stat1,surv1)))) {tag<-c(tag,index);status<-c(status,stat1);surv<-c(surv,surv1);height<-c(height,height1);ninfl<-c(ninfl,ninfl1)}
  if(!any(is.na(c(stat2,surv2)))) {tag<-c(tag,index);status<-c(status,stat2);surv<-c(surv,surv2);height<-c(height,height2);ninfl<-c(ninfl,ninfl2)}
  
  if(!(length(tag)==length(status))) {print(index)}
}

surv.data<-data.frame(tag=tag,status=status,surv=surv,height=height,ninfl=ninfl)
surv.data<-droplevels(surv.data)
plot(jitter(as.numeric(surv.data$status)),jitter(surv.data$surv),xlab="status",ylab="",axes=F)
axis(1,at=c(1,2),labels=c("diseased","healthy"))
axis(2,at=c(0,1),labels = c("died","survived"))

mod<-glm(surv~status,data=surv.data)
summary(mod)

plot(surv.data$ninfl,surv.data$surv)

mod1<-glm(surv~status*ninfl*height,data=surv.data)
mod2<-glm(surv~status+ninfl*height,data=surv.data)
mod3<-glm(surv~status*ninfl+height,data=surv.data)
mod4<-glm(surv~status*height+ninfl,data=surv.data)
mod5<-glm(surv~status+ninfl+height,data=surv.data)
mod6<-glm(surv~ninfl*height,data=surv.data)
mod7<-glm(surv~ninfl+height,data=surv.data)
mod8<-glm(surv~status*height,data=surv.data)
mod9<-glm(surv~status+height,data=surv.data)
mod10<-glm(surv~status*ninfl,data=surv.data)
mod11<-glm(surv~status+ninfl,data=surv.data)
mod12<-glm(surv~status,data=surv.data)
mod13<-glm(surv~ninfl,data=surv.data)
mod14<-glm(surv~height,data=surv.data)

AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14)

summary(mod)
