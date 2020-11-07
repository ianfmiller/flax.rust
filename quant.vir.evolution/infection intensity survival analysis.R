library(lubridate)
library(viridis)

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

#### analyze total tissue infected, including healhty plant data. Can't look at anything besides total infection metric because final stem#s for healthy plants weren't tracked.

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
  if(length(index)>0) {death.metrics<-c(death.metrics,sub.surv.data[index,"death"])} else {death.metrics<-c(death.metrics,NA)}
}

final.data2<-data.frame(predictor.data,"death"=death.metrics,"p.stems.inf"=predictor.data$N.D.Stems/predictor.data$N.Stems)

### fit models
## N.D.stems as predictor
mod1<-glm(death~N.D.Stems*N.Stems*max.height+site,data=final.data2)
mod2<-glm(death~N.D.Stems*N.Stems+max.height+site,data=final.data2)
mod3<-glm(death~N.D.Stems*N.Stems+site,data=final.data2)
mod4<-glm(death~N.D.Stems*max.height+N.Stems+site,data=final.data2)
mod5<-glm(death~N.D.Stems*max.height+site,data=final.data2)
mod6<-glm(death~N.D.Stems+N.Stems+max.height+site,data=final.data2)
mod7<-glm(death~N.D.Stems+max.height+site,data=final.data2)
mod8<-glm(death~N.D.Stems+N.Stems+site,data=final.data2)
mod9<-glm(death~N.D.Stems+site,data=final.data2)

mod10<-glm(death~N.D.Stems*N.Stems*max.height,data=final.data2)
mod11<-glm(death~N.D.Stems*N.Stems+max.height,data=final.data2)
mod12<-glm(death~N.D.Stems*N.Stems,data=final.data2)
mod13<-glm(death~N.D.Stems*max.height+N.Stems,data=final.data2)
mod14<-glm(death~N.D.Stems*max.height,data=final.data2)
mod15<-glm(death~N.D.Stems+N.Stems+max.height,data=final.data2)
mod16<-glm(death~N.D.Stems+max.height,data=final.data2)
mod17<-glm(death~N.D.Stems+N.Stems,data=final.data2)
mod18<-glm(death~N.D.Stems,data=final.data2)

AICs<-AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18)

final.stem.mod<-mod5

layout(matrix(c(1,1,1,1,1,2),1,6,byrow = T))
par(mar=c(5.1,4.1,4.1,2.1))
plot(final.data2$N.D.Stems,jitter(final.data2$death,amount=.05),xlab="N diseased stems",ylab="",axes=F)
axis(1)
axis(2,at=c(0,1),labels = c("survived","died"))
box()
ilogit<-function(x) {exp(x)/(exp(x)+1)}
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.05,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.05,na.rm = T)+coef(final.stem.mod)[4]),add=T,col=viridis(5,direction=-1)[1],lwd=2,lty=1)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.25,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.25,na.rm = T)+coef(final.stem.mod)[4]),add=T,col=viridis(5,direction=-1)[2],lwd=2,lty=1)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.5,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.5,na.rm = T)+coef(final.stem.mod)[4]),add=T,col=viridis(5,direction=-1)[3],lwd=2,lty=1)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.75,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.75,na.rm = T)+coef(final.stem.mod)[4]),add=T,col=viridis(5,direction=-1)[4],lwd=2,lty=1)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.95,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.95,na.rm = T)+coef(final.stem.mod)[4]),add=T,col=viridis(5,direction=-1)[5],lwd=2,lty=1)

curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.05,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.05,na.rm = T)),add=T,col=viridis(5,direction=-1)[1],lwd=2,lty=2)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.25,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.25,na.rm = T)),add=T,col=viridis(5,direction=-1)[2],lwd=2,lty=2)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.5,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.5,na.rm = T)),add=T,col=viridis(5,direction=-1)[3],lwd=2,lty=2)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.75,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.75,na.rm = T)),add=T,col=viridis(5,direction=-1)[4],lwd=2,lty=2)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.95,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.95,na.rm = T)),add=T,col=viridis(5,direction=-1)[5],lwd=2,lty=2)

curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.05,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.05,na.rm = T)+coef(final.stem.mod)[5]),add=T,col=viridis(5,direction=-1)[1],lwd=2,lty=3)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.25,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.25,na.rm = T)+coef(final.stem.mod)[5]),add=T,col=viridis(5,direction=-1)[2],lwd=2,lty=3)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.5,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.5,na.rm = T)+coef(final.stem.mod)[5]),add=T,col=viridis(5,direction=-1)[3],lwd=2,lty=3)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.75,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.75,na.rm = T)+coef(final.stem.mod)[5]),add=T,col=viridis(5,direction=-1)[4],lwd=2,lty=3)
curve(ilogit(coef(final.stem.mod)[1]+coef(final.stem.mod)[2]*x+coef(final.stem.mod)[3]*quantile(final.data2$max.height,.95,na.rm = T)+coef(final.stem.mod)[6]*x*quantile(final.data2$max.height,.95,na.rm = T)+coef(final.stem.mod)[5]),add=T,col=viridis(5,direction=-1)[5],lwd=2,lty=3)

par(mar=c(1,0,1,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",axes = F,xlab="",ylab="")
legend(.5,.4,legend=c("5%","25%","50%","75%","95%"),lwd=2,col=viridis(5,direction=-1),title="CC\nquantile height",xjust = .5,yjust = .5,cex=1,bty="n")
legend(.5,.6,legend=c("CC","BT","GM"),lwd=2,lty=c(1,2,3),col="black",title="Site",xjust = .5,yjust = .5,cex=1,bty="n")


## p.stems.inf as predictor
mod1<-glm(death~p.stems.inf*N.Stems*max.height+site,data=final.data2)
mod2<-glm(death~p.stems.inf*N.Stems+max.height+site,data=final.data2)
mod3<-glm(death~p.stems.inf*N.Stems+site,data=final.data2)
mod4<-glm(death~p.stems.inf*max.height+N.Stems+site,data=final.data2)
mod5<-glm(death~p.stems.inf*max.height+site,data=final.data2)
mod6<-glm(death~p.stems.inf+N.Stems+max.height+site,data=final.data2)
mod7<-glm(death~p.stems.inf+max.height+site,data=final.data2)
mod8<-glm(death~p.stems.inf+N.Stems+site,data=final.data2)
mod9<-glm(death~p.stems.inf+site,data=final.data2)

mod10<-glm(death~p.stems.inf*N.Stems*max.height,data=final.data2)
mod11<-glm(death~p.stems.inf*N.Stems+max.height,data=final.data2)
mod12<-glm(death~p.stems.inf*N.Stems,data=final.data2)
mod13<-glm(death~p.stems.inf*max.height+N.Stems,data=final.data2)
mod14<-glm(death~p.stems.inf*max.height,data=final.data2)
mod15<-glm(death~p.stems.inf+N.Stems+max.height,data=final.data2)
mod16<-glm(death~p.stems.inf+max.height,data=final.data2)
mod17<-glm(death~p.stems.inf+N.Stems,data=final.data2)
mod18<-glm(death~p.stems.inf,data=final.data2)

AICs<-AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18)

final.p.stem.mod<-mod7

## tot.inf.metric as predictor
mod1<-glm(death~tot.inf.metric*N.Stems*max.height+site,data=final.data2)
mod2<-glm(death~tot.inf.metric*N.Stems+max.height+site,data=final.data2)
mod3<-glm(death~tot.inf.metric*N.Stems+site,data=final.data2)
mod4<-glm(death~tot.inf.metric*max.height+N.Stems+site,data=final.data2)
mod5<-glm(death~tot.inf.metric*max.height+site,data=final.data2)
mod6<-glm(death~tot.inf.metric+N.Stems+max.height+site,data=final.data2)
mod7<-glm(death~tot.inf.metric+max.height+site,data=final.data2)
mod8<-glm(death~tot.inf.metric+N.Stems+site,data=final.data2)
mod9<-glm(death~tot.inf.metric+site,data=final.data2)

mod10<-glm(death~tot.inf.metric*N.Stems*max.height,data=final.data2)
mod11<-glm(death~tot.inf.metric*N.Stems+max.height,data=final.data2)
mod12<-glm(death~tot.inf.metric*N.Stems,data=final.data2)
mod13<-glm(death~tot.inf.metric*max.height+N.Stems,data=final.data2)
mod14<-glm(death~tot.inf.metric*max.height,data=final.data2)
mod15<-glm(death~tot.inf.metric+N.Stems+max.height,data=final.data2)
mod16<-glm(death~tot.inf.metric+max.height,data=final.data2)
mod17<-glm(death~tot.inf.metric+N.Stems,data=final.data2)
mod18<-glm(death~tot.inf.metric,data=final.data2)

AICs<-AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18)

final.metric.mod<-mod14

layout(matrix(c(1,1,1,1,1,2),1,6,byrow = T))
par(mar=c(5.1,4.1,4.1,2.1))
plot(final.data2$tot.inf.metric,jitter(final.data2$death,amount=.05),xlab="total infection load metric",ylab="",axes = F)
axis(1)
axis(2,at=c(0,1),labels = c("survived","died"))
box()
ilogit<-function(x) {exp(x)/(exp(x)+1)}
curve(ilogit(coef(final.metric.mod)[1]+coef(final.metric.mod)[2]*x+coef(final.metric.mod)[3]*quantile(final.data2$max.height,.05,na.rm = T)+coef(final.metric.mod)[4]*x*quantile(final.data2$max.height,.05,na.rm = T)),add=T,col=viridis(5,direction=-1)[1],lwd=2)
curve(ilogit(coef(final.metric.mod)[1]+coef(final.metric.mod)[2]*x+coef(final.metric.mod)[3]*quantile(final.data2$max.height,.25,na.rm = T)+coef(final.metric.mod)[4]*x*quantile(final.data2$max.height,.25,na.rm = T)),add=T,col=viridis(5,direction=-1)[2],lwd=2)
curve(ilogit(coef(final.metric.mod)[1]+coef(final.metric.mod)[2]*x+coef(final.metric.mod)[3]*quantile(final.data2$max.height,.5,na.rm = T)+coef(final.metric.mod)[4]*x*quantile(final.data2$max.height,.5,na.rm = T)),add=T,col=viridis(5,direction=-1)[3],lwd=2)
curve(ilogit(coef(final.metric.mod)[1]+coef(final.metric.mod)[2]*x+coef(final.metric.mod)[3]*quantile(final.data2$max.height,.75,na.rm = T)+coef(final.metric.mod)[4]*x*quantile(final.data2$max.height,.75,na.rm = T)),add=T,col=viridis(5,direction=-1)[4],lwd=2)
curve(ilogit(coef(final.metric.mod)[1]+coef(final.metric.mod)[2]*x+coef(final.metric.mod)[3]*quantile(final.data2$max.height,.95,na.rm = T)+coef(final.metric.mod)[4]*x*quantile(final.data2$max.height,.95,na.rm = T)),add=T,col=viridis(5,direction=-1)[5],lwd=2)
par(mar=c(1,0,1,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",axes = F,xlab="",ylab="")
legend(.5,.5,legend=c("5%","25%","50%","75%","95%"),lwd=2,col=viridis(5,direction=-1),title="quantile height",xjust = .5,yjust = .5,cex=1)

