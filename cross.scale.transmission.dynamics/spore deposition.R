spores<-read.csv("~/Documents/GitHub/flax.rust/data/spore counts.csv")
spores$Distance..cm.<-as.numeric(spores$Distance..cm.)
spores<-spores[-which(spores$Tag %in% c("?","8,11")),]
spores<-spores[-which(spores$Direction=="?"),]

source('~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant data prep.R')

get.x.direction<-function(x) {
  if(x=="D") {out<-0}
  if(x=="R") {out<-1}
  if(x=="U") {out<-0}
  if(x=="L") {out<--1}
  out
}

get.y.direction<-function(x) {
  if(x=="D") {out<- -1}
  if(x=="R") {out<-0}
  if(x=="U") {out<-1}
  if(x=="L") {out<-0}
  out
}

mean.zero.depos<-c()
tags<-c()
dates<-c()
n.stems<-c()
n.d.stems<-c()
max.heights<-c()
plant.inf.intens<-c()

par(mfrow=c(5,5))
for(tag in unique(spores$Tag))
{
  sub.spores.1<-spores[which(spores$Tag==tag),]
  for(date in unique(sub.spores.1$Date.collected))
  {
    sub.spores.2<-sub.spores.1[which(sub.spores.1$Date.collected==date),]
    plot(0,0,type="n",xlim=c(-100,100),ylim=c(-100,100),xlab="",ylab="")
    points(unlist(lapply(sub.spores.2$Direction,get.x.direction))*sub.spores.2$Distance..cm.,unlist(lapply(sub.spores.2$Direction,get.y.direction))*sub.spores.2$Distance..cm.,cex=sub.spores.2$Pustules/100)
    
    zeros<-sub.spores.2[which(sub.spores.2$Distance..cm.==0),]
    mean.zero.depo<-mean(zeros$Pustules/zeros$X..squares.counted)
    
    sub.plants<-plants[which(plants$Tag==tag),]
    sub.plants<-sub.plants[which(sub.plants$Date==as.Date(date,tryFormats = "%m/%d/%y")),]
    
    if(dim(sub.plants)[1]>0)
    {
      new.n.stems<-sub.plants$N.Stems
      new.n.d.stems<-sub.plants$N.D.Stems
      new.max.height<-sub.plants$max.height
      new.plant.inf.intens<-sub.plants$plant.inf.intens
    } else {
      new.n.stems<-NA
      new.n.d.stems<-NA
      new.max.height<-NA
      new.plant.inf.intens<-NA
    }

    mean.zero.depos<-c(mean.zero.depos,mean.zero.depo)
    tags<-c(tags,tag)
    dates<-c(dates,date)
    n.stems<-c(n.stems,new.n.stems)
    n.d.stems<-c(n.d.stems,new.n.d.stems)
    max.heights<-c(max.heights,new.max.height)
    plant.inf.intens<-c(plant.inf.intens,new.plant.inf.intens)
   }
}

spore.dep<-data.frame(tag=tags,date=dates,mean.zero.depo=mean.zero.depos,n.stems=n.stems,n.d.stems=n.d.stems,max.height=max.heights,plant.inf.intens=plant.inf.intens)
spore.dep$date<-as.Date(spore.dep$date,tryFormats = "%m/%d/%y")

par(mfrow=c(1,2))
plot(spore.dep$n.d.stems,spore.dep$mean.zero.depo)
plot(spore.dep$plant.inf.intens,spore.dep$mean.zero.depo)
mod1<-lmer(mean.zero.depo~n.d.stems+(1|tag),data = spore.dep,REML = F)
mod2<-lmer(mean.zero.depo~plant.inf.intens+(1|tag),data = spore.dep,REML = F)

plot(spore.dep$n.d.stems,log10(spore.dep$mean.zero.depo*100+1))
plot(spore.dep$plant.inf.intens,log10(spore.dep$mean.zero.depo*100+1))
mod3<-lmer(log10(mean.zero.depo*100+1)~n.d.stems+(1|tag),data=spore.dep,REML=F)
mod4<-lmer(log10(mean.zero.depo*100+1)~plant.inf.intens+(1|tag),data=spore.dep,REML=F)


