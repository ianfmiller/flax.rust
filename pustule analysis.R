library(mgcv)

source("prep.enviro.data.R")

# clean data

## load data
pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv")

## add in area
area<-pi*(pustules$max.diam..mm./4)*(pustules$min.diam..mm./4)
pustules$date<-as.Date(pustules$date,tryFormats = "%m/%d/%y")
pustules$area<-area

## clean data
pustules<-pustules[which(pustules$pustule.ID.confidence=="Yes"),]
pustules<-pustules[which(pustules$area>0),]

## make new data object for change in pustule size
temp.rh.sub.func<-function(x,lower.bound,upper.bound) {out<-subset(x,temp.c>=lower.bound); out<-subset(out,temp.c<=upper.bound); out}

tags<-c()
stem.iters<-c()
leaf.iters<-c()
pustule.nums<-c()
start.vals<-c()
end.vals<-c()
days<-c()
temp.days.16.22<-c()
temp.days.7.30<-c()
temp.days<-c()
mean.temp<-c()
dew.point.days<-c()
mean.dew.point<-c()
temp.dew.point.days<-c()
temp.16.22.dew.point.days<-c()
temp.7.30.dew.point.days<-c()
mean.wetnesss<-c()
tot.rains<-c()
mean.solars<-c()
mean.wind.speeds<-c()
mean.gust.speeds<-c()



for (tag in unique(pustules$tag)) #916, 917, 920
{
  sub.pustules1<-pustules[which(pustules$tag==tag),]
  
  for (color in unique(sub.pustules1$color))
  {
    sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
    {
      sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
      
      for(pustule.number in unique(sub.pustules3$pustule.number))
      {
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        if(dim(sub.pustules4)[1]>=2)
        {
          if(length(unique(sub.pustules4$date))<dim(sub.pustules4)[1]) #get rid of duplicate data resulting from two pics of same pustule on same date
          {
            for(date in unique(sub.pustules4$date))
            {
              date.indicies<-which(sub.pustules4$date==date)
              if(length(date.indicies)>1) {sub.pustules4<-sub.pustules4[-date.indicies[2:length(date.indicies)],]}
            }
          }
          for(i in 1:(dim(sub.pustules4)[1]-1))
          {
            #pull reference data
            date0<-sub.pustules4[i,"date"]
            date1<-sub.pustules4[i+1,"date"]
            site<-sub.pustules4[i,"site"]
            
            #subset temp data to relevant window
            temp.rh.sub<-subset(all.temp.rh,site==site) #pull out temp data for site
            temp.rh.sub<-subset(temp.rh.sub,date.time<=date1) #pull out relevant data
            temp.rh.sub<-subset(temp.rh.sub,date.time>=date0) #pull out relevant data
            temp.rh.sub<-subset(temp.rh.sub,!is.na(temp.c)) #throw out NAs
            temp.rh.sub<-cbind(temp.rh.sub,interval.length=c(diff(as.numeric(temp.rh.sub$date.time))/(60*60*24),NA)) #add interval length in days
            
            #subset weather data to relevant window
            weath.sub<-subset(all.weath,site==site)
            weath.sub<-subset(weath.sub,date<=date1) #pull out relevant data
            weath.sub<-subset(weath.sub,date>=date0) #pull out relevant data
            
            #calculate environmental variable metrics
            new.temp.days.16.22<-sum(temp.rh.sub.func(temp.rh.sub,16,22)$temp.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T) #temperature days for temp between 16 and 22 celsius
            new.temp.days.7.30<-sum(temp.rh.sub.func(temp.rh.sub,7,30)$temp.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T) #temperature days for temp between 7 and 30 celsius
            new.temp.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$interval.length,na.rm = T) #temperature days
            new.mean.temp<-mean(temp.rh.sub$temp.c,na.rm = T) #temperature days
            new.dew.point.days<-mean(temp.rh.sub$dew.pt.c,na.rm = T) #Dew point days
            new.mean.dew.point<-mean(temp.rh.sub$dew.pt.c,na.rm = T) #temperature days

            #calculate joint environmental variable metrics--accounts for temporal co-occurence of environmental variables
            new.temp.dew.point.days<-sum(temp.rh.sub$temp.c*temp.rh.sub$dew.pt.c*temp.rh.sub$interval.length,na.rm = T)
            new.temp.16.22.dew.point.days<-sum(temp.rh.sub.func(temp.rh.sub,16,22)$temp.c*temp.rh.sub.func(temp.rh.sub,16,22)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,16,22)$interval.length,na.rm = T)
            new.temp.7.30.dew.point.days<-sum(temp.rh.sub.func(temp.rh.sub,7,30)$temp.c*temp.rh.sub.func(temp.rh.sub,7,30)$dew.pt.c*temp.rh.sub.func(temp.rh.sub,7,30)$interval.length,na.rm = T)
            
            #calculate weather metrics
            new.mean.wetness<-mean(weath.sub$wetness)
            new.tot.rain<-sum(weath.sub$rain)
            new.mean.solar<-mean(weath.sub$solar.radiation)
            new.mean.wind.speed<-mean(weath.sub$wind.speed)
            new.mean.gust.speed<-mean(weath.sub$gust.speed)
            
            #pull out core predictors
            start.val<-sub.pustules4[i,"area"]
            end.val<-sub.pustules4[i+1,"area"]
            delta.days<-date1-date0
            
            #store values
            tags<-c(tags,tag)
            stem.iters<-c(stem.iters,color)
            leaf.iters<-c(leaf.iters,leaf.iteration)
            pustule.nums<-c(pustule.nums,pustule.number)
            
            start.vals<-c(start.vals,start.val)
            end.vals<-c(end.vals,end.val)
            days<-c(days,delta.days)
            
            temp.days.16.22<-c(temp.days.16.22,new.temp.days.16.22)
            temp.days.7.30<-c(temp.days.7.30,new.temp.days.7.30)
            temp.days<-c(temp.days,new.temp.days)
            mean.temp<-c(mean.temp,new.mean.temp)
            dew.point.days<-c(dew.point.days,new.dew.point.days)
            mean.dew.point<-c(mean.dew.point,new.mean.dew.point)
            temp.dew.point.days<-c(temp.dew.point.days,new.temp.dew.point.days)
            temp.16.22.dew.point.days<-c(temp.16.22.dew.point.days,new.temp.16.22.dew.point.days)
            temp.7.30.dew.point.days<-c(temp.7.30.dew.point.days,new.temp.7.30.dew.point.days)
            
            mean.wetnesss<-c(mean.wetnesss,new.mean.wetness)
            tot.rains<-c(tot.rains,new.tot.rain)
            mean.solars<-c(mean.solars,new.mean.solar)
            mean.wind.speeds<-c(mean.wind.speeds,new.mean.wind.speed)
            mean.gust.speeds<-c(mean.gust.speeds,new.mean.gust.speed)
            
            
          } 
        }
      }
    }
  }
}

delta.pustules<-data.frame(tag=factor(tags),stem.iter=stem.iters,leaf.iter=leaf.iters,pustule.num=pustule.nums,area=start.vals,area.next=end.vals,time=days,
                           temp.days.16.22=temp.days.16.22,temp.days.7.30=temp.days.7.30,temp.days=temp.days,mean.temp=mean.temp,dew.point.days=dew.point.days,mean.dew.point=mean.dew.point,temp.16.22.dew.point.days=temp.16.22.dew.point.days,temp.7.30.dew.point.days=temp.7.30.dew.point.days,temp.dew.point.days=temp.dew.point.days,
                           mean.wetness=mean.wetnesss,tot.rain=tot.rains,mean.solar=mean.solars,mean.wind.speed=mean.wind.speeds,mean.gust.speed=mean.gust.speeds)

# visualize data

## size histogram
hist(pustules$area,main="pustule area",breaks=100)

## plot trajectories
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="date",ylab="pustule area")
i<-0

plot.cols<-sample(rainbow(2042))

for (tag in unique(pustules$tag))
{
  sub.pustules1<-pustules[which(pustules$tag==tag),]
  
  for (color in unique(sub.pustules1$color))
  {
    sub.pustules2<-sub.pustules1[which(sub.pustules1$color==color),]
    
    for(leaf.iteration in unique(sub.pustules2$leaf.iteration)) 
    {
      sub.pustules3<-sub.pustules2[which(sub.pustules2$leaf.iteration==leaf.iteration),]
      
      for(pustule.number in unique(sub.pustules3$pustule.number))
      {
        i<-i+1
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        points(sub.pustules4$date,sub.pustules4$area,col=plot.cols[i],type="l",lwd=.5)
      }
    }
  }
}

## plot change
plot(delta.pustules$area,delta.pustules$area.next,col="grey")
abline(0,1)

#plot relationships
par(mfrow=c(3,2))
plot(delta.pustules$time,delta.pustules$area.next)
plot(delta.pustules$time,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$temp.days.16.22,delta.pustules$area.next)
plot(delta.pustules$temp.days.16.22,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$temp.days.7.30,delta.pustules$area.next)
plot(delta.pustules$temp.days.7.30,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$mean.temp,delta.pustules$area.next)
plot(delta.pustules$mean.temp,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$dew.point.days,delta.pustules$area.next)
plot(delta.pustules$dew.point.days,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$mean.dew.point ,delta.pustules$area.next)
plot(delta.pustules$mean.dew.point ,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$temp.16.22.dew.point.days,delta.pustules$area.next)
plot(delta.pustules$temp.16.22.dew.point.days,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$temp.7.30.dew.point.days,delta.pustules$area.next)
plot(delta.pustules$temp.7.30.dew.point.days,delta.pustules$area-delta.pustules$area.next)

plot(delta.pustules$temp.dew.point.days,delta.pustules$area.next)
plot(delta.pustules$temp.dew.point.days,delta.pustules$area-delta.pustules$area.next)

# analyze data
mod<-gam(area.next~s(area)+s(time)+s(temp.dew.point.days),data = delta.pustules)
