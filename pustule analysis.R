pustules<-read.csv("~/Documents/GitHub/flax.rust/data/pustule measurements.csv")
area<-pi*(pustules$max.diam..mm./4)*(pustules$min.diam..mm./4)
pustules$date<-as.Date(pustules$date,tryFormats = "%m/%d/%y")
pustules$area<-area

pustules<-pustules[which(pustules$pustule.ID.confidence=="Yes"),]
pustules<-pustules[which(pustules$area>0),]
hist(pustules$area)

### plot trajectories
plot(c(min(pustules$date),max(pustules$date)),c(0,max(pustules$area)),type="n",xlab="date",ylab="pustule area")

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
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        points(sub.pustules4$date,sub.pustules4$area,col="grey60",type="l",lwd=.5)
      }
    }
  }
}

### plot change
start.vals<-c()
change<-c()

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
        sub.pustules4<-sub.pustules3[which(sub.pustules3$pustule.number==pustule.number),]
        sub.pustules4<-sub.pustules4[order(sub.pustules4$date),]
        for(i in 1:(dim(sub.pustules4)[1]-1))
        {
          start.vals<-c(start.vals,sub.pustules4[i,"area"])
          change<-c(change,sub.pustules4[i+1,"area"]-sub.pustules4[i,"area"])
        }
      }
    }
  }
}

plot(start.vals,change,col="grey")
abline(0,1)

