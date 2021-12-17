# gaussian plume functions

library(rlist)
library(parallel)

## parameters

### I is quantity of spores--need to change to function of wind speed and tot. infection intensity
### k is constant realting infection intensity to spore availibility at source
### H is height from which spores are disperesed--assume 0.5*plant height
### s is wind speed
### x is distance in direction of wind
### y is distance orthoganal to direction of wind
### sigma is sd of dispersal in y and z directions


## tilted gaussian plume


tilted.plume<-function(I,H,k,Ws,A,s,x,y)
{
  if(s==0) {s<-.33/2}
  sigma<-2*A*x/s
  out<-((I*k*Ws)/(2*pi*s*sigma^2))*exp(((-y^2)/(2*sigma^2))-((H-Ws*x/s)^2/(2*sigma^2)))
  out
}

## function to find x and y coordinates to plug into gaussian plume for a given wind direction and direction and distance of spore trap (assuming N/S/E/W cord system), assume counterclockwise is positive
get.plume.xy<-function(degree,xorigin,yorigin,xtarget,ytarget,plot=F)
{
  # find point on wind vector that when connected with x,y forms a right angle with wind vector
  # https://stackoverflow.com/questions/6630596/find-a-line-intersecting-a-known-line-at-right-angle-given-a-point
  x1<-xorigin
  x2<-x1+sin(degree)
  y1<-yorigin
  y2<-y1+cos(degree)
  t<-((xtarget-x1)*(x2-x1)+(ytarget-y1)*(y2-y1))/((x2-x1)^2+(y2-y1)^2)
  xp<-x1+t*(x2-x1)
  yp<-y1+t*(y2-y1)
  
  # line connecting x1,y1 and x2,y2 is y=((y2-y1)/(x2-x1))*x-((y2-y1)/(x2-x1))*x1+y1
  # line passing throughh x1,y1 and perpendicular to line connecting x1,y1 and x2,y2 is s=-1*((y2-y1)/(x2-x1))^-1; y=s*x+y1+s*x1
  
  s=-1*((y2-y1)/(x2-x1))^-1 #slope of perpendicular line
  
  
  if(y2>s*x2+y1+s*x1)
  {
    if(yp<(s*xp+y1+s*x1))
    {
      newx<--1*((xp-x1)^2+(yp-y1)^2)^.5
      c<-((xtarget-x1)^2+(ytarget-y1)^2)^.5
      newy<-((c^2-newx^2))^.5
    } else
    {
      newx<-((xp-x1)^2+(yp-y1)^2)^.5
      c<-((xtarget-x1)^2+(ytarget-y1)^2)^.5
      newy<-((c^2-newx^2))^.5
    }
  } else {
    if(yp>(s*xp+y1+s*x1))
    {
      newx<--1*((xp-x1)^2+(yp-y1)^2)^.5
      c<-((xtarget-x1)^2+(ytarget-y1)^2)^.5
      newy<-((c^2-newx^2))^.5
    } else
    {
      newx<-((xp-x1)^2+(yp-y1)^2)^.5
      c<-((xtarget-x1)^2+(ytarget-y1)^2)^.5
      newy<-((c^2-newx^2))^.5
    }
  }
  
  
  if(plot)
  {
    plot(0,0,type="n",xlim=c(min(x1,x2,xp,xtarget)-.5,max(x1,x2,xp,xtarget)+.5),ylim=c(min(y1,y2,yp,ytarget,-1.5),max(y1,y2,yp,ytarget,1.5)),asp=1,xlab="plot X-axis",ylab="plot Y-axis")
    points(x1,y1,col="darkgreen",pch=16,cex=4)
    arrows(x1,y1,x2,y2)
    curve(((y2-y1)/(x2-x1))*x-((y2-y1)/(x2-x1))*x1+y1,add=T,lty=2)
    curve(s*x+y1+s*x1,col="red",lty=2,add=T)
    points(xtarget,ytarget,col="lightgreen",pch=16,cex=4)
    #points(xp,yp,col="lightgreen",cex=4)
    arrows(x1,y1,xp,yp,col="lightblue")
    arrows(xp,yp,xtarget,ytarget,col="darkblue")
    legend("topright",legend=c("source plant","target plant","wind","perpendicular","x","y"),pch=c(16,16,NA,NA,NA,NA),lty=c(NA,NA,1,2,1,1),col=c("darkgreen","lightgreen","black","red","lightblue","darkblue"))
  }
  return(c(newx,newy))
}

## function to correct wind direction to flip from direction of origin to direction of flow; correct for transect direction

correct.wind.degree.func<-function(x,site="blank")
{
  ### correct wind direction so that 360/0 corresponds to plot "UP" direction
  if(site=="GM")
  {
    newx<-x-309 #measured in google earth pro
    if(newx<0) {newx<-360+newx}
  }
  
  if(site=="HM")
  {
    newx<-x-40
    if(newx<0) {newx<-360+newx}
  }
  if(site=="BT")
  {
    newx<-x-75
    if(newx<0) {newx<-360+newx}
  }
  
  if(site=="CC")
  {
    newx<-x-332
    if(newx<0) {newx<-360+newx}
  }
  
  if(site=="blank") {newx<-x}
  
  newx<-newx-180
  if(newx<0) {newx<-360+newx}
  newx
}

correct.wind.degree<-function(x,site="blank")
{
  sapply(x,correct.wind.degree.func,site=site)
}

#plot(0,0,type="n",xlim=c(-max(wind.data$wind.speed),max(wind.data$wind.speed)),ylim=c(-max(wind.data$wind.speed),max(wind.data$wind.speed)))
#arrows(0,0,wind.data$wind.speed*cos(2*pi*wind.data$wind.direction/360),wind.data$wind.speed*sin(2*pi*wind.data$wind.direction/360))
#arrows(0,0,wind.data$wind.speed*cos(2*pi*correct.wind.degree(wind.data$wind.direction,site=site)/360),wind.data$wind.speed*sin(2*pi*correct.wind.degree(wind.data$wind.direction,site=site)/360),col="blue")
#points(0,0,col="red",pch=15)

predict.kernel.tilted.plume.inst<-function(i,I,H,k,Ws,A,xtarget,ytarget,wind.data,site)
{
  delta.t<-wind.data[i+1,"date"]-wind.data[i,"date"]
  cords<-mapply(get.plume.xy,2*pi*correct.wind.degree(wind.data[i,"wind.direction"],site = site)/360,MoreArgs=list(xorigin=0,yorigin=0,xtarget=xtarget,ytarget=ytarget))
  tilted.plume(I=I,H=H,k=k,Ws=Ws,A=A,s=wind.data[i,"wind.speed"],x=cords[1,],y=cords[2,])
}

predict.kernel.tilted.plume<-function(I,H,k,Ws,A,xtarget,ytarget,wind.data)
{
  suppressWarnings(predict.kernel.tilted.plume.inst(1:(dim(wind.data)[1]-1),I=I,H=H,k=k,Ws=Ws,A=A,xtarget=xtarget,ytarget=ytarget,wind.data=wind.data,site=wind.data[1,"site"]))->tot.dep
  sum(tot.dep,na.rm = T)
}

param.search.optim.tilted.plume<-function(x,return.out=F)
{
  
  out<-list.rbind(mcmapply(param.search.optim.tilted.plume.tag,unique(spore.deposition$Tag),MoreArgs = list(kval=x[1],Wsval=x[2],Aval=x[3]),SIMPLIFY=F,mc.cores = 6))
  if(return.out==T) {out} else {sum(out$val)}
  
}

param.search.optim.tilted.plume.tag<-function(tag,kval,Wsval,Aval)
{
  tags<-c()
  dists<-c()
  directions<-c()
  Is<-c()
  Hs<-c()
  val<-c()
  preds<-c()
  obs<-c()
  
  site<-spore.deposition[which(spore.deposition$Tag==tag)[1],"Site"]
  plantx<-demog[which(demog$tag==tag),"X"]+demog[which(demog$tag==tag),"x"]
  planty<-demog[which(demog$tag==tag),"Y"]+demog[which(demog$tag==tag),"y"]
  sub.1.spore.deposition<-spore.deposition[which(spore.deposition$Tag==tag),]
  
  for(date in unique(sub.1.spore.deposition$Date.collected))
  {
    sub.2.spore.deposition<-sub.1.spore.deposition[which(sub.1.spore.deposition$Date.collected==date),]
    deploy.date<-sub.2.spore.deposition[1,"Date.deployed"]
    I<-sub.2.spore.deposition[1,"plant.inf.intens"]
    
    wind.data<-all.weath[which(all.weath$site==site),]
    wind.data<-wind.data[which(wind.data$date>(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),]
    wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC")+60*60*24*1)),] ### fit to one day post spore trap deploy
    #wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(deploy.date,"%m/%d/%y")," 12:00:00"),tz="UTC")+60*60*24*2)),] ### fit to two days post spore trap deploy
    #wind.data<-wind.data[which(wind.data$date<=(as.POSIXct(paste0(as.Date(date,"%m/%d/%y")," 12:00:00"),tz="UTC"))),] ### fit to full spore trap deploy period
    
    for(j in 1:dim(sub.2.spore.deposition)[1])
    {
      xtarget<-0
      ytarget<-0
      
      H<-.5*sub.2.spore.deposition[1,"Height.cm"]/100
      
      if(sub.2.spore.deposition[j,"Direction"]=="U") {ytarget<-as.numeric(sub.2.spore.deposition[j,"Distance.cm"])/100}
      if(sub.2.spore.deposition[j,"Direction"]=="R") {xtarget<-as.numeric(sub.2.spore.deposition[j,"Distance.cm"])/100}
      if(sub.2.spore.deposition[j,"Direction"]=="D") {ytarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance.cm"])/100}
      if(sub.2.spore.deposition[j,"Direction"]=="L") {xtarget<-(-1)*as.numeric(sub.2.spore.deposition[j,"Distance.cm"])/100}
      
      new.pred<-predict.kernel.tilted.plume(I=I,H=H,k=kval,Ws=Wsval,A=Aval,xtarget=xtarget,ytarget=ytarget,wind.data=wind.data)
      if(Wsval>5) {new.pred<- -888}
      new.obs<-sub.2.spore.deposition[j,"spores.per.square.mm"]
      tags<-c(tags,tag)
      dists<-c(dists,as.numeric(sub.2.spore.deposition[j,"Distance.cm"]))
      directions<-c(directions,sub.2.spore.deposition[j,"Direction"])
      Is<-c(Is,I)
      Hs<-c(Hs,H)
      val<-c(val,(new.obs-new.pred)^2)
      preds<-c(preds,new.pred)
      obs<-c(obs,new.obs)
    }
  }
  return(data.frame("tag"=tags,"dist"=dists,"direction"=directions,"plant.inf.intens"=Is,"height.cm"=Hs,"val"=val,"pred"=preds,"obs"=obs))
}
