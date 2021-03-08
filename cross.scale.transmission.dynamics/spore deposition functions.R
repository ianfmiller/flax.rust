# gaussian plume functions

vis<-F # set to F so plots aren't generated on source

## parameters

### q is quantity of spores--need to change to function of wind speed and tot. infection intensity
### k is constant realting infection intensity to spore availibility at source
### s is wind speed
### x is distance in direction of wind
### y is distance orthoganal to direction of wind
### alphay is sd of dispersal in y direction
### c is coefficient used to calculate decay constant along x as a function of s


## two D gaussian plume with decay along x a function of wind speed
decay.plume<-function(q,k,s,x,y,alphay,c)
{
  ifelse(x>=0,1,0)*q*k*exp(-( x^2/(2*(c*s)^2) + y^2/(2*alphay^2)))
}

if(vis)
{
  test.mat<-data.frame(x=rep(seq(-1,1,.1),each=21),y=rep(seq(-1,1,.1),times=21))
  out<-mapply(decay.plume, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=266.4167,k=.001,s=1,alphay=0.15906877,c=0.03382694))
  res.mat<-matrix(out,21,21,byrow = T)
  filled.contour(res.mat)
}

## function to find x and y coordinates to plug into tilted gaussian plume for a given wind direction and direction and distance of spore trap (assuming N/S/E/W cord system), assume counterclockwise is positive
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

correct.wind.degree<-function(x,site="blank")
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
