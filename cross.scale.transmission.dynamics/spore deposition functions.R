# gaussian plume functions

vis<-F # set to F so plots aren't generated on source

## parameters

### q is quantity of spores--need to change to function of wind speed and tot. infection intensity
### H is height from which spores are disperesed--assume 0.5*plant height
### s is wind speed
### x is distance in direction of wind
### y is distance orthoganal to direction of wind
### alphay is sd of dispersal in y direction
### alphaz is sd of dispersal in z direction (vertical)
### Ws is falling velocity of spores
### c is coefficient used to calculate decay constant along x as a function of s

## tilted gaussian plume https://www.jstor.org/stable/pdf/1937537.pdf?casa_token=wXe8JWI0TrIAAAAA:g6XJxzyOWPEPZd_yw_lMiam3aO9QRO6Usa_WYPKY7zNdRu3PGeXnCoeGv2tazU3CBimc0zfXEiHI39brv7GhMlG7zuUYAZGF2DrFWnH4GRIHul6m4Trcodel
tilted.plume<-function(q,H,s,x,y,alphay,alphaz)
{
  Ws<-100 #terminal downward velocity of spores
  ifelse(x>=0,1,0)*((q*Ws)/(2*pi*s*alphay*alphaz))*
    exp(
    (-(y^2)/(2*alphay^2))+
      (-((0-(H-Ws*x/s))^2)/(2*alphaz^2))
    )
}

## visualize
if(vis)
{
  test.mat<-data.frame(x=rep(seq(-100,100,1),each=201),y=rep(seq(-100,100,1),times=201))
  out<-mapply(D, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=1,H=5,s=500,alphay=10,alphaz=10))
  res.mat<-matrix(out,201,201,byrow = T)
  filled.contour(res.mat)
}

## two D gaussian plume with decay along x a function of wind speed
decay.plume<-function(q,s,x,y,alphay,c)
{
  ifelse(x>=0,1,0)*q*exp(-( x^2/(2*(c*s)^2) + y^2/(2*alphay^2)))
}

if(vis)
{
  test.mat<-data.frame(x=rep(seq(-100,100,1),each=201),y=rep(seq(-100,100,1),times=201))
  out<-mapply(two.D.gauss, x = test.mat[,1],y=test.mat[,2], MoreArgs = list(q=1,s=10,alphay=10,c=4))
  res.mat<-matrix(out,201,201,byrow = T)
  filled.contour(res.mat)
}

## options moving forward: fit tilted gaussian, estimate 4 params while including height; fit two D gaussian with 3 params

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
    plot(0,0,type="n",xlim=c(min(x1,x2,xp,xtarget)-.5,max(x1,x2,xp,xtarget)+.5),ylim=c(min(y1,y2,yp,ytarget,-1.5),max(y1,y2,yp,ytarget,1.5)),asp=1)
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
