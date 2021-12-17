par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(0,0,xlim=c(-2.5,2.5),ylim=c(-2.5,2),type="n",asp=1,xlab="",ylab="")
xorigin<-.25
yorigin<-.25
degree<-45
xs<-xorigin
xw<-xs+sin(degree)
ys<-yorigin
yw<-ys+cos(degree)
xt<--1
yt<--1
c<-((xt-xs)*(xw-xs)+(yt-ys)*(yw-ys))/((xw-xs)^2+(yw-ys)^2)
xp<-xs+c*(xw-xs)
yp<-ys+c*(yw-ys)

arrows(xs,ys,xs+1.5*(xw-xs),ys+1.5*(yw-ys),lwd=4)
segments(xs,ys,xp,yp,lty=2,col="blue",lwd=4)
segments(xp,yp,xt,yt,lty=2,col="blue",lwd=4)
segments(xs,ys,xt,yt,lty=2,col="blue",lwd=4)

points(xp,yp,cex=3,pch=16)
text(xp,yp,expression((x[P]*','*y[P])),pos=4,cex=1.5)
points(xs,ys,cex=3,pch=16)
text(xs,ys,expression((x[S]*','*y[S])),pos=2,cex=1.5)
points(xw,yw,cex=3,pch=16)
text(xw,yw,expression((x[W]*','*y[W])),pos=4,cex=1.5)
points(xt,yt,cex=3,pch=16)
text(xt,yt,expression((x[T]*','*y[T])),pos=4,cex=1.5)
text(xs+1.5*(xw-xs),ys+1.5*(yw-ys),"W",pos=3,cex=1.5)

s=-1*((yw-ys)/(xw-xs))^-1
curve(ys+s*(x-xs),add=T)

points(xp,yp)
text(xp,yp,expression((x[P]*','*y[P])),pos=4)
segments(xs,ys,xp,yp,lty=2,col="blue")
segments(xp,yp,xtarget,ytarget,lty=2,col="blue")
segments(xs,ys,xtarget,ytarget,lty=2,col="blue")

s=-1*((yw-ys)/(xw-xs))^-1 #slope of perpendicular line

curve(((yw-ys)/(xw-xs))*x-((yw-ys)/(xw-xs))*xs+ys,add=T,lty=2,col="red")
curve(s*x+ys+s*xs,add=T,col="red",lty=2)
