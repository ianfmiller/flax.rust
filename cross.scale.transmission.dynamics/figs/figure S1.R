par(mfrow=c(1,2),mar=c(2,2,2,2))

plot(0,0,xlim=c(-2.5,2),ylim=c(-2.5,2),type="n",asp=1,xlab="",ylab="",main="X < 0",cex.main=2)
mtext("A",adj=1,font=2)
xorigin<-0
yorigin<-0
degree<-2*pi*(45/360)
xs<-xorigin
xw<-xs+sin(degree)
ys<-yorigin
yw<-ys+cos(degree)
xt<--1.5
yt<--3
c<-((xt-xs)*(xw-xs)+(yt-ys)*(yw-ys))/((xw-xs)^2+(yw-ys)^2)
xp<-xs+c*(xw-xs)
yp<-ys+c*(yw-ys)

arrows(xs,ys,xs+2*(xw-xs),ys+2*(yw-ys),lwd=4,col="grey")
segments(xs,ys,xp,yp,lty=2,col="red",lwd=4)
segments(xp,yp,xt,yt,lty=2,col="blue",lwd=4)

points(xp,yp,cex=2,pch=16,col="purple")
text(xp,yp,expression((x[P]*','*y[P])),pos=4,cex=1.5,offset=1,col="purple")
points(xs,ys,cex=3,pch=16)
text(xs,ys,expression((x[S]*','*y[S])),pos=2,cex=1.5)
points(xw,yw,cex=2,pch=16,col="grey")
text(xw,yw,expression((x[W]*','*y[W])),pos=4,cex=1.5,offset=1,col="grey")
points(xt,yt,cex=3,pch=16,col="orange")
text(xt,yt,expression((x[T]*','*y[T])),pos=4,cex=1.5,col="orange")
text(xs+2*(xw-xs),ys+2*(yw-ys),"W",pos=3,cex=1.5,offset=1,col="grey")
text(-1.5,-1,"X",cex=1.5,col="red")
text(-2.1,-2.9,"Y",cex=1.5,col="blue")

plot(0,0,xlim=c(-1,3.5),ylim=c(-1,3.5),type="n",asp=1,xlab="",ylab="",main="X > 0",cex.main=2)
mtext("B",adj=1,font=2)
xorigin<-0
yorigin<-0
degree<-2*pi*(45/360)
xs<-xorigin
xw<-xs+sin(degree)
ys<-yorigin
yw<-ys+cos(degree)
xt<-1
yt<-4
c<-((xt-xs)*(xw-xs)+(yt-ys)*(yw-ys))/((xw-xs)^2+(yw-ys)^2)
xp<-xs+c*(xw-xs)
yp<-ys+c*(yw-ys)

arrows(xs,ys,xs+2*(xw-xs),ys+2*(yw-ys),lwd=4,col="grey")
segments(xs,ys,xp,yp,lty=2,col="red",lwd=4)
segments(xp,yp,xt,yt,lty=2,col="blue",lwd=4)

points(xp,yp,cex=2,pch=16,col="purple")
text(xp,yp,expression((x[P]*','*y[P])),pos=4,cex=1.5,offset=1,col="purple")
points(xs,ys,cex=3,pch=16)
text(xs,ys,expression((x[S]*','*y[S])),pos=2,cex=1.5)
points(xw,yw,cex=2,pch=16,col="grey")
text(xw,yw,expression((x[W]*','*y[W])),pos=4,cex=1.5,offset=1,col="grey")
points(xt,yt,cex=3,pch=16,col="orange")
text(xt,yt,expression((x[T]*','*y[T])),pos=4,cex=1.5,col="orange")
text(xs+2*(xw-xs),ys+2*(yw-ys),"W",pos=3,cex=1.5,offset=1,col="grey")
text(2.3,2,"X",cex=1.5,col="red")
text(2.1,3.25,"Y",cex=1.5,col="blue")

