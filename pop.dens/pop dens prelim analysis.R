pop.dens<-read.csv("~/Documents/GitHub/flax.rust/data/density.transects.csv")

n.h1<-pop.dens$first.pop.n.healthy
n.d1<-pop.dens$first.pop.n.diseased
elev1<-pop.dens$first.pop.elevation
n.h2<-pop.dens$dis.pop.n.healthy
n.d2<-pop.dens$dis.pop.n.diseased
elev2<-pop.dens$dis.pop.elevation

xx<-data.frame(n.h=c(n.h1,n.h2),n.d=c(n.d1,n.d2),elev=c(elev1,elev2))

plot(xx$n.h+xx$n.d,xx$n.d/(xx$n.d+xx$n.h),xlab="density",ylab="prevalence",main="How does pop. density affect prevalence?")
plot(xx$n.h+xx$n.d,(xx$n.d)>0,xlab="density",ylab="incidence",main="How does pop. density affect incidence?")

plot(xx$elev,xx$n.h+xx$n.d,xlab="elevation",ylab="pop. density",main="How does elevation affect pop. density?")
plot(xx$elev,xx$n.d/(xx$n.d+xx$n.h),xlab="elevation",ylab="prevalence",main="How does elevation affect prevalence?")
plot(xx$elev,xx$n.d>0,xlab="elevation",ylab="incidence",main="How does elevation affect incidence?")
