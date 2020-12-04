### temp data visualization ###
library(lubridate)

source("~/Documents/Github/flax.rust/prep.enviro.data.R")

make.transparent<-function(col,alpha)
{
  rgb<-col2rgb(col)
  r<-rgb[1]/255
  g<-rgb[2]/255
  b<-rgb[3]/255
  rgb(r,g,b,alpha)
}

### plot raw data

par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(5.1,4.1,4.1,2.1))


plot.func<-function(x,title)
{
  plot(x[,1],x[,2],type="l",col="red",main=title,ylim=c(-5,100),xlim=c(min(all.temp.rh[,1],na.rm=T),max(all.temp.rh[,1],na.rm = T)),xlab="date",ylab="")
  mtext("temp",2,cex=par()$cex,line=2.25,col="red")
  lines(x[,1],x[,3],type="l",col="blue")
  mtext("RH",2,cex=par()$cex,line=3,col="blue",las=0)
}

plot.func(hm,"high meadow")
plot.func(gm,"gothic mountain")
plot.func(bt,"bus turnaround")
plot.func(cc,"cement creek")

mtext("raw temp + humidity data",outer=T,cex=1.75)

par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(2.1,4.1,2.1,2.1))

plot(cc.weath$date,cc.weath$rain,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,4),ylab="rainfall.temp.rh (mm)",xlab="")
lines(bt.weath$date,bt.weath$rain,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$rain,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$rain,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("rainfall.temp.rh")

plot(cc.weath$date,cc.weath$wetness,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,100),ylab="leaf wetness (%)",xlab="")
lines(bt.weath$date,bt.weath$wetness,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$wetness,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$wetness,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("wetness")

plot(cc.weath$date,cc.weath$soil.moisture,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,.3),ylab=expression("soil moisture ("*m^3/m^3*")"),xlab="")
lines(bt.weath$date,bt.weath$soil.moisture,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$soil.moisture,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$soil.moisture,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("soil moisture")

par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(2.1,4.1,2.1,2.1))

plot(cc.weath$date,cc.weath$wind.speed,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,8),ylab="wind speed (m/s)",xlab="")
lines(bt.weath$date,bt.weath$wind.speed,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$wind.speed,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$wind.speed,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("wind speed")

plot(cc.weath$date,cc.weath$gust.speed,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,15),ylab="gust speed (m/s)",xlab="")
lines(bt.weath$date,bt.weath$gust.speed,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$gust.speed,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$gust.speed,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("gust speed")

plot(cc.weath$date,cc.weath$wind.direction,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,360),ylab="wind direction (degrees)",xlab="")
lines(bt.weath$date,bt.weath$wind.direction,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$wind.direction,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$wind.direction,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("wind direction")

par(mfrow=c(1,1))
plot(cc.weath$date,cc.weath$solar.radiation,type="l",col=make.transparent("yellow2",alpha = .5),lwd=2,xlim=c(1592577358,1596028858),ylim=c(0,1300),ylab=expression("solar radiation ("*W/m^2*")"),xlab="")
lines(bt.weath$date,bt.weath$solar.radiation,col=make.transparent("orange",alpha = .5),lwd=2)
lines(gm.weath$date,gm.weath$solar.radiation,col=make.transparent("red",alpha = .5),lwd=2)
lines(hm.weath$date,hm.weath$solar.radiation,col=make.transparent("purple",alpha = .5),lwd=2)
mtext("solar radiation")

### plot daily max/min data

pull.daily.data<-function(input,output.name)
{

output<-as.data.frame(matrix(NA,length(unique(as.Date(input[,1]))),10))
output[,1]<-unique(as.Date(input[,1]))
colnames(output)<-c("date","max.temp","min.temp","mean.temp","max.RH","min.RH","mean.RH","max.dewpoint","min.dewpoint","mean.dewpoint")

for (x in 1:(dim(output)[1]))
{
  index<-which(as.Date(input[,1])==output[x,1])
  output[x,"max.temp"]<-max(input[index,2],na.rm = T)
  output[x,"min.temp"]<-min(input[index,2],na.rm = T)
  output[x,"mean.temp"]<-mean(input[index,2],na.rm = T)
  output[x,"max.RH"]<-max(input[index,3],na.rm = T)
  output[x,"min.RH"]<-min(input[index,3],na.rm = T)
  output[x,"mean.RH"]<-mean(input[index,3],na.rm = T)
  output[x,"max.dewpoint"]<-max(input[index,4],na.rm = T)
  output[x,"min.dewpoint"]<-min(input[index,4],na.rm = T)
  output[x,"mean.dewpoint"]<-min(input[index,4],na.rm = T)
}

assign(output.name,output,envir=.GlobalEnv)

}

pull.daily.data(hm,"hm.day")
pull.daily.data(gm,"gm.day")
pull.daily.data(bt,"bt.day")
pull.daily.data(cc,"cc.day")

plot.func2<-function(x,title)
{
  plot(x[,1],x[,"max.temp"],type="l",col="red",main=title,ylim=c(-5,100),xlim=as.Date(c(min(all.temp.rh[,1],na.rm=T),max(all.temp.rh[,1],na.rm = T))),xlab="date",ylab="")
  lines(x[,1],x[,"min.temp"],type="l",col="red",lty=2)
  mtext("temp",2,cex=par()$cex,line=2.25,col="red")
  
  lines(x[,1],x[,"max.RH"],type="l",col="blue")
  lines(x[,1],x[,"min.RH"],type="l",col="blue",lty=2)
  
  mtext("RH",2,cex=par()$cex,line=3,col="blue",las=0)
}

par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(5.1,4.1,4.1,2.1))

plot.func2(hm.day,"high meadow")
plot.func2(gm.day,"gothic mountain")
plot.func2(bt.day,"bus turnaround")
plot.func2(cc.day,"cement creek")

mtext("daily max/min temp + humidity",outer=T,cex=1.75)


### site comparison ###

par(mfrow=c(1,3),mar=c(3.1,3.1,3.1,1.1),oma=c(0,0,2,0))

plot(cc.day[,1],cc.day[,"max.temp"],type = "l",col="yellow2",xlab="date",ylab="temp",ylim=c(15,36),xlim=c(18424,18472),main="daily max",lwd=2)
points(bt.day[,1],bt.day[,"max.temp"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"max.temp"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"max.temp"],type = "l",col="purple",lwd=2)

plot(cc.day[,1],cc.day[,"min.temp"],type = "l",col="yellow2",xlab="date",ylab="temp",ylim=c(-6,12),xlim=c(18424,18472),main="daily min",lwd=2)
points(bt.day[,1],bt.day[,"min.temp"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"min.temp"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"min.temp"],type = "l",col="purple",lwd=2)

plot(cc.day[,1],cc.day[,"mean.temp"],type = "l",col="red",xlab="date",ylab="temp",ylim=c(7,21),xlim=c(18424,18472),main="daily mean",lwd=2)
points(bt.day[,1],bt.day[,"mean.temp"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"mean.temp"],type = "l",col="yellow2",lwd=2)
points(hm.day[,1],hm.day[,"mean.temp"],type = "l",col="purple",lwd=2)

legend("bottomright",c("CC","BT","GM","HM"),col=c("yellow2","orange","red","purple"),lty=1,lwd=2)

mtext("Temperature",3,outer=T,cex=1.5)

par(mfrow=c(1,3),mar=c(3.1,3.1,3.1,1.1),oma=c(0,0,2,0))

plot(cc.day[,1],cc.day[,"max.RH"],type = "l",col="yellow2",xlab="date",ylab="RH",ylim=c(25,100),xlim=c(18424,18472),main="daily max",lwd=2)
points(bt.day[,1],bt.day[,"max.RH"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"max.RH"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"max.RH"],type = "l",col="purple",lwd=2)

plot(cc.day[,1],cc.day[,"min.RH"],type = "l",col="yellow2",xlab="date",ylab="RH",ylim=c(5,65),xlim=c(18424,18472),main="daily min",lwd=2)
points(bt.day[,1],bt.day[,"min.RH"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"min.RH"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"min.RH"],type = "l",col="purple",lwd=2)

plot(cc.day[,1],cc.day[,"mean.RH"],type = "l",col="yellow2",xlab="date",ylab="RH",ylim=c(18,90),xlim=c(18424,18472),main="daily mean",lwd=2)
points(bt.day[,1],bt.day[,"mean.RH"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"mean.RH"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"mean.RH"],type = "l",col="purple",lwd=2)

legend("bottomright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),lty=1)

mtext("Relative Humidity",3,outer=T,cex=1.5)

par(mfrow=c(1,3),mar=c(3.1,3.1,3.1,1.1),oma=c(0,0,2,0))

plot(cc.day[,1],cc.day[,"max.dewpoint"],type = "l",col="yellow2",xlab="date",ylab="dewpoint",xlim=c(18424,18472),ylim=c(0,17),main="daily max",lwd=2)
points(bt.day[,1],bt.day[,"max.dewpoint"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"max.dewpoint"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"max.dewpoint"],type = "l",col="purple",lwd=2)

plot(cc.day[,1],cc.day[,"min.dewpoint"],type = "l",col="yellow2",xlab="date",ylab="dewpoint",xlim=c(18424,18472),ylim=c(-17,10),main="daily min",lwd=2)
points(bt.day[,1],bt.day[,"min.dewpoint"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"min.dewpoint"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"min.dewpoint"],type = "l",col="purple",lwd=2)

plot(cc.day[,1],cc.day[,"mean.dewpoint"],type = "l",col="yellow2",xlab="date",ylab="dewpoint",xlim=c(18424,18472),ylim=c(-17,10),main="daily mean",lwd=2)
points(bt.day[,1],bt.day[,"mean.dewpoint"],type = "l",col="orange",lwd=2)
points(gm.day[,1],gm.day[,"mean.dewpoint"],type = "l",col="red",lwd=2)
points(hm.day[,1],hm.day[,"mean.dewpoint"],type = "l",col="purple",lwd=2)

legend("bottomright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),lty=1)

mtext("Dewpoint",3,outer=T,cex=1.5)

################################################################################################
################################ distribution comparisons ######################################

par(mfrow=c(2,3))
hist(cc.day[,"min.temp"],col=make.transparent("yellow2",.5),xlim=c(-7,12),ylim=c(0,21),main="daily min temp",xlab=expression(''*~degree*C*''))
hist(bt.day[,"min.temp"],add=T,col=make.transparent("orange",.5))
hist(gm.day[,"min.temp"],add=T,col=make.transparent("red",.5))
hist(hm.day[,"min.temp"],add=T,col=make.transparent("purple",.5))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.day[,"max.temp"],col=make.transparent("yellow2",.5),xlim=c(15,36),ylim=c(0,17),main="daily max temp",xlab=expression(''*~degree*C*''))
hist(bt.day[,"max.temp"],add=T,col=make.transparent("orange",.5))
hist(gm.day[,"max.temp"],add=T,col=make.transparent("red",.5))
hist(hm.day[,"max.temp"],add=T,col=make.transparent("purple",.5))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.day[,"mean.temp"],col=make.transparent("yellow2",.5),xlim=c(6,22),ylim=c(0,22),main="daily mean temp",xlab=expression(''*~degree*C*''))
hist(bt.day[,"mean.temp"],add=T,col=make.transparent("orange",.5))
hist(gm.day[,"mean.temp"],add=T,col=make.transparent("red",.5))
hist(hm.day[,"mean.temp"],add=T,col=make.transparent("purple",.5))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

plot(density(cc.day[,"min.temp"]),col="yellow2",xlim=c(-7,12),ylim=c(0,.25),main = "daily min temp")
lines(density(bt.day[,"min.temp"]),col="orange")
lines(density(gm.day[,"min.temp"]),col="red")
lines(density(hm.day[,"min.temp"]),col="purple")

plot(density(cc.day[,"max.temp"]),col="yellow2",xlim=c(15,36),ylim=c(0,.175),main = "daily max temp")
lines(density(bt.day[,"max.temp"]),col="orange")
lines(density(gm.day[,"max.temp"]),col="red")
lines(density(hm.day[,"max.temp"]),col="purple")

plot(density(cc.day[,"mean.temp"]),col="yellow2",xlim=c(6,22),main = "daily mean temp")
lines(density(bt.day[,"mean.temp"]),col="orange")
lines(density(gm.day[,"mean.temp"]),col="red")
lines(density(hm.day[,"mean.temp"]),col="purple")

par(mfrow=c(2,3))
hist(cc.day[,"min.dewpoint"],col=make.transparent("yellow2",.5),xlim=c(-17,11),breaks=20,main="daily min dewpoint",xlab=expression(''*~degree*C*''))
hist(bt.day[,"min.dewpoint"],add=T,col=make.transparent("orange",.5),breaks=20)
hist(gm.day[,"min.dewpoint"],add=T,col=make.transparent("red",.5),breaks=20)
hist(hm.day[,"min.dewpoint"],add=T,col=make.transparent("purple",.5),breaks=20)
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.day[,"max.dewpoint"],col=make.transparent("yellow2",.5),xlim=c(0,17),breaks=20,main="daily max dewpoint",xlab=expression(''*~degree*C*''))
hist(bt.day[,"max.dewpoint"],add=T,col=make.transparent("orange",.5),breaks=20)
hist(gm.day[,"max.dewpoint"],add=T,col=make.transparent("red",.5),breaks=20)
hist(hm.day[,"max.dewpoint"],add=T,col=make.transparent("purple",.5),breaks=20)
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.day[,"mean.dewpoint"],col=make.transparent("yellow2",.5),xlim=c(-17,11),breaks=20,main="daily mean dewpoint",xlab=expression(''*~degree*C*''))
hist(bt.day[,"mean.dewpoint"],add=T,col=make.transparent("orange",.5),breaks=20)
hist(gm.day[,"mean.dewpoint"],add=T,col=make.transparent("red",.5),breaks=20)
hist(hm.day[,"mean.dewpoint"],add=T,col=make.transparent("purple",.5),breaks=20)
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

plot(density(cc.day[,"min.dewpoint"]),col="yellow2",xlim=c(-22,17),main = "daily min dewpoint")
lines(density(bt.day[,"min.dewpoint"]),col="orange")
lines(density(gm.day[,"min.dewpoint"]),col="red")
lines(density(hm.day[,"min.dewpoint"]),col="purple")

plot(density(cc.day[,"max.dewpoint"]),col="yellow2",xlim=c(-5,22),ylim=c(0,.175),main = "daily max dewpoint")
lines(density(bt.day[,"max.dewpoint"]),col="orange")
lines(density(gm.day[,"max.dewpoint"]),col="red")
lines(density(hm.day[,"max.dewpoint"]),col="purple")

plot(density(cc.day[,"mean.dewpoint"]),col="yellow2",xlim=c(-25,14),main = "daily mean dewpoint")
lines(density(bt.day[,"mean.dewpoint"]),col="orange")
lines(density(gm.day[,"mean.dewpoint"]),col="red")
lines(density(hm.day[,"mean.dewpoint"]),col="purple")

par(mfrow=c(2,3))
hist(cc.weath[,"rain"],col=make.transparent("yellow2",.5),xlim=c(0,4),ylim=c(0,20000),breaks=seq(0,4,length.out = 20),main="rainfall.temp.rh",xlab="mm")
hist(bt.weath[,"rain"],add=T,col=make.transparent("orange",.5),breaks=seq(0,4,length.out = 20))
hist(gm.weath[,"rain"],add=T,col=make.transparent("red",.5),breaks=seq(0,4,length.out = 20))
hist(hm.weath[,"rain"],add=T,col=make.transparent("purple",.5),breaks=seq(0,4,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.weath[,"wetness"],col=make.transparent("yellow2",.5),xlim=c(0,100),breaks=seq(0,100,length.out = 20),main="leaf wetness",xlab="%")
hist(bt.weath[,"wetness"],add=T,col=make.transparent("orange",.5),breaks=seq(0,100,length.out = 20))
hist(gm.weath[,"wetness"],add=T,col=make.transparent("red",.5),breaks=seq(0,100,length.out = 20))
hist(hm.weath[,"wetness"],add=T,col=make.transparent("purple",.5),breaks=seq(0,100,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.weath[,"soil.moisture"],col=make.transparent("yellow2",.5),xlim=c(0,.3),breaks=seq(0,.3,length.out = 20),main="soil moisture",xlab=expression("soil moisture ("*m^3/m^3*")"))
hist(bt.weath[,"soil.moisture"],add=T,col=make.transparent("orange",.5),breaks=seq(0,.3,length.out = 20))
hist(gm.weath[,"soil.moisture"],add=T,col=make.transparent("red",.5),breaks=seq(-0.01578947,.3,length.out = 21))
hist(hm.weath[,"soil.moisture"],add=T,col=make.transparent("purple",.5),breaks=seq(0,.3,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

plot(density(cc.weath[,"rain"],na.rm=T),col="yellow2",xlim=c(0,4),main = "rainfall.temp.rh")
lines(density(bt.weath[,"rain"],na.rm=T),col="orange")
lines(density(gm.weath[,"rain"],na.rm=T),col="red")
lines(density(hm.weath[,"rain"],na.rm=T),col="purple")

plot(density(cc.weath[,"wetness"],na.rm=T),col="yellow2",xlim=c(0,100),ylim=c(0,.5),main = "leaf wetness")
lines(density(bt.weath[,"wetness"],na.rm=T),col="orange")
lines(density(gm.weath[,"wetness"],na.rm=T),col="red")
lines(density(hm.weath[,"wetness"],na.rm=T),col="purple")

plot(density(cc.weath[,"soil.moisture"],na.rm=T),col="yellow2",xlim=c(0,.3),main = "soil moisture")
lines(density(bt.weath[,"soil.moisture"],na.rm=T),col="orange")
lines(density(gm.weath[,"soil.moisture"],na.rm=T),col="red")
lines(density(hm.weath[,"soil.moisture"],na.rm=T),col="purple")

par(mfrow=c(2,3))
hist(cc.weath[,"wind.speed"],col=make.transparent("yellow2",.5),xlim=c(0,8),ylim=c(0,7500),breaks=seq(0,8,length.out = 20),main="wind speed",xlab="m/s")
hist(bt.weath[,"wind.speed"],add=T,col=make.transparent("orange",.5),breaks=seq(0,8,length.out = 20))
hist(gm.weath[,"wind.speed"],add=T,col=make.transparent("red",.5),breaks=seq(0,8,length.out = 20))
hist(hm.weath[,"wind.speed"],add=T,col=make.transparent("purple",.5),breaks=seq(0,8,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.weath[,"gust.speed"],col=make.transparent("yellow2",.5),xlim=c(0,15),ylim=c(0,5000),breaks=seq(0,15,length.out = 20),main="gust speed",xlab="m/s")
hist(bt.weath[,"gust.speed"],add=T,col=make.transparent("orange",.5),breaks=seq(0,15,length.out = 20))
hist(gm.weath[,"gust.speed"],add=T,col=make.transparent("red",.5),breaks=seq(0,15,length.out = 20))
hist(hm.weath[,"gust.speed"],add=T,col=make.transparent("purple",.5),breaks=seq(0,15,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

hist(cc.weath[,"wind.direction"],col=make.transparent("yellow2",.5),xlim=c(0,360),ylim=c(0,4000),breaks=seq(0,360,length.out = 20),main="wind direction",xlab="degrees")
hist(bt.weath[,"wind.direction"],add=T,col=make.transparent("orange",.5),breaks=seq(0,360,length.out = 20))
hist(gm.weath[,"wind.direction"],add=T,col=make.transparent("red",.5),breaks=seq(0,360,length.out = 20))
hist(hm.weath[,"wind.direction"],add=T,col=make.transparent("purple",.5),breaks=seq(0,360,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

plot(density(cc.weath[,"wind.speed"],na.rm=T),col="yellow2",xlim=c(0,8),main = "wind speed")
lines(density(bt.weath[,"wind.speed"],na.rm=T),col="orange")
lines(density(gm.weath[,"wind.speed"],na.rm=T),col="red")
lines(density(hm.weath[,"wind.speed"],na.rm=T),col="purple")

plot(density(cc.weath[,"gust.speed"],na.rm=T),col="yellow2",xlim=c(0,15),ylim=c(0,.5),main = "gust speed")
lines(density(bt.weath[,"gust.speed"],na.rm=T),col="orange")
lines(density(gm.weath[,"gust.speed"],na.rm=T),col="red")
lines(density(hm.weath[,"gust.speed"],na.rm=T),col="purple")

plot(density(cc.weath[,"wind.direction"],na.rm=T),col="yellow2",xlim=c(0,360),main = "wind direction")
lines(density(bt.weath[,"wind.direction"],na.rm=T),col="orange")
lines(density(gm.weath[,"wind.direction"],na.rm=T),col="red")
lines(density(hm.weath[,"wind.direction"],na.rm=T),col="purple")

par(mfrow=c(2,1))
hist(cc.weath[,"solar.radiation"],col=make.transparent("yellow2",.5),xlim=c(0,1300),ylim=c(0,7500),breaks=seq(0,1300,length.out = 20),main="solar radiation",xlab=expression(""*W/m^2*""))
hist(bt.weath[,"solar.radiation"],add=T,col=make.transparent("orange",.5),breaks=seq(0,1300,length.out = 20))
hist(gm.weath[,"solar.radiation"],add=T,col=make.transparent("red",.5),breaks=seq(0,1300,length.out = 20))
hist(hm.weath[,"solar.radiation"],add=T,col=make.transparent("purple",.5),breaks=seq(0,1300,length.out = 20))
legend("topright",c("HM","GM","BT","CC"),col=c("purple","red","orange","yellow2"),pch=15)

plot(density(cc.weath[,"solar.radiation"],na.rm=T),col="yellow2",xlim=c(0,1300),main = "solar radiation")
lines(density(bt.weath[,"solar.radiation"],na.rm=T),col="orange")
lines(density(gm.weath[,"solar.radiation"],na.rm=T),col="red")
lines(density(hm.weath[,"solar.radiation"],na.rm=T),col="purple")

