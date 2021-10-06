plot<-F

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant data prep.R")
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant loc dataset building.R")

epi.obs.dates<-list("CC"=c("2020-06-22","2020-06-29","2020-07-06","2020-07-13","2020-07-20","2020-07-27"),"BT"=c("2020-06-24","2020-07-01","2020-07-08"),"GM"=c("2020-06-23","2020-06-30","2020-07-02","2020-07-07","2020-07-09","2020-07-15"),"HM"=c("2020-06-25","2020-07-02","2020-07-07","2020-07-09","2020-07-15"))    

data.vec<-c() 
for(tag in unique(plants$Tag)) # loop to pull out plant inf intens of newly diseased plants
{
  if((tag %in% corrected.epi$Tag))
  {
    inf.date<-corrected.epi[which(corrected.epi$Tag==tag),"Date.First.Observed.Diseased"]
    sub.plant.data<-plants[which(plants$Tag==tag),]
    site<-sub.plant.data[1,"Site"]
    if(site=="BT") {min.date<-"2020-06-24"}
    if(site=="CC") {min.date<-"2020-06-24"}
    if(site=="GM") {min.date<-"2020-06-23"}
    if(site=="HM") {min.date<-"2020-06-25"}
    if(inf.date>=min.date) {
      index<-which(sub.plant.data$Date==inf.date)
      data.vec<-c(data.vec,sub.plant.data[index,"plant.inf.intens"])
    }
  }
}

data.vec<-data.vec[which(data.vec<=100)] #subset out observations likely to not be truly new infections
model.data<-hist(data.vec,breaks=100,plot=F)
y<-model.data$density
x<-model.data$mids
mod<-nls(y~c*exp(-c*x),start=list(c=.5))
lambda<-summary(mod)$coefficients[1]

starting.plant.inf.intens.mod<-function(x) {log(1-x)/-2.103}

start.inf.intens.mod<-function(x) {1-exp(-2.103*x)}

if(plot)
{
plot(model.data$mids,model.data$density,xlab="plant inf. intens.",ylab="prob")
curve(2.103*exp(-2.103*x),add=T)
legend("topright",legend=c("observed","fitted"),pch=c(1,NA),lty=c(NA,1))
curve(1-exp(-2.103*x),ylim=c(0,1),xlim=c(0,5),xlab="plant inf. intens.",ylab="fitted cumulative prob")
}

      