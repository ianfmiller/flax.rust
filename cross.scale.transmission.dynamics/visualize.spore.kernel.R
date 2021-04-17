source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/prep.enviro.data.R")
plant.inf.intens<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/summarized data/plants.RDS")
site<-"BT"
tag<-86
date0<-"2020-06-24"
date1<-"2020-07-01"

wind.data<-all.weath[which(all.weath$site==site),]
wind.data<-wind.data[which(wind.data$date>as.POSIXct(paste0(date0," 12:00:00"),tz="UTC")),]
wind.data<-wind.data[which(wind.data$date<=as.POSIXct(paste0(date1," 12:00:00"),tz="UTC")),]

plant.inf.intens.index<-intersect(which(plant.inf.intens$Date==date0),which(plant.inf.intens$Tag==tag))
q<-plant.inf.intens[plant.inf.intens.index,"plant.inf.intens"]
half.height<-.5*plant.inf.intens[plant.inf.intens.index,"max.height"]/100

preds<-c()
xxvals<-seq(-2,2,.2)
yyvals<-seq(-2,2,.2)
for(xx in xxvals)
{
  for(yy in yyvals)
  {
    preds<-c(preds,predict.kernel.tilted.plume(q=q,H=half.height,k=5.803369e-07,alphaz=1.596314e-01,Ws=1.100707e+00,xtarget=xx,ytarget=yy,wind.data=wind.data))
  }
}
