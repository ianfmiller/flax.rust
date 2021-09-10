# constructing sets of predictor variables--offset, random effects not included
################################################################################################################################

### construct all combinations of predictors
predictors<-c("time",
               "mean.temp","max.temp","min.temp",
               "mean.abs.hum","max.abs.hum","min.abs.hum",
               "mean.temp,mean.abs.hum",
               "mean.vpd","max.vpd","min.vpd",
               "mean.wetness",
               "tot.rain",
               "mean.solar"
               
)

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("mean.vpd",c("mean.temp","max.temp","min.temp","mean.abs.hum","max.abs.hum","min.abs.hum","mean.temp,mean.abs.hum")),
  list("max.vpd",c("mean.temp","max.temp","min.temp","mean.abs.hum","max.abs.hum","min.abs.hum","mean.temp,mean.abs.hum")),
  list("min.vpd",c("mean.temp","max.temp","min.temp","mean.abs.hum","max.abs.hum","min.abs.hum","mean.temp,mean.abs.hum")),
  list("mean.temp",c("mean.vpd","max.vpd","min.vpd","mean.temp,mean.abs.hum")),
  list("max.temp",c("mean.vpd","max.vpd","min.vpd","mean.temp,mean.abs.hum")),
  list("min.temp",c("mean.vpd","max.vpd","min.vpd","mean.temp,mean.abs.hum")),
  list("mean.abs.hum",c("mean.vpd","max.vpd","min.vpd","mean.temp,mean.abs.hum")),
  list("max.abs.hum",c("mean.vpd","max.vpd","min.vpd","mean.temp,mean.abs.hum")),
  list("min.abs.hum",c("mean.vpd","max.vpd","min.vpd","mean.temp,mean.abs.hum")),
  list("mean.temp,mean.abs.hum",c(c("mean.temp","max.temp","min.temp","mean.abs.hum","max.abs.hum","min.abs.hum","mean.vpd","max.vpd","min.vpd")))
)

### turn variables off
for (i in 1:length(switch.vars))
{
  on.var<-unlist(switch.vars[[i]][1])
  off.vars<-unlist(switch.vars[[i]][2])
  off.rows<-which(pred.mat[,on.var]==T)
  pred.mat[off.rows,off.vars]<-F
}



### subset to only unique combinations
pred.mat<-unique(pred.mat)

