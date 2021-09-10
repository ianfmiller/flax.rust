# constructing sets of predictor variables--offset, random effects not included
################################################################################################################################

### construct all combinations of predictors
predictors<-c("time",
               "mean.temp","max.temp","min.temp",
               "mean.dew.point","max.dew.point","min.dew.point",
               "mean.wetness",
               "tot.rain",
               "mean.solar"
               
)

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("vpd",c("mean.temp","max.temp","min.temp","mean.dew.point","max.dew.point","min.dew.point")),
  list("mean.temp",c("vpd")),
  list("max.temp",c("vpd")),
  list("min.temp",c("vpd")),
  list("mean.dew.point",c("vpd")),
  list("max.dew.point",c("vpd")),
  list("min.dew.point",c("vpd")),
)

### turn variables off
for (i in 1:length(switch.vars))
{
  on.var<-unlist(switch.vars[[i]][1])
  off.vars<-unlist(switch.vars[[i]][2])
  off.rows<-which(pred.mat[,on.var]==T)
  pred.mat[off.rows,off.vars]<-F
}

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### subset to only unique combinations
pred.mat<-unique(pred.mat)

