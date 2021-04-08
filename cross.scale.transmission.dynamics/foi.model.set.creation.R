# constructing sets of predictor variables--offset, random effects not included
################################################################################################################################

### construct all combinations of predictors
predictors<-c("mean.temp.days","mean.temp.days.16.22","mean.temp.days.7.30",
              "mean.dew.point.days",
              "mean.temp.dew.point.days","mean.temp.16.22.dew.point.days","mean.temp.7.30.dew.point.days",
              "mean.wetness.days",
              "mean.temp.wetness.days","mean.temp.16.22.wetness.days","mean.temp.7.30.wetness.days",
              "mean.tot.rain","mean.pred.pustule.diam.growth","mean.pred.pustule.num.increase"
)

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("mean.temp.days",c("mean.temp.days.16.22","mean.temp.days.7.30")),
  list("mean.temp.days.16.22",c("mean.temp.days","mean.temp.days.7.30")),
  list("mean.temp.days.7.30",c("mean.temp.days","mean.temp.days.16.22")),
  list("mean.temp.dew.point.days",c("mean.temp.16.22.dew.point.days","mean.temp.7.30.dew.point.days")),
  list("mean.temp.16.22.dew.point.days",c("mean.temp.dew.point.days","mean.temp.7.30.dew.point.days")),
  list("mean.temp.7.30.dew.point.days",c("mean.temp.dew.point.days","mean.temp.16.22.dew.point.days")),
  list("mean.temp.wetness.days",c("mean.temp.16.22.wetness.days","mean.temp.7.30.wetness.days")),
  list("mean.temp.16.22.wetness.days",c("mean.temp.wetness.days","mean.temp.7.30.wetness.days")),
  list("mean.temp.7.30.wetness.days",c("mean.temp.wetness.days","mean.temp.16.22.wetness.days"))
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
