# constructing sets of predictor variables--offset, random effects not included
################################################################################################################################

### construct all combinations of predictors
predictors<-c("time",
               "temp.days","temp.days.16.22","temp.days.7.30",
               "dew.point.days",
               "temp.dew.point.days","temp.16.22.dew.point.days","temp.7.30.dew.point.days",
               "wetness.days",
               "temp.wetness.days","temp.16.22.wetness.days","temp.7.30.wetness.days",
               "tot.rain","solar.days","wind.speed.days","gust.speed.days"
)

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("temp.days",c("temp.days.16.22","temp.days.7.30")),
  list("temp.days.16.22",c("temp.days","temp.days.7.30")),
  list("temp.days.7.30",c("temp.days","temp.days.16.22")),
  list("temp.dew.point.days",c("temp.16.22.dew.point.days","temp.7.30.dew.point.days")),
  list("temp.16.22.dew.point.days",c("temp.dew.point.days","temp.7.30.dew.point.days")),
  list("temp.7.30.dew.point.days",c("temp.dew.point.days","temp.16.22.dew.point.days")),
  list("temp.wetness.days",c("temp.16.22.wetness.days","temp.7.30.wetness.days")),
  list("temp.16.22.wetness.days",c("temp.wetness.days","temp.7.30.wetness.days")),
  list("temp.7.30.wetness.days",c("temp.wetness.days","temp.16.22.wetness.days"))
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

