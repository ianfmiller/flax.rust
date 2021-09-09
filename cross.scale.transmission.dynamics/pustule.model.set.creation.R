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

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### subset to only unique combinations
pred.mat<-unique(pred.mat)

