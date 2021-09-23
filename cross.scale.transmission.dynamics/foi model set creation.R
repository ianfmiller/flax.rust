# constructing sets of predictor variables--offset, random effects not included
################################################################################################################################

### construct all combinations of predictors
predictors<-c(
  "height.cm",
  "mean.temp","max.temp","min.temp",
  "mean.abs.hum","max.abs.hum","min.abs.hum",
  #"mean.vpd","max.vpd","min.vpd",
  "tot.rain",
  "mean.solar"
  #"pred.pustule.diam.growth","pred.pustule.num.increase","pred.plant.inf.intens.increase"
)

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### subset to only unique combinations
pred.mat<-unique(pred.mat)
