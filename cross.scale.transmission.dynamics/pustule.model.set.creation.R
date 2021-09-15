# constructing sets of predictor variables--offset, random effects not included
################################################################################################################################

### construct all combinations of predictors
predictors<-c( "mean.temp*time","s(max.temp,k=4,by=time)",
               "max.temp*time","s(max.temp,k=4,by=time)",
               "min.temp*time","s(min.temp,k=4,by=time)",
               "mean.abs.hum*time","s(mean.abs.hum,k=4,by=time)",
               "max.abs.hum*time","s(max.abs.hum,k=4,by=time)",
               "min.abs.hum*time","s(min.abs.hum,k=4,by=time)",
               "mean.wetness*time","s(mean.wetness,k=4,by=time)",
               "tot.rain*time","s(tot.rain,k=4,by=time)",
               "mean.solar*time","s(mean.solar,k=4,by=time)"
)

pred.mat<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat)<-predictors

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("mean.temp*time",c("s(mean.temp,k=4,by=time)")),
  list("max.temp*time",c("s(max.temp,k=4,by=time)")),
  list("min.temp*time",c("s(min.temp,k=4,by=time)")),
  list("mean.abs.hum*time",c("s(mean.abs.hum,k=4,by=time)")),
  list("max.abs.hum*time",c("s(max.abs.hum,k=4,by=time)")),
  list("min.abs.hum*time",c("s(min.abs.hum,k=4,by=time)")),
  list("mean.wetness*time",c("s(mean.wetness,k=4,by=time)")),
  list("tot.rain*time",c("s(tot.rain,k=4,by=time)")),
  list("mean.solar*time",c("s(mean.solar,k=4,by=time)"))
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

