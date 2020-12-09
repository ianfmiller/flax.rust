# constructing sets of predictor variables--offset, random effects not included

################################################################################################################################

## no temp. subsetting

### construct all combinations of predictors
predictors1<-c("area","time","mean.temp","mean.dew.point","mean.wetness","mean.solar","mean.wind.speed","tot.rain", 
              "mean.temp:time","temp.days", 
              "mean.dew.point:time","dew.point.days",
              "mean.wetness:time", 
              "mean.solar:time",
              "mean.wind.speed:time",
              "mean.temp:mean.wetness",
              "mean.temp:mean.wetness:time",
              "mean.temp:mean.dew.point", 
              "mean.temp:mean.dew.point:time","temp.dew.point.days"
              )

pred.mat1<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat1)<-predictors1

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("mean.temp:time",c("temp.days")),
  list("temp.days",c("mean.temp:time")),
  list("mean.dew.point:time",c("dew.point.days")),
  list("dew.point.days",c("mean.dew.point:time")),
  list("mean.temp:mean.dew.point:time",c("temp.dew.point.days")),
  list("temp.dew.point.days",c("mean.temp:mean.dew.point:time"))
)

### define sets of variables that need to occur together
force.include.vars<-list(
  list("mean.temp:time",c("mean.temp","time")),
  list("mean.dew.point:time",c("mean.dew.point","time")),
  list("mean.wetness:time",c("mean.wetness","time")),
  list("mean.solar:time",c("mean.solar","time")),
  list("mean.wind.speed:time",c("mean.wind.speed","time")),
  list("mean.temp:mean.dew.point",c("mean.temp","mean.dew.point")),
  list("mean.temp:mean.wetness",c("mean.temp","mean.wetness")),
  list("mean.temp:mean.dew.point:time",c("mean.temp","mean.dew.point","time")),
  list("mean.temp:mean.wetness:time",c("mean.temp","mean.wetness","time"))
)

### turn variables off
for (i in 1:length(switch.vars))
{
  on.var<-unlist(switch.vars[[i]][1])
  off.vars<-unlist(switch.vars[[i]][2])
  off.rows<-which(pred.mat1[,on.var]==T)
  pred.mat1[off.rows,off.vars]<-F
}

### turn variables on
for (i in 1:length(force.include.vars))
{
  on.var<-unlist(force.include.vars[[i]][1])
  force.on.vars<-unlist(force.include.vars[[i]][2])
  on.rows<-which(pred.mat1[,on.var]==T)
  pred.mat1[on.rows,force.on.vars]<-T
}

### subset to only unique combinations
pred.mat1<-unique(pred.mat1)

################################################################################################################################

## subset temp between 16 and 22

### construct all combinations of predictors
predictors2<-c("area","time","temp.days.16.22","mean.dew.point","mean.wetness","mean.solar","mean.wind.speed","tot.rain", 
              "mean.dew.point:time","dew.point.days",
              "mean.wetness:time", 
              "mean.solar:time",
              "mean.wind.speed:time",
              "temp.days.16.22:mean.wetness",
              "temp.days.16.22:mean.dew.point"
)

pred.mat2<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat2)<-predictors2

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("mean.dew.point:time",c("dew.point.days")),
  list("dew.point.days",c("mean.dew.point:time"))
)

### define sets of variables that need to occur together
force.include.vars<-list(
  list("mean.dew.point:time",c("mean.dew.point","time")),
  list("mean.wetness:time",c("mean.wetness","time")),
  list("mean.solar:time",c("mean.solar","time")),
  list("mean.wind.speed:time",c("mean.wind.speed","time")),
  list("temp.days.16.22:mean.wetness",c("temp.days.16.22","mean.wetness")),
  list("temp.days.16.22:mean.dew.point",c("temp.days.16.22","mean.dew.point"))
)

### turn variables off
for (i in 1:length(switch.vars))
{
  on.var<-unlist(switch.vars[[i]][1])
  off.vars<-unlist(switch.vars[[i]][2])
  off.rows<-which(pred.mat2[,on.var]==T)
  pred.mat2[off.rows,off.vars]<-F
}

### turn variables on
for (i in 1:length(force.include.vars))
{
  on.var<-unlist(force.include.vars[[i]][1])
  force.on.vars<-unlist(force.include.vars[[i]][2])
  on.rows<-which(pred.mat2[,on.var]==T)
  pred.mat2[on.rows,force.on.vars]<-T
}

### subset to only unique combinations
pred.mat2<-unique(pred.mat2)

################################################################################################################################

## subset temp between 7 and 30

### construct all combinations of predictors
predictors3<-c("area","time","temp.days.7.30","mean.dew.point","mean.wetness","mean.solar","mean.wind.speed","tot.rain", 
              "mean.dew.point:time","dew.point.days",
              "mean.wetness:time", 
              "mean.solar:time",
              "mean.wind.speed:time",
              "temp.days.7.30:mean.wetness",
              "temp.days.7.30:mean.dew.point"
)

pred.mat3<-expand.grid(c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F),c(T,F))
names(pred.mat3)<-predictors3

### define sets of variables that explain the same thing
switch.vars<-list( 
  list("mean.dew.point:time",c("dew.point.days")),
  list("dew.point.days",c("mean.dew.point:time"))
)

### define sets of variables that need to occur together
force.include.vars<-list(
  list("mean.dew.point:time",c("mean.dew.point","time")),
  list("mean.wetness:time",c("mean.wetness","time")),
  list("mean.solar:time",c("mean.solar","time")),
  list("mean.wind.speed:time",c("mean.wind.speed","time")),
  list("temp.days.7.30:mean.wetness",c("temp.days.7.30","mean.wetness")),
  list("temp.days.7.30:mean.dew.point",c("temp.days.7.30","mean.dew.point"))
)

### turn variables off
for (i in 1:length(switch.vars))
{
  on.var<-unlist(switch.vars[[i]][1])
  off.vars<-unlist(switch.vars[[i]][2])
  off.rows<-which(pred.mat3[,on.var]==T)
  pred.mat3[off.rows,off.vars]<-F
}

### turn variables on
for (i in 1:length(force.include.vars))
{
  on.var<-unlist(force.include.vars[[i]][1])
  force.on.vars<-unlist(force.include.vars[[i]][2])
  on.rows<-which(pred.mat3[,on.var]==T)
  pred.mat3[on.rows,force.on.vars]<-T
}

### subset to only unique combinations
pred.mat3<-unique(pred.mat3)
