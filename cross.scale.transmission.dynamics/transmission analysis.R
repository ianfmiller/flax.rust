library(mgcv)
library(viridis)

# load data

source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/transmission data set building.R")
transmission.data<-transmission.data[which(transmission.data$status==0),]
transmission.data$site<-as.factor(transmission.data$site)
for(i in 1:nrow(transmission.data))
{
  if(is.na(transmission.data[i,"tag"])) {transmission.data[i,"tag"]<-paste0(transmission.data[i,"site"],"X",transmission.data[i,"X"]+transmission.data[i,"x"],"Y",transmission.data[i,"Y"]+transmission.data[i,"y"])}
}
transmission.data$tag<-as.factor(transmission.data$tag)

# analyze data

## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS"))
{
  set.seed(9843525)
  mod<-bam(status.next~
             te(log10(tot.spore.deposition),height.cm,sp=c(1,-1,-1))+ #When fitting the tensor, we fixed the smoothing parameter for the marginal smooth for log10 total spore deposition to 1e10. This effectively sets the shape of this marginal smooth to be linear, preventing any overfitting involving non-monotonic effects of spore deposition.
             s(mean.temp)+
             s(max.temp)+
             s(min.temp)+
             s(mean.abs.hum)+
             s(mean.daily.rain)+
             s(tag,bs="re")+
             s(site,bs="re"),
           family=binomial(link="logit"),
           select = T,
           method="fREML",
           data=transmission.data,
           control = list(nthreads=4),
           discrete = T)
  saveRDS(mod,file="~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")
}

## load model

transmission.model<-readRDS("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/transmission.model.RDS")

## model checking
par(mfrow=c(2,2))
set.seed(9843525)
transmission.data$log.10.tot.spore.deposition<-log10(transmission.data$tot.spore.deposition)
dummy.mod<-bam(status.next~
           te(log.10.tot.spore.deposition,height.cm,sp=c(1,-1,-1))+ #Fit dummy mod with hard coded log.10.tot.spore.deposition as predictor because the log10 function throws off gam.check
           s(mean.temp)+
           s(max.temp)+
           s(min.temp)+
           s(mean.abs.hum)+
           s(mean.daily.rain)+
           s(tag,bs="re")+
           s(site,bs="re"),
         family=binomial(link="logit"),
         select = T,
         method="fREML",
         data=transmission.data,
         control = list(nthreads=4),
         discrete = T)
gam.check(dummy.mod) #No issues
concurvity(transmission.model,full=F) #Concurvity is expected between all weather variables. The non-pessimistic estimations of concurvity (under $estimate and $observed) indicate that it is not a serious issue in the model fit.
