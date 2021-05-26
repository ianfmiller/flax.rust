library(mgcv)
library(lme4)
library(lmerTest)
library(progress)

# load data
source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant.growth.R")

delta.height<-subset(delta.height,time<=7)

# visualize data
par(mfrow =c(1,1))
colors <- ifelse(delta.height$height.next > delta.height$height, "green", "red")
plot(height.next~height, data = delta.height, col = colors)


## fit models--only if not already fit

if(!file.exists("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS"))
{
  ### construct all combinations of predictors
  source("~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/plant.model.set.creation.R")
  
  ### create all sets of models
  model.set <-apply(pred.mat, 1, function(x) as.formula( paste(c("height.next ~ height",predictors[x],'(1|site)'),collapse=" + ")))
  names(model.set)<-seq(1,length(model.set),1)
  
  ## run to search for best  model
  all.fit.models<-c()
  AIC.benchmark<- AIC(lmer(height.next~height+(1|site),data=delta.height,REML = F)) #cutoff to limit memory usage
  pb <- progress_bar$new(total = length(model.set),format = " fitting models [:bar] :percent eta: :eta")
  for (i in 1:length(model.set))
  {
    suppressMessages(new.mod<-lmer(model.set[[i]],data=delta.height,REML = F))
    AIC.new.mod<-AIC(new.mod)
    if(AIC.new.mod<=(AIC.benchmark)) {all.fit.models<-append(all.fit.models,list(new.mod))}
    pb$tick()
  }
  
  
  ## compare between models
  AICs<-unlist(lapply(all.fit.models,AIC))
  delta.AICs<-AICs-min(AICs)
  index<-1
  best.model<-all.fit.models[[order(AICs)[index]]]
  saveRDS(best.model,file="~/Dropbox/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")
}

## load best model

plant.growth.model<-readRDS("~/Dropbox/flax.rust/cross.scale.transmission.dynamics/models/plant.growth.model.RDS")

